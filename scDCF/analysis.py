"""
Ultra-optimized Monte Carlo analysis with Quick Wins optimizations.
This implementation follows the current scDCF workflow while keeping the
optimized computation path.
"""

import os
import logging
import time
import numpy as np
import pandas as pd
from scipy.stats import ttest_1samp, ttest_ind
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import contextlib
from scipy.sparse import issparse

# Constants
DEFAULT_N_NEAREST_CELLS = 1000
DEFAULT_N_SAMPLES_PER_ITER = 100
DEFAULT_N_CONTROL_GENES = 10
BATCH_SIZE = 500

@contextlib.contextmanager
def nullcontext():
    yield


def _to_one_sided_pvalue(t_stat, p_two_sided):
    """Convert a two-sided t-test p-value to the one-sided alternative mean>0."""
    if np.isnan(t_stat) or np.isnan(p_two_sided):
        return np.nan
    return p_two_sided / 2 if t_stat > 0 else 1 - (p_two_sided / 2)

def _validate_monte_carlo_inputs(adata, cell_type, cell_type_column, 
                                  significant_genes_df, disease_marker, 
                                  rna_count_column):
    """Validate inputs - same as before"""
    import anndata as ad
    
    if not isinstance(adata, ad.AnnData):
        raise TypeError(f"adata must be AnnData object, got {type(adata)}")
    
    required_cols = [cell_type_column, disease_marker, rna_count_column]
    missing_cols = [col for col in required_cols if col not in adata.obs.columns]
    
    if missing_cols:
        available = list(adata.obs.columns[:10])
        raise ValueError(
            f"Missing required columns in adata.obs: {missing_cols}\n"
            f"Available columns (first 10): {available}"
        )
    
    if cell_type not in adata.obs[cell_type_column].unique():
        available_types = list(adata.obs[cell_type_column].unique())
        raise ValueError(
            f"Cell type '{cell_type}' not found in column '{cell_type_column}'.\n"
            f"Available cell types: {available_types}"
        )
    
    if 'gene_name' not in significant_genes_df.columns:
        raise ValueError(
            f"significant_genes_df must have 'gene_name' column.\n"
            f"Found columns: {list(significant_genes_df.columns)}"
        )
    
    genes_in_adata = set(adata.var_names)
    genes_in_df = set(significant_genes_df['gene_name'])
    overlap = genes_in_adata & genes_in_df
    
    if len(overlap) == 0:
        raise ValueError(
            f"No gene overlap! Check gene naming conventions.\n"
            f"AnnData genes (first 10): {list(genes_in_adata)[:10]}\n"
            f"Input genes (first 10): {list(genes_in_df)[:10]}"
        )
    
    if len(overlap) < 10:
        logging.warning(f"Only {len(overlap)} genes overlap. Results may be unreliable.")
    
    ct_mask = adata.obs[cell_type_column] == cell_type
    n_cells = ct_mask.sum()
    
    if n_cells < 100:
        raise ValueError(
            f"Insufficient cells for {cell_type}: only {n_cells} found. "
            "Need at least 100 cells for reliable analysis."
        )
    
    logging.info("✅ Input validation passed")

def get_nearest_cells(target_cells, reference_cells, rna_count_column,
                     n_pool=DEFAULT_N_NEAREST_CELLS, exclude_self=False):
    """
    Build per-target pools of nearest reference cells by library size.

    The returned pools are reused across Monte Carlo iterations, while the
    actual reference subset is resampled from each pool at every iteration.
    """
    if rna_count_column not in target_cells.obs.columns:
        raise ValueError(f"Column '{rna_count_column}' not found in target_cells.obs")
    if rna_count_column not in reference_cells.obs.columns:
        raise ValueError(f"Column '{rna_count_column}' not found in reference_cells.obs")

    target_counts = pd.to_numeric(target_cells.obs[rna_count_column], errors='coerce').to_numpy()
    reference_counts = pd.to_numeric(reference_cells.obs[rna_count_column], errors='coerce').to_numpy()
    reference_cell_ids = reference_cells.obs_names.values

    if np.isnan(target_counts).any() or np.isnan(reference_counts).any():
        raise ValueError(
            f"Column '{rna_count_column}' contains non-numeric or missing values."
        )

    results = {}
    max_pool = min(int(n_pool), len(reference_counts))

    for i, target_cell_id in enumerate(target_cells.obs_names):
        eligible_mask = np.ones(len(reference_counts), dtype=bool)
        if exclude_self:
            eligible_mask &= (reference_cell_ids != target_cell_id)

        eligible_indices = np.flatnonzero(eligible_mask)
        if eligible_indices.size == 0:
            raise ValueError(
                f"No eligible reference cells remain for target cell '{target_cell_id}'."
            )

        if eligible_indices.size <= max_pool:
            nearest_indices = eligible_indices
        else:
            eligible_differences = np.abs(reference_counts[eligible_indices] - target_counts[i])
            kth = max_pool - 1
            nearest_local = np.argpartition(eligible_differences, kth)[:max_pool]
            nearest_indices = eligible_indices[nearest_local]

        results[str(target_cell_id)] = nearest_indices.astype(int, copy=False)

    pool_sizes = np.fromiter((pool.size for pool in results.values()), dtype=int)
    logging.info(
        f"Reference pools determined for {len(results)} target cells "
        f"(pool size range {pool_sizes.min()}-{pool_sizes.max()}, exclude_self={exclude_self})."
    )
    return results

# ============================================================================
# QUICK WIN OPTIMIZATION 1: Pre-extract Expression Matrices
# ============================================================================
def _preextract_expression_matrices(target_cells, reference_cells, valid_genes, 
                                    control_genes_filtered, adata_var_names):
    """
    OPTIMIZATION: Extract all expression data ONCE before iterations.
    Instead of extracting 4.5M times, extract once and reuse.
    
    GUARANTEE: Same data, just extracted more efficiently.
    """
    logging.info("Pre-extracting expression matrices (one-time cost)...")
    import time
    start = time.time()
    
    # Get gene indices for valid genes
    gene_indices = np.array([adata_var_names.get_loc(gene) for gene in valid_genes])
    
    # Extract target cell expressions for all valid genes
    if issparse(target_cells.X):
        target_expr = target_cells.X[:, gene_indices].tocsr()
    else:
        target_expr = np.asarray(target_cells.X[:, gene_indices])
    
    # Extract reference cell expressions
    if issparse(reference_cells.X):
        ref_expr = reference_cells.X[:, gene_indices].tocsr()
    else:
        ref_expr = np.asarray(reference_cells.X[:, gene_indices])
    
    # Pre-extract control gene expressions
    all_ctrl_genes = set()
    for gene in valid_genes:
        if gene in control_genes_filtered:
            all_ctrl_genes.update(control_genes_filtered[gene])
    
    all_ctrl_genes = sorted(all_ctrl_genes)
    ctrl_gene_indices = np.array([adata_var_names.get_loc(g) for g in all_ctrl_genes if g in adata_var_names])
    
    if len(ctrl_gene_indices) > 0:
        if issparse(target_cells.X):
            target_ctrl_expr = target_cells.X[:, ctrl_gene_indices].tocsr()
            ref_ctrl_expr = reference_cells.X[:, ctrl_gene_indices].tocsr()
        else:
            target_ctrl_expr = np.asarray(target_cells.X[:, ctrl_gene_indices])
            ref_ctrl_expr = np.asarray(reference_cells.X[:, ctrl_gene_indices])
    else:
        target_ctrl_expr = None
        ref_ctrl_expr = None
    
    ctrl_gene_to_col = {gene: i for i, gene in enumerate(all_ctrl_genes)}
    
    elapsed = time.time() - start
    logging.info(f"✅ Expression matrices pre-extracted in {elapsed:.2f}s")
    logging.info(f"   Target shape: {target_expr.shape}")
    logging.info(f"   Reference shape: {ref_expr.shape}")
    
    return (target_expr, ref_expr, target_ctrl_expr, ref_ctrl_expr,
            ctrl_gene_to_col, all_ctrl_genes)

def monte_carlo_comparison_optimized(adata, cell_type, cell_type_column, significant_genes_df, 
                                    disease_control_genes=None, healthy_control_genes=None, 
                                    output_dir=".", rna_count_column='nCount_RNA', 
                                    iterations=10, target_group="disease", 
                                    disease_marker='disease_numeric', 
                                    disease_value=1, healthy_value=0, show_progress=False,
                                    batch_size=BATCH_SIZE, random_seed=None,
                                    write_iteration_files=False,
                                    write_combined_output=False):
    """
    Optimized Monte Carlo comparison for per-cell scDCF testing.

    The computation is vectorized where possible, but the statistical target
    remains the current manuscript definition: weighted gene-level
    difference-of-differences tested against zero within each iteration.
    """
    import time
    total_start = time.time()
    
    try:
        # Validation (same as before)
        _validate_monte_carlo_inputs(
            adata, cell_type, cell_type_column, 
            significant_genes_df, disease_marker, rna_count_column
        )

        if target_group not in {"disease", "healthy"}:
            raise ValueError("target_group must be either 'disease' or 'healthy'")
        
        logging.info(f"Starting OPTIMIZED Monte Carlo for {cell_type}, {target_group}, {iterations} iterations")
        rng = np.random.default_rng(random_seed)

        # String-safe comparison (same as before)
        def value_equals(a, b):
            if isinstance(a, bool) or isinstance(b, bool):
                return bool(a) == bool(b)
            return str(a) == str(b)
        
        # Filter cell type (same as before)
        if cell_type and cell_type_column in adata.obs.columns:
            adata_subset = adata[adata.obs[cell_type_column] == cell_type].copy()
        else:
            adata_subset = adata.copy()
        
        # Normalize columns (same as before)
        significant_genes_df = significant_genes_df.copy()
        significant_genes_df.columns = significant_genes_df.columns.str.lower().str.strip()
        
        if 'zstat' not in significant_genes_df.columns or 'gene_name' not in significant_genes_df.columns:
            raise ValueError(
                "significant_genes_df must have 'zstat' and 'gene_name' columns.\n"
                f"Found: {list(significant_genes_df.columns)}"
            )
        
        # Split groups (same as before)
        exclude_self_reference = False
        reference_group = "healthy"
        if target_group == "disease":
            target_mask = [value_equals(v, disease_value) for v in adata_subset.obs[disease_marker]]
            target_cells = adata_subset[target_mask].copy()
            ref_mask = [value_equals(v, healthy_value) for v in adata_subset.obs[disease_marker]]
            reference_cells = adata_subset[ref_mask].copy()
            control_genes = disease_control_genes
        else:
            target_mask = [value_equals(v, healthy_value) for v in adata_subset.obs[disease_marker]]
            target_cells = adata_subset[target_mask].copy()
            ref_mask = [value_equals(v, healthy_value) for v in adata_subset.obs[disease_marker]]
            reference_cells = adata_subset[ref_mask].copy()
            control_genes = healthy_control_genes
            exclude_self_reference = True

        logging.info(
            f"Target group '{target_group}': {len(target_cells)} target cells; "
            f"reference group '{reference_group}': {len(reference_cells)} cells "
            f"(self-exclusion={exclude_self_reference})"
        )

        if len(target_cells) == 0:
            raise ValueError(f"No target cells found for {cell_type}, {target_group}")
        if len(reference_cells) == 0:
            raise ValueError(f"No reference cells found for {cell_type}, {target_group}")
        if exclude_self_reference and len(reference_cells) < 2:
            raise ValueError(
                f"Healthy target analysis for {cell_type} requires at least 2 healthy cells "
                "so each target cell has a non-self healthy reference."
            )
        
        # Create output directory (same as before)
        cell_type_dir = os.path.join(output_dir, cell_type)
        os.makedirs(cell_type_dir, exist_ok=True)
        
        # Get nearest cells (same algorithm)
        matched_indices = get_nearest_cells(
            target_cells, reference_cells, rna_count_column,
            n_pool=DEFAULT_N_NEAREST_CELLS,
            exclude_self=exclude_self_reference
        )
        
        # Filter valid genes (same as before)
        valid_genes = [
            gene for gene in significant_genes_df['gene_name'] 
            if gene in adata.var_names
        ]
        
        if not valid_genes:
            raise ValueError("No valid genes found in dataset after filtering")
        
        gene_weights = np.maximum(
            significant_genes_df[
                significant_genes_df['gene_name'].isin(valid_genes)
            ]['zstat'].values,
            0.0
        )
        if np.sum(gene_weights) > 0:
            gene_weights = gene_weights / np.sum(gene_weights)
        else:
            logging.warning(
                "All z-stat weights are non-positive after filtering; falling back to uniform weights."
            )
            gene_weights = np.ones(len(valid_genes), dtype=float) / len(valid_genes)
        
        target_idx_map = {cell_id: i for i, cell_id in enumerate(target_cells.obs_names)}
        
        # Filter control genes (same as before)
        control_genes_filtered = {}
        if control_genes:
            for gene in valid_genes:
                if gene in control_genes:
                    control_genes_filtered[gene] = [
                        ctrl for ctrl in control_genes[gene] 
                        if ctrl in adata.var_names
                    ]
        
        # ===== OPTIMIZATION 1: Pre-extract expression matrices =====
        (target_expr_all, ref_expr_all, target_ctrl_expr, ref_ctrl_expr,
         ctrl_gene_to_col, all_ctrl_genes) = _preextract_expression_matrices(
            target_cells, reference_cells, valid_genes, 
            control_genes_filtered, adata.var_names
        )
        
        # ===== OPTIMIZATION 4: Pre-compute control gene mappings =====
        # Map each significant gene to its control gene columns
        gene_ctrl_map = {}
        for gene_idx, gene in enumerate(valid_genes):
                if gene in control_genes_filtered and control_genes_filtered[gene]:
                    # Map to column indices in ctrl expression matrix
                    mapped = [
                        ctrl_gene_to_col[ctrl_gene]
                        for ctrl_gene in control_genes_filtered[gene]
                        if ctrl_gene in ctrl_gene_to_col
                    ]
                    if mapped:
                        gene_ctrl_map[gene_idx] = np.array(mapped, dtype=int)

        if not gene_ctrl_map:
            raise ValueError(
                f"No usable control-gene matches remain for {cell_type} ({target_group}) after filtering."
            )
        
        logging.info(f"✅ Pre-computation complete. Starting iterations...")
        
        # Run iterations
        iteration_files = []
        all_iteration_frames = []
        prog_context = tqdm(total=iterations, desc="Iterations") if show_progress else nullcontext()
        
        with prog_context as prog:
            for iteration in range(iterations):
                if show_progress:
                    prog.update(1)
                
                iter_start = time.time()
                logging.info(f"Iteration {iteration + 1}/{iterations}")
                all_results = []
                
                # Process cells in batches
                cell_ids = list(matched_indices.keys())
                
                for batch_start in range(0, len(cell_ids), batch_size):
                    batch_end = min(batch_start + batch_size, len(cell_ids))
                    batch_cells = cell_ids[batch_start:batch_end]
                    
                    for idx in batch_cells:
                        reference_pool = matched_indices[idx]
                        if reference_pool.size <= DEFAULT_N_SAMPLES_PER_ITER:
                            ref_indices = reference_pool
                        else:
                            ref_indices = rng.choice(
                                reference_pool,
                                DEFAULT_N_SAMPLES_PER_ITER,
                                replace=False
                            )

                        target_idx = target_idx_map[idx]
                        
                        # ===== OPTIMIZATION 1 & 2: Use pre-extracted matrices =====
                        # Get expressions from pre-extracted arrays (FAST!)
                        if issparse(target_expr_all):
                            target_expr = target_expr_all.getrow(target_idx).toarray().ravel()
                        else:
                            target_expr = target_expr_all[target_idx]
                        
                        if issparse(ref_expr_all):
                            ref_expr_mean = np.asarray(ref_expr_all[ref_indices].mean(axis=0)).ravel()
                        else:
                            ref_expr_mean = ref_expr_all[ref_indices].mean(axis=0)
                        
                        # Compute weighted gene-level deltas:
                        # w_g * (|trait deviation| - |control deviation|)
                        sig_abs_diffs = np.abs(target_expr - ref_expr_mean)
                        delta_diffs = []
                        ctrl_abs_diffs = []
                        if gene_ctrl_map and target_ctrl_expr is not None and ref_ctrl_expr is not None:
                            if issparse(target_ctrl_expr):
                                target_ctrl_row = target_ctrl_expr.getrow(target_idx).toarray().ravel()
                            else:
                                target_ctrl_row = target_ctrl_expr[target_idx]
                            
                            if issparse(ref_ctrl_expr):
                                ref_ctrl_mean_vec = np.asarray(ref_ctrl_expr[ref_indices].mean(axis=0)).ravel()
                            else:
                                ref_ctrl_mean_vec = ref_ctrl_expr[ref_indices].mean(axis=0)
                            for gene_idx, ctrl_cols in gene_ctrl_map.items():
                                if ctrl_cols.size == 0:
                                    continue
                                ctrl_col_idx = rng.choice(ctrl_cols)
                                target_ctrl_val = target_ctrl_row[ctrl_col_idx]
                                ref_ctrl_val = ref_ctrl_mean_vec[ctrl_col_idx]
                                trait_component = gene_weights[gene_idx] * sig_abs_diffs[gene_idx]
                                ctrl_component = gene_weights[gene_idx] * abs(target_ctrl_val - ref_ctrl_val)
                                delta_diffs.append(trait_component - ctrl_component)
                                ctrl_abs_diffs.append(ctrl_component)

                        if len(delta_diffs) < 2:
                            continue

                        delta_diffs = np.array(delta_diffs, dtype=float)
                        ctrl_abs_diffs = np.array(ctrl_abs_diffs, dtype=float)
                        delta_diffs = delta_diffs[~np.isnan(delta_diffs)]

                        if len(delta_diffs) < 2:
                            continue

                        t_stat, p_val_two_sided = ttest_1samp(delta_diffs, popmean=0.0, nan_policy='omit')
                        p_val_one_tailed = _to_one_sided_pvalue(t_stat, p_val_two_sided)
                        if np.isnan(p_val_one_tailed):
                            continue

                        total_sig_diff = float(np.sum(gene_weights * sig_abs_diffs))
                        total_ctrl_diff = float(np.sum(ctrl_abs_diffs))
                        mean_delta_diff = float(np.mean(delta_diffs))

                        all_results.append({
                            'cell_id': idx,
                            't_stat': t_stat,
                            'p_value': p_val_one_tailed,
                            'sig_diff': total_sig_diff,
                            'ctrl_diff': total_ctrl_diff,
                            'delta_mean': mean_delta_diff,
                            'n_gene_pairs_tested': len(delta_diffs),
                            'significant': p_val_one_tailed < 0.05,
                            'reference_group': reference_group,
                            'reference_pool_size': int(reference_pool.size),
                            'reference_sample_size': int(len(ref_indices)),
                            'reference_self_excluded': bool(exclude_self_reference),
                            'iteration': iteration + 1,
                            'target_group': target_group
                        })
                
                if not all_results:
                    logging.warning(f"No results for iteration {iteration + 1}")
                    continue
                
                # FDR correction (same as before)
                results_df = pd.DataFrame(all_results)
                results_df['p_value_adj'] = multipletests(
                    results_df['p_value'], method='fdr_bh'
                )[1]
                all_iteration_frames.append(results_df.copy())
                
                # Save iteration (same as before)
                if write_iteration_files:
                    iteration_file = os.path.join(
                        cell_type_dir, 
                        f"{cell_type}_{target_group}_monte_carlo_results_iteration{iteration + 1}.csv"
                    )
                    results_df.to_csv(iteration_file, index=False)
                    iteration_files.append(iteration_file)
                
                iter_time = time.time() - iter_start
                logging.info(f"Iteration {iteration + 1} complete in {iter_time:.1f}s")
                
                del results_df, all_results
        
        # Combine results (same as before)
        if all_iteration_frames:
            if iteration_files:
                logging.info(f"Combining {len(iteration_files)} iteration files...")
            else:
                logging.info(f"Combining {len(all_iteration_frames)} iteration results in memory...")

            combined_results = pd.concat(all_iteration_frames, ignore_index=True)
            if write_combined_output:
                combined_output_file = os.path.join(
                    cell_type_dir, 
                    f"{cell_type}_{target_group}_monte_carlo_results.csv"
                )
                combined_results.to_csv(combined_output_file, index=False)
            
            total_time = time.time() - total_start
            logging.info(f"✅ Analysis complete in {total_time/60:.2f} minutes")
            if write_combined_output:
                logging.info(f"Combined results saved to {combined_output_file}")
            return combined_results
        else:
            logging.warning(f"No results generated for {cell_type} ({target_group})")
            return pd.DataFrame()
    
    except ValueError as e:
        logging.error(f"Invalid input: {e}")
        raise
    except MemoryError as e:
        logging.error(f"Out of memory: {e}\nTry reducing batch_size")
        raise
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        import traceback
        logging.error(traceback.format_exc())
        raise

def compare_groups(disease_df, healthy_df):
    """Compare groups - same as before"""
    logging.info("Comparing results between disease and healthy groups.")

    if disease_df.empty or healthy_df.empty:
        logging.warning("One of the result DataFrames is empty. Cannot perform comparison.")
        return {}

    t_stat_pval, t_pval_pval = ttest_ind(
        disease_df['p_value'], healthy_df['p_value'], 
        equal_var=False, nan_policy='omit'
    )

    comparison_results = {
        't_stat_pval': t_stat_pval,
        't_pval_pval': t_pval_pval,
    }

    logging.info(f"Comparison completed: {comparison_results}")
    return comparison_results

# Alias for backward compatibility
monte_carlo_comparison = monte_carlo_comparison_optimized
