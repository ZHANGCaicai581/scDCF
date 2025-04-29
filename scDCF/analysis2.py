# scDCF/analysis.py

import os
import logging
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm  # For progress bars
import anndata
import tempfile
from functools import lru_cache  # For caching function results
from scipy import sparse
from scipy.sparse import issparse
import numba  # For just-in-time compilation for numeric operations

# Configure Numba for performance
numba.config.THREADING_LAYER = 'threadsafe'

# Constants for optimization thresholds
NUMBA_SIZE_THRESHOLD = 5000  # Only use Numba for datasets larger than this
CACHE_SIZE_SMALL = 32  # Cache size for smaller functions
CACHE_SIZE_MEDIUM = 128  # Default cache size
CACHE_SIZE_LARGE = 256  # Cache size for frequently called functions

# Create a decorator for caching with size-based key and conditional application
def smart_cache(maxsize=CACHE_SIZE_MEDIUM, min_size_for_caching=1000):
    """Cache decorator with conditional application based on data size"""
    def decorator(func):
        # Create the cache
        cached_func = lru_cache(maxsize=maxsize)(func)
        
        def wrapper(*args, **kwargs):
            # Check if any argument is a large array or sparse matrix
            use_cache = True
            for arg in args:
                if isinstance(arg, np.ndarray):
                    if arg.size < min_size_for_caching:
                        use_cache = False  # Skip caching for small arrays
                elif isinstance(arg, sparse.spmatrix):
                    if arg.nnz < min_size_for_caching:
                        use_cache = False  # Skip caching for small sparse matrices
            
            if use_cache:
                return cached_func(*args, **kwargs)
            else:
                return func(*args, **kwargs)
        
        # Preserve the function's metadata
        wrapper.__name__ = func.__name__
        wrapper.__doc__ = func.__doc__
        
        return wrapper
    return decorator

# Define both Numba and non-Numba versions of core functions
@numba.jit(nopython=True, parallel=True)
def _calculate_abs_diff_numba(target_counts, reference_counts):
    """
    Numba-accelerated calculation of absolute differences.
    For large datasets only.
    """
    n_targets = len(target_counts)
    n_refs = len(reference_counts)
    differences_matrix = np.zeros((n_targets, n_refs), dtype=np.float32)
    
    for i in numba.prange(n_targets):
        target_value = target_counts[i]
        for j in range(n_refs):
            differences_matrix[i, j] = abs(target_value - reference_counts[j])
            
    return differences_matrix

def _calculate_abs_diff_numpy(target_counts, reference_counts):
    """
    Vectorized NumPy implementation of absolute differences.
    More efficient for smaller datasets.
    """
    return np.abs(target_counts[:, np.newaxis] - reference_counts)

def calculate_absolute_differences(target_counts, reference_counts):
    """
    Intelligently choose the best implementation based on data size.
    """
    # Use Numba only for larger datasets where it provides real benefit
    if len(target_counts) * len(reference_counts) > NUMBA_SIZE_THRESHOLD:
        return _calculate_abs_diff_numba(target_counts, reference_counts)
    else:
        return _calculate_abs_diff_numpy(target_counts, reference_counts)

def get_nearest_cells(target_cells, reference_cells, rna_count_column, n_samples=100):
    """
    Find the nearest reference cells for each target cell based on RNA counts.
    Optimized with vectorized operations and conditional JIT compilation.
    """
    if rna_count_column not in target_cells.obs.columns or rna_count_column not in reference_cells.obs.columns:
        raise ValueError(f"Column '{rna_count_column}' not found in AnnData object's obs.")

    # Extract counts as NumPy arrays for faster computation
    target_counts = target_cells.obs[rna_count_column].values.astype(np.float32)
    reference_counts = reference_cells.obs[rna_count_column].values.astype(np.float32)
    reference_cell_ids = reference_cells.obs_names.values

    # Pre-allocate results dictionary
    results = {}
    
    # Calculate differences matrix using the appropriate method
    differences_matrix = calculate_absolute_differences(target_counts, reference_counts)
    
    # Process all target cells
    for i, target_cell_id in enumerate(target_cells.obs_names):
        # Get distances for this target cell
        differences = differences_matrix[i]
        
        # Get indices of smallest differences based on dataset size
        if len(reference_counts) <= 1000:
            nearest_indices = np.argsort(differences)
        else:
            # For large arrays, use partial sort which is much faster
            nearest_indices = np.argpartition(differences, min(1000, len(differences)-1))[:1000]
            # Sort just these indices
            nearest_indices = nearest_indices[np.argsort(differences[nearest_indices])]
        
        # Sample indices without replacement
        if n_samples >= len(nearest_indices):
            sampled_indices = nearest_indices
        else:
            sampled_indices = np.random.choice(nearest_indices, n_samples, replace=False)

        # Convert to cell IDs
        results[target_cell_id] = reference_cell_ids[sampled_indices]

    logging.info(f"Nearest cells determined for {len(results)} target cells.")
    return results

# Cache the preprocessing step (high computational impact)
@smart_cache(maxsize=CACHE_SIZE_LARGE)
def preprocess_gene_data(adata_var_names, significant_genes, control_genes=None):
    """
    Pre-filter and prepare gene data for processing.
    This function is cached because it's repeatedly called with the same arguments.
    
    Args:
        adata_var_names: Gene names in the AnnData object
        significant_genes: List of significant gene names
        control_genes: Dictionary of control genes
        
    Returns:
        Tuple of (valid_genes, valid_weights, control_genes_filtered)
    """
    # Convert gene names to set for O(1) lookup
    adata_genes_set = set(adata_var_names)
    
    # Filter significant genes that exist in the dataset
    valid_genes = []
    valid_weights = []
    
    for gene, weight in significant_genes:
        if gene in adata_genes_set:
            valid_genes.append(gene)
            valid_weights.append(weight)
    
    # Filter control genes that exist in the dataset
    control_genes_filtered = {}
    if control_genes:
        for gene in valid_genes:
            if gene in control_genes:
                # Get control genes that exist in the dataset
                valid_controls = [g for g in control_genes[gene] if g in adata_genes_set]
                if valid_controls:
                    control_genes_filtered[gene] = valid_controls
    
    return valid_genes, valid_weights, control_genes_filtered

# Optimized function for extracting expression values
def extract_expression_values(matrix, indices, is_sparse=None):
    """
    Extract values from expression matrix efficiently for both sparse and dense formats.
    
    Args:
        matrix: Expression matrix (X) from AnnData
        indices: Indices to extract
        is_sparse: Whether matrix is sparse (if None, will detect)
        
    Returns:
        Array of expression values
    """
    # Auto-detect if not specified
    if is_sparse is None:
        is_sparse = issparse(matrix)
    
    if is_sparse:
        # For sparse matrices, efficiently extract values
        if len(indices) == 1:
            # Single index access
            result = matrix[indices[0]].toarray().flatten()
        else:
            # Multiple indices - use efficient sparse slicing
            result = matrix[indices].toarray()
    else:
        # For dense matrices
        if len(indices) == 1:
            result = matrix[indices[0]]
        else:
            result = matrix[indices]
    
    return result

def monte_carlo_comparison(
    adata,
    cell_type,
    significant_genes_df,
    disease_control_genes,
    healthy_control_genes,
    output_dir,
    rna_count_column='nCount_RNA',
    iterations=100,
    target_group="disease",
    disease_marker="disease_numeric",
    disease_value=1,
    healthy_value=0,
    show_progress=False,
    chunk_size=None,
    max_memory=None
):
    """
    Perform Monte Carlo comparison to find cells with significant differences.
    
    Args:
        adata: AnnData object with expression data
        cell_type: Cell type to analyze
        significant_genes_df: DataFrame with significant genes and weights
        disease_control_genes: Dictionary of control genes for disease
        healthy_control_genes: Dictionary of control genes for healthy
        output_dir: Directory to save output
        rna_count_column: Column name with RNA count
        iterations: Number of iterations to run
        target_group: Which group to target ('disease' or 'healthy')
        disease_marker: Column in obs with disease status
        disease_value: Value indicating disease
        healthy_value: Value indicating healthy
        show_progress: Whether to show progress bars
        chunk_size: Size of chunks for processing (for memory management)
        max_memory: Maximum memory usage in GB (None for no limit)
    
    Returns:
        DataFrame with results
    """
    try:
        logging.info(f"Starting Monte Carlo comparison for cell type: {cell_type}, target group: {target_group}")

        # Check if required columns exist
        required_columns = [disease_marker, rna_count_column]
        if not all(col in adata.obs.columns for col in required_columns):
            missing = [col for col in required_columns if col not in adata.obs.columns]
            raise ValueError(f"Missing required columns in AnnData object: {missing}")
        
        # Normalize gene dataframe 
        if 'gene_name' not in significant_genes_df.columns:
            # Try to find a column containing gene names
            gene_cols = [col for col in significant_genes_df.columns if 'gene' in col.lower()]
            if gene_cols:
                significant_genes_df = significant_genes_df.rename(columns={gene_cols[0]: 'gene_name'})
            else:
                raise ValueError("Cannot find gene name column in significant_genes_df")
        
        # Normalize zstat column
        weight_cols = [col for col in significant_genes_df.columns 
                  if col != 'gene_name' and significant_genes_df[col].dtype == float]
        
        if not weight_cols:
            raise ValueError("Cannot find weight column in significant_genes_df")
        
        # Use the first weight column
        weight_col = weight_cols[0]
        if weight_col != 'zstat':
            significant_genes_df = significant_genes_df.rename(columns={weight_col: 'zstat'})
        
        logging.info(f"Columns in significant_genes_df after normalization: {list(significant_genes_df.columns)}")
        logging.info(f"First few rows of significant_genes_df:\n{significant_genes_df.head()}")

        # Create cell type directory
        cell_type_dir = os.path.join(output_dir, cell_type)
        os.makedirs(cell_type_dir, exist_ok=True)
        
        # Select the appropriate control genes based on target group
        control_genes = disease_control_genes if target_group == "disease" else healthy_control_genes
        
        # Filter anndata object for the specific cell type
        if cell_type not in adata.obs['celltype_major'].unique():
            raise ValueError(f"Cell type '{cell_type}' not found in adata")
        
        adata_ct = adata[adata.obs['celltype_major'] == cell_type].copy()
        
        # Identify target and reference cells based on condition
        if target_group == "disease":
            target_cells = adata_ct[adata_ct.obs[disease_marker] == disease_value].copy()
            reference_cells = adata_ct[adata_ct.obs[disease_marker] == healthy_value].copy()
        else:  # target_group == "healthy"
            target_cells = adata_ct[adata_ct.obs[disease_marker] == healthy_value].copy()
            reference_cells = adata_ct[adata_ct.obs[disease_marker] == disease_value].copy()

        logging.info(f"{len(target_cells)} target cells and {len(reference_cells)} reference cells identified for {target_group} group.")

        if len(target_cells) == 0 or len(reference_cells) == 0:
            logging.warning(f"No {'target' if len(target_cells) == 0 else 'reference'} cells found for {cell_type}, {target_group}.")
            return pd.DataFrame()  # Return empty dataframe
        
        # Determine nearest reference cells for each target cell
        nearest_cells = get_nearest_cells(
            target_cells=target_cells,
            reference_cells=reference_cells,
            rna_count_column=rna_count_column,
            n_samples=100  # Number of nearest reference cells to sample
        )
        
        # Prepare for iterations
        all_iterations_results = []

        # Check if expression matrix is sparse to optimize operations
        is_sparse_matrix = issparse(adata.X)
        
        # Create gene pairs with weights
        significant_genes_pairs = list(zip(
            significant_genes_df['gene_name'], 
            significant_genes_df['zstat']
        ))
        
        # Pre-filter gene data once to avoid repeated filtering
        valid_genes, valid_weights, control_genes_filtered = preprocess_gene_data(
            adata.var_names, 
            significant_genes_pairs,
            control_genes
        )
        
        # Get indices for faster access
        valid_gene_indices = [adata.var_names.get_loc(gene) for gene in valid_genes]
        target_cell_indices = {cell: i for i, cell in enumerate(target_cells.obs_names)}
        
        # Progress bar context manager
        prog_context = tqdm(total=iterations, desc="Iterations") if show_progress else nullcontext()
        
        # Run iterations
        with prog_context as prog:
            for iteration in range(iterations):
                if show_progress:
                    prog.update(1)
                    
                all_results = []
                
                # Calculate how many cells to process in each chunk
                if chunk_size:
                    cell_chunks = [list(nearest_cells.keys())[i:i+chunk_size] 
                                  for i in range(0, len(nearest_cells), chunk_size)]
                else:
                    cell_chunks = [list(nearest_cells.keys())]
                
                # Process cells in chunks for memory efficiency
                for chunk in cell_chunks:
                    # Process each target cell
                    for idx in chunk:
                        # Get the nearest reference cells
                        nearest_reference_cells = nearest_cells[idx]
                        
                        # Sample a random reference cell from the nearest ones
                        reference_cell = np.random.choice(nearest_reference_cells)
                        reference_idx = np.where(reference_cells.obs_names == reference_cell)[0][0]
                        
                        # Initialize arrays for differences
                        sig_diffs = np.zeros(len(valid_genes), dtype=np.float32)
                        ctrl_diffs = []
                        
                        # Get target and reference expression values for all genes at once
                        # This is much more efficient than accessing genes one by one
                        target_cell_expression = {}
                        reference_cell_expression = {}
                        
                        # Loop through each gene once
                        for i, gene in enumerate(valid_genes):
                            # Get gene indices for faster access
                            gene_idx = valid_gene_indices[i]
                            weight = valid_weights[i]
                            
                            # Get expression values
                            if gene not in target_cell_expression:
                                if is_sparse_matrix:
                                    val = target_cells.X[target_cell_indices[idx], gene_idx]
                                    target_cell_expression[gene] = val.toarray()[0, 0] if issparse(val) else val
                                else:
                                    target_cell_expression[gene] = target_cells.X[target_cell_indices[idx], gene_idx]
                                    
                            if gene not in reference_cell_expression:
                                if is_sparse_matrix:
                                    val = reference_cells.X[reference_idx, gene_idx]
                                    reference_cell_expression[gene] = val.toarray()[0, 0] if issparse(val) else val
                                else:
                                    reference_cell_expression[gene] = reference_cells.X[reference_idx, gene_idx]
                            
                            # Calculate the weighted difference for the target gene
                            target_value = target_cell_expression[gene]
                            reference_value = reference_cell_expression[gene]
                            sig_diffs[i] = weight * abs(target_value - reference_value)

                            # For control genes, use control genes for this target gene
                            if gene in control_genes_filtered:
                                controls = control_genes_filtered[gene]
                                if controls:
                                    # Random choice from valid controls
                                    ctrl_gene = np.random.choice(controls)
                                    # Check if control gene indices exists
                                    ctrl_idx = adata.var_names.get_loc(ctrl_gene) if ctrl_gene in adata.var_names else None
                                    if ctrl_idx is not None:
                                        # Get control gene expression values
                                        if ctrl_gene not in target_cell_expression:
                                            if is_sparse_matrix:
                                                val = target_cells.X[target_cell_indices[idx], ctrl_idx]
                                                target_cell_expression[ctrl_gene] = val.toarray()[0, 0] if issparse(val) else val
                                            else:
                                                target_cell_expression[ctrl_gene] = target_cells.X[target_cell_indices[idx], ctrl_idx]
                                                
                                        if ctrl_gene not in reference_cell_expression:
                                            if is_sparse_matrix:
                                                val = reference_cells.X[reference_idx, ctrl_idx]
                                                reference_cell_expression[ctrl_gene] = val.toarray()[0, 0] if issparse(val) else val
                                            else:
                                                reference_cell_expression[ctrl_gene] = reference_cells.X[reference_idx, ctrl_idx]
                                        
                                        ctrl_target_value = target_cell_expression[ctrl_gene]
                                        ctrl_reference_value = reference_cell_expression[ctrl_gene]
                                        
                            ctrl_diff = weight * abs(ctrl_target_value - ctrl_reference_value)
                            ctrl_diffs.append(ctrl_diff)

                        # Skip cells with insufficient data
                        if len(sig_diffs) == 0 or len(ctrl_diffs) == 0:
                            continue

                        # Convert control differences to array for faster calculations
                        ctrl_diffs = np.array(ctrl_diffs, dtype=np.float32)

                        # Remove zeros and NaNs for accurate t-test
                        sig_diffs = sig_diffs[~np.isnan(sig_diffs) & (sig_diffs != 0)]
                        ctrl_diffs = ctrl_diffs[~np.isnan(ctrl_diffs) & (ctrl_diffs != 0)]
                        
                        if len(sig_diffs) == 0 or len(ctrl_diffs) == 0:
                            continue

                        # Perform t-test efficiently
                        t_stat, p_val = ttest_ind(sig_diffs, ctrl_diffs, equal_var=False)
                        
                        # Calculate totals using numpy sum for better performance
                        total_sig_diff = np.sum(sig_diffs)
                        total_ctrl_diff = np.sum(ctrl_diffs)

                        # Calculate one-tailed p-value
                        if total_sig_diff > total_ctrl_diff:
                            p_val_one_tailed = p_val / 2
                        else:
                            p_val_one_tailed = 1 - (p_val / 2)

                        # Store results
                        all_results.append({
                            'cell_id': idx,
                            't_stat': t_stat,
                            'p_value': p_val_one_tailed,
                            'sig_diff': total_sig_diff,
                            'ctrl_diff': total_ctrl_diff,
                            'significant': p_val_one_tailed < 0.05,
                            'iteration': iteration + 1,
                            'target_group': target_group
                        })

                if not all_results:
                    logging.warning(f"No results generated for iteration {iteration + 1} in {cell_type} ({target_group}).")
                    continue

                # Create DataFrame and adjust p-values
                results_df = pd.DataFrame(all_results)
                
                # Apply FDR correction
                results_df['p_value_adj'] = multipletests(results_df['p_value'], method='fdr_bh')[1]

                # Save iteration file
                iteration_file = os.path.join(cell_type_dir, f"{cell_type}_{target_group}_monte_carlo_results_iteration{iteration + 1}.csv")
                results_df.to_csv(iteration_file, index=False)
                logging.info(f"Monte Carlo comparison results saved for iteration {iteration + 1}")

                all_iterations_results.append(results_df)
        
        # Combine results from all iterations
        if all_iterations_results:
            combined_results = pd.concat(all_iterations_results, ignore_index=True)
            combined_output_file = os.path.join(cell_type_dir, f"{cell_type}_{target_group}_monte_carlo_results.csv")
            combined_results.to_csv(combined_output_file, index=False)
            logging.info(f"Combined results saved to {combined_output_file}")
            return combined_results
        else:
            logging.warning(f"No results generated for {cell_type} ({target_group}) in any iteration.")
            return pd.DataFrame()
            
    except Exception as e:
        logging.error(f"Unexpected error in monte_carlo_comparison: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return pd.DataFrame()

# Simple context manager for progress bars
class nullcontext:
    def __enter__(self): return None
    def __exit__(self, *excinfo): pass

def compare_groups(disease_df, healthy_df):
    logging.info("Comparing results between disease and healthy groups.")

    if disease_df.empty or healthy_df.empty:
        logging.warning("One of the result DataFrames is empty. Cannot perform comparison.")
        return {}

    t_stat_pval, t_pval_pval = ttest_ind(disease_df['p_value'], healthy_df['p_value'], equal_var=False, nan_policy='omit')

    comparison_results = {
        't_stat_pval': t_stat_pval,
        't_pval_pval': t_pval_pval,
    }

    logging.info(f"Comparison completed: {comparison_results}")
    return comparison_results
