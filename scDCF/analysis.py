import os
import logging
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import contextlib
from scipy.sparse import issparse

# Simple context manager for when tqdm is not used
@contextlib.contextmanager
def nullcontext():
    yield

def get_nearest_cells(target_cells, reference_cells, rna_count_column, n_samples=100):
    """
    Find the nearest reference cells for each target cell based on RNA counts.

    Args:
        target_cells (AnnData): AnnData object containing target cells.
        reference_cells (AnnData): AnnData object containing reference cells.
        rna_count_column (str): Column name for RNA counts in the AnnData object's obs.
        n_samples (int): Number of nearest reference cells to find for each target cell.

    Returns:
        dict: Dictionary mapping each target cell ID to a list of nearest reference cell IDs.
    """
    if rna_count_column not in target_cells.obs.columns or rna_count_column not in reference_cells.obs.columns:
        raise ValueError(f"Column '{rna_count_column}' not found in AnnData object's obs.")

    results = {}
    target_counts = target_cells.obs[rna_count_column].values
    reference_counts = reference_cells.obs[rna_count_column].values
    reference_cell_ids = reference_cells.obs_names.values

    for i, target_cell_id in enumerate(target_cells.obs_names):
        differences = np.abs(reference_counts - target_counts[i])
        nearest_indices = np.argsort(differences)[:min(len(reference_counts), 1000)]

        if n_samples >= len(nearest_indices):
            sampled_indices = nearest_indices
        else:
            sampled_indices = np.random.choice(nearest_indices, n_samples, replace=False)

        results[target_cell_id] = reference_cell_ids[sampled_indices]

    logging.info(f"Nearest cells determined for {len(results)} target cells.")
    return results

def monte_carlo_comparison(adata, cell_type, cell_type_column, significant_genes_df, disease_control_genes=None, 
                          healthy_control_genes=None, output_dir=".", rna_count_column='nCount_RNA', 
                          iterations=100, target_group="disease", disease_marker='disease_numeric', 
                          disease_value=1, healthy_value=0, show_progress=False):
    """
    Perform Monte Carlo comparison for differential correlation analysis.

    Args:
        adata: AnnData object containing single-cell data
        cell_type: Cell type to analyze
        cell_type_column: Column in adata.obs containing cell type information
        significant_genes_df: DataFrame containing significant genes with gene_name and zstat columns
        disease_control_genes: Dictionary mapping significant genes to control genes for disease cells
        healthy_control_genes: Dictionary mapping significant genes to control genes for healthy cells
        output_dir: Directory to save output files
        rna_count_column: Column in adata.obs containing RNA count information
        iterations: Number of iterations for Monte Carlo simulation
        target_group: Target group for analysis ('disease' or 'healthy')
        disease_marker: Column in adata.obs containing disease status
        disease_value: Value in disease_marker column indicating disease (can be numeric or string)
        healthy_value: Value in disease_marker column indicating healthy (can be numeric or string)
        show_progress: Whether to show progress bars
        
    Returns:
        DataFrame containing results
    """
    try:
        logging.info(f"Starting Monte Carlo comparison for cell type: {cell_type}, target group: {target_group}")

        # Create a string-safe comparison function to handle both numeric and string types
        def value_equals(a, b):
            """Compare values that could be strings, numbers, or bool"""
            # Handle boolean special case
            if isinstance(a, bool) or isinstance(b, bool):
                return bool(a) == bool(b)
            # Try string comparison (works for both numbers and strings)
            return str(a) == str(b)
        
        # Filter cell type data
        if cell_type and cell_type_column in adata.obs.columns:
            adata_subset = adata[adata.obs[cell_type_column] == cell_type].copy()
        else:
            adata_subset = adata.copy()
            if cell_type:
                logging.warning(f"Column '{cell_type_column}' not found, using all data")
        
        # Normalize column names and ensure they exist
        significant_genes_df.columns = significant_genes_df.columns.str.lower().str.strip()
        
        if 'zstat' not in significant_genes_df.columns or 'gene_name' not in significant_genes_df.columns:
            logging.error("Columns 'zstat' or 'gene_name' not found in significant_genes_df.")
            return pd.DataFrame()
        
        # Based on target_group, filter for disease or healthy cells
        if target_group == "disease":
            # Use string-safe comparison
            mask = [value_equals(v, disease_value) for v in adata_subset.obs[disease_marker]]
            target_cells = adata_subset[mask].copy()
            mask = [value_equals(v, healthy_value) for v in adata_subset.obs[disease_marker]]
            reference_cells = adata_subset[mask].copy()
            control_genes = disease_control_genes
        else:  # target_group == "healthy"
            # Use string-safe comparison
            mask = [value_equals(v, healthy_value) for v in adata_subset.obs[disease_marker]]
            target_cells = adata_subset[mask].copy()
            mask = [value_equals(v, disease_value) for v in adata_subset.obs[disease_marker]]
            reference_cells = adata_subset[mask].copy()
            control_genes = healthy_control_genes

        # Log info about the selected cells
        logging.info(f"{len(target_cells)} target cells and {len(reference_cells)} reference cells identified for {target_group} group.")

        if len(target_cells) == 0 or len(reference_cells) == 0:
            logging.warning(f"No {'target' if len(target_cells) == 0 else 'reference'} cells found for {cell_type}, {target_group}.")
            return pd.DataFrame()  # Return empty dataframe
        
        # Create cell type directory
        cell_type_dir = os.path.join(output_dir, cell_type)
        os.makedirs(cell_type_dir, exist_ok=True)
        
        # Get nearest cells based on RNA counts
        matched_indices = get_nearest_cells(target_cells, reference_cells, rna_count_column, n_samples=100)
        
        # Prepare for iterations
        all_iterations_results = []
        
        # Check if expression matrix is sparse to optimize operations
        is_sparse_matrix = issparse(adata.X)
        
        # Prepare gene weights
        valid_genes = [gene for gene in significant_genes_df['gene_name'] if gene in adata.var_names]
        gene_weights = np.abs(significant_genes_df[significant_genes_df['gene_name'].isin(valid_genes)]['zstat'].values)
        gene_weights = gene_weights / np.sum(gene_weights) if np.sum(gene_weights) > 0 else gene_weights
        
        if not valid_genes:
            logging.error("No valid genes found in the dataset")
            return pd.DataFrame()
        
        # Preprocess control genes
        control_genes_filtered = {}
        if control_genes:
            for gene in valid_genes:
                if gene in control_genes:
                    # Keep only control genes that exist in the dataset
                    control_genes_filtered[gene] = [ctrl for ctrl in control_genes[gene] if ctrl in adata.var_names]
        
        # Progress bar context manager
        prog_context = tqdm(total=iterations, desc="Iterations") if show_progress else nullcontext()
        
        # Run iterations
        with prog_context as prog:
            for iteration in range(iterations):
                if show_progress:
                    prog.update(1)
                
                logging.info(f"Iteration {iteration + 1} of {iterations}")
                all_results = []
                
                # Process each target cell with its matched reference cells
                for idx, reference_cell_ids in matched_indices.items():
                    if is_sparse_matrix:
                        # For sparse matrices, manually convert to DataFrame
                        idx_pos = np.where(target_cells.obs_names == idx)[0][0]
                        target_row = target_cells[idx_pos].X
                        if issparse(target_row):
                            target_row = target_row.toarray()[0]
                        target_cell_expression = pd.Series(target_row, index=target_cells.var_names)
                        
                        # Get reference cell expressions
                        reference_rows = []
                        for ref_id in reference_cell_ids:
                            ref_pos = np.where(reference_cells.obs_names == ref_id)[0][0]
                            ref_row = reference_cells[ref_pos].X
                            if issparse(ref_row):
                                ref_row = ref_row.toarray()[0]
                            reference_rows.append(ref_row)
                        
                        reference_expression = pd.DataFrame(reference_rows, columns=reference_cells.var_names)
                        sampled_reference_cells_expression = reference_expression.mean()
                    else:
                        # For dense matrices, use to_df()
                        idx_pos = np.where(target_cells.obs_names == idx)[0][0]
                        target_cell_expression = target_cells[idx_pos].to_df().iloc[0]
                        
                        # Get reference cells and calculate mean expression
                        ref_positions = [np.where(reference_cells.obs_names == ref_id)[0][0] for ref_id in reference_cell_ids]
                        reference_cells_subset = reference_cells[ref_positions]
                        sampled_reference_cells_expression = reference_cells_subset.to_df().mean()
                    
                    # Initialize arrays for differences
                    sig_diffs = []
                    ctrl_diffs = []
                    
                    # Process each significant gene
                    for i, gene in enumerate(valid_genes):
                        weight = gene_weights[i]
                        
                        # Calculate the difference for the target gene
                        target_value = target_cell_expression.get(gene, 0)
                        reference_value = sampled_reference_cells_expression.get(gene, 0)
                        sig_diff = weight * abs(target_value - reference_value)
                        
                        # Ensure sig_diff is numeric before adding
                        if np.isscalar(sig_diff):
                            sig_diffs.append(sig_diff)
                        
                        # For control genes, use the list of control genes for this target gene
                        if gene in control_genes_filtered and control_genes_filtered[gene]:
                            ctrl_gene = np.random.choice(control_genes_filtered[gene])
                            
                            # Get control gene expression values
                            ctrl_target_value = target_cell_expression.get(ctrl_gene, 0)
                            ctrl_reference_value = sampled_reference_cells_expression.get(ctrl_gene, 0)
                            
                            # Calculate control difference
                            ctrl_diff = weight * abs(ctrl_target_value - ctrl_reference_value)
                            
                            # Ensure ctrl_diff is numeric before adding
                            if np.isscalar(ctrl_diff):
                                ctrl_diffs.append(ctrl_diff)
                    
                    # Skip cells with insufficient data
                    if len(sig_diffs) == 0 or len(ctrl_diffs) == 0:
                        continue
                    
                    # Convert to numpy arrays and ensure float type
                    sig_diffs = np.array([float(x) for x in sig_diffs], dtype=float)
                    ctrl_diffs = np.array([float(x) for x in ctrl_diffs], dtype=float)
                    
                    # Remove zeros and NaNs
                    sig_diffs = sig_diffs[~np.isnan(sig_diffs) & (sig_diffs != 0)]
                    ctrl_diffs = ctrl_diffs[~np.isnan(ctrl_diffs) & (ctrl_diffs != 0)]
                    
                    if len(sig_diffs) == 0 or len(ctrl_diffs) == 0:
                        continue
                    
                    # Perform t-test
                    t_stat, p_val = ttest_ind(sig_diffs, ctrl_diffs, equal_var=False)
                    
                    # Calculate totals
                    total_sig_diff = np.sum(sig_diffs)
                    total_ctrl_diff = np.sum(ctrl_diffs)
                    
                    # Calculate one-tailed p-value based on direction
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

def compare_groups(disease_df, healthy_df):
    """
    Compare results between disease and healthy groups.

    Args:
        disease_df: DataFrame containing disease group results
        healthy_df: DataFrame containing healthy group results
        
    Returns:
        dict: Dictionary containing comparison statistics
    """
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