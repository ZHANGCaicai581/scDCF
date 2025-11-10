#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import logging
from scipy.stats import combine_pvalues, norm

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def load_monte_carlo_results(results_file):
    """
    Load results from Monte Carlo iterations.
    
    Args:
        results_file: Path to the CSV file containing Monte Carlo results
        
    Returns:
        DataFrame containing Monte Carlo results
    """
    if not os.path.exists(results_file):
        logging.error(f"Results file not found: {results_file}")
        return pd.DataFrame()
    
    try:
        results = pd.read_csv(results_file)
        logging.info(f"Loaded {len(results)} results from {results_file}")
        return results
    except Exception as e:
        logging.error(f"Error loading results from {results_file}: {e}")
        return pd.DataFrame()

def combine_p_values_across_iterations(combined_results, output_dir, cell_type, target_group):
    """
    Combine p-values across Monte Carlo iterations.
    
    Args:
        combined_results: DataFrame containing Monte Carlo results
        output_dir: Directory to save the combined p-values
        cell_type: Cell type identifier
        target_group: Target group identifier (disease or healthy)
        
    Returns:
        DataFrame with combined p-values
    """
    logging.info(f"Combining p-values for {cell_type}, {target_group} group")
    
    # Debug: Print column names to help diagnose issues
    logging.info(f"Available columns in results: {combined_results.columns.tolist()}")
    
    # Determine gene column name - try different possible names
    gene_column = None
    possible_gene_columns = ['cell_id', 'gene', 'gene_name', 'significant_gene', 'sig_gene']
    
    for col in possible_gene_columns:
        if col in combined_results.columns:
            gene_column = col
            logging.info(f"Using column '{gene_column}' for gene identifiers")
            break
    
    # If no gene column found, try to infer it
    if gene_column is None:
        # Look for columns that might contain gene names (containing 'gene' in the name)
        gene_like_columns = [col for col in combined_results.columns if 'gene' in col.lower()]
        if gene_like_columns:
            gene_column = gene_like_columns[0]
            logging.info(f"Inferred gene column: '{gene_column}'")
        else:
            # If we still can't find it, use the first column as a fallback
            gene_column = combined_results.columns[0]
            logging.warning(f"Could not identify gene column, using first column: '{gene_column}'")
    
    # Get unique genes
    genes = combined_results[gene_column].unique()
    logging.info(f"Found {len(genes)} unique genes")
    
    # Create output directory for cell type
    output_dir_cell_type = os.path.join(output_dir, cell_type)
    os.makedirs(output_dir_cell_type, exist_ok=True)
    
    # Initialize list to store combined p-values
    combined_p_values = []
    
    # Determine iteration column 
    iteration_column = 'iteration'
    if iteration_column not in combined_results.columns:
        # Try to find an alternative
        if 'iter' in combined_results.columns:
            iteration_column = 'iter'
        elif 'monte_carlo_iteration' in combined_results.columns:
            iteration_column = 'monte_carlo_iteration'
        else:
            # If no iteration column, assume all are from the same iteration
            combined_results['iteration'] = 1
            iteration_column = 'iteration'
            logging.warning("No iteration column found, assuming all results are from a single iteration")
    
    # Combine p-values for each gene
    for gene in genes:
        gene_results = combined_results[combined_results[gene_column] == gene]
        iterations = gene_results[iteration_column].unique()
        
        # Extract p-values across iterations
        p_values = []
        
        # Determine p-value column
        p_value_column = None
        possible_p_columns = ['p_value', 'pvalue', 'p-value', 'p']
        for col in possible_p_columns:
            if col in gene_results.columns:
                p_value_column = col
                break
        
        if p_value_column is None:
            # Try to infer p-value column
            p_like_columns = [col for col in gene_results.columns if 'p' in col.lower() and 'value' in col.lower()]
            if p_like_columns:
                p_value_column = p_like_columns[0]
            else:
                logging.warning(f"Could not find p-value column for gene {gene}, skipping")
                continue
        
        for iteration in iterations:
            iter_results = gene_results[gene_results[iteration_column] == iteration]
            if not iter_results.empty and p_value_column in iter_results.columns:
                p_value = iter_results[p_value_column].iloc[0]
                if not np.isnan(p_value):
                    p_values.append(p_value)
        
        if len(p_values) > 0:
            # Combine p-values using Fisher's method
            try:
                combined_pvalue_fisher = combine_pvalues(p_values, method='fisher')[1]
                # Combine p-values using Stouffer's method
                combined_pvalue_stouffer = combine_pvalues(p_values, method='stouffer')[1]
                
                combined_p_values.append({
                    'gene': gene,
                    'num_iterations': len(p_values),
                    'combined_p_value_fisher': combined_pvalue_fisher,
                    'combined_p_value_stouffer': combined_pvalue_stouffer,
                    'min_p_value': min(p_values),
                    'max_p_value': max(p_values),
                    'mean_p_value': np.mean(p_values)
                })
            except Exception as e:
                logging.warning(f"Error combining p-values for gene {gene}: {e}")
    
    # Create DataFrame from combined p-values
    combined_p_values_df = pd.DataFrame(combined_p_values)
    if not combined_p_values_df.empty and 'gene' in combined_p_values_df.columns and 'cell_id' not in combined_p_values_df.columns:
        combined_p_values_df = combined_p_values_df.rename(columns={'gene': 'cell_id'})
    
    # Save combined p-values
    combined_p_values_file = os.path.join(output_dir_cell_type, f"{cell_type}_{target_group}_combined.csv")
    combined_p_values_df.to_csv(combined_p_values_file, index=False)
    logging.info(f"Combined p-values saved to {combined_p_values_file}")

    return combined_p_values_df

def visualize_combined_p_values(disease_combined, healthy_combined, cell_type, output_dir='.'):
    """Placeholder visualization hook (disabled)."""
    logging.info(f"Visualization skipped for {cell_type} (plotting disabled)")
    return


def export_final_celltype_summary(cell_type, disease_combined, healthy_combined,
                                  output_dir,
                                  include_metadata=True, adata=None,
                                  metadata_columns=None):
    """
    Create a final per-cell-type summary CSV with key statistics.

    The summary includes combined p-values, derived statistics, and (by default)
    AnnData metadata so researchers have a single file to inspect.
    """
    cell_dir = os.path.join(output_dir, cell_type)
    os.makedirs(cell_dir, exist_ok=True)

    frames = []
    for group_name, df in [('disease', disease_combined), ('healthy', healthy_combined)]:
        if df is None or df.empty:
            continue
        group_df = df.copy()
        if 'gene' in group_df.columns and 'cell_id' not in group_df.columns:
            group_df = group_df.rename(columns={'gene': 'cell_id'})
        group_df['cell_id'] = group_df['cell_id'].astype(str)
        group_df['target_group'] = group_df.get('target_group', group_name)
        group_df['group'] = group_name
        # Derived statistics
        if 'combined_p_value_fisher' in group_df.columns:
            group_df['combined_minus_log10_p'] = -np.log10(group_df['combined_p_value_fisher'].clip(lower=1e-300))
            group_df['combined_fisher_z'] = norm.isf(group_df['combined_p_value_fisher'].clip(lower=1e-300) / 2)
        if 'combined_p_value_stouffer' in group_df.columns:
            group_df['combined_stouffer_z'] = norm.isf(group_df['combined_p_value_stouffer'].clip(lower=1e-300) / 2)
        frames.append(group_df)

    if not frames:
        logging.warning(f"No combined results available for {cell_type}; final summary not generated.")
        return pd.DataFrame()

    final_df = pd.concat(frames, ignore_index=True)
    final_df['original_cell_id'] = final_df['cell_id']

    if include_metadata and adata is not None:
        try:
            from .cell_metadata import add_cell_metadata
            final_df = add_cell_metadata(
                final_df,
                adata=adata,
                cell_id_column='cell_id',
                metadata_columns=metadata_columns
            )
        except Exception as exc:
            logging.warning(f"Unable to append metadata for {cell_type}: {exc}")

    columns_order = [
        'cell_id', 'original_cell_id', 'group',
        'combined_p_value_fisher', 'combined_p_value_stouffer',
        'combined_minus_log10_p', 'combined_fisher_z', 'combined_stouffer_z',
        'num_iterations', 'min_p_value', 'max_p_value', 'mean_p_value'
    ]
    existing_order = [col for col in columns_order if col in final_df.columns]
    remaining_cols = [col for col in final_df.columns if col not in existing_order]
    final_df = final_df[existing_order + remaining_cols]

    final_path = os.path.join(cell_dir, f"{cell_type}_final_summary.csv")
    final_df.to_csv(final_path, index=False)
    logging.info(f"Final summary saved to {final_path} ({len(final_df)} rows).")

    return final_df

def examine_results_format(results_file):
    """
    Examine and log the format of a Monte Carlo results file
    
    Args:
        results_file: Path to the CSV file to examine
    """
    if not os.path.exists(results_file):
        logging.error(f"Results file not found: {results_file}")
        return
    
    try:
        results = pd.read_csv(results_file)
        logging.info(f"File: {results_file}")
        logging.info(f"  Columns: {results.columns.tolist()}")
        logging.info(f"  Shape: {results.shape}")
        
        # Print first few rows
        logging.info(f"  First rows:\n{results.head(2)}")
        
        return results
    except Exception as e:
        logging.error(f"Error examining results file {results_file}: {e}")
        return None

def organize_results(source_dir, dest_dir="output", cell_types=None):
    """
    Organize the analysis results into a cell-type centric structure.
    
    Args:
        source_dir: Directory containing the raw analysis results
        dest_dir: Directory to create the organized structure in
        cell_types: List of cell types to organize. If None, will autodetect from source_dir.
        
    Returns:
        str: Path to the organized output directory
    """
    import os
    import pandas as pd
    import shutil
    import glob
    import logging
    
    logging.info(f"Organizing results from {source_dir} to {dest_dir}")
    
    # Auto-detect cell types if not provided
    if cell_types is None:
        cell_types = []
        # Detect cell types from directories in source_dir
        if os.path.exists(source_dir):
            for item in os.listdir(source_dir):
                item_path = os.path.join(source_dir, item)
                if os.path.isdir(item_path):
                    # Check if this looks like a cell type directory
                    monte_carlo_files = glob.glob(os.path.join(item_path, f"{item}_*_monte_carlo_results.csv"))
                    if monte_carlo_files:
                        cell_types.append(item)
        
        if not cell_types:
            raise ValueError(f"No cell types could be detected in {source_dir}. Please specify cell_types manually.")
        else:
            logging.info(f"Detected cell types: {cell_types}")
    
    # Create base output directory
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    
    # Process each cell type
    for cell_type in cell_types:
        # Create directory structure for this cell type
        cell_dir = os.path.join(dest_dir, cell_type)
        supporting_dir = os.path.join(cell_dir, "supporting_data")
        iterations_dir = os.path.join(supporting_dir, "monte_carlo_iterations")
        control_genes_dir = os.path.join(supporting_dir, "control_genes")
        
        for dir_path in [cell_dir, supporting_dir, iterations_dir, control_genes_dir]:
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
        
        # 1. Copy Monte Carlo iteration files to supporting_data/monte_carlo_iterations
        iteration_pattern = os.path.join(source_dir, cell_type, f"{cell_type}_*_monte_carlo_results_iteration*.csv")
        iteration_files = glob.glob(iteration_pattern)
        for file in iteration_files:
            filename = os.path.basename(file)
            dest_file = os.path.join(iterations_dir, filename)
            shutil.copy2(file, dest_file)
        
        # 2. Copy combined files and KS test results to supporting_data
        key_files = [
            f"{cell_type}_disease_combined.csv",
            f"{cell_type}_healthy_combined.csv",
            f"{cell_type}_disease_monte_carlo_results.csv",
            f"{cell_type}_healthy_monte_carlo_results.csv",
            f"{cell_type}_trait_scores.csv",
            "trait_association_scores.csv",
            "disease_combined_p_values.csv",
            "healthy_combined_p_values.csv"
        ]
        
        for filename in key_files:
            src_file = os.path.join(source_dir, cell_type, filename)
            if os.path.exists(src_file):
                dest_file = os.path.join(supporting_dir, filename)
                shutil.copy2(src_file, dest_file)
        
        # 3. Create comprehensive cell_metrics.csv file
        cell_data = None
        monte_carlo_file = os.path.join(source_dir, cell_type, f"{cell_type}_disease_monte_carlo_results.csv")
        
        if os.path.exists(monte_carlo_file):
            try:
                cell_data = pd.read_csv(monte_carlo_file)
                # Extract essential columns
                essential_cols = ['cell_id', 't_stat', 'p_value']
                if 'p_value_adj' in cell_data.columns:
                    essential_cols.append('p_value_adj')
                
                cell_data = cell_data[essential_cols].copy()
                cell_data['cell_type'] = cell_type
            except Exception as e:
                logging.error(f"Error processing Monte Carlo results for {cell_type}: {e}")
        
        # Add enhanced scores (z-scores and -log10 p-values)
        if cell_data is not None:
            enhanced_file = os.path.join(source_dir, cell_type, "trait_association_scores.csv")
            if os.path.exists(enhanced_file):
                try:
                    enhanced = pd.read_csv(enhanced_file)
                    # Find enhanced score columns
                    score_cols = [col for col in enhanced.columns 
                                if 'z_score' in col or 'neg_log10' in col]
                    
                    if 'cell_id' in enhanced.columns and score_cols:
                        # Merge with cell data
                        cell_data = pd.merge(
                            cell_data, 
                            enhanced[['cell_id'] + score_cols],
                            on='cell_id', how='left'
                        )
                except Exception as e:
                    logging.error(f"Error adding enhanced scores for {cell_type}: {e}")
        
            # Add trait correlations
            trait_file = os.path.join(source_dir, cell_type, f"{cell_type}_trait_scores.csv")
            if os.path.exists(trait_file):
                try:
                    trait_data = pd.read_csv(trait_file)
                    for _, row in trait_data.iterrows():
                        trait = row.get('trait')
                        if trait:
                            cell_data[f'trait_{trait}_score'] = row.get('correlation')
                except Exception as e:
                    logging.error(f"Error adding trait correlations for {cell_type}: {e}")
        
            # Save the comprehensive cell_metrics.csv
            if cell_data is not None:
                try:
                    cell_data.to_csv(os.path.join(cell_dir, "cell_metrics.csv"), index=False)
                except Exception as e:
                    logging.error(f"Error saving cell_metrics.csv for {cell_type}: {e}")
    
    # Copy metadata files
    for filename in ["simulated_trait_data.csv"]:
        src_file = os.path.join(source_dir, filename)
        if os.path.exists(src_file):
            dest_file = os.path.join(dest_dir, filename)
            shutil.copy2(src_file, dest_file)
    
    # Copy control gene files to cell-specific supporting_data/control_genes folders
    for cell_type in cell_types:
        control_file = os.path.join(source_dir, f"{cell_type}_control_genes.json")
        if os.path.exists(control_file):
            # Cell type specific location in the control_genes subfolder
            dest_file = os.path.join(dest_dir, cell_type, "supporting_data", "control_genes", f"{cell_type}_control_genes.json")
            shutil.copy2(control_file, dest_file)
            
            # Also maintain a copy at the root level for backward compatibility
            root_dest_file = os.path.join(dest_dir, f"{cell_type}_control_genes.json")
            shutil.copy2(control_file, root_dest_file)
    
    logging.info(f"Results organized in {dest_dir}/")
    return dest_dir