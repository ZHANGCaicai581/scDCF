#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import logging
from scipy.stats import combine_pvalues, norm, fisher_exact
from statsmodels.stats.multitest import multipletests

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

CELL_Q_ALPHA = 0.05
CELLTYPE_Q_ALPHA = 0.05


def _sanitize_p_values(p_values):
    """Clip p-values into the open interval (0, 1] and drop invalid entries."""
    clean = pd.to_numeric(pd.Series(p_values), errors='coerce').dropna()
    if clean.empty:
        return np.array([], dtype=float)
    clean = clean.clip(lower=np.finfo(float).tiny, upper=1.0)
    return clean.to_numpy(dtype=float)


def _add_fdr_column(df, p_column, q_column):
    """Append BH-adjusted q-values for the supplied p-value column."""
    out = df.copy()
    out[q_column] = np.nan

    if out.empty or p_column not in out.columns:
        return out

    valid_mask = out[p_column].notna()
    if not valid_mask.any():
        return out

    out.loc[valid_mask, q_column] = multipletests(
        out.loc[valid_mask, p_column].astype(float),
        method='fdr_bh'
    )[1]
    return out

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

def combine_p_values_across_iterations(combined_results, output_dir, cell_type, target_group, write_output=False):
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
    
    # Get unique cells/genes
    unit_ids = combined_results[gene_column].dropna().unique()
    logging.info(f"Found {len(unit_ids)} unique entries to aggregate")
    
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
    for unit_id in unit_ids:
        unit_results = combined_results[combined_results[gene_column] == unit_id]
        iterations = unit_results[iteration_column].dropna().unique()
        
        # Extract p-values across iterations
        p_values = []
        
        # Determine p-value column
        p_value_column = None
        possible_p_columns = ['p_value', 'pvalue', 'p-value', 'p']
        for col in possible_p_columns:
            if col in unit_results.columns:
                p_value_column = col
                break

        if p_value_column is None:
            # Try to infer p-value column
            p_like_columns = [col for col in unit_results.columns if 'p' in col.lower() and 'value' in col.lower()]
            if p_like_columns:
                p_value_column = p_like_columns[0]
            else:
                logging.warning(f"Could not find p-value column for '{unit_id}', skipping")
                continue

        for iteration in iterations:
            iter_results = unit_results[unit_results[iteration_column] == iteration]
            if not iter_results.empty and p_value_column in iter_results.columns:
                p_value = iter_results[p_value_column].iloc[0]
                if not np.isnan(p_value):
                    p_values.append(p_value)

        if len(p_values) > 0:
            # Combine p-values using Fisher's method
            try:
                clean_p_values = _sanitize_p_values(p_values)
                if clean_p_values.size == 0:
                    continue

                combined_pvalue_fisher = combine_pvalues(clean_p_values, method='fisher')[1]
                # Combine p-values using Stouffer's method
                combined_pvalue_stouffer = combine_pvalues(clean_p_values, method='stouffer')[1]

                combined_p_values.append({
                    'cell_id': str(unit_id),
                    'num_iterations': int(clean_p_values.size),
                    'combined_p_value_fisher': combined_pvalue_fisher,
                    'combined_p_value_stouffer': combined_pvalue_stouffer,
                    'min_p_value': float(np.min(clean_p_values)),
                    'max_p_value': float(np.max(clean_p_values)),
                    'mean_p_value': float(np.mean(clean_p_values)),
                    'target_group': target_group
                })
            except Exception as e:
                logging.warning(f"Error combining p-values for '{unit_id}': {e}")

    # Create DataFrame from combined p-values
    combined_p_values_df = pd.DataFrame(combined_p_values)
    if combined_p_values_df.empty:
        combined_p_values_df = pd.DataFrame(columns=[
            'cell_id',
            'num_iterations',
            'combined_p_value_fisher',
            'combined_p_value_stouffer',
            'min_p_value',
            'max_p_value',
            'mean_p_value',
            'target_group'
        ])
    else:
        combined_p_values_df = _add_fdr_column(
            combined_p_values_df,
            p_column='combined_p_value_fisher',
            q_column='q_cell'
        )
        combined_p_values_df['scdcf_significant'] = combined_p_values_df['q_cell'] <= CELL_Q_ALPHA

    if write_output:
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
        group_df['analyzed_cell_type'] = cell_type
        # Derived statistics
        if 'combined_p_value_fisher' in group_df.columns:
            group_df['combined_minus_log10_p'] = -np.log10(group_df['combined_p_value_fisher'].clip(lower=1e-300))
            group_df['combined_fisher_z'] = norm.isf(group_df['combined_p_value_fisher'].clip(lower=1e-300))
        if 'combined_p_value_stouffer' in group_df.columns:
            group_df['combined_stouffer_z'] = norm.isf(group_df['combined_p_value_stouffer'].clip(lower=1e-300))
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
        'q_cell', 'scdcf_significant',
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


def load_final_summaries(output_dir, cell_types=None):
    """Load per-cell-type final summaries from disk."""
    summaries = {}

    if cell_types is None:
        if not os.path.exists(output_dir):
            return summaries
        cell_types = [
            entry for entry in os.listdir(output_dir)
            if os.path.isdir(os.path.join(output_dir, entry))
        ]

    for cell_type in cell_types:
        summary_path = os.path.join(output_dir, cell_type, f"{cell_type}_final_summary.csv")
        if not os.path.exists(summary_path):
            logging.warning(f"Final summary not found for cell type '{cell_type}': {summary_path}")
            continue
        try:
            summaries[cell_type] = pd.read_csv(summary_path)
        except Exception as exc:
            logging.warning(f"Unable to load final summary for '{cell_type}': {exc}")

    return summaries


def apply_dataset_level_cell_fdr(combined_results_by_cell_type, output_dir='.', alpha=CELL_Q_ALPHA, write_output=False):
    """
    Recompute q_cell across all tested cells within each target-group analysis.

    The input is a mapping of:
        {cell_type: {'disease': disease_df, 'healthy': healthy_df}}
    """
    adjusted = {}

    for target_group in ('disease', 'healthy'):
        group_frames = []

        for cell_type, group_map in combined_results_by_cell_type.items():
            group_df = group_map.get(target_group)
            if group_df is None or group_df.empty:
                continue

            tmp = group_df.copy()
            tmp['cell_type'] = cell_type
            group_frames.append(tmp)

        if not group_frames:
            continue

        all_group_df = pd.concat(group_frames, ignore_index=True)
        all_group_df = _add_fdr_column(
            all_group_df,
            p_column='combined_p_value_fisher',
            q_column='q_cell'
        )
        all_group_df['scdcf_significant'] = all_group_df['q_cell'] <= alpha

        for cell_type, cell_df in all_group_df.groupby('cell_type', sort=False):
            split_df = cell_df.drop(columns=['cell_type']).reset_index(drop=True)
            adjusted.setdefault(cell_type, {})[target_group] = split_df

            if write_output:
                combined_path = os.path.join(output_dir, cell_type, f"{cell_type}_{target_group}_combined.csv")
                split_df.to_csv(combined_path, index=False)
                logging.info(
                    f"Dataset-level q_cell updated for {cell_type} ({target_group}) and saved to {combined_path}"
                )

    return adjusted


def compute_celltype_enrichment(final_summaries=None, output_dir='.', cell_types=None,
                                alpha=CELL_Q_ALPHA, celltype_alpha=CELLTYPE_Q_ALPHA):
    """
    Compare significant-cell proportions between disease and healthy groups per cell type.

        Fisher's exact test is applied independently within each analyzed cell type with the
        one-sided alternative that the disease group has a higher significant-cell proportion,
        then the resulting p-values are BH-adjusted across the tested cell types.
    """
    if final_summaries is None:
        final_summaries = load_final_summaries(output_dir, cell_types=cell_types)

    records = []

    for cell_type, df in final_summaries.items():
        if df is None or df.empty:
            continue

        if 'group' not in df.columns:
            logging.warning(f"Skipping cell-type enrichment for '{cell_type}': missing 'group' column.")
            continue

        if 'q_cell' not in df.columns:
            logging.warning(f"Skipping cell-type enrichment for '{cell_type}': missing 'q_cell' column.")
            continue

        group_series = df['group'].astype(str)
        disease_mask = group_series == 'disease'
        healthy_mask = group_series == 'healthy'

        n_disease = int(disease_mask.sum())
        n_healthy = int(healthy_mask.sum())

        if n_disease == 0 or n_healthy == 0:
            logging.warning(
                f"Skipping cell-type enrichment for '{cell_type}': "
                f"disease cells={n_disease}, healthy cells={n_healthy}."
            )
            continue

        disease_sig = int((df.loc[disease_mask, 'q_cell'] <= alpha).sum())
        healthy_sig = int((df.loc[healthy_mask, 'q_cell'] <= alpha).sum())
        disease_non_sig = n_disease - disease_sig
        healthy_non_sig = n_healthy - healthy_sig

        table = [[disease_sig, disease_non_sig], [healthy_sig, healthy_non_sig]]
        odds_ratio, p_value = fisher_exact(table, alternative='greater')

        disease_prop = disease_sig / n_disease
        healthy_prop = healthy_sig / n_healthy

        records.append({
            'cell_type': cell_type,
            'disease_significant_cells': disease_sig,
            'disease_non_significant_cells': disease_non_sig,
            'healthy_significant_cells': healthy_sig,
            'healthy_non_significant_cells': healthy_non_sig,
            'disease_total_cells': n_disease,
            'healthy_total_cells': n_healthy,
            'disease_prop_significant': disease_prop,
            'healthy_prop_significant': healthy_prop,
            'odds_ratio': odds_ratio,
            'p_value': p_value,
            'fisher_alternative': 'greater',
            'direction': 'disease_higher' if disease_prop > healthy_prop else 'healthy_higher_or_equal'
        })

    enrichment_df = pd.DataFrame(records)
    if enrichment_df.empty:
        logging.warning("No cell-type enrichment results were generated.")
        return enrichment_df

    enrichment_df = _add_fdr_column(enrichment_df, p_column='p_value', q_column='q_type')
    enrichment_df['disease_enriched'] = (
        (enrichment_df['disease_prop_significant'] > enrichment_df['healthy_prop_significant']) &
        (enrichment_df['q_type'] <= celltype_alpha)
    )
    enrichment_df = enrichment_df.sort_values(
        ['disease_enriched', 'q_type', 'p_value', 'cell_type'],
        ascending=[False, True, True, True]
    ).reset_index(drop=True)

    output_path = os.path.join(output_dir, "celltype_enrichment_summary.csv")
    enrichment_df.to_csv(output_path, index=False)
    logging.info(f"Cell-type enrichment summary saved to {output_path} ({len(enrichment_df)} cell types).")

    return enrichment_df

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
