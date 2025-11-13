"""
Control gene selection for scDCF package.
"""

import pandas as pd
import numpy as np
import logging
import os
import json
import scipy.sparse as sp

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------
# Caches to reuse AnnData slices and per-group statistics across calls.
# Keys include the id of the AnnData object to avoid cross-dataset reuse.
# -------------------------------------------------------------------------
_ADATA_SUBSET_CACHE = {}
_GROUP_STATS_CACHE = {}


def _value_equals(a, b):
    """Robust equality check for values that may be strings, numbers, or bools."""
    if isinstance(a, bool) or isinstance(b, bool):
        return bool(a) == bool(b)
    return str(a) == str(b)


def _get_subset_indices(adata, cell_type, cell_type_column):
    """Return cached indices for the requested cell type subset."""
    cache_key = (id(adata), str(cell_type) if cell_type is not None else "__all__", cell_type_column)
    cached = _ADATA_SUBSET_CACHE.get(cache_key)
    if cached is not None:
        return cached

    if cell_type is None or (isinstance(cell_type, str) and cell_type == ""):
        indices = np.arange(adata.n_obs, dtype=int)
    else:
        if cell_type_column not in adata.obs.columns:
            raise KeyError(f"Column '{cell_type_column}' not found in AnnData.obs.")
        obs_vals = adata.obs[cell_type_column].astype(str).to_numpy()
        mask = obs_vals == str(cell_type)
        indices = np.where(mask)[0]

    _ADATA_SUBSET_CACHE[cache_key] = indices
    return indices


def _get_group_stats(adata, subset_key, subset_indices, disease_marker, target_value):
    """
    Return cached per-gene statistics (mean/var) for the specified group.

    Returns:
        stats_df (pd.DataFrame): index=gene, columns=['mean', 'var']
        n_cells (int): number of cells contributing to the statistics
    """
    stats_key = (id(adata), subset_key, disease_marker, str(target_value))
    cached = _GROUP_STATS_CACHE.get(stats_key)
    if cached is not None:
        return cached

    if disease_marker not in adata.obs.columns:
        raise KeyError(f"Column '{disease_marker}' not found in AnnData.obs.")

    obs_values = adata.obs.iloc[subset_indices][disease_marker].to_numpy()
    mask = np.fromiter((_value_equals(v, target_value) for v in obs_values), dtype=bool, count=len(obs_values))
    group_indices = subset_indices[mask]
    n_cells = int(group_indices.size)

    if n_cells == 0:
        stats_df = pd.DataFrame(index=adata.var_names, data={'mean': np.nan, 'var': np.nan})
        _GROUP_STATS_CACHE[stats_key] = (stats_df, n_cells)
        return stats_df, n_cells

    group_matrix = adata.X[group_indices, :]
    if sp.issparse(group_matrix):
        means = np.asarray(group_matrix.mean(axis=0)).ravel()
        squared_means = np.asarray(group_matrix.power(2).mean(axis=0)).ravel()
        vars_ = squared_means - means ** 2
    else:
        means = np.asarray(group_matrix.mean(axis=0)).ravel()
        vars_ = np.asarray(group_matrix.var(axis=0)).ravel()

    vars_ = np.maximum(vars_, 0.0)

    stats_df = pd.DataFrame({'mean': means, 'var': vars_}, index=adata.var_names)
    stats_df.index.name = 'gene'

    _GROUP_STATS_CACHE[stats_key] = (stats_df, n_cells)
    return stats_df, n_cells

def generate_control_genes(adata, significant_genes_df, cell_type, cell_type_column,
                           n_control_genes=10, disease_marker='disease_numeric',
                           disease_value=1, healthy_value=0, output_dir=None):
    """
    Generate control genes for differential correlation analysis.
    
    Args:
        adata: AnnData object containing single-cell data
        significant_genes_df: DataFrame containing significant genes with gene_name and zstat columns
        cell_type: Cell type to analyze
        cell_type_column: Column in adata.obs containing cell type information
        n_control_genes: Number of control genes to generate for each significant gene
        disease_marker: Column in adata.obs containing disease status
        disease_value: Value in disease_marker column indicating disease
        healthy_value: Value in disease_marker column indicating healthy
        output_dir: Directory to save control genes JSON file
        
    Returns:
        Tuple of (disease_control_genes, healthy_control_genes)
    """
    logger.info(f"Generating control genes for cell type: {cell_type}")

    subset_label = str(cell_type) if isinstance(cell_type, str) and cell_type else "__all__"
    subset_key = (subset_label, cell_type_column)
    subset_indices = _get_subset_indices(adata, cell_type, cell_type_column)

    if subset_indices.size == 0:
        logger.warning(f"No cells found for cell type '{cell_type}'. Cannot generate control genes.")
        return {}, {}

    # Get list of significant genes that exist in the dataset
    sig_genes = [gene for gene in significant_genes_df['gene_name'] if gene in adata.var_names]

    if len(sig_genes) == 0:
        logger.warning("No significant genes found in the dataset. Cannot generate control genes.")
        return {}, {}

    logger.info(f"Found {len(sig_genes)} significant genes in the dataset")

    disease_stats, disease_cells = _get_group_stats(
        adata, subset_key, subset_indices, disease_marker, disease_value
    )
    healthy_stats, healthy_cells = _get_group_stats(
        adata, subset_key, subset_indices, disease_marker, healthy_value
    )

    logger.info(f"Split data into {disease_cells} disease cells and {healthy_cells} healthy cells")

    # Clean statistics (remove NaNs/Infs) and precompute pools excluding significant genes
    disease_stats = disease_stats.replace([np.inf, -np.inf], np.nan).dropna(subset=['mean', 'var'])
    healthy_stats = healthy_stats.replace([np.inf, -np.inf], np.nan).dropna(subset=['mean', 'var'])
    disease_pool = disease_stats.drop(sig_genes, errors='ignore')
    healthy_pool = healthy_stats.drop(sig_genes, errors='ignore')

    # Generate control genes
    disease_control_genes = {}
    healthy_control_genes = {}

    for gene in sig_genes:
        # Find control genes for disease group
        if gene in disease_stats.index:
            gene_stats = disease_stats.loc[gene]
            if gene_stats[['mean', 'var']].isnull().any():
                disease_control_genes[gene] = []
            elif not disease_pool.empty:
                filtered_stats = disease_pool
                mean_diff = (filtered_stats['mean'].values - gene_stats['mean']) / (gene_stats['mean'] + 1e-6)
                var_diff = (filtered_stats['var'].values - gene_stats['var']) / (gene_stats['var'] + 1e-6)
                distances = np.sqrt(mean_diff**2 + var_diff**2)
                ranked = np.argsort(distances)
                selected = filtered_stats.index.values[ranked[:n_control_genes]].tolist()
                disease_control_genes[gene] = selected
            else:
                disease_control_genes[gene] = []
        else:
            disease_control_genes[gene] = []

        # Find control genes for healthy group
        if gene in healthy_stats.index:
            gene_stats = healthy_stats.loc[gene]
            if gene_stats[['mean', 'var']].isnull().any():
                healthy_control_genes[gene] = []
            elif not healthy_pool.empty:
                filtered_stats = healthy_pool
                mean_diff = (filtered_stats['mean'].values - gene_stats['mean']) / (gene_stats['mean'] + 1e-6)
                var_diff = (filtered_stats['var'].values - gene_stats['var']) / (gene_stats['var'] + 1e-6)
                distances = np.sqrt(mean_diff**2 + var_diff**2)
                ranked = np.argsort(distances)
                selected = filtered_stats.index.values[ranked[:n_control_genes]].tolist()
                healthy_control_genes[gene] = selected
            else:
                healthy_control_genes[gene] = []
        else:
            healthy_control_genes[gene] = []

    # Save control genes to a file if output_dir is provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        control_genes_dict = {
            'disease_control_genes': disease_control_genes,
            'healthy_control_genes': healthy_control_genes
        }
        
        # Save as JSON
        cell_type_safe = cell_type.replace('/', '_').replace(' ', '_') if isinstance(cell_type, str) else 'all'
        control_genes_file = os.path.join(output_dir, f"{cell_type_safe}_control_genes.json")
        with open(control_genes_file, 'w') as f:
            json.dump(control_genes_dict, f, indent=4)
        
        logger.info(f"Control genes saved to {control_genes_file}")
    
    return disease_control_genes, healthy_control_genes