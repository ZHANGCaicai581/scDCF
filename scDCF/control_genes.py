"""
Control gene selection for scDCF package.
"""

import pandas as pd
import numpy as np
import scanpy as sc
import logging
import os
import json
from scipy.stats import pearsonr
import scipy.sparse as sp

logger = logging.getLogger(__name__)

def generate_control_genes(adata, significant_genes_df, cell_type, cell_type_column, 
                          n_control_genes=5, disease_marker='disease_numeric', 
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
    
    # Create string-safe comparison function to handle different types
    def value_equals(a, b):
        """Compare values that could be strings, numbers, or bool"""
        # Handle boolean special case
        if isinstance(a, bool) or isinstance(b, bool):
            return bool(a) == bool(b)
        # Try string comparison (works for both numbers and strings)
        return str(a) == str(b)
    
    # Filter data by cell type
    if isinstance(cell_type, str) and cell_type:
        adata_subset = adata[adata.obs[cell_type_column] == cell_type].copy()
    else:
        adata_subset = adata.copy()
    
    # Get list of significant genes that exist in the dataset
    sig_genes = [gene for gene in significant_genes_df['gene_name'] if gene in adata_subset.var_names]
    
    if len(sig_genes) == 0:
        logger.warning("No significant genes found in the dataset. Cannot generate control genes.")
        return {}, {}
    
    logger.info(f"Found {len(sig_genes)} significant genes in the dataset")
    
    # Split data into disease and healthy
    disease_mask = [value_equals(v, disease_value) for v in adata_subset.obs[disease_marker]]
    disease_adata = adata_subset[disease_mask].copy()
    
    healthy_mask = [value_equals(v, healthy_value) for v in adata_subset.obs[disease_marker]]
    healthy_adata = adata_subset[healthy_mask].copy()
    
    logger.info(f"Split data into {disease_adata.n_obs} disease cells and {healthy_adata.n_obs} healthy cells")
    
    # Calculate mean expression for disease and healthy cells
    if sp.issparse(disease_adata.X):
        disease_means = disease_adata.X.mean(axis=0).A1
        disease_vars = disease_adata.X.var(axis=0).A1
    else:
        disease_means = disease_adata.X.mean(axis=0)
        disease_vars = disease_adata.X.var(axis=0)
    
    if sp.issparse(healthy_adata.X):
        healthy_means = healthy_adata.X.mean(axis=0).A1
        healthy_vars = healthy_adata.X.var(axis=0).A1
    else:
        healthy_means = healthy_adata.X.mean(axis=0)
        healthy_vars = healthy_adata.X.var(axis=0)
    
    # Create DataFrames with gene statistics
    disease_stats = pd.DataFrame({
        'gene': adata_subset.var_names,
        'mean': disease_means,
        'var': disease_vars,
        'cv': disease_vars / (disease_means + 1e-6)  # Coefficient of variation
    })
    
    healthy_stats = pd.DataFrame({
        'gene': adata_subset.var_names,
        'mean': healthy_means,
        'var': healthy_vars,
        'cv': healthy_vars / (healthy_means + 1e-6)  # Coefficient of variation
    })
    
    # Generate control genes
    disease_control_genes = {}
    healthy_control_genes = {}
    
    for gene in sig_genes:
        # Find control genes for disease group
        if gene in disease_stats['gene'].values:
            gene_idx = disease_stats.index[disease_stats['gene'] == gene].tolist()[0]
            gene_stats = disease_stats.iloc[gene_idx]
            
            # Filter to exclude the gene itself and other significant genes
            filtered_stats = disease_stats[(disease_stats['gene'] != gene) & 
                                          (~disease_stats['gene'].isin(sig_genes))]
            
            if not filtered_stats.empty:
                # Find genes with similar statistics (mean, variance, CV)
                distances = np.sqrt(
                    ((filtered_stats['mean'] - gene_stats['mean']) / (gene_stats['mean'] + 1e-6))**2 +
                    ((filtered_stats['var'] - gene_stats['var']) / (gene_stats['var'] + 1e-6))**2 +
                    ((filtered_stats['cv'] - gene_stats['cv']) / (gene_stats['cv'] + 1e-6))**2
                )
                
                # Get the n closest genes
                sorted_indices = np.argsort(distances)
                control_genes = filtered_stats.iloc[sorted_indices[:n_control_genes]]['gene'].tolist()
                disease_control_genes[gene] = control_genes
            else:
                disease_control_genes[gene] = []
        
        # Find control genes for healthy group
        if gene in healthy_stats['gene'].values:
            gene_idx = healthy_stats.index[healthy_stats['gene'] == gene].tolist()[0]
            gene_stats = healthy_stats.iloc[gene_idx]
            
            # Filter to exclude the gene itself and other significant genes
            filtered_stats = healthy_stats[(healthy_stats['gene'] != gene) & 
                                          (~healthy_stats['gene'].isin(sig_genes))]
            
            if not filtered_stats.empty:
                # Find genes with similar statistics (mean, variance, CV)
                distances = np.sqrt(
                    ((filtered_stats['mean'] - gene_stats['mean']) / (gene_stats['mean'] + 1e-6))**2 +
                    ((filtered_stats['var'] - gene_stats['var']) / (gene_stats['var'] + 1e-6))**2 +
                    ((filtered_stats['cv'] - gene_stats['cv']) / (gene_stats['cv'] + 1e-6))**2
                )
                
                # Get the n closest genes
                sorted_indices = np.argsort(distances)
                control_genes = filtered_stats.iloc[sorted_indices[:n_control_genes]]['gene'].tolist()
                healthy_control_genes[gene] = control_genes
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