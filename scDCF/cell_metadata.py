"""
Helper functions for adding cell metadata to scDCF results.

This module makes it easy for researchers to enhance results with cell information
from their original AnnData object.
"""

import pandas as pd
import logging

def add_cell_metadata(results_df, adata, cell_id_column='cell_id', metadata_columns=None):
    """
    Add cell metadata from AnnData to scDCF results.
    
    This function enriches scDCF results with cell-level metadata from the original
    AnnData object, making it easier to interpret and filter results.
    
    Args:
        results_df (pd.DataFrame): Results from monte_carlo_comparison
        adata (AnnData): Original AnnData object used in analysis
        cell_id_column (str): Column name containing cell IDs in results_df
        metadata_columns (list): Specific metadata columns to add (None = all from adata.obs)
        
    Returns:
        pd.DataFrame: Enhanced results with cell metadata
        
    Example:
        >>> import scDCF
        >>> import scanpy as sc
        >>> 
        >>> # Run analysis
        >>> adata = sc.read_h5ad("data.h5ad")
        >>> results = scDCF.monte_carlo_comparison(...)
        >>> 
        >>> # Add metadata
        >>> enhanced = scDCF.add_cell_metadata(results, adata)
        >>> 
        >>> # Now has: cell_id, p_value, celltype, sample, batch, etc.
        >>> enhanced.to_csv("results_with_metadata.csv")
        >>> 
        >>> # Filter by metadata
        >>> significant_in_sample1 = enhanced[
        ...     (enhanced['p_value'] < 0.05) & 
        ...     (enhanced['sample'] == 'Sample1')
        ... ]
    """
    logging.info("Adding cell metadata to results...")
    
    # Make a copy to avoid modifying original
    results_enhanced = results_df.copy()
    
    # Get metadata from adata
    if metadata_columns is None:
        # Use all columns from adata.obs
        metadata = adata.obs.copy()
    else:
        # Use only specified columns
        metadata = adata.obs[metadata_columns].copy()
    
    # Add the obs_names as a column for merging
    metadata['original_cell_name'] = adata.obs_names
    metadata['cell_id_for_merge'] = adata.obs_names
    
    # Ensure cell_id types match for merging
    results_enhanced['cell_id_for_merge'] = results_enhanced[cell_id_column].astype(str)
    metadata['cell_id_for_merge'] = metadata['cell_id_for_merge'].astype(str)
    
    # Merge
    enhanced = pd.merge(
        results_enhanced,
        metadata,
        on='cell_id_for_merge',
        how='left',
        suffixes=('', '_from_adata')
    )
    
    # Clean up merge column
    enhanced = enhanced.drop('cell_id_for_merge', axis=1)
    
    # Report
    n_matched = enhanced['original_cell_name'].notna().sum()
    logging.info(f"✅ Metadata added for {n_matched}/{len(enhanced)} cells")
    
    if metadata_columns:
        logging.info(f"   Columns added: {metadata_columns}")
    else:
        new_cols = [c for c in enhanced.columns if c not in results_df.columns]
        logging.info(f"   Columns added: {len(new_cols)} metadata fields")
    
    return enhanced


def create_cell_id_mapping(adata, cell_type=None, cell_type_column=None,
                           disease_marker=None, disease_value=None):
    """
    Create a complete cell ID mapping file for reference.
    
    This creates a standalone file mapping cell IDs to all available metadata,
    useful for downstream analysis and result interpretation.
    
    Args:
        adata (AnnData): AnnData object
        cell_type (str): Specific cell type to map (None = all cells)
        cell_type_column (str): Column name for cell types
        disease_marker (str): Column name for disease status
        disease_value: Value indicating disease cells
        
    Returns:
        pd.DataFrame: Cell ID mapping with metadata
        
    Example:
        >>> import scDCF
        >>> import scanpy as sc
        >>> 
        >>> adata = sc.read_h5ad("data.h5ad")
        >>> 
        >>> # Create mapping for all cells
        >>> mapping = scDCF.create_cell_id_mapping(adata)
        >>> mapping.to_csv("cell_id_mapping.csv", index=False)
        >>> 
        >>> # Or for specific cell type
        >>> mapping = scDCF.create_cell_id_mapping(
        ...     adata, 
        ...     cell_type="T_cell",
        ...     cell_type_column="celltype"
        ... )
    """
    logging.info("Creating cell ID mapping...")
    
    # Filter to specific cell type if requested
    if cell_type and cell_type_column:
        adata_subset = adata[adata.obs[cell_type_column] == cell_type]
    else:
        adata_subset = adata
    
    # Further filter by disease status if requested
    if disease_marker and disease_value is not None:
        adata_subset = adata_subset[adata_subset.obs[disease_marker] == disease_value]
    
    # Create mapping DataFrame
    mapping = adata_subset.obs.copy()
    mapping['cell_id'] = adata_subset.obs_names
    mapping['cell_index_position'] = range(len(adata_subset))
    
    # Reorder columns (cell_id first)
    cols = ['cell_id'] + [c for c in mapping.columns if c != 'cell_id']
    mapping = mapping[cols]
    
    logging.info(f"✅ Cell ID mapping created for {len(mapping)} cells")
    logging.info(f"   Columns: {len(mapping.columns)} metadata fields")
    
    return mapping


def validate_cell_ids(results_df, adata, cell_id_column='cell_id'):
    """
    Validate that cell IDs in results match cells in AnnData.
    
    This is a quality control function to ensure results can be traced
    back to the original data.
    
    Args:
        results_df (pd.DataFrame): Results from scDCF
        adata (AnnData): Original AnnData object
        cell_id_column (str): Column with cell IDs in results
        
    Returns:
        dict: Validation report
        
    Example:
        >>> report = scDCF.validate_cell_ids(results, adata)
        >>> print(report['all_matched'])  # True if all IDs valid
    """
    result_ids = set(results_df[cell_id_column].astype(str).unique())
    adata_ids = set(adata.obs_names.astype(str))
    
    matched = result_ids & adata_ids
    missing = result_ids - adata_ids
    
    report = {
        'total_result_ids': len(result_ids),
        'total_adata_ids': len(adata_ids),
        'matched': len(matched),
        'missing': len(missing),
        'match_percentage': 100 * len(matched) / len(result_ids) if result_ids else 0,
        'all_matched': len(missing) == 0,
        'missing_ids': list(missing) if missing else []
    }
    
    if report['all_matched']:
        logging.info(f"✅ All {len(result_ids)} cell IDs validated successfully")
    else:
        logging.warning(f"⚠️  {len(missing)} cell IDs not found in AnnData:")
        logging.warning(f"   Missing IDs: {list(missing)[:10]}...")
    
    return report


def get_cell_info_by_id(adata, cell_ids):
    """
    Get complete information for specific cells.
    
    Args:
        adata (AnnData): AnnData object
        cell_ids (list): List of cell IDs to look up
        
    Returns:
        pd.DataFrame: Cell information
        
    Example:
        >>> # Get info for significant cells
        >>> sig_cells = results[results['p_value'] < 0.05]['cell_id']
        >>> cell_info = scDCF.get_cell_info_by_id(adata, sig_cells)
        >>> print(cell_info[['celltype', 'sample', 'batch']])
    """
    cell_ids_str = [str(cid) for cid in cell_ids]
    
    # Filter adata to these cells
    mask = adata.obs_names.isin(cell_ids_str)
    cells_subset = adata[mask]
    
    # Get metadata
    info = cells_subset.obs.copy()
    info['cell_id'] = cells_subset.obs_names
    
    logging.info(f"Retrieved information for {len(info)} cells")
    
    return info

