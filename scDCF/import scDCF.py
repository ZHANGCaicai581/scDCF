#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
import matplotlib.pyplot as plt
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, 
                   format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

def generate_test_data(n_cells=500, n_genes=200, n_cell_types=3, save_dir='test_data'):
    """Generate simulated scRNA-seq data for testing scDCF."""
    os.makedirs(save_dir, exist_ok=True)
    
    # Generate random gene expression data
    np.random.seed(42)  # For reproducibility
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes))
    
    # Create gene names and cell barcodes
    gene_names = [f"gene_{i}" for i in range(n_genes)]
    cell_barcodes = [f"cell_{i}" for i in range(n_cells)]
    
    # Create AnnData object
    adata = AnnData(X, var=pd.DataFrame(index=gene_names), obs=pd.DataFrame(index=cell_barcodes))
    
    # Add cell type annotation
    cell_types = ['T_cell', 'B_cell', 'NK_cell'][:n_cell_types]
    adata.obs['cell_type'] = np.random.choice(cell_types, size=n_cells)
    
    # Add disease status (0 = healthy, 1 = disease)
    adata.obs['disease_numeric'] = np.random.choice([0, 1], size=n_cells, p=[0.4, 0.6])
    
    # Add RNA count information (required by scDCF)
    adata.obs['nCount_RNA'] = np.random.randint(1000, 5000, size=n_cells)
    
    # Save the AnnData object
    adata_file = os.path.join(save_dir, 'test_adata.h5ad')
    adata.write_h5ad(adata_file)
    logging.info(f"AnnData saved to {adata_file}")
    
    # Create a dataframe of significant genes with z-statistics
    sig_genes = np.random.choice(gene_names, size=20, replace=False)
    gene_df = pd.DataFrame({
        'gene_name': sig_genes,
        'zstat': np.random.normal(2.5, 0.5, size=len(sig_genes))
    })
    
    # Save the gene list
    gene_file = os.path.join(save_dir, 'significant_genes.csv')
    gene_df.to_csv(gene_file, index=False)
    logging.info(f"Significant genes saved to {gene_file}")
    
    return adata, gene_df, adata_file, gene_file

def run_scDCF_test(adata_file, gene_file, output_dir='test_output'):
    """Run scDCF on the test data"""
    try:
        import scDCF
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        logging.info("Running scDCF through Python API...")
        
        # Load data
        adata = sc.read_h5ad(adata_file)
        significant_genes_df = scDCF.utils.read_gene_symbols(gene_file)
        
        # Generate control genes
        cell_types = adata.obs['cell_type'].unique().tolist()
        
        for cell_type in cell_types:
            # Generate control genes
            logging.info(f"Generating control genes for {cell_type}...")
            
            # Fixed parameter order: Now passing cell_type correctly
            disease_control_genes, healthy_control_genes = scDCF.control_genes.generate_control_genes(
                adata=adata, 
                significant_genes_df=significant_genes_df,
                cell_type=cell_type,  # Passing a string, not DataFrame
                cell_type_column='cell_type',
                output_dir=output_dir
            )
            
            # Run monte carlo comparison
            logging.info(f"Running monte carlo comparison for {cell_type}...")
            disease_results = scDCF.analysis.monte_carlo_comparison(
                adata=adata,
                cell_type=cell_type,
                cell_type_column='cell_type',
                significant_genes_df=significant_genes_df,
                disease_control_genes=disease_control_genes,
                healthy_control_genes=healthy_control_genes,
                output_dir=output_dir,
                rna_count_column='nCount_RNA',
                iterations=2,  # Use a small number for testing
                disease_marker='disease_numeric',
                target_group='disease',
                show_progress=True
            )
            
            healthy_results = scDCF.analysis.monte_carlo_comparison(
                adata=adata,
                cell_type=cell_type,
                cell_type_column='cell_type',
                significant_genes_df=significant_genes_df,
                disease_control_genes=disease_control_genes,
                healthy_control_genes=healthy_control_genes,
                output_dir=output_dir,
                rna_count_column='nCount_RNA',
                iterations=2,  # Use a small number for testing
                disease_marker='disease_numeric',
                target_group='healthy',
                show_progress=True
            )
            
            # Compare results
            if not disease_results.empty and not healthy_results.empty:
                comparison = scDCF.analysis.compare_groups(disease_results, healthy_results)
                logging.info(f"Comparison results: {comparison}")
                
                # Post-analysis
                scDCF.post_analysis.visualize_combined_p_values(
                    disease_results,
                    healthy_results,
                    cell_type,
                    output_dir=os.path.join(output_dir, cell_type)
                )
        
        logging.info("scDCF test completed successfully!")
        return True
        
    except Exception as e:
        logging.error(f"Error running scDCF test: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return False

if __name__ == "__main__":
    # Step 1: Generate test data with smaller size for quicker testing
    adata, gene_df, adata_file, gene_file = generate_test_data(n_cells=200, n_genes=100)
    
    # Step 2: Run scDCF test
    success = run_scDCF_test(adata_file, gene_file)
    
    if success:
        logging.info("✅ Test completed successfully!")
    else:
        logging.error("❌ Test failed!")


