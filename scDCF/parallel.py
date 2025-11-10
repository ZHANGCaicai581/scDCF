"""
Parallel processing utilities for scDCF Monte Carlo analysis.

This module provides parallel execution of Monte Carlo iterations for improved performance.
Typical speedup: 4-8x on multi-core systems.
"""

import os
import logging
import multiprocessing as mp
import numpy as np
import pandas as pd
from functools import partial
import traceback

from .analysis import monte_carlo_comparison as _monte_carlo_single

def _worker_single_iteration(iteration_num, adata_path, genes_df, disease_ctrl, healthy_ctrl,
                             output_dir, cell_type, target_group, config):
    """
    Worker function to process a single Monte Carlo iteration.
    Runs in a separate process.
    """
    try:
        import scanpy as sc
        from scDCF.utils import read_gene_symbols
        
        # Set random seed for reproducibility
        np.random.seed(42 + iteration_num)
        
        # Each worker loads its own data copy
        adata = sc.read_h5ad(adata_path)
        
        # Convert genes_df back to DataFrame if it's a dict
        if isinstance(genes_df, dict):
            significant_genes_df = pd.DataFrame(genes_df)
        else:
            significant_genes_df = genes_df
        
        # Run single iteration
        result = _monte_carlo_single(
            adata=adata,
            cell_type=cell_type,
            cell_type_column=config['cell_type_column'],
            significant_genes_df=significant_genes_df,
            disease_control_genes=disease_ctrl,
            healthy_control_genes=healthy_ctrl,
            output_dir=output_dir,
            disease_marker=config['disease_marker'],
            disease_value=config['disease_value'],
            healthy_value=config['healthy_value'],
            rna_count_column=config['rna_count_column'],
            iterations=1,  # Single iteration per worker
            target_group=target_group,
            batch_size=config.get('batch_size', 500),
            show_progress=False
        )
        
        if result.empty:
            return {'success': False, 'iteration': iteration_num, 'error': 'Empty result'}
        
        # Update iteration number
        result['iteration'] = iteration_num
        
        # Save iteration file
        cell_type_dir = os.path.join(output_dir, cell_type)
        os.makedirs(cell_type_dir, exist_ok=True)
        
        iter_file = os.path.join(
            cell_type_dir,
            f"{cell_type}_{target_group}_monte_carlo_results_iteration{iteration_num}.csv"
        )
        result.to_csv(iter_file, index=False)
        
        return {
            'success': True,
            'iteration': iteration_num,
            'n_cells': len(result),
            'n_significant': (result['p_value'] < 0.05).sum(),
            'result': result
        }
        
    except Exception as e:
        return {
            'success': False,
            'iteration': iteration_num,
            'error': str(e),
            'traceback': traceback.format_exc()
        }

def parallel_monte_carlo_comparison(adata, cell_type, cell_type_column, significant_genes_df,
                                   disease_control_genes=None, healthy_control_genes=None,
                                   output_dir=".", rna_count_column='nCount_RNA',
                                   iterations=10, target_group="disease",
                                   disease_marker='disease_numeric',
                                   disease_value=1, healthy_value=0,
                                   n_workers=None, batch_size=500, show_progress=True):
    """
    Parallel Monte Carlo comparison with automatic multi-core utilization.
    
    This function provides 4-8x speedup over serial execution by running
    Monte Carlo iterations in parallel across multiple CPU cores.
    
    Args:
        adata: AnnData object containing single-cell data
        cell_type: Cell type to analyze
        cell_type_column: Column in adata.obs with cell type labels
        significant_genes_df: DataFrame with 'gene_name' and 'zstat' columns
        disease_control_genes: Dict mapping genes to control genes (disease)
        healthy_control_genes: Dict mapping genes to control genes (healthy)
        output_dir: Directory to save results
        rna_count_column: Column in adata.obs with RNA counts
        iterations: Number of Monte Carlo iterations
        target_group: 'disease' or 'healthy'
        disease_marker: Column in adata.obs with disease status
        disease_value: Value indicating disease cells
        healthy_value: Value indicating healthy cells
        n_workers: Number of parallel workers (None = auto-detect)
        batch_size: Cells to process per batch
        show_progress: Show progress updates
        
    Returns:
        DataFrame: Combined results from all iterations
        
    Example:
        >>> import scDCF
        >>> import scanpy as sc
        >>> 
        >>> adata = sc.read_h5ad("data.h5ad")
        >>> genes_df = scDCF.read_gene_symbols("genes.txt")
        >>> 
        >>> # Generate control genes first
        >>> disease_ctrl, healthy_ctrl = scDCF.generate_control_genes(
        ...     adata, genes_df, "T_cell", "celltype"
        ... )
        >>> 
        >>> # Run parallel Monte Carlo (4-8x faster!)
        >>> results = scDCF.parallel_monte_carlo_comparison(
        ...     adata=adata,
        ...     cell_type="T_cell",
        ...     cell_type_column="celltype",
        ...     significant_genes_df=genes_df,
        ...     disease_control_genes=disease_ctrl,
        ...     healthy_control_genes=healthy_ctrl,
        ...     output_dir="results/",
        ...     iterations=100,
        ...     n_workers=8  # or None for auto
        ... )
    """
    import time
    import tempfile
    
    logging.info(f"Starting parallel Monte Carlo analysis ({iterations} iterations)")
    
    # Auto-detect workers (leave at least one core and cap for fairness)
    total_cpus = mp.cpu_count() or 1
    if n_workers is None or n_workers <= 0:
        if total_cpus <= 2:
            auto_workers = 1
        else:
            auto_workers = max(1, min(total_cpus - 1, 8))
        n_workers = auto_workers
        logging.info(
            f"Auto-selected {n_workers} workers (total CPUs: {total_cpus}, cap=8, reserve=1)"
        )
    else:
        n_workers = max(1, min(n_workers, total_cpus))
        logging.info(
            f"Using user-specified worker count: {n_workers} (total CPUs: {total_cpus})"
        )

    n_workers = min(n_workers, iterations)  # Don't use more workers than iterations
    if n_workers == 0:
        n_workers = 1
    
    logging.info(f"Parallel worker pool size: {n_workers}")
    logging.info(f"Expected speedup: up to ~{n_workers}x")
    
    # Save adata to temp file for workers
    temp_h5ad = tempfile.NamedTemporaryFile(suffix='.h5ad', delete=False)
    temp_h5ad_path = temp_h5ad.name
    temp_h5ad.close()
    
    try:
        # Save adata for workers to load
        adata.write_h5ad(temp_h5ad_path)
        
        # Convert genes_df to dict for pickling
        genes_dict = significant_genes_df.to_dict('list')
        
        # Configuration dict
        config = {
            'cell_type_column': cell_type_column,
            'disease_marker': disease_marker,
            'disease_value': disease_value,
            'healthy_value': healthy_value,
            'rna_count_column': rna_count_column,
            'batch_size': batch_size
        }
        
        # Prepare worker arguments
        worker_args = [
            (i+1, temp_h5ad_path, genes_dict, disease_control_genes, healthy_control_genes,
             output_dir, cell_type, target_group, config)
            for i in range(iterations)
        ]
        
        # Run parallel processing
        start_time = time.time()
        results = []
        failed = []
        
        with mp.Pool(processes=n_workers) as pool:
            for i, result_dict in enumerate(pool.starmap(_worker_single_iteration, worker_args)):
                if result_dict['success']:
                    results.append(result_dict)
                    if show_progress:
                        logging.info(
                            f"✅ Iteration {result_dict['iteration']}/{iterations} complete "
                            f"({i+1}/{iterations} finished)"
                        )
                else:
                    failed.append(result_dict)
                    logging.error(f"❌ Iteration {result_dict['iteration']} failed: {result_dict['error']}")
        
        elapsed = time.time() - start_time
        
        logging.info(f"\n{'='*80}")
        logging.info(f"Parallel processing complete: {target_group}")
        logging.info(f"Time: {elapsed/60:.2f} minutes")
        logging.info(f"Successful: {len(results)}/{iterations}")
        logging.info(f"Failed: {len(failed)}/{iterations}")
        
        if not results:
            logging.error("All iterations failed!")
            return pd.DataFrame()
        
        # Combine results
        all_dfs = [r['result'] for r in results if r['result'] is not None]
        combined = pd.concat(all_dfs, ignore_index=True)
        
        # Save combined
        cell_type_dir = os.path.join(output_dir, cell_type)
        combined_file = os.path.join(cell_type_dir, f"{cell_type}_{target_group}_monte_carlo_results.csv")
        combined.to_csv(combined_file, index=False)
        
        logging.info(f"Combined results saved: {combined_file}")
        logging.info(f"Total cells: {combined['cell_id'].nunique()}")
        
        # Statistics
        sig_count = (combined['p_value'] < 0.05).sum()
        fdr_count = (combined['p_value_adj'] < 0.05).sum()
        logging.info(f"Nominal p<0.05: {sig_count} ({100*sig_count/len(combined):.2f}%)")
        logging.info(f"FDR<0.05: {fdr_count} ({100*fdr_count/len(combined):.2f}%)")
        
        return combined
        
    finally:
        # Clean up temp file
        try:
            os.unlink(temp_h5ad_path)
        except:
            pass

# Convenience function
def auto_monte_carlo(adata, cell_type, cell_type_column, significant_genes_df,
                    disease_control_genes, healthy_control_genes, output_dir,
                    iterations=10, use_parallel=None, **kwargs):
    """
    Automatically choose between serial and parallel execution.
    
    Defaults to serial execution for predictable resource usage.
    Parallel mode is enabled only when explicitly requested via
    `use_parallel=True` or by providing `n_workers`.
    
    Args:
        use_parallel: Force parallel (True) or serial (False), or None for auto
        **kwargs: Additional arguments passed to monte_carlo_comparison
        
    Returns:
        DataFrame: Combined results
    """
    # Auto-detect
    total_cpus = mp.cpu_count() or 1
    n_workers = kwargs.pop('n_workers', None)

    if use_parallel is None:
        # If user specified worker count we respect it and enable parallel,
        # otherwise stay in serial mode by default.
        use_parallel = n_workers is not None
    
    if use_parallel and iterations >= 2:
        logging.info("Using parallel processing (auto-monte-carlo)")
        return parallel_monte_carlo_comparison(
            adata, cell_type, cell_type_column, significant_genes_df,
            disease_control_genes, healthy_control_genes, output_dir,
            iterations=iterations, n_workers=n_workers, **kwargs
        )

    logging.info("Using serial processing (auto-monte-carlo)")
    return _monte_carlo_single(
        adata, cell_type, cell_type_column, significant_genes_df,
        disease_control_genes, healthy_control_genes, output_dir,
        iterations=iterations, **kwargs
    )

