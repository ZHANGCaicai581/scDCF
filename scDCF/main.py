# scDCF/main.py

import os
import sys
import argparse
import logging
import warnings
import scanpy as sc
import pandas as pd
import numpy as np
import anndata

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Import functions from your package modules
from scDCF.dependencies import check_and_install_dependencies
from scDCF.utils import (
    read_gene_symbols,
    filter_valid_genes,
    load_control_genes
)
from scDCF.control_genes import generate_control_genes
try:
    from scDCF.parallel import auto_monte_carlo as monte_carlo_comparison
except ImportError:
    from scDCF.analysis import monte_carlo_comparison
from scDCF.post_analysis import (
    load_monte_carlo_results,
    combine_p_values_across_iterations,
    visualize_combined_p_values,
    organize_results,
    export_final_celltype_summary
)

def setup_logging(log_file=None):
    """Set up logging configuration"""
    log_format = '%(asctime)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_format)
    
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter(log_format))
        logging.getLogger().addHandler(file_handler)

def read_gene_list(file_path):
    """Read gene list from a file (supports both CSV and TXT formats)"""
    # Use the read_gene_symbols utility function instead of duplicating logic
    return read_gene_symbols(file_path)

def main():
    """
    Main function to execute the scDCF analysis workflow.
    """
    # Check and install missing dependencies
    check_and_install_dependencies()

    parser = argparse.ArgumentParser(description='scDCF Analysis')
    parser.add_argument('--csv_file', 
                        help='Path to the file containing gene symbols (CSV, TXT, or TSV). Required if --gene_list_file is not provided.')
    parser.add_argument('--gene_list_file', 
                        help='Path to a text file with gene names (one per line). Required if --csv_file is not provided.')
    parser.add_argument('--h5ad_file', required=True, 
                        help='Path to the AnnData h5ad file.')
    parser.add_argument('--output_dir', required=True, 
                        help='Directory for output files.')
    parser.add_argument('--celltype_column', default='celltype_major', 
                        help='Column name containing cell type labels in AnnData obs.')
    parser.add_argument('--cell_types', nargs='+', 
                        help='One or more cell types to analyze.')
    parser.add_argument('--disease_marker', default='disease_numeric', 
                        help='Column name containing disease status in AnnData obs.')
    parser.add_argument('--disease_value', 
                        help='Value indicating disease cells. Can be number or text.')
    parser.add_argument('--healthy_value', 
                        help='Value indicating healthy cells. Can be number or text.')
    parser.add_argument('--rna_count_column', default='nCount_RNA', 
                        help='Column name containing RNA counts in AnnData obs.')
    parser.add_argument('--iterations', type=int, default=10, 
                        help='Number of Monte Carlo iterations.')
    parser.add_argument('--log_file', 
                        help='Path to log file. If not provided, logs will only be written to stdout.')
    parser.add_argument('--show_progress', action='store_true', 
                        help='Show progress bars.')
    parser.add_argument('--control_genes_file', 
                        help='Path to existing control genes JSON file. If not provided, control genes will be generated.')
    parser.add_argument('--control_genes_dir', 
                        help='Directory to save generated control genes. Required if generating new control genes.')
    parser.add_argument('--parallel', action='store_true',
                        help='Enable parallel execution (auto-select worker pool unless overridden).')
    parser.add_argument('--parallel_workers', type=int,
                        help='Number of worker processes for parallel Monte Carlo (implies --parallel).')
    parser.add_argument('--serial', action='store_true',
                        help='Force serial execution even if parallel processing is available.')
    parser.add_argument('--top_n', type=int, default=1000, help='Number of top genes to select based on ZSTAT. Default is 1000.')
    parser.add_argument('--step', default='all', choices=['all', 'monte_carlo', 'post_analysis'],
                        help='Analysis step to run. Default is "all".')
    parser.add_argument('--no_metadata', action='store_true',
                        help='Skip adding AnnData.obs metadata to final summaries (default includes all columns).')
    parser.add_argument('--metadata_columns', nargs='+',
                        help='Specific metadata columns from AnnData.obs to include (default: all columns unless --no_metadata).')

    args = parser.parse_args()

    # Set up logging
    setup_logging(args.log_file)

    if args.parallel and args.serial:
        parser.error("Cannot use --parallel and --serial at the same time.")

    # Check required arguments
    if not args.csv_file and not args.gene_list_file:
        parser.error("Either --csv_file or --gene_list_file must be provided.")

    # Load gene list using the utility function
    if args.csv_file:
        significant_genes_df = read_gene_list(args.csv_file)
    else:
        significant_genes_df = read_gene_list(args.gene_list_file)

    # Normalize column names for consistency
    significant_genes_df = significant_genes_df.copy()
    significant_genes_df.columns = significant_genes_df.columns.str.lower().str.strip()

    if 'gene_name' not in significant_genes_df.columns:
        parser.error("Input gene file must contain a 'gene_name' column.")

    if 'zstat' not in significant_genes_df.columns:
        logging.warning("No 'zstat' column found. Defaulting all z-stat values to 1.0.")
        significant_genes_df['zstat'] = 1.0
    else:
        significant_genes_df['zstat'] = pd.to_numeric(significant_genes_df['zstat'], errors='coerce').fillna(0.0)

    significant_genes_df['gene_name'] = significant_genes_df['gene_name'].astype(str).str.strip()
    significant_genes_df = significant_genes_df[significant_genes_df['gene_name'] != ""]
    significant_genes_df = significant_genes_df.drop_duplicates(subset='gene_name')

    if args.top_n and args.top_n > 0 and len(significant_genes_df) > args.top_n:
        significant_genes_df = (
            significant_genes_df
            .assign(_abs_z=np.abs(significant_genes_df['zstat']))
            .sort_values('_abs_z', ascending=False)
            .head(args.top_n)
            .drop(columns='_abs_z')
            .reset_index(drop=True)
        )
        logging.info(f"Selected top {len(significant_genes_df)} genes by |zstat| (requested top_n={args.top_n}).")
    else:
        logging.info(f"Using {len(significant_genes_df)} genes from input list.")

    # Load data
    logging.info(f"Loading AnnData from {args.h5ad_file}")
    try:
        adata = sc.read_h5ad(args.h5ad_file)
    except Exception as e:
        logging.error(f"Failed to read AnnData file: {e}")
        sys.exit(1)

    # Ensure required columns are in adata.obs
    required_columns = [args.disease_marker, args.celltype_column, args.rna_count_column]
    missing_columns = [col for col in required_columns if col not in adata.obs.columns]
    if missing_columns:
        logging.error(f"Columns {missing_columns} not found in the AnnData object's obs.")
        sys.exit(1)

    # If cell_types is not provided, use unique values from the celltype_column
    if args.cell_types is None:
        args.cell_types = adata.obs[args.celltype_column].unique().tolist()
        logging.info(f"Cell types not provided. Using unique values from '{args.celltype_column}': {args.cell_types}")

    # Ensure all cell_types are strings and filter out NaNs
    args.cell_types = [str(cell_type) for cell_type in args.cell_types if pd.notna(cell_type)]

    # Check if disease_value and healthy_value need to be converted to numeric
    if args.disease_value is not None:
        try:
            disease_value = int(args.disease_value)
        except ValueError:
            try:
                disease_value = float(args.disease_value)
            except ValueError:
                # Keep as string if not numeric
                disease_value = args.disease_value
    else:
        disease_value = 1  # Default
    
    if args.healthy_value is not None:
        try:
            healthy_value = int(args.healthy_value)
        except ValueError:
            try:
                healthy_value = float(args.healthy_value)
            except ValueError:
                # Keep as string if not numeric
                healthy_value = args.healthy_value
    else:
        healthy_value = 0  # Default

    # Prepare control genes configuration
    preloaded_control_genes = None
    if args.control_genes_file:
        preloaded_control_genes = load_control_genes(args.control_genes_file)

    control_genes_dir = args.control_genes_dir
    if not control_genes_dir and not args.control_genes_file:
        control_genes_dir = os.path.join(args.output_dir, "control_genes")

    control_genes_cache = {}

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Run analysis for each cell type
    logging.info(f"Starting analysis for {len(args.cell_types)} cell types")
    
    for cell_type in args.cell_types:
        logging.info(f"Processing cell type: {cell_type}")
        
        # Resolve control genes for this cell type
        if preloaded_control_genes is not None:
            disease_control_genes, healthy_control_genes = preloaded_control_genes
        else:
            disease_control_genes, healthy_control_genes = control_genes_cache.get(cell_type, (None, None))

            if disease_control_genes is None or healthy_control_genes is None:
                cell_type_safe = str(cell_type).replace('/', '_').replace(' ', '_')
                control_genes_file = None

                if control_genes_dir:
                    control_genes_file = os.path.join(control_genes_dir, f"{cell_type_safe}_control_genes.json")
                    if os.path.exists(control_genes_file):
                        disease_control_genes, healthy_control_genes = load_control_genes(control_genes_file)

                if not disease_control_genes or not healthy_control_genes:
                    if control_genes_dir:
                        os.makedirs(control_genes_dir, exist_ok=True)
                        logging.info(f"Generating control genes for {cell_type} (saving to {control_genes_dir}).")
                    else:
                        logging.info(f"Generating control genes for {cell_type} (in-memory).")

                    disease_control_genes, healthy_control_genes = generate_control_genes(
                        adata=adata,
                        significant_genes_df=significant_genes_df,
                        cell_type=cell_type,
                        cell_type_column=args.celltype_column,
                        disease_marker=args.disease_marker,
                        disease_value=disease_value,
                        healthy_value=healthy_value,
                        output_dir=control_genes_dir
                    )

                control_genes_cache[cell_type] = (disease_control_genes, healthy_control_genes)

        # Step: Monte Carlo Analysis
        if args.step in ['all', 'monte_carlo']:
            logging.info(f"Running Monte Carlo comparison for {cell_type}.")
            
            additional_kwargs = {}
            if args.serial:
                additional_kwargs['use_parallel'] = False
            elif args.parallel or args.parallel_workers is not None:
                additional_kwargs['use_parallel'] = True
            
            if args.parallel_workers is not None:
                additional_kwargs['n_workers'] = args.parallel_workers

            # Run for disease group
            disease_results = monte_carlo_comparison(
                adata=adata,
                cell_type=cell_type,
                cell_type_column=args.celltype_column,
                significant_genes_df=significant_genes_df,
                disease_control_genes=disease_control_genes,
                healthy_control_genes=healthy_control_genes,
                output_dir=args.output_dir,
                rna_count_column=args.rna_count_column,
                iterations=args.iterations,
                target_group="disease",
                disease_marker=args.disease_marker,
                disease_value=disease_value,
                healthy_value=healthy_value,
                show_progress=args.show_progress,
                **additional_kwargs
            )
            
            # Run for healthy group
            healthy_results = monte_carlo_comparison(
                adata=adata,
                cell_type=cell_type,
                cell_type_column=args.celltype_column,
                significant_genes_df=significant_genes_df,
                disease_control_genes=disease_control_genes,
                healthy_control_genes=healthy_control_genes,
                output_dir=args.output_dir,
                rna_count_column=args.rna_count_column,
                iterations=args.iterations,
                target_group="healthy",
                disease_marker=args.disease_marker,
                disease_value=disease_value,
                healthy_value=healthy_value,
                show_progress=args.show_progress,
                **additional_kwargs
            )

        # Step: Post-Analysis
        if args.step in ['all', 'post_analysis']:
            logging.info(f"Post-analysis for {cell_type}.")

            # Load Monte Carlo results for both disease and healthy groups (by file path)
            disease_file = os.path.join(args.output_dir, cell_type, f"{cell_type}_disease_monte_carlo_results.csv")
            healthy_file = os.path.join(args.output_dir, cell_type, f"{cell_type}_healthy_monte_carlo_results.csv")
            disease_results = load_monte_carlo_results(disease_file)
            healthy_results = load_monte_carlo_results(healthy_file)

            # Check if either dataset is empty before proceeding
            if disease_results.empty or healthy_results.empty:
                logging.warning(f"One of the result DataFrames is empty for {cell_type}. Cannot perform post-analysis.")
                # Skip to the next cell type instead of using continue which would exit the function
                continue

            # Combine p-values across iterations
            disease_combined = combine_p_values_across_iterations(disease_results, args.output_dir, cell_type, 'disease')
            healthy_combined = combine_p_values_across_iterations(healthy_results, args.output_dir, cell_type, 'healthy')

            # Visualize combined p-values and significant cell counts
            visualize_combined_p_values(disease_combined, healthy_combined, cell_type, args.output_dir)

            export_final_celltype_summary(
                cell_type=cell_type,
                disease_combined=disease_combined,
                healthy_combined=healthy_combined,
                output_dir=args.output_dir,
                include_metadata=not args.no_metadata if args.metadata_columns is None else True,
                adata=adata if not args.no_metadata else None,
                metadata_columns=args.metadata_columns
            )

    logging.info("Analysis complete.")

def organize_output(source_dir, dest_dir="organized_output", cell_types=None):
    """
    Organize analysis results into a clean, cell-type-centric structure.
    
    This function takes the raw output from scDCF analysis and organizes it into a
    well-structured directory format that's easy to navigate and interpret.
    
    Args:
        source_dir (str): Directory containing the raw analysis results
        dest_dir (str): Directory to create the organized structure in
        cell_types (list, optional): List of cell types to organize. If None, will autodetect.
    
    Returns:
        str: Path to the organized output directory
    
    Example:
        >>> import scDCF
        >>> # After running your analysis with raw output in "results_dir"
        >>> scDCF.organize_output("results_dir", "clean_results") 
    """
    return organize_results(source_dir, dest_dir, cell_types)

if __name__ == "__main__":
    main()
