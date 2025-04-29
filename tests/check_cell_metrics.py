#!/usr/bin/env python
import pandas as pd
import os
import logging
import glob
import sys
import platform
import argparse

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('root')

def check_cell_metrics_files(output_dir="output"):
    """Check the structure and content of cell_metrics.csv files"""
    logger.info(f"Checking cell_metrics.csv files in {output_dir}")
    logger.info(f"Running on: {platform.system()} {platform.release()}")
    
    # Enhanced directory check
    if not os.path.exists(output_dir):
        logger.error(f"❌ Output directory '{output_dir}' not found!")
        logger.error(f"Please make sure you've run organize_final_output() first.")
        logger.error(f"You can specify a different directory with --output_dir argument.")
        print(f"\nERROR: Directory '{output_dir}' not found. Please check your path.")
        return False
    
    # Detect available cell types by listing subdirectories
    available_cell_types = []
    default_cell_types = ['T_cell', 'B_cell', 'NK_cell']
    
    for item in os.listdir(output_dir):
        item_path = os.path.join(output_dir, item)
        if os.path.isdir(item_path):
            # Check if this directory contains cell_metrics.csv
            if os.path.exists(os.path.join(item_path, "cell_metrics.csv")):
                available_cell_types.append(item)
    
    # If no cell types found, check for default ones and warn
    if not available_cell_types:
        logger.warning(f"⚠️ No directories with cell_metrics.csv found in '{output_dir}'.")
        
        # Check if default cell type directories exist
        any_default_dir = False
        for cell_type in default_cell_types:
            if os.path.exists(os.path.join(output_dir, cell_type)):
                any_default_dir = True
                logger.warning(f"  Directory '{cell_type}' exists but doesn't contain cell_metrics.csv")
        
        if any_default_dir:
            logger.warning(f"Some default cell type directories exist but don't have cell_metrics.csv files.")
        else:
            logger.warning(f"Expected structure: {output_dir}/[CELL_TYPE]/cell_metrics.csv")
        
        print(f"\nWARNING: No cell metrics files found in '{output_dir}'. Is this the correct directory?")
        return False
    
    logger.info(f"Found {len(available_cell_types)} cell types: {', '.join(available_cell_types)}")
    
    # Process each available cell type
    for cell_type in available_cell_types:
        # Use os.path.join for cross-platform path handling
        cell_metrics_file = os.path.join(output_dir, cell_type, "cell_metrics.csv")
        
        if not os.path.exists(cell_metrics_file):
            logger.error(f"❌ cell_metrics.csv not found for {cell_type}")
            continue
            
        try:
            # Load the file
            df = pd.read_csv(cell_metrics_file)
            
            # Check basic info
            logger.info(f"\n{cell_type} cell_metrics.csv:")
            logger.info(f"  Contains {len(df)} cells with {len(df.columns)} metrics")
            logger.info(f"  Columns: {df.columns.tolist()}")
            
            # Check for duplicate cell IDs - IMPORTANT CHECK
            if 'cell_id' in df.columns:
                unique_ids = df['cell_id'].nunique()
                total_rows = len(df)
                if unique_ids < total_rows:
                    logger.error(f"  ❌ DUPLICATE CELL IDs FOUND: {total_rows - unique_ids} duplicates detected")
                    # Show the duplicated IDs
                    duplicated = df[df.duplicated(subset=['cell_id'], keep=False)]
                    logger.error(f"  Duplicated cell IDs: {duplicated['cell_id'].unique().tolist()}")
                else:
                    logger.info(f"  ✅ All cell IDs are unique ({unique_ids} unique IDs)")
            
            # Check source files as well
            supporting_dir = os.path.join(output_dir, cell_type, "supporting_data")
            disease_file = os.path.join(supporting_dir, f"{cell_type}_disease_monte_carlo_results.csv")
            healthy_file = os.path.join(supporting_dir, f"{cell_type}_healthy_monte_carlo_results.csv")
            
            if os.path.exists(disease_file):
                try:
                    disease_df = pd.read_csv(disease_file)
                    if 'cell_id' in disease_df.columns:
                        unique_disease_ids = disease_df['cell_id'].nunique()
                        disease_rows = len(disease_df)
                        if unique_disease_ids < disease_rows:
                            logger.error(f"  ❌ DUPLICATE CELL IDs IN DISEASE FILE: {disease_rows - unique_disease_ids} duplicates")
                        else:
                            logger.info(f"  ✅ All disease cell IDs are unique ({unique_disease_ids} IDs)")
                except Exception as e:
                    logger.error(f"  Error reading disease file: {e}")
            
            if os.path.exists(healthy_file):
                try:
                    healthy_df = pd.read_csv(healthy_file)
                    if 'cell_id' in healthy_df.columns:
                        unique_healthy_ids = healthy_df['cell_id'].nunique()
                        healthy_rows = len(healthy_df)
                        if unique_healthy_ids < healthy_rows:
                            logger.error(f"  ❌ DUPLICATE CELL IDs IN HEALTHY FILE: {healthy_rows - unique_healthy_ids} duplicates")
                        else:
                            logger.info(f"  ✅ All healthy cell IDs are unique ({unique_healthy_ids} IDs)")
                except Exception as e:
                    logger.error(f"  Error reading healthy file: {e}")
            
            # Check required columns
            essential_cols = ['cell_id', 't_stat', 'p_value']
            missing = [col for col in essential_cols if col not in df.columns]
            if missing:
                logger.warning(f"  ❌ Missing essential columns: {missing}")
            else:
                logger.info(f"  ✅ Contains all essential columns")
                
            # Check for enhanced scores
            z_score_cols = [col for col in df.columns if 'z_score' in col.lower()]
            log_cols = [col for col in df.columns if 'log' in col.lower()]
            
            if z_score_cols:
                logger.info(f"  ✅ Contains z-score columns: {z_score_cols}")
            else:
                logger.warning(f"  ❌ No z-score columns found")
                
            if log_cols:
                logger.info(f"  ✅ Contains -log10 columns: {log_cols}")
            else:
                logger.warning(f"  ❌ No -log10 columns found")
            
            # Check for trait scores
            trait_cols = [col for col in df.columns if 'trait_' in col.lower()]
            if trait_cols:
                logger.info(f"  ✅ Contains trait score columns: {trait_cols}")
                
                # Verify trait scores are consistent for all cells
                for col in trait_cols:
                    unique_values = df[col].nunique()
                    if unique_values == 1:
                        logger.info(f"  ✅ {col} has consistent value across all cells (as expected)")
                    else:
                        logger.warning(f"  ⚠️ {col} has {unique_values} different values across cells")
            else:
                logger.warning(f"  ❌ No trait score columns found")
            
            # Sample data
            if len(df) > 0:
                logger.info(f"\n  Sample rows from {cell_type} cell_metrics.csv:")
                logger.info(f"{df.head(3)}")
            
            # Check supporting files
            if os.path.exists(supporting_dir):
                # Use cross-platform pattern for glob
                supporting_files = glob.glob(os.path.join(supporting_dir, "*.csv"))
                logger.info(f"\n  Supporting files ({len(supporting_files)} found):")
                for file in supporting_files[:5]:  # Show first 5
                    # Skip mentioning KS test files if they happen to be there
                    if "ks_test" not in file.lower():
                        logger.info(f"    - {os.path.basename(file)}")
                if len(supporting_files) > 5:
                    logger.info(f"    - ... and {len(supporting_files)-5} more")
                
                # Check Monte Carlo iterations
                mc_dir = os.path.join(supporting_dir, "monte_carlo_iterations")
                if os.path.exists(mc_dir):
                    mc_files = glob.glob(os.path.join(mc_dir, "*.csv"))
                    logger.info(f"  Monte Carlo iterations: {len(mc_files)} files")
                else:
                    logger.warning(f"  ❌ Monte Carlo iterations folder not found")
            
            # Check for control genes
            control_genes_dir = os.path.join(supporting_dir, "control_genes")
            if os.path.exists(control_genes_dir):
                control_files = glob.glob(os.path.join(control_genes_dir, "*.json"))
                if control_files:
                    logger.info(f"  ✅ Found {len(control_files)} control gene files in control_genes/")
                    for file in control_files:
                        logger.info(f"    - {os.path.basename(file)}")
                else:
                    logger.warning(f"  ⚠️ Control genes directory exists but no JSON files found")
            else:
                # Try to find control genes at root level (backward compatibility)
                root_control_file = os.path.join(output_dir, f"{cell_type}_control_genes.json")
                if os.path.exists(root_control_file):
                    logger.info(f"  ✅ Found control genes file at root level (legacy format)")
                    logger.warning(f"  ⚠️ Consider moving it to {cell_type}/supporting_data/control_genes/")
                else:
                    logger.warning(f"  ❌ No control genes folder or files found")
                
        except pd.errors.EmptyDataError:
            logger.error(f"  ❌ The file {cell_metrics_file} is empty")
        except pd.errors.ParserError:
            logger.error(f"  ❌ Parsing error in {cell_metrics_file} - not a valid CSV file")
        except Exception as e:
            logger.error(f"Error analyzing cell_metrics.csv for {cell_type}: {e}")
    
    return True

def main():
    """Main function with command line argument parsing"""
    parser = argparse.ArgumentParser(description="Check cell_metrics.csv files in the output directory")
    parser.add_argument("--output_dir", "-o", default="output", 
                       help="Directory containing the organized output (default: 'output')")
    parser.add_argument("--create_example", action="store_true",
                       help="Create an example of the expected directory structure")
    
    args = parser.parse_args()
    
    if args.create_example:
        create_example_structure()
        return
    
    # Show usage reminder
    print(f"Checking cell metrics files in: {args.output_dir}")
    print(f"You can specify a different directory with --output_dir YOUR_PATH")
    
    # Run the check
    success = check_cell_metrics_files(args.output_dir)
    
    if not success:
        sys.exit(1)  # Exit with error code for scripts

def create_example_structure():
    """Create an example directory structure to show users what's expected"""
    example_dir = "example_output"
    
    if os.path.exists(example_dir):
        print(f"Example directory {example_dir} already exists.")
        return
    
    # Create directories for sample cell types
    for cell_type in ['B_cell', 'T_cell', 'Macrophages']:  # Added a non-default cell type
        os.makedirs(os.path.join(example_dir, cell_type, "supporting_data", "monte_carlo_iterations"), exist_ok=True)
        
        # Create a minimal example cell_metrics.csv
        example_df = pd.DataFrame({
            'cell_id': [f'{cell_type}_1', f'{cell_type}_2', f'{cell_type}_3'],
            'cell_type': [cell_type, cell_type, cell_type],
            't_stat': [2.5, 1.8, 3.1],
            'p_value': [0.01, 0.05, 0.001],
            'p_value_adj': [0.03, 0.12, 0.005],
            'disease_z_score': [0.8, 0.5, 0.9],
            'trait_age_score': [0.3, 0.3, 0.3]
        })
        
        example_df.to_csv(os.path.join(example_dir, cell_type, "cell_metrics.csv"), index=False)
    
    print(f"\nCreated example directory structure in '{example_dir}'")
    print(f"Run the check on this example with: python check_cell_metrics.py --output_dir {example_dir}")
    print(f"This example includes a custom cell type 'Macrophages' to demonstrate flexibility.")

if __name__ == "__main__":
    main()