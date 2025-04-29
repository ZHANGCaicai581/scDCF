#!/usr/bin/env python
"""
Test for reorganizing existing output with control genes support.
"""
import os
import sys
import shutil

def test_with_existing():
    """Reorganize existing output with control genes in subfolder structure"""
    print("Reorganizing existing output data...")
    
    # Check if test_output directory exists
    if not os.path.exists("test_output"):
        print("Error: test_output directory not found. Run test_organization.py first.")
        return False
    
    # Use the function directly from post_trait_test
    try:
        from post_trait_test import organize_final_output
        # Reorganize the existing output
        organize_final_output(source_dir="test_output", dest_dir="test_output_reorganized")
    except ImportError:
        # Fallback to using the package
        try:
            from scDCF import organize_output
            organize_output(source_dir="test_output", dest_dir="test_output_reorganized")
        except ImportError:
            print("ERROR: Neither post_trait_test.organize_final_output nor scDCF.organize_output could be imported.")
            return False
    
    # Check the new structure
    print("Running check_cell_metrics.py on reorganized output...")
    result = os.system("python check_cell_metrics.py --output_dir test_output_reorganized")
    
    return result == 0

if __name__ == "__main__":
    # Clean up any previous test output
    if os.path.exists("test_output_reorganized"):
        shutil.rmtree("test_output_reorganized")
    
    success = test_with_existing()
    print("Test completed successfully" if success else "Test failed")
    sys.exit(0 if success else 1)

#!/usr/bin/env python
import os
import json
import shutil
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def test_package_integration():
    """Test that the package correctly handles the reorganized structure"""
    # Set up the test directory
    test_dir = "integration_test"
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir)
    
    # Create test directory structure
    for cell_type in ['T_cell', 'B_cell']:
        # Create both types of structure - old and new
        # New structure (with subfolder)
        cell_dir = os.path.join(test_dir, cell_type)
        support_dir = os.path.join(cell_dir, "supporting_data")
        control_dir = os.path.join(support_dir, "control_genes")
        os.makedirs(control_dir, exist_ok=True)
        
        # Create control genes file in new location
        control_genes = {
            "disease": ["GENE1", "GENE2"],
            "healthy": ["GENE3", "GENE4"]
        }
        with open(os.path.join(control_dir, f"{cell_type}_control_genes.json"), 'w') as f:
            json.dump(control_genes, f, indent=2)
            
        # Create test data files
        with open(os.path.join(cell_dir, "cell_metrics.csv"), 'w') as f:
            f.write("cell_id,cell_type,t_stat,p_value\n")
            f.write(f"cell1,{cell_type},2.5,0.01\n")
    
    # Test loading control genes from new location
    try:
        from scDCF.utils import load_control_genes
        for cell_type in ['T_cell', 'B_cell']:
            # Try to load from the new location
            control_genes_path = os.path.join(test_dir, cell_type, "supporting_data", "control_genes", f"{cell_type}_control_genes.json")
            controls = load_control_genes(control_genes_path)
            logging.info(f"Successfully loaded control genes for {cell_type} from new location")
            logging.info(f"Control genes: {controls}")
    except ImportError:
        logging.warning("Could not import scDCF.utils - skipping this test")
    except Exception as e:
        logging.error(f"Error loading control genes: {e}")
    
    # Run check_cell_metrics.py on the test directory
    logging.info("Running check_cell_metrics.py on test directory")
    os.system(f"python check_cell_metrics.py --output_dir {test_dir}")
    
    logging.info("Integration test completed!")

if __name__ == "__main__":
    test_package_integration()

# Example: Organizing output with the enhanced structure
from scDCF import organize_output

# Original analysis output directory
source_dir = "my_analysis_results"

# Create organized output with enhanced structure
organized_dir = organize_output(source_dir, "organized_results")

# The control genes files will be in:
# - organized_results/[CELL_TYPE]/supporting_data/control_genes/[CELL_TYPE]_control_genes.json
# - organized_results/[CELL_TYPE]_control_genes.json (for backward compatibility)

## Key Files and Directories

# Basic usage for organizing output
from scDCF import organize_output

# Organize an existing analysis output directory
organize_output(source_dir="my_analysis_results", dest_dir="organized_results")

# Other key functions from the package
from scDCF import read_gene_symbols, generate_control_genes
from scDCF import monte_carlo_comparison, load_monte_carlo_results
from scDCF import get_trait_association_scores

