#!/usr/bin/env python
"""
Simple test to verify control genes organization with existing data.
"""
import os
import glob
from post_trait_test import organize_final_output

# Check if test_output directory exists
if not os.path.exists("test_output"):
    print("Error: test_output directory not found. Run the main tests first.")
    exit(1)

# Create a simplified reorganization
organize_final_output(source_dir="test_output", dest_dir="simple_test_output")

# Verify the control genes were copied to both locations
all_ok = True
for cell_type in ['T_cell', 'B_cell', 'NK_cell']:
    # Check root level file
    root_file = os.path.join("simple_test_output", f"{cell_type}_control_genes.json")
    if os.path.exists(root_file):
        print(f"✅ {cell_type} control genes found at root level")
    else:
        print(f"❌ {cell_type} control genes missing at root level")
        all_ok = False
    
    # Check subfolder file
    subfolder_file = os.path.join(
        "simple_test_output", cell_type, "supporting_data", "control_genes", 
        f"{cell_type}_control_genes.json"
    )
    if os.path.exists(subfolder_file):
        print(f"✅ {cell_type} control genes found in subfolder")
    else:
        print(f"❌ {cell_type} control genes missing in subfolder")
        all_ok = False

# Check for control_genes directories
control_dirs = glob.glob("simple_test_output/*/supporting_data/control_genes")
print(f"Found {len(control_dirs)} control_genes directories")

if all_ok:
    print("\n✅ All control genes are correctly organized!")
    exit(0)
else:
    print("\n❌ Some control genes are missing!")
    exit(1) 