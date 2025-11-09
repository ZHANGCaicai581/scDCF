#!/usr/bin/env python3
"""
Example: Working with Cell IDs in scDCF

This example shows how to:
1. Understand what cell IDs represent
2. Add metadata to results
3. Trace results back to original cells
4. Create enhanced result files
"""

import scDCF
import scanpy as sc
import pandas as pd

print("=" * 80)
print("🔖 Cell ID Workflow Example")
print("=" * 80)

# ============================================================================
# Understanding Cell IDs
# ============================================================================
print("\n📚 Understanding Cell IDs in scDCF")
print("-" * 80)

print("""
Cell IDs in scDCF results are the ORIGINAL cell identifiers from your AnnData.

If your adata.obs_names are: ['AAACCTGA', 'AAACGGG', 'TTTCCCAA', ...]
Then results will have: cell_id = 'AAACCTGA', 'AAACGGG', 'TTTCCCAA', ...

If your adata.obs_names are: [0, 1, 2, 3, ...]
Then results will have: cell_id = 0, 1, 2, 3, ...

They're always consistent! ✅
""")

# ============================================================================
# Example 1: Basic Workflow
# ============================================================================
print("\n" + "=" * 80)
print("Example 1: Basic Analysis with Metadata")
print("=" * 80)

# Load data
adata = sc.read_h5ad("path/to/data.h5ad")
genes_df = scDCF.read_gene_symbols("genes.txt")

print(f"\n📊 Your data:")
print(f"   Cell names: {list(adata.obs_names[:5])}...")
print(f"   Cell name type: {type(adata.obs_names[0])}")

# Generate control genes
disease_ctrl, healthy_ctrl = scDCF.generate_control_genes(
    adata, genes_df, "T_cell", "celltype"
)

# Run analysis
results = scDCF.monte_carlo_comparison(
    adata, "T_cell", "celltype", genes_df,
    disease_ctrl, healthy_ctrl, "results/",
    iterations=10
)

print(f"\n✅ Analysis complete: {len(results)} result rows")
print(f"   Unique cells: {results['cell_id'].nunique()}")
print(f"   Cell IDs: {list(results['cell_id'].unique()[:5])}...")

# ============================================================================
# Example 2: Add Metadata to Results (NEW!)
# ============================================================================
print("\n" + "=" * 80)
print("Example 2: Enhance Results with Cell Metadata")
print("=" * 80)

# Add all metadata
enhanced_results = scDCF.add_cell_metadata(results, adata)

print(f"\n✅ Metadata added!")
print(f"   Original columns: {len(results.columns)}")
print(f"   Enhanced columns: {len(enhanced_results.columns)}")
print(f"   New metadata: {len(enhanced_results.columns) - len(results.columns)} fields")

# Show example
print(f"\nExample enhanced results:")
cols_to_show = ['cell_id', 'p_value', 'significant']
metadata_cols = [c for c in enhanced_results.columns if c not in results.columns][:3]
print(enhanced_results[cols_to_show + metadata_cols].head())

# Save enhanced results
enhanced_results.to_csv("results_with_metadata.csv", index=False)
print(f"\n💾 Saved to: results_with_metadata.csv")

# ============================================================================
# Example 3: Add Specific Metadata Only
# ============================================================================
print("\n" + "=" * 80)
print("Example 3: Add Only Specific Metadata Columns")
print("=" * 80)

# Add only columns you need
enhanced_selected = scDCF.add_cell_metadata(
    results, 
    adata,
    metadata_columns=['celltype', 'sample', 'batch', 'patient_id']
)

print(f"✅ Added only: celltype, sample, batch, patient_id")
print(enhanced_selected[['cell_id', 'p_value', 'celltype', 'sample']].head())

# ============================================================================
# Example 4: Create Cell ID Reference File
# ============================================================================
print("\n" + "=" * 80)
print("Example 4: Create Cell ID Mapping File")
print("=" * 80)

# Create mapping for all cells
mapping = scDCF.create_cell_id_mapping(adata)
mapping.to_csv("cell_id_mapping.csv", index=False)

print(f"✅ Cell ID mapping created: {len(mapping)} cells")
print(f"   Columns: {list(mapping.columns)}")
print(f"\nNow you can always look up cell information:")
print(mapping.head())

# ============================================================================
# Example 5: Validate Cell IDs
# ============================================================================
print("\n" + "=" * 80)
print("Example 5: Validate Cell IDs")
print("=" * 80)

# Check that results match your data
validation = scDCF.validate_cell_ids(results, adata)

print(f"Validation Report:")
print(f"   Total cell IDs in results: {validation['total_result_ids']}")
print(f"   Matched with AnnData: {validation['matched']}")
print(f"   Match percentage: {validation['match_percentage']:.1f}%")
print(f"   All matched: {validation['all_matched']}")

if validation['all_matched']:
    print(f"\n✅ Perfect! All cell IDs trace back to your AnnData")
else:
    print(f"\n⚠️  Some cell IDs not found:")
    print(f"   Missing: {validation['missing_ids']}")

# ============================================================================
# Example 6: Filter and Analyze by Metadata
# ============================================================================
print("\n" + "=" * 80)
print("Example 6: Filter Results by Cell Metadata")
print("=" * 80)

# Get enhanced results
enhanced = scDCF.add_cell_metadata(results, adata)

# Filter significant cells from specific sample
significant_sample1 = enhanced[
    (enhanced['p_value'] < 0.05) &
    (enhanced['sample'] == 'Sample1')
]

print(f"Significant cells in Sample1: {len(significant_sample1)}")

# Group by batch
batch_summary = enhanced.groupby('batch').agg({
    'cell_id': 'count',
    'p_value': 'mean',
    'significant': 'sum'
})
batch_summary.columns = ['n_cells', 'mean_p_value', 'n_significant']

print(f"\nSummary by batch:")
print(batch_summary)

# ============================================================================
# Example 7: Get Info for Significant Cells
# ============================================================================
print("\n" + "=" * 80)
print("Example 7: Look Up Significant Cell Information")
print("=" * 80)

# Get cells with p < 0.01
highly_significant = results[results['p_value'] < 0.01]
sig_cell_ids = highly_significant['cell_id'].tolist()

print(f"Found {len(sig_cell_ids)} highly significant cells")

# Get their full information
sig_cell_info = scDCF.get_cell_info_by_id(adata, sig_cell_ids)

print(f"\nInformation for highly significant cells:")
print(sig_cell_info[['celltype', 'sample', 'batch']].head())

# Save for downstream analysis
sig_cell_info.to_csv("significant_cells_info.csv")
print(f"\n💾 Saved to: significant_cells_info.csv")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 80)
print("✅ Complete Cell ID Workflow")
print("=" * 80)

print("""
Key takeaways:

1. Cell IDs in results = adata.obs_names (always consistent!)

2. Use add_cell_metadata() to enhance results:
   enhanced = scDCF.add_cell_metadata(results, adata)

3. Use create_cell_id_mapping() for reference:
   mapping = scDCF.create_cell_id_mapping(adata)

4. Use validate_cell_ids() to verify:
   report = scDCF.validate_cell_ids(results, adata)

5. Use get_cell_info_by_id() to look up specific cells:
   info = scDCF.get_cell_info_by_id(adata, cell_ids)

All functions work seamlessly with your data! ✅
""")

print("\n🎉 You can now easily trace any result back to your original cells!")

