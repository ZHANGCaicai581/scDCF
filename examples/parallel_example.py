#!/usr/bin/env python3
"""
Example: Using scDCF with Parallel Processing (v0.1.12+)

This example demonstrates the new parallel processing features for 4-8x speedup.
"""

import scDCF
import scanpy as sc

# ============================================================================
# Example 1: Default Serial Mode (Out-of-the-box behavior)
# ============================================================================
print("=" * 80)
print("Example 1: Default Serial Mode (Easiest)")
print("=" * 80)

# Load your data
adata = sc.read_h5ad("path/to/your/data.h5ad")
genes_df = scDCF.read_gene_symbols("path/to/genes.txt")

# Generate control genes
disease_ctrl, healthy_ctrl = scDCF.generate_control_genes(
    adata=adata,
    significant_genes_df=genes_df,
    cell_type="T_cell",
    cell_type_column="celltype_major",
    n_control_genes=10
)

results = scDCF.auto_monte_carlo(
    adata=adata,
    cell_type="T_cell",
    cell_type_column="celltype_major",
    significant_genes_df=genes_df,
    disease_control_genes=disease_ctrl,
    healthy_control_genes=healthy_ctrl,
    output_dir="results/",
    iterations=10  # Serial by default
)

print(f"✅ Analysis complete! Analyzed {len(results)} cells")

# ============================================================================
# Example 2: Explicit Parallel (Advanced Users)
# ============================================================================
print("\n" + "=" * 80)
print("Example 2: Enable Parallel Processing")
print("=" * 80)

# Run with auto_monte_carlo in parallel mode (auto-select workers)
parallel_results = scDCF.auto_monte_carlo(
    adata=adata,
    cell_type="T_cell",
    cell_type_column="celltype_major",
    significant_genes_df=genes_df,
    disease_control_genes=disease_ctrl,
    healthy_control_genes=healthy_ctrl,
    output_dir="results/",
    iterations=100,
    use_parallel=True,      # enable parallel processing
    target_group="disease"
)

# Run with a specific worker count
healthy_results = scDCF.parallel_monte_carlo_comparison(
    adata=adata,
    cell_type="T_cell",
    cell_type_column="celltype_major",
    significant_genes_df=genes_df,
    disease_control_genes=disease_ctrl,
    healthy_control_genes=healthy_ctrl,
    output_dir="results/",
    iterations=100,
    target_group="healthy",
    n_workers=4,            # use 4 cores explicitly
    show_progress=True
)

print(f"✅ Disease: {len(parallel_results)} cells analyzed")
print(f"✅ Healthy: {len(healthy_results)} cells analyzed")

# ============================================================================
# Example 3: Original Serial (Backward Compatibility)
# ============================================================================
print("\n" + "=" * 80)
print("Example 3: Original Serial Method (Still Supported)")
print("=" * 80)

# Original method still works (useful for single iterations or debugging)
serial_results = scDCF.monte_carlo_comparison(
    adata=adata,
    cell_type="T_cell",
    cell_type_column="celltype_major",
    significant_genes_df=genes_df,
    disease_control_genes=disease_ctrl,
    healthy_control_genes=healthy_ctrl,
    output_dir="results/",
    iterations=10  # Works fine with serial for small iteration counts
)

print(f"✅ Serial analysis complete: {len(serial_results)} cells")

# ============================================================================
# Example 4: Performance Comparison
# ============================================================================
print("\n" + "=" * 80)
print("Example 4: Compare Serial vs Parallel Performance")
print("=" * 80)

import time

# Time serial execution
start = time.time()
serial_results = scDCF.monte_carlo_comparison(
    adata=adata,
    cell_type="T_cell",
    cell_type_column="celltype_major",
    significant_genes_df=genes_df,
    disease_control_genes=disease_ctrl,
    healthy_control_genes=healthy_ctrl,
    output_dir="results/serial/",
    iterations=3
)
serial_time = time.time() - start

# Time parallel execution
start = time.time()
parallel_results = scDCF.parallel_monte_carlo_comparison(
    adata=adata,
    cell_type="T_cell",
    cell_type_column="celltype_major",
    significant_genes_df=genes_df,
    disease_control_genes=disease_ctrl,
    healthy_control_genes=healthy_ctrl,
    output_dir="results/parallel/",
    iterations=3,
    n_workers=None  # Auto-detect (≤ min(CPUs-1, 8))
)
parallel_time = time.time() - start

print(f"\n📊 Performance Comparison:")
print(f"   Serial:   {serial_time/60:.2f} minutes")
print(f"   Parallel: {parallel_time/60:.2f} minutes")
print(f"   Speedup:  {serial_time/parallel_time:.1f}x")

# ============================================================================
# Example 5: Complete Workflow with Parallel Processing
# ============================================================================
print("\n" + "=" * 80)
print("Example 5: Complete Workflow")
print("=" * 80)

# 1. Load data
adata = sc.read_h5ad("data.h5ad")
genes_df = scDCF.read_gene_symbols("genes.txt")

# 2. Generate control genes
disease_ctrl, healthy_ctrl = scDCF.generate_control_genes(
    adata, genes_df, "T_cell", "celltype_major"
)

# 3. Run parallel analysis (both groups)
print("Running disease group (parallel)...")
disease_results = scDCF.auto_monte_carlo(
    adata, "T_cell", "celltype_major", genes_df,
    disease_ctrl, healthy_ctrl, "results/",
    iterations=100, target_group="disease", use_parallel=True
)

print("Running healthy group (parallel)...")
healthy_results = scDCF.auto_monte_carlo(
    adata, "T_cell", "celltype_major", genes_df,
    disease_ctrl, healthy_ctrl, "results/",
    iterations=100, target_group="healthy", use_parallel=True
)

# 4. Post-analysis
disease_combined = scDCF.combine_p_values_across_iterations(
    disease_results, "results/", "T_cell", "disease"
)
healthy_combined = scDCF.combine_p_values_across_iterations(
    healthy_results, "results/", "T_cell", "healthy"
)

# 5. Final per-cell summary (writes T_cell_final_summary.csv)
final_summary = scDCF.export_final_celltype_summary(
    cell_type="T_cell",
    disease_combined=disease_combined,
    healthy_combined=healthy_combined,
    output_dir="results/"
)

print("\n✅ Complete workflow finished!")
print(f"Final summary saved to results/T_cell/T_cell_final_summary.csv with {len(final_summary)} rows")

# ============================================================================
# Tips for Users
# ============================================================================
print("\n" + "=" * 80)
print("💡 Tips for Best Performance")
print("=" * 80)

print("""
1. Use auto_monte_carlo() with `use_parallel=True` for multi-core speedups.
   - Without the flag it stays single-core (predictable resource usage).

2. For 100+ iterations, enabling parallel mode is strongly recommended.
   - 5-10x faster on typical machines
   - Saves hours of computation time

3. Check your system:
   import multiprocessing as mp
   print(f"CPU cores: {mp.cpu_count()}")

4. Adjust workers for your hardware:
   n_workers=None  # Auto (capped at ≤ min(CPUs-1, 8))
   n_workers=4     # Specific number
   n_workers=mp.cpu_count()  # Use all cores (consider reserving 1-2 cores)

5. Trade-offs:
   - More workers = faster (up to # of cores)
   - More workers = more memory (500MB per worker)
   - Diminishing returns beyond # of cores

6. For servers/clusters:
   - Use n_workers=24 or more
   - Can achieve 15-20x speedup!
""")

print("\n✅ Happy analyzing! 🚀")

