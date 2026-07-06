# scDCF Method Summary

This note summarizes the current public implementation of `scDCF`.

## Inputs

`scDCF` expects:

- an AnnData `.h5ad` object with cell metadata in `adata.obs`
- a prioritized gene list with at least a `gene_name` column
- optionally a `zstat` column used as non-negative gene weights

If `zstat` is missing, the package assigns uniform weights.

## Control Genes

For each analyzed cell type, `scDCF` constructs disease-group and healthy-group
control-gene dictionaries by matching each prioritized gene to genes with
similar mean and variance in the corresponding cell population.

## Reference Matching

For each target cell, `scDCF` builds a pool of up to 1,000 healthy reference
cells matched by RNA count.

- Disease-target analysis:
  disease cells are evaluated against healthy reference cells.
- Healthy-target analysis:
  healthy cells are evaluated against healthy reference cells, with the target
  cell excluded from its own reference pool.

At each Monte Carlo iteration, up to 100 cells are sampled from the matched
healthy reference pool for each target cell.

## Cell-Level Test

Within each iteration, `scDCF`:

1. computes absolute target-versus-reference deviations for prioritized genes
2. computes matched control-gene deviations
3. forms weighted gene-level difference-of-differences values
4. applies a one-sided one-sample t-test against zero for each target cell

The one-sided alternative tests whether prioritized-gene deviations exceed the
matched-control baseline on average.

## Aggregation and FDR

Across Monte Carlo iterations, per-cell p-values are combined with Fisher's
method. The resulting per-cell p-values are then corrected with
Benjamini-Hochberg across all tested cells within each target-group analysis,
yielding `q_cell`.

## Cell-Type Enrichment

After defining scDCF-significant cells from `q_cell`, the package compares
significant versus non-significant cells between disease and healthy groups
within each analyzed cell type using Fisher's exact test. These p-values are
then Benjamini-Hochberg corrected across tested cell types to produce `q_type`.

## Default Outputs

By default, the package writes only user-facing outputs:

- one `*_final_summary.csv` per analyzed cell type
- one `celltype_enrichment_summary.csv` per run

Intermediate Monte Carlo and combined tables can be exported with
`--export_intermediate`.
