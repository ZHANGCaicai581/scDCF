# Output Structure

## Default Outputs

For a standard full run (`python -m scDCF ...` with the default `--step all`),
the package writes only the user-facing results:

- `<output_dir>/<cell_type>/<cell_type>_final_summary.csv`
- `<output_dir>/celltype_enrichment_summary.csv`

## Per-Cell Final Summary

Each `*_final_summary.csv` contains the combined cell-level results for both the
disease and healthy target analyses, including:

- `cell_id`
- `group`
- `combined_p_value_fisher`
- `combined_p_value_stouffer`
- `q_cell`
- `scdcf_significant`

If metadata export is enabled, columns from `adata.obs` are also merged into the
final summary.

## Cell-Type Enrichment Summary

`celltype_enrichment_summary.csv` reports cell-type-level enrichment statistics,
including:

- `cell_type`
- disease and healthy significant-cell counts
- disease and healthy total-cell counts
- `odds_ratio`
- `p_value`
- `q_type`
- `disease_enriched`

## Optional Intermediate Outputs

If `--export_intermediate` is supplied, the package additionally writes:

- `<cell_type>_disease_monte_carlo_results.csv`
- `<cell_type>_healthy_monte_carlo_results.csv`
- `<cell_type>_disease_combined.csv`
- `<cell_type>_healthy_combined.csv`

These files are mainly intended for debugging, auditing, and method development.
