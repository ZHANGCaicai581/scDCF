# Matched directionality simulations

Reproduces the analysis workflow for **Additional file 2, Tables S7–S8**:
40%/60% down-regulation with exact reversal of intended effects, and the
additional 10%/20% settings with nested down-regulated gene assignments.

**Simulation is implemented in R with Splatter.** One generator handles all
four settings and reuses unchanged helper functions from the author's original
R simulator. Python calls the R scripts, runs scDCF/scDRS, aligns cell IDs and
calculates metrics. There is no second simulation implementation in Python.

## Files

| File | Purpose |
|---|---|
| `generate_scenarios.R` | Generate paired counts, labels and perturbation plans |
| `score_reference_methods.R` | Seurat, reference AUCell, full-rank enrichment and VAM scores |
| `run_benchmark.py` | Run stages and calculate AUROC/AP with fixed score direction |
| `audit_matched.py` | Check the 40%/60% count differences and effect reversal |
| `plot_directionality.R` | Optional score-distribution and perturbation-count plots |
| `test_audit_matched.py` | Count-support and table-definition checks |
| `reference_metrics.csv` | Saved seed-46 values underlying S7–S8 |
| `provenance.json` | Original R source identity, software versions and scDCF source hashes |
| `VALIDATION.txt` | Checks performed when preparing this code package |

For the exact source used in the saved analyses, download and extract the
versioned Zenodo archive below, then run commands from its repository root.
The archive includes the scDCF source snapshot; it need not be copied into
this folder. GitHub `main` can contain later changes to the scDCF package.
License: the repository's MIT license.

## Archive identity and publication

Analysis archive version: **directionality-2026.11.03**.
This identifies the S7–S8 code snapshot, separately from scDCF release v1.0.0
and the included package metadata version 0.1.15.

Rebuild the ZIP from the repository root with
`python3 analysis/build_directionality_archive.py`. The builder checks the
recorded scDCF source hashes and includes only the listed code, metadata and
small reference files. The ZIP and its `.sha256` checksum are in `analysis/`.
`archive_metadata.json` contains the prepared Zenodo title, creator, scope,
version and MIT license; it is also included as `.zenodo.json` in the ZIP.

Published archive DOI: **[10.5281/zenodo.23918448](https://doi.org/10.5281/zenodo.23918448)**.
The ZIP is the frozen analysis snapshot. Its README records the DOI reservation
at packaging time; publication was completed on 2026-11-03. This GitHub README
updates the publication status without changing the archived analysis code.
The existing DOI `10.5281/zenodo.18409083` identifies the earlier scDCF v1.0.0
release and must not be cited as containing these additional analyses.

Suggested citation: Zhang, Caicai (2026). *scDCF: matched directionality
analysis code for Supplementary Tables S7–S8*, version
directionality-2026.11.03. Zenodo. https://doi.org/10.5281/zenodo.23918448.

The scDCF v1.0.0 Zenodo record's license metadata was corrected to MIT on
2026-11-03. This analysis archive also uses MIT.

## Environment

The saved simulations used **R 4.2.2**, Splatter 1.22.1,
SingleCellExperiment 1.20.1, Matrix 1.6-5 and jsonlite 2.0.0. Scoring used
Seurat 4.4.0, SeuratObject 5.2.0, AUCell 1.20.2, VAM 1.1.0 and matrixStats
1.5.0. The R scorer retains the process-local `colRanks(useNames=NA)`
compatibility adjustment needed by that AUCell/matrixStats combination.
Install these R packages in a compatible R/Bioconductor environment.
Optional plotting also uses ggplot2 and patchwork.

Use Python 3.12 for the scDCF/evaluation environment:

```sh
python3 -m venv .venv-analysis
.venv-analysis/bin/python -m pip install -r analysis/simulation/requirements.txt
.venv-analysis/bin/python -m pip install --no-deps -e .
```

The saved scDRS run used a separate environment:

```sh
python3 -m venv .venv-scdrs
.venv-scdrs/bin/python -m pip install -r analysis/simulation/requirements-scdrs.txt
```

`Rscript` must be on PATH. Pass the scDRS executable explicitly when using a
separate environment. The requirements record the packages used here; they
are not a lock of every transitive dependency. Exact numerical reproduction
should be checked against the saved metrics and recorded R versions.

## Generate data directly in R

```sh
Rscript analysis/simulation/generate_scenarios.R results/directionality_10_20 46 10,20 nested
Rscript analysis/simulation/generate_scenarios.R results/directionality_40_60 46 40,60 sign_reverse
```

Each dataset contains 12,000 genes and 2,400 cells, with 800 cells per type,
50% disease labels within each type, and 40 positive disease CellTypeC cells.
The input list contains 500 genes; 350 are selected per positive cell.
The down-regulated counts are 35, 70, 140 and 210 for the four fractions.
Magnitudes are sampled in log1p-count space from `N(log(1.7), 0.25^2)`, clipped
to `[log(1.1), log(4)]`. Perturbed counts are clipped and rounded. Intended
sign reversal therefore need not yield perfectly mirrored observed counts.
Splatter dropout is disabled for these dense matched settings.

The two paired runs share the same seed and baseline-generation procedure.
Use separate output directories for different parameter configurations.
Completed generation is reused only after configuration and checksum checks.

## Score and evaluate

```sh
.venv-analysis/bin/python analysis/simulation/run_benchmark.py \
  --output-root results/directionality_10_20 --fractions 10,20 --pair-mode nested \
  --scdrs-executable .venv-scdrs/bin/scdrs --check-reference

.venv-analysis/bin/python analysis/simulation/run_benchmark.py \
  --output-root results/directionality_40_60 --fractions 40,60 --pair-mode sign_reverse \
  --scdrs-executable .venv-scdrs/bin/scdrs --check-reference
```

Alternatively add `--generate` to either command to call the same R generator
first. Defaults are seed 46, 10 scDCF iterations and 1,000 scDRS control sets.
Use `--stages summarize --check-reference` to rebuild metrics from saved scores,
or `--generate --stages validate` to generate and check counts without scoring.
Methods can take substantially longer than generation; no full-run time estimate
is asserted. A failed or incomplete requested scoring stage exits unsuccessfully.

Scores are aligned by cell ID. Higher scores always indicate predicted
positives; labels never determine score orientation. AP is calculated with
`average_precision_score`; integral AUPRC from the older benchmark is a
different summary. A single seed is one matched realization, not a replicate
estimate of performance variability.

Outputs include:

- `table_metrics.csv`: the five statistics reported in S7–S8, including
  **Full-rank enrichment** (`LegacyFullRank`).
- `primary_results.csv` and `results_table.tex`: reference-method results,
  including the actual **AUCell package**. These are distinct from S7–S8.
- `metrics_fixed_direction.csv`: all metrics and evaluation compartments.
- `paired_design_audit.json`, `cell_alignment_audit.csv`, per-scenario manifests
  and logs: design and execution checks.

The table exporter keeps full-rank enrichment separate from default AUCell.
`--check-reference` checks all reported method/scenario values at four decimal
places against `reference_metrics.csv`; it does not alter results to match.

## Validate and optionally plot

```sh
.venv-analysis/bin/python -m unittest discover -s analysis/simulation -p 'test_*.py' -v
.venv-analysis/bin/python analysis/simulation/audit_matched.py --output-root results/directionality_40_60
Rscript analysis/simulation/plot_directionality.R results/directionality_40_60 46
```

Plots are optional and are written as PDF plus PNG, with source data.
The count audit verifies that differences stay within the planned gene–cell
pairs and agree with intended directions, while permitting ties after rounding.

## Scope and provenance

This is the compact **S7–S8 analysis package**, not the historical Groups A–H
grid or the real-data analyses. The original R file's exploratory run loop is
excluded: one old `10down` name actually used `mixed_frac_up=0.70` (30% down).
The current R generator sets and checks all four down fractions explicitly.
The original source hash and the four reused helpers are recorded in
`provenance.json`; no file under a user's Downloads directory is required.
Large matrices, cached scores, manuscript drafts and duplicate figure exports
are regenerated outputs and are not included in the code archive.
