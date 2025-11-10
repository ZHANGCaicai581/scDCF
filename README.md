# scDCF

*A Framework for Detecting Disease-associated Cells in Single-cell RNA-seq  
Leveraging Healthy Reference Panels and GWAS Findings*

[![PyPI version](https://img.shields.io/pypi/v/scDCF.svg)](https://pypi.org/project/scDCF/)
[![Python versions](https://img.shields.io/pypi/pyversions/scDCF.svg)](https://pypi.org/project/scDCF/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

![scDCF workflow](scDCF/docs/scDCF_workflow.png)

---

## Table of Contents
1. [Introduction](#1-introduction)  
2. [Key Features](#2-key-features)  
3. [Installation](#3-installation)  
4. [Quick Start](#4-quick-start)  
5. [Datasets and Methods](#5-datasets-and-methods)  
6. [Data Sources](#6-data-sources)  
7. [Contact](#7-contact)  
8. [License](#8-license)

## 1. Introduction
Genome-wide association studies (GWAS) have uncovered thousands of risk loci, but the cell types through which these variants act remain unclear. **scDCF (single-cell Disease Cell Finder)** integrates GWAS-derived gene sets with single-cell RNA-seq data, using a library-size-matched healthy reference panel, control-gene matching, and Monte-Carlo statistics to pinpoint cells whose expression profiles are genuinely perturbed by inherited risk.

## 2. Key Features
| Capability | Summary |
|------------|---------|
| **MAGMA/TWAS integration** | Transforms GWAS SNP statistics to gene-level Z-scores; intersects with expressed genes (G* = G ∩ E). |
| **Library-size-matched reference pools** | Constructs 1,000-cell healthy reference pools per target cell; samples 100 cells per Monte Carlo iteration. |
| **Cell-type-specific control matching** | Assigns 10 control genes per prioritized gene, matched on mean/variance within cell type and disease status. |
| **Difference-of-differences framework** | Isolates disease signal via δ_target - δ_control, weighted by MAGMA Z-scores and averaged across genes. |
| **Fisher meta-analysis** | Aggregates iteration-level p-values using Fisher's method; applies Benjamini-Hochberg FDR correction. |
| **Cell-type enrichment testing** | Fisher's exact test on 2×2 contingency tables of significant cells vs. disease status per cell type. |
| **Scalable implementation** | Python ≥ 3.9; optimized for sparse matrices; supports custom gene lists and flexible annotations. |

## 3. Installation
```bash
# Install from PyPI (recommended)
pip install scDCF

# Or install latest from GitHub
pip install git+https://github.com/ZHANGCaicai581/scDCF.git

# Verify installation
python -c "import scDCF; print(f'scDCF version: {scDCF.__version__}')"
```

**Requirements**: Python ≥ 3.9

## 4. Quick Start

### Python API

```python
import scDCF
import scanpy as sc

# 1. Load your preprocessed scRNA-seq data
adata = sc.read_h5ad("path/to/data.h5ad")

# 2. Load GWAS/MAGMA prioritized genes
significant_genes_df = scDCF.read_gene_symbols("genes.txt")  # One gene per line

# 3. Generate control genes (10 per significant gene, matched on expression)
disease_ctrl, healthy_ctrl = scDCF.generate_control_genes(
    adata=adata,                          # Your AnnData object
    significant_genes_df=significant_genes_df,  # GWAS genes
    cell_type="T_cell",                   # Cell type to analyze
    cell_type_column="celltype_major"     # Column with cell type labels
)

# 4. Run Monte Carlo analysis
results = scDCF.monte_carlo_comparison(
    adata=adata,                          # Your AnnData object
    cell_type="T_cell",                   # Cell type to analyze
    cell_type_column="celltype_major",    # Column with cell type labels
    significant_genes_df=significant_genes_df,  # GWAS genes
    disease_control_genes=disease_ctrl,   # From step 3
    healthy_control_genes=healthy_ctrl,   # From step 3
    output_dir="results/",                # Where to save results
    iterations=10                         # Number of iterations (default: 10)
)

# Optional: Add cell metadata to results (get celltype, sample, batch, etc.)
enhanced = scDCF.add_cell_metadata(results, adata)
enhanced.to_csv("results_with_metadata.csv")
```

**For faster analysis** (recommended for 100+ iterations):
```python
# Enable parallel processing (4-8x speedup on multi-core systems)
results = scDCF.auto_monte_carlo(
    adata=adata,
    cell_type="T_cell",
    cell_type_column="celltype_major",
    significant_genes_df=significant_genes_df,
    disease_control_genes=disease_ctrl,
    healthy_control_genes=healthy_ctrl,
    output_dir="results/",
    iterations=100,
    use_parallel=True  # Set True to allow multi-core execution
)
```

For detailed examples, see the [examples directory](examples/). Also see the methods summary in [scDCF/docs/methods.md](scDCF/docs/methods.md).

### Command Line Usage

**Basic usage** (replace with your file paths):
```bash
python -m scDCF \
  --h5ad_file YOUR_DATA.h5ad \
  --gene_list_file YOUR_GENES.txt \
  --output_dir results/ \
  --celltype_column YOUR_CELLTYPE_COLUMN \
  --disease_marker YOUR_DISEASE_COLUMN \
  --rna_count_column YOUR_RNA_COUNT_COLUMN
```

**Example with real data** (uses default 10 iterations):
```bash
python -m scDCF \
  --h5ad_file pbmc_data.h5ad \
  --gene_list_file sle_genes.txt \
  --output_dir results/ \
  --celltype_column celltype_major \
  --disease_marker disease_status \
  --disease_value "SLE" \
  --healthy_value "Control" \
  --rna_count_column nCount_RNA
```

**Quick test** (bundled synthetic data, completes in ~5 min):
```bash
python -m scDCF \
  --h5ad_file data/test/sim_adata.h5ad \
  --gene_list_file data/test/genes.txt \
  --output_dir test_results/ \
  --celltype_column cell_type \
  --disease_marker disease_numeric \
  --rna_count_column nCount_RNA \
  --iterations 2
```

> **Note:** scDCF now runs Monte Carlo iterations serially by default (single core). Enable parallel mode with `--parallel` (auto-selects a capped worker pool, `≤ min(total CPUs - 1, 8)`) or specify `--parallel_workers N`. Use `--serial` to force single-core behavior explicitly.


### Methods at a glance

For a concise overview, see the detailed methodology in `scDCF/docs/methods.md`. The README intentionally stays brief to focus on usage.

### Command-line parameters

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `--csv_file` | path | None | Path to CSV/TSV file containing prioritized genes (must include gene name and preferably Z-stat). |
| `--gene_list_file` | path | None | Path to a plain-text file with one gene per line. |
| `--h5ad_file` | path | required | Path to AnnData `.h5ad` file. |
| `--output_dir` | path | required | Output directory for results. |
| `--celltype_column` | str | `celltype_major` | Column in `adata.obs` with cell type labels. |
| `--cell_types` | list[str] | None | Subset of cell types to analyze; defaults to all in `celltype_column`. |
| `--disease_marker` | str | `disease_numeric` | Column in `adata.obs` indicating disease status. |
| `--disease_value` | (str|int|float) | `1` | Value indicating disease cells. |
| `--healthy_value` | (str|int|float) | `0` | Value indicating healthy cells. |
| `--rna_count_column` | str | `nCount_RNA` | Column in `adata.obs` for library size / RNA counts. |
| `--iterations` | int | `10` | Number of Monte Carlo iterations. |
| `--show_progress` | flag | `False` | Show per-iteration progress bar. |
| `--log_file` | path | None | Optional log file path. |
| `--control_genes_file` | path | None | JSON file with precomputed control genes. |
| `--control_genes_dir` | path | None | Directory to save newly generated control genes. |
| `--step` | {`all`,`monte_carlo`,`post_analysis`} | `all` | Run full pipeline or a specific step only. |
| `--parallel` | flag | `False` | Enable parallel execution with auto-selected worker pool. |
| `--parallel_workers` | int | auto (≤ min(CPUs-1, 8)) | Limit worker processes for Monte Carlo iterations. |
| `--serial` | flag | `False` | Force single-core execution (disables parallel pool). |

For the methodological details, see [scDCF/docs/methods.md](scDCF/docs/methods.md).

#### Advanced CLI examples

```bash
# Use CSV gene list with custom columns
python -m scDCF --csv_file magma_genes.csv --h5ad_file data.h5ad --output_dir results/

# Enable parallel processing with 4 workers
python -m scDCF --gene_list_file genes.txt --h5ad_file data.h5ad \
                --cell_types T_cell B_cell --iterations 100 --output_dir results/ \
                --parallel --parallel_workers 4

# Reuse precomputed control genes
python -m scDCF --csv_file genes.csv --h5ad_file data.h5ad \
                --control_genes_file control_genes.json --output_dir results/

# Run only post-analysis step
python -m scDCF --gene_list_file genes.txt --h5ad_file data.h5ad \
                --step post_analysis --output_dir results/
```

#### Quick test with bundled synthetic data

```bash
python -m scDCF \
  --h5ad_file data/test/sim_adata.h5ad \
  --gene_list_file data/test/genes.txt \
  --control_genes_file data/test/control_genes.json \
  --output_dir quick_test \
  --celltype_column cell_type \
  --disease_marker disease_numeric \
  --rna_count_column nCount_RNA \
  --cell_types T_cell B_cell \
  --iterations 2 \
  --show_progress
```

## 5. Datasets and Methods

### GWAS Gene Selection
scDCF accepts MAGMA- or TWAS-derived gene sets as input. Readers should define and apply their own study-specific selection criteria (e.g., p-value thresholds, top-N rules) appropriate to their dataset and statistical power.

### scRNA-seq Requirements
The framework works with standard scRNA-seq datasets, but performs best with:

- At least 1,000 cells per condition
- Clear cell type annotations
- Matched healthy controls

### Statistical Approach
scDCF implements a rigorous statistical framework:

1. **Library-size matching**: Each target cell matched to 1,000 nearest healthy cells by RNA count; 100 sampled per Monte Carlo iteration
2. **Control gene selection**: 10 control genes per prioritized gene, matched on mean/variance/CV within cell type and disease status
3. **Difference-of-differences**: Target-reference differences minus control-reference differences, weighted by MAGMA Z-scores
4. **Fisher meta-analysis**: Iteration-level p-values combined via Fisher's method; Benjamini-Hochberg FDR correction across cells
5. **Cell-type enrichment**: Fisher's exact test on disease-associated cell proportions between patient and control groups

## 6. Data Sources

See [data/DATA_SOURCES.md](data/DATA_SOURCES.md) for information about the datasets used in scDCF analyses, including SLE, SJS, and CKD datasets with download links.

## 7. Contact

For questions or further information, please contact Caicai Zhang at u3009162@connect.hku.hk.

## 8. License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
