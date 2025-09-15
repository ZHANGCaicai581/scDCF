open# scDCF

*A Framework for Detecting Disease-associated Cells in Single-cell RNA-seq  
Leveraging Healthy Reference Panels and GWAS Findings*

[![PyPI version](https://img.shields.io/pypi/v/scDCF.svg)](https://pypi.org/project/scDCF/)
[![Python versions](https://img.shields.io/pypi/pyversions/scDCF.svg)](https://pypi.org/project/scDCF/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

![scDCF workflow](scDCF/docs/scDCF_workflow.png)

> **Figure 1  –  scDCF analytical workflow.**  
> GWAS-prioritised genes are paired with matched control genes and tested against a 1 000-cell healthy reference panel via Monte-Carlo difference-of-differences statistics, yielding disease-associated cells and enriched cell types.

---

## Table of Contents
1. [Introduction](#1-introduction)  
2. [Key Features](#2-key-features)  
3. [Installation](#3-installation)  
4. [Quick Start](#4-quick-start)  
5. [Datasets and Methods](#5-datasets-and-methods)  
6. [Reproducing the Paper Results](#6-reproducing-the-paper-results)  
7. [Citation](#7-citation)  
8. [Contact](#8-contact)  
9. [License](#9-license)

## 1. Introduction
Genome-wide association studies (GWAS) have uncovered thousands of risk loci, but the cell types through which these variants act remain unclear. **scDCF (single-cell Disease Cell Finder)** integrates GWAS-derived gene sets with single-cell RNA-seq data, using a library-size-matched healthy reference panel, control-gene matching, and Monte-Carlo statistics to pinpoint cells whose expression profiles are genuinely perturbed by inherited risk.

## 2. Key Features
| Capability | Summary |
|------------|---------|
| **Healthy-panel normalisation** | Controls intra-type transcriptional variability by conditioning on 1 000 matched healthy cells. |
| **Control-gene matching** | Pairs each GWAS gene with 10 expression-profile-matched background genes. |
| **Monte-Carlo difference-of-differences test** | Iteratively samples reference cells + control genes; one-tailed *t*-test → Fisher aggregation. |
| **Cell-type enrichment** | Two-tailed Fisher's exact test on proportions of disease-associated cells. |
| **Scalable & interpretable** | Python ≥ 3.9; Scanpy / Pandas stack; outputs tidy tables + UMAP & density plots. |
| **Flexible gene sets** | Accepts MAGMA, TWAS, or any user-supplied list. |

## 3. Installation
```bash
# Stable release
pip install scDCF

# Development (latest) version
pip install git+https://github.com/ZHANGCaicai581/scDCF.git
```

## 4. Quick Start

```python
import scDCF
import scanpy as sc

# Load preprocessed data
adata = sc.read_h5ad("path/to/data.h5ad")

# Prepare GWAS genes
gwas_genes = ["GENE1", "GENE2", "GENE3", ...]  # or load from file
```

For detailed examples, see the [examples directory](examples/). Also see the methods summary in [scDCF/docs/methods.md](scDCF/docs/methods.md).

### Command Line Usage

```bash
# Run scDCF with basic parameters
# You can now invoke the package directly:
python -m scDCF --h5ad_file data/test/sim_adata.h5ad --gene_list_file data/test/genes.txt --output_dir results \
                --celltype_column cell_type --disease_marker disease_numeric --rna_count_column nCount_RNA --iterations 2

# Or with your dataset:
python -m scDCF --h5ad_file data.h5ad --gene_list_file genes.txt --output_dir results/

# Run with additional options
python -m scDCF --h5ad_file data.h5ad \
                --csv_file magma_genes.csv \
                --iterations 1000 \
                --output_dir results/
```

#### Typical CLI recipes

```bash
# 1) Use a CSV/TSV gene list file (auto-detect column)
python -m scDCF \
  --csv_file magma_genes.csv \
  --h5ad_file data.h5ad \
  --output_dir results/

# 2) Use a plain-text gene list (one gene per line)
python -m scDCF \
  --gene_list_file genes.txt \
  --h5ad_file data.h5ad \
  --output_dir results/

# 3) Specify disease/healthy labels and cell type column
python -m scDCF \
  --csv_file magma_genes.csv \
  --h5ad_file data.h5ad \
  --celltype_column celltype_major \
  --disease_marker disease_numeric \
  --disease_value 1 \
  --healthy_value 0 \
  --output_dir results/

# 4) Limit to selected cell types and increase iterations
python -m scDCF \
  --csv_file magma_genes.csv \
  --h5ad_file data.h5ad \
  --cell_types T_cell B_cell NK_cell \
  --iterations 1000 \
  --output_dir results/

# 5) Reuse precomputed control genes (JSON) or save newly generated ones
python -m scDCF \
  --csv_file magma_genes.csv \
  --h5ad_file data.h5ad \
  --control_genes_file control_genes.json \
  --output_dir results/

# or generate control genes and save them to a directory
python -m scDCF \
  --csv_file magma_genes.csv \
  --h5ad_file data.h5ad \
  --control_genes_dir control_genes/ \
  --output_dir results/

# 6) Run a specific step only (e.g., post_analysis on existing results)
python -m scDCF \
  --csv_file magma_genes.csv \
  --h5ad_file data.h5ad \
  --step post_analysis \
  --iterations 1000 \
  --output_dir results/

# 7) Enable logs and progress bars
python -m scDCF \
  --csv_file magma_genes.csv \
  --h5ad_file data.h5ad \
  --log_file scDCF.log \
  --show_progress \
  --output_dir results/
```

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

For the methodological details, see [scDCF/docs/methods.md](scDCF/docs/methods.md).

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

1. **Control gene matching**: Each GWAS gene is matched to 10 control genes with similar expression properties
2. **Monte Carlo sampling**: Repeated sampling from healthy reference panel
3. **Difference-of-differences test**: Compares disease vs. healthy differential expression
4. **Multiple testing correction**: FDR control for cell-type enrichment

## 6. Reproducing the Paper Results

To reproduce the results from the paper, follow these steps:

1. Ensure you have the necessary dependencies installed.
2. Download the dataset from the specified GEO accession.
3. Run the preprocessing pipeline to clean and annotate the data.
4. Apply the scDCF framework to identify disease-associated cells.
5. Interpret the results and create visualizations.

### Data Sources

See [data/DATA_SOURCES.md](data/DATA_SOURCES.md) for information about the datasets used in scDCF analyses, including SLE, SJS, and CKD datasets with download links.

## 7. Contact

For questions or further information, please contact Caicai Zhang at u3009162@connect.hku.hk.

## 8. License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
