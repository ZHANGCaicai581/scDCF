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
pip install git+https://github.com/YourUserName/scDCF.git
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

For detailed examples, see the [examples directory](examples/).

### Command Line Usage

```bash
# Run scDCF with basic parameters
python -m scDCF --input data.h5ad --gwas-genes genes.txt --output results/

# Run with additional options
python -m scDCF --input data.h5ad \
                --gwas-genes genes.txt \
                --control-count 10 \
                --iterations 1000 \
                --output results/
```

## 5. Datasets and Methods

### GWAS Gene Selection
scDCF uses MAGMA or TWAS-derived gene sets as input. For optimal performance:

- Select a statistical threshold appropriate for your dataset (typically p < 0.05)
- Consider the top 300-1000 ranked genes depending on the statistical power of your GWAS

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

## 7. Citation

If you use scDCF in your research, please cite our preprint:

```
Zhang, C. (2023). scDCF: A Framework for Detecting Disease-associated Cells in Single-cell RNA-seq Leveraging Healthy Reference Panels and GWAS Findings.
```

## 8. Contact

For questions or further information, please contact Caicai Zhang at u3009162@connect.hku.hk.

## 9. License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
