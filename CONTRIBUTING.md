# Contributing to scDCF

Thank you for your interest in contributing to scDCF! We welcome contributions from the community.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally
3. Create a new branch for your changes
4. Make your changes and test them
5. Submit a pull request

## Development Setup

```bash
git clone https://github.com/YourUsername/scDCF.git
cd scDCF
pip install -e .
```

## Testing

Run the test suite with synthetic data:

```bash
python -m scDCF \
  --h5ad_file data/test/sim_adata.h5ad \
  --gene_list_file data/test/genes.txt \
  --control_genes_file data/test/control_genes.json \
  --output_dir test_run \
  --celltype_column cell_type \
  --disease_marker disease_numeric \
  --rna_count_column nCount_RNA \
  --cell_types T_cell B_cell \
  --iterations 2
```

## Reporting Issues

Please use the GitHub issue tracker to report bugs or request features. Include:

- A clear description of the problem
- Steps to reproduce the issue
- Your environment (Python version, OS, etc.)
- Sample data or code if applicable

## Code Style

- Follow PEP 8 for Python code
- Use meaningful variable and function names
- Add docstrings for new functions
- Keep functions focused and modular

## Pull Request Guidelines

- Include a clear description of your changes
- Add tests for new functionality
- Update documentation if needed
- Ensure all tests pass
- Keep commits focused and atomic

## Questions?

Feel free to open an issue for questions about contributing or reach out to the maintainers.
