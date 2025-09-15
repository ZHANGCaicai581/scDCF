#!/usr/bin/env python

from setuptools import setup, find_packages
import os

# Read the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Read version from __init__.py
def get_version():
    with open("scDCF/__init__.py", "r") as f:
        for line in f:
            if line.startswith("__version__"):
                return line.split("=")[1].strip().strip("'\"")
    return "0.1.11"

setup(
    name="scDCF",
    version=get_version(),
    author="Caicai Zhang",
    author_email="u3009162@connect.hku.hk",
    description="A Framework for Detecting Disease-associated Cells in Single-cell RNA-seq Leveraging Healthy Reference Panels and GWAS Findings",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ZHANGCaicai581/scDCF",
    project_urls={
        "Bug Tracker": "https://github.com/ZHANGCaicai581/scDCF/issues",
        "Repository": "https://github.com/ZHANGCaicai581/scDCF",
    },
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    install_requires=[
        "scanpy>=1.9.0",
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "tqdm>=4.60.0",
        "statsmodels>=0.12.0",
        "anndata>=0.8.0",
    ],
    entry_points={
        "console_scripts": [
            "scDCF=scDCF.main:main",
        ],
    },
    package_data={
        "scDCF": ["docs/*"],
    },
    include_package_data=True,
    keywords="single-cell genomics GWAS bioinformatics disease cell-type",
    license="MIT",
)