#!/usr/bin/env python
"""
Minimal end-to-end smoke test for the public scDCF pipeline.

This test exercises the packaged CLI on the bundled synthetic dataset and
verifies that the default user-facing outputs are final summaries only.
"""

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import scanpy as sc
REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = Path(tempfile.gettempdir()) / "scdcf_pipeline_smoke"
INTERMEDIATE_OUTPUT_DIR = Path(tempfile.gettempdir()) / "scdcf_pipeline_smoke_intermediate"


def _assert_file(path):
    if not path.exists():
        raise AssertionError(f"Expected file not found: {path}")


def main():
    if OUTPUT_DIR.exists():
        shutil.rmtree(OUTPUT_DIR)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if INTERMEDIATE_OUTPUT_DIR.exists():
        shutil.rmtree(INTERMEDIATE_OUTPUT_DIR)
    INTERMEDIATE_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mpl"))
    env.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "cache"))
    env.setdefault("NUMBA_CACHE_DIR", str(Path(tempfile.gettempdir()) / "numba"))

    cmd = [
        sys.executable,
        "-m",
        "scDCF",
        "--h5ad_file", str(REPO_ROOT / "data" / "test" / "sim_adata.h5ad"),
        "--gene_list_file", str(REPO_ROOT / "data" / "test" / "genes.txt"),
        "--control_genes_file", str(REPO_ROOT / "data" / "test" / "control_genes.json"),
        "--output_dir", str(OUTPUT_DIR),
        "--celltype_column", "cell_type",
        "--cell_types", "T_cell", "B_cell",
        "--disease_marker", "disease_numeric",
        "--rna_count_column", "nCount_RNA",
        "--iterations", "1",
        "--random_seed", "123"
    ]

    result = subprocess.run(cmd, cwd=REPO_ROOT, env=env, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise SystemExit(result.returncode)

    enrichment_path = OUTPUT_DIR / "celltype_enrichment_summary.csv"
    _assert_file(enrichment_path)
    enrichment_df = pd.read_csv(enrichment_path)
    required_enrichment_cols = {"cell_type", "p_value", "q_type", "disease_enriched"}
    if not required_enrichment_cols.issubset(enrichment_df.columns):
        raise AssertionError(
            f"Missing enrichment columns: {required_enrichment_cols - set(enrichment_df.columns)}"
        )

    for cell_type in ["T_cell", "B_cell"]:
        disease_mc_path = OUTPUT_DIR / cell_type / f"{cell_type}_disease_monte_carlo_results.csv"
        healthy_mc_path = OUTPUT_DIR / cell_type / f"{cell_type}_healthy_monte_carlo_results.csv"
        combined_path = OUTPUT_DIR / cell_type / f"{cell_type}_disease_combined.csv"
        summary_path = OUTPUT_DIR / cell_type / f"{cell_type}_final_summary.csv"

        _assert_file(summary_path)

        summary_df = pd.read_csv(summary_path)

        if disease_mc_path.exists():
            raise AssertionError(f"Intermediate Monte Carlo output should not exist by default: {disease_mc_path}")
        if healthy_mc_path.exists():
            raise AssertionError(f"Intermediate Monte Carlo output should not exist by default: {healthy_mc_path}")
        if combined_path.exists():
            raise AssertionError(f"Intermediate combined output should not exist by default: {combined_path}")
        if "q_cell" not in summary_df.columns:
            raise AssertionError(f"'q_cell' missing from {summary_path}")
        if "scdcf_significant" not in summary_df.columns:
            raise AssertionError(f"'scdcf_significant' missing from {summary_path}")

    intermediate_cmd = [
        sys.executable,
        "-m",
        "scDCF",
        "--h5ad_file", str(REPO_ROOT / "data" / "test" / "sim_adata.h5ad"),
        "--gene_list_file", str(REPO_ROOT / "data" / "test" / "genes.txt"),
        "--control_genes_file", str(REPO_ROOT / "data" / "test" / "control_genes.json"),
        "--output_dir", str(INTERMEDIATE_OUTPUT_DIR),
        "--celltype_column", "cell_type",
        "--cell_types", "T_cell",
        "--disease_marker", "disease_numeric",
        "--rna_count_column", "nCount_RNA",
        "--iterations", "1",
        "--random_seed", "123",
        "--export_intermediate"
    ]

    result = subprocess.run(intermediate_cmd, cwd=REPO_ROOT, env=env, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise SystemExit(result.returncode)

    adata = sc.read_h5ad(REPO_ROOT / "data" / "test" / "sim_adata.h5ad")
    disease_mc_path = INTERMEDIATE_OUTPUT_DIR / "T_cell" / "T_cell_disease_monte_carlo_results.csv"
    healthy_mc_path = INTERMEDIATE_OUTPUT_DIR / "T_cell" / "T_cell_healthy_monte_carlo_results.csv"
    combined_path = INTERMEDIATE_OUTPUT_DIR / "T_cell" / "T_cell_disease_combined.csv"

    _assert_file(disease_mc_path)
    _assert_file(healthy_mc_path)
    _assert_file(combined_path)

    disease_mc_df = pd.read_csv(disease_mc_path)
    healthy_mc_df = pd.read_csv(healthy_mc_path)
    combined_df = pd.read_csv(combined_path)

    required_mc_cols = {
        "reference_group",
        "reference_pool_size",
        "reference_sample_size",
        "reference_self_excluded"
    }
    if not required_mc_cols.issubset(disease_mc_df.columns):
        raise AssertionError(
            f"Missing Monte Carlo columns in {disease_mc_path}: "
            f"{required_mc_cols - set(disease_mc_df.columns)}"
        )
    if not required_mc_cols.issubset(healthy_mc_df.columns):
        raise AssertionError(
            f"Missing Monte Carlo columns in {healthy_mc_path}: "
            f"{required_mc_cols - set(healthy_mc_df.columns)}"
        )

    ct_mask = adata.obs["cell_type"].astype(str) == "T_cell"
    healthy_count = int((ct_mask & (adata.obs["disease_numeric"] == 0)).sum())
    expected_disease_pool = healthy_count
    expected_healthy_pool = healthy_count - 1

    if not disease_mc_df["reference_group"].astype(str).eq("healthy").all():
        raise AssertionError("Disease Monte Carlo did not use healthy references for T_cell")
    if not healthy_mc_df["reference_group"].astype(str).eq("healthy").all():
        raise AssertionError("Healthy Monte Carlo did not use healthy references for T_cell")
    if disease_mc_df["reference_self_excluded"].astype(bool).any():
        raise AssertionError("Disease Monte Carlo should not self-exclude references for T_cell")
    if not healthy_mc_df["reference_self_excluded"].astype(bool).all():
        raise AssertionError("Healthy Monte Carlo should self-exclude references for T_cell")
    if not disease_mc_df["reference_pool_size"].eq(expected_disease_pool).all():
        raise AssertionError("Unexpected disease reference pool size for T_cell")
    if not healthy_mc_df["reference_pool_size"].eq(expected_healthy_pool).all():
        raise AssertionError("Unexpected healthy reference pool size for T_cell")

    expected_disease_sample = min(100, expected_disease_pool)
    expected_healthy_sample = min(100, expected_healthy_pool)
    if not disease_mc_df["reference_sample_size"].eq(expected_disease_sample).all():
        raise AssertionError("Unexpected disease reference sample size for T_cell")
    if not healthy_mc_df["reference_sample_size"].eq(expected_healthy_sample).all():
        raise AssertionError("Unexpected healthy reference sample size for T_cell")
    if "q_cell" not in combined_df.columns or "scdcf_significant" not in combined_df.columns:
        raise AssertionError(f"Combined output missing q_cell/scdcf_significant: {combined_path}")

    print("scDCF pipeline smoke test passed.")


if __name__ == "__main__":
    main()
