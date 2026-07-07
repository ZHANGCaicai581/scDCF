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

REPO_ROOT = Path(__file__).resolve().parents[1]


def _assert_file(path):
    if not path.exists():
        raise AssertionError(f"Expected file not found: {path}")


def _runtime_env():
    tmp_root = Path(tempfile.gettempdir())
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", str(tmp_root / "mpl"))
    env.setdefault("XDG_CACHE_HOME", str(tmp_root / "cache"))
    env.setdefault("NUMBA_CACHE_DIR", str(tmp_root / "numba"))
    return env


def _run_cli(output_dir, extra_args=None):
    cmd = [
        sys.executable,
        "-m",
        "scDCF",
        "--h5ad_file", str(REPO_ROOT / "data" / "test" / "sim_adata.h5ad"),
        "--gene_list_file", str(REPO_ROOT / "data" / "test" / "genes.txt"),
        "--control_genes_file", str(REPO_ROOT / "data" / "test" / "control_genes.json"),
        "--output_dir", str(output_dir),
        "--celltype_column", "cell_type",
        "--disease_marker", "disease_numeric",
        "--rna_count_column", "nCount_RNA",
        "--iterations", "1",
        "--random_seed", "123",
    ]
    if extra_args:
        cmd.extend(extra_args)

    result = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        env=_runtime_env(),
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise AssertionError(
            "scDCF CLI failed\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def run_smoke_test(base_dir=None):
    base_path = Path(base_dir) if base_dir is not None else Path(tempfile.gettempdir())
    output_dir = base_path / "scdcf_pipeline_smoke"
    intermediate_output_dir = base_path / "scdcf_pipeline_smoke_intermediate"

    if output_dir.exists():
        shutil.rmtree(output_dir)
    if intermediate_output_dir.exists():
        shutil.rmtree(intermediate_output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    intermediate_output_dir.mkdir(parents=True, exist_ok=True)

    _run_cli(output_dir, extra_args=["--cell_types", "T_cell", "B_cell"])

    enrichment_path = output_dir / "celltype_enrichment_summary.csv"
    _assert_file(enrichment_path)
    enrichment_df = pd.read_csv(enrichment_path)
    required_enrichment_cols = {"cell_type", "p_value", "q_type", "disease_enriched"}
    missing_enrichment_cols = required_enrichment_cols - set(enrichment_df.columns)
    if missing_enrichment_cols:
        raise AssertionError(f"Missing enrichment columns: {missing_enrichment_cols}")

    for cell_type in ["T_cell", "B_cell"]:
        disease_mc_path = output_dir / cell_type / f"{cell_type}_disease_monte_carlo_results.csv"
        healthy_mc_path = output_dir / cell_type / f"{cell_type}_healthy_monte_carlo_results.csv"
        combined_path = output_dir / cell_type / f"{cell_type}_disease_combined.csv"
        summary_path = output_dir / cell_type / f"{cell_type}_final_summary.csv"

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

    _run_cli(
        intermediate_output_dir,
        extra_args=["--cell_types", "T_cell", "--export_intermediate"],
    )

    from anndata import read_h5ad

    adata = read_h5ad(REPO_ROOT / "data" / "test" / "sim_adata.h5ad")
    disease_mc_path = intermediate_output_dir / "T_cell" / "T_cell_disease_monte_carlo_results.csv"
    healthy_mc_path = intermediate_output_dir / "T_cell" / "T_cell_healthy_monte_carlo_results.csv"
    combined_path = intermediate_output_dir / "T_cell" / "T_cell_disease_combined.csv"

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
        "reference_self_excluded",
    }
    missing_disease_cols = required_mc_cols - set(disease_mc_df.columns)
    missing_healthy_cols = required_mc_cols - set(healthy_mc_df.columns)
    if missing_disease_cols:
        raise AssertionError(f"Missing Monte Carlo columns in {disease_mc_path}: {missing_disease_cols}")
    if missing_healthy_cols:
        raise AssertionError(f"Missing Monte Carlo columns in {healthy_mc_path}: {missing_healthy_cols}")

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


def test_pipeline_smoke(tmp_path):
    run_smoke_test(base_dir=tmp_path)


def main():
    run_smoke_test()
    print("scDCF pipeline smoke test passed.")


if __name__ == "__main__":
    main()
