#!/usr/bin/env python
"""
Regression tests for the public organize_output helper.
"""

import json
import sys
import tempfile
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import scDCF
from scDCF.utils import load_control_genes


def _write_csv(path, rows):
    pd.DataFrame(rows).to_csv(path, index=False)


def _build_source_tree(base_dir):
    source_dir = Path(base_dir) / "raw_results"
    cell_type = "T_cell"
    cell_dir = source_dir / cell_type
    cell_dir.mkdir(parents=True, exist_ok=True)

    disease_rows = [
        {"cell_id": "cell_1", "t_stat": 2.5, "p_value": 0.01},
        {"cell_id": "cell_2", "t_stat": 1.8, "p_value": 0.03},
    ]
    healthy_rows = [
        {"cell_id": "healthy_1", "t_stat": 0.7, "p_value": 0.42},
        {"cell_id": "healthy_2", "t_stat": 0.9, "p_value": 0.33},
    ]
    trait_assoc_rows = [
        {"cell_id": "cell_1", "disease_z_score": 1.2, "neg_log10_p_disease": 2.0},
        {"cell_id": "cell_2", "disease_z_score": 0.8, "neg_log10_p_disease": 1.5},
    ]
    trait_score_rows = [{"trait": "age", "correlation": 0.7}]

    _write_csv(cell_dir / f"{cell_type}_disease_monte_carlo_results.csv", disease_rows)
    _write_csv(cell_dir / f"{cell_type}_healthy_monte_carlo_results.csv", healthy_rows)
    _write_csv(cell_dir / f"{cell_type}_disease_monte_carlo_results_iteration1.csv", disease_rows)
    _write_csv(cell_dir / f"{cell_type}_healthy_monte_carlo_results_iteration1.csv", healthy_rows)
    _write_csv(cell_dir / "trait_association_scores.csv", trait_assoc_rows)
    _write_csv(cell_dir / f"{cell_type}_trait_scores.csv", trait_score_rows)

    control_genes = {
        "disease_control_genes": {"GENE1": ["CTRL1", "CTRL2"]},
        "healthy_control_genes": {"GENE1": ["CTRL3", "CTRL4"]},
    }
    with open(source_dir / f"{cell_type}_control_genes.json", "w", encoding="utf-8") as handle:
        json.dump(control_genes, handle, indent=2)

    return source_dir, cell_type, control_genes


def run_organization_test(base_dir):
    source_dir, cell_type, control_genes = _build_source_tree(base_dir)
    dest_dir = Path(base_dir) / "organized_output"

    returned_dir = Path(scDCF.organize_output(str(source_dir), str(dest_dir)))
    if returned_dir != dest_dir:
        raise AssertionError(f"Unexpected return value from organize_output: {returned_dir}")

    cell_dir = dest_dir / cell_type
    supporting_dir = cell_dir / "supporting_data"
    iterations_dir = supporting_dir / "monte_carlo_iterations"
    control_dir = supporting_dir / "control_genes"
    cell_metrics_path = cell_dir / "cell_metrics.csv"

    expected_paths = [
        cell_metrics_path,
        supporting_dir / f"{cell_type}_disease_monte_carlo_results.csv",
        supporting_dir / f"{cell_type}_healthy_monte_carlo_results.csv",
        iterations_dir / f"{cell_type}_disease_monte_carlo_results_iteration1.csv",
        iterations_dir / f"{cell_type}_healthy_monte_carlo_results_iteration1.csv",
        control_dir / f"{cell_type}_control_genes.json",
        dest_dir / f"{cell_type}_control_genes.json",
    ]
    for path in expected_paths:
        if not path.exists():
            raise AssertionError(f"Expected organized file not found: {path}")

    cell_metrics = pd.read_csv(cell_metrics_path)
    expected_columns = {
        "cell_id",
        "t_stat",
        "p_value",
        "cell_type",
        "disease_z_score",
        "neg_log10_p_disease",
        "trait_age_score",
    }
    missing_columns = expected_columns - set(cell_metrics.columns)
    if missing_columns:
        raise AssertionError(f"cell_metrics.csv missing columns: {missing_columns}")
    if not cell_metrics["cell_type"].eq(cell_type).all():
        raise AssertionError("cell_metrics.csv contains an unexpected cell_type value")
    if not cell_metrics["trait_age_score"].eq(0.7).all():
        raise AssertionError("Trait scores were not propagated into cell_metrics.csv")

    disease_ctrl, healthy_ctrl = load_control_genes(control_dir / f"{cell_type}_control_genes.json")
    if disease_ctrl != control_genes["disease_control_genes"]:
        raise AssertionError("Disease control genes were not preserved during organization")
    if healthy_ctrl != control_genes["healthy_control_genes"]:
        raise AssertionError("Healthy control genes were not preserved during organization")


def test_organize_output(tmp_path):
    run_organization_test(tmp_path)


def main():
    with tempfile.TemporaryDirectory(prefix="scdcf_organize_") as tmp_dir:
        run_organization_test(tmp_dir)
    print("scDCF organization test passed.")


if __name__ == "__main__":
    main()
