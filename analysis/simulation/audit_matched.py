#!/usr/bin/env python3
"""Audit the matched sign-reversal pair and tabulate fixed-direction results."""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread


HERE = Path(__file__).resolve().parent
METHODS = ["scDCF", "scDRS", "Seurat", "AUCell", "VAM"]


def audit_count_support(first, second, plan: pd.DataFrame) -> dict:
    """Check count changes without requiring equality after clipping/rounding."""
    assert first.shape == second.shape
    rows = plan.gene_index.to_numpy(dtype=int) - 1
    cols = plan.cell_index.to_numpy(dtype=int) - 1
    keys = rows * first.shape[1] + cols
    assert len(np.unique(keys)) == len(keys), "Duplicated gene-cell effects"
    assert ((rows >= 0) & (rows < first.shape[0])).all()
    assert ((cols >= 0) & (cols < first.shape[1])).all()
    difference = (first - second).tocoo()
    difference.eliminate_zeros()
    changed_keys = difference.row * first.shape[1] + difference.col
    assert np.isin(changed_keys, keys).all(), "Counts changed outside the perturbation plan"
    first_values = np.asarray(first.tocsr()[rows, cols]).ravel()
    second_values = np.asarray(second.tocsr()[rows, cols]).ravel()
    ordered = np.sign(plan.signed_delta.to_numpy()) * (first_values - second_values)
    assert (ordered >= 0).all(), "Observed counts contradict intended sign reversal"
    return {
        "matrix_shape_genes_by_cells": list(first.shape),
        "planned_gene_cell_pairs": int(len(plan)),
        "changed_gene_cell_pairs": int(difference.nnz),
        "unchanged_planned_pairs_after_clipping_rounding": int((ordered == 0).sum()),
        "counts_identical_outside_perturbation_plan": True,
        "observed_pairwise_count_order_matches_intended_sign": True,
    }


def audit(root: Path, seed: int):
    scenarios = [root / f"T2_mixed{fraction}down_clean_dense_geneDelta_seed{seed}"
                 for fraction in [40, 60]]
    runner_path = HERE / "run_benchmark.py"
    spec = importlib.util.spec_from_file_location("directionality_benchmark", runner_path)
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    for scenario in scenarios:
        manifest, _ = runner.validate_scenario(scenario)
        assert manifest["pair_mode"] == "sign_reverse"
    runner.validate_pairs(root, scenarios)
    status = json.loads((root / "run_status.json").read_text())
    assert status["complete"], "Do not report an incomplete benchmark"
    method_markers = [json.loads((scenario / "scDCF_output_all/benchmark_complete.json").read_text())
                      for scenario in scenarios]
    assert method_markers[0]["source_sha256"] == method_markers[1]["source_sha256"], (
        "scDCF source changed between scenarios")
    data_a, data_b = [scenario / "data" for scenario in scenarios]
    plan = pd.read_csv(data_a / "perturbation_plan.csv")
    counts_a = mmread(data_a / "counts.mtx").tocsr()
    counts_b = mmread(data_b / "counts.mtx").tocsr()
    report = audit_count_support(counts_a, counts_b, plan)
    report.update({"seed": seed, "fractions_percent": [40, 60],
                   "exact_reversal_of_intended_log_effects": True,
                   "same_genes_cells_labels_and_effect_magnitudes": True,
                   "same_scdcf_source_hashes_across_scenarios": True,
                   "n_replicates": 1,
                   "input_sha256": {
                       scenario.name: {name: runner.digest(scenario / "data" / name)
                                       for name in ["counts.mtx", "metadata.csv",
                                                    "magma_genes_with_z.csv", "perturbation_plan.csv"]}
                       for scenario in scenarios}})
    (root / "matched_count_audit.json").write_text(json.dumps(report, indent=2) + "\n")

    direction_rows = []
    for fraction, scenario in zip([40, 60], scenarios):
        directions = pd.read_csv(scenario / "data/direction_audit.csv")
        row = {"down_percent": fraction, "n_perturbed_cells": len(directions)}
        row.update(directions.drop(columns="cell_id").sum().to_dict())
        direction_rows.append(row)
    pd.DataFrame(direction_rows).to_csv(root / "observed_direction_totals.csv", index=False)

    metrics = pd.read_csv(root / "metrics_fixed_direction.csv")
    metrics = metrics[(metrics.seed == seed) & metrics.method.isin(METHODS)]
    columns = ["compartment", "method", "n_cells", "n_positive", "auroc",
               "auprc_average_precision", "mean_score_positive", "mean_score_negative"]
    first = metrics[np.isclose(metrics.down_fraction, 0.4)][columns]
    second = metrics[np.isclose(metrics.down_fraction, 0.6)][columns]
    comparison = first.merge(second, on=["compartment", "method", "n_cells", "n_positive"],
                             suffixes=("_40", "_60"), validate="one_to_one")
    assert len(comparison) == len(METHODS) * 3
    for fraction in [40, 60]:
        comparison[f"mean_score_gap_{fraction}"] = (
            comparison[f"mean_score_positive_{fraction}"] -
            comparison[f"mean_score_negative_{fraction}"])
    comparison["auroc_change_60_minus_40"] = comparison.auroc_60 - comparison.auroc_40
    comparison.to_csv(root / "matched_method_comparison.csv", index=False)
    lines = [r"\begin{tabular}{@{}lrrrr@{}}", r"\toprule",
             r" & \multicolumn{2}{c}{40\% down} & \multicolumn{2}{c}{60\% down} \\",
             r"\cmidrule(lr){2-3}\cmidrule(lr){4-5}",
             r"Method & AUROC & AP & AUROC & AP \\", r"\midrule"]
    for method in METHODS:
        row = comparison[(comparison.compartment == "all_cells") &
                         (comparison.method == method)].iloc[0]
        lines.append(f"{method} & {row.auroc_40:.3f} & {row.auprc_average_precision_40:.3f} "
                     f"& {row.auroc_60:.3f} & {row.auprc_average_precision_60:.3f} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (root / "matched_comparison_table.tex").write_text("\n".join(lines) + "\n")
    print(json.dumps(report, indent=2))
    print(comparison.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=HERE / "output")
    parser.add_argument("--seed", type=int, default=46)
    args = parser.parse_args()
    audit(args.output_root.resolve(), args.seed)


if __name__ == "__main__":
    main()
