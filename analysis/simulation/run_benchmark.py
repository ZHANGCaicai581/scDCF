#!/usr/bin/env python3
"""Generate and evaluate paired down-regulation simulations."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd
from scipy.io import mmread
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
METHODS = ["scDCF", "scDRS", "Seurat", "AUCell", "VAM"]
SCORE_FILES = {
    "Seurat": "reference_scores/Seurat.csv",
    "AUCell": "reference_scores/AUCell.csv",
    "VAM": "reference_scores/VAM.csv",
    "VAM_distance": "reference_scores/VAM_distance.csv",
    "LegacyFullRank": "reference_scores/LegacyFullRank.csv",
}


TABLE_METHODS = {"scDCF": "scDCF", "scDRS": "scDRS", "Seurat": "Seurat",
                 "Full-rank enrichment": "LegacyFullRank", "VAM": "VAM"}


def table_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    selected = frame[(frame.compartment == "all_cells") &
                     frame.method.isin(TABLE_METHODS.values())].copy()
    selected["source_score"] = selected.method
    selected["method"] = selected.method.map({v: k for k, v in TABLE_METHODS.items()})
    selected = selected.rename(columns={"auprc_average_precision": "average_precision"})
    return selected[["seed", "down_fraction", "method", "source_score", "auroc", "average_precision"]].sort_values(
        ["seed", "down_fraction", "method"])


def check_reference(actual: pd.DataFrame):
    if actual.empty:
        raise ValueError("Reference check requires completed table metrics")
    expected = pd.read_csv(HERE / "reference_metrics.csv")
    keys = ["seed", "down_fraction", "method", "source_score"]
    wanted = expected[expected.seed.isin(actual.seed) & expected.down_fraction.isin(actual.down_fraction)]
    joined = wanted.merge(actual, on=keys, how="outer", suffixes=("_expected", "_actual"),
                          indicator=True, validate="one_to_one")
    if len(wanted) != len(actual) or not joined._merge.eq("both").all():
        raise ValueError("Reference check needs all five table methods for saved seed-46 fractions")
    for metric in ["auroc", "average_precision"]:
        if not (joined[metric + "_expected"].round(4) == joined[metric + "_actual"].round(4)).all():
            raise ValueError(f"{metric} differs from the saved S7/S8 values at four decimal places")
    print(f"Reference check passed: {len(actual)} method/scenario rows", flush=True)


def save_json(path: Path, data):
    path.write_text(json.dumps(data, indent=2) + "\n")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def execution_env():
    env = os.environ.copy()
    cache = Path(tempfile.gettempdir()) / "scdcf_directionality"
    cache.mkdir(exist_ok=True)
    for key, directory in {
        "MPLCONFIGDIR": "matplotlib", "NUMBA_CACHE_DIR": "numba", "XDG_CACHE_HOME": "cache"
    }.items():
        location = cache / directory
        location.mkdir(exist_ok=True)
        env[key] = str(location)
    env["MPLBACKEND"] = "Agg"
    env["PYTHONHASHSEED"] = "0"
    for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
        env[key] = "1"
    return env


def execute(command: list[str], log_path: Path):
    print(f"Running {log_path.parent.name}/{log_path.stem}", flush=True)
    start = time.monotonic()
    with log_path.open("w") as handle:
        handle.write(json.dumps(command) + "\n")
        handle.flush()
        result = subprocess.run(command, cwd=REPO, env=execution_env(), stdout=handle,
                                stderr=subprocess.STDOUT, timeout=3600)
    if result.returncode:
        tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-25:])
        raise RuntimeError(f"Command failed (see {log_path}):\n{tail}")
    print(f"Completed {log_path.stem} in {time.monotonic() - start:.1f}s", flush=True)


def validate_scenario(scenario: Path):
    data = scenario / "data"
    manifest = json.loads((data / "scenario_manifest.json").read_text())
    assert (data / "generation_complete.json").exists(), "Generation is incomplete"
    f = manifest["tp"]["mixed_frac_down"]
    assert 0 < f < 1
    assert np.isclose(manifest["tp"]["mixed_frac_up"], 1 - f)
    meta = pd.read_csv(data / "metadata.csv")
    assert len(meta) == 2400 and meta.Cell_ID.is_unique
    assert meta.TruePositive.eq("TP").sum() == 40
    assert (meta.loc[meta.TruePositive.eq("TP"), "Cell_Type"] == "CellTypeC").all()
    assert (meta.loc[meta.TruePositive.eq("TP"), "Condition"] == "Disease").all()
    audit = pd.read_csv(data / "direction_audit.csv")
    assert len(audit) == 40
    assert np.allclose(audit.intended_down / audit.selected_genes, f)
    return manifest, meta


def prepare(scenario: Path):
    import anndata as ad

    data = scenario / "data"
    h5ad = data / "scdrs_input_raw.h5ad"
    gene_list = data / "gene_list_scdcf.csv"
    gs = data / "magma_genes.gs"
    if all(p.exists() for p in [h5ad, gene_list, gs]):
        return
    manifest, meta = validate_scenario(scenario)
    genes = pd.read_csv(data / "genes.csv").Gene_name.astype(str)
    counts = mmread(data / "counts.mtx").T.tocsr().astype(np.float32)
    assert counts.shape == (len(meta), len(genes))
    obs = meta.set_index("Cell_ID")
    obs["nCount_RNA"] = np.asarray(counts.sum(axis=1)).ravel()
    var = pd.DataFrame(index=pd.Index(genes, name="gene_name"))
    ad.AnnData(X=counts, obs=obs, var=var).write_h5ad(h5ad, compression="gzip")
    magma = pd.read_csv(data / "magma_genes_with_z.csv")
    assert len(magma) == manifest["magma"]["actual_size"]
    magma[["Gene_name", "ZSTAT"]].rename(columns={"Gene_name": "gene_name", "ZSTAT": "zstat"}).to_csv(gene_list, index=False)
    gs.write_text("TRAIT\tGENESET\ntrait1\t" + ",".join(magma.Gene_name) + "\n")
    print(f"Prepared inputs: {scenario.name}", flush=True)


def validate_pairs(root: Path, scenarios: list[Path]):
    rows = []
    for first, second in zip(scenarios[::2], scenarios[1::2]):
        data_a, data_b = first / "data", second / "data"
        for relative in ["metadata.csv", "genes.csv", "causal_genes.csv", "magma_genes_with_z.csv"]:
            assert digest(data_a / relative) == digest(data_b / relative), relative
        a = pd.read_csv(data_a / "perturbation_plan.csv")
        b = pd.read_csv(data_b / "perturbation_plan.csv")
        columns = ["cell_id", "gene_name", "direction_order", "delta_magnitude"]
        pd.testing.assert_frame_equal(a[columns], b[columns])
        manifest_a = json.loads((data_a / "scenario_manifest.json").read_text())
        manifest_b = json.loads((data_b / "scenario_manifest.json").read_text())
        mode = manifest_a.get("pair_mode", "nested")
        assert mode == manifest_b.get("pair_mode", "nested")
        if mode == "sign_reverse":
            assert (a.direction != b.direction).all()
            assert np.allclose(a.signed_delta, -b.signed_delta, rtol=0, atol=1e-14)
            assert np.isclose(manifest_a["tp"]["mixed_frac_down"] + manifest_b["tp"]["mixed_frac_down"], 1)
        else:
            assert ((a.direction != "down") | (b.direction == "down")).all()
        rows.append({"seed": manifest_a["seed"],
                     "same_cells": True, "same_gene_set": True,
                     "same_gene_cell_effect_magnitudes": True,
                     "pair_mode": mode, "exact_sign_reversal": mode == "sign_reverse",
                     "perturbed_gene_cell_pairs": len(a)})
    save_json(root / "paired_design_audit.json", rows)


def run_reference(scenario: Path):
    marker = scenario / "reference_scores/complete.json"
    if marker.exists() and all((scenario / relative).exists() for relative in SCORE_FILES.values()):
        return
    execute(["Rscript", str(HERE / "score_reference_methods.R"), str(scenario)],
            scenario / "reference.log")


def run_scdrs(scenario: Path, n_ctrl: int, executable: str):
    score_file = scenario / "scdrs_out/trait1.score.gz"
    marker = scenario / "scdrs_out/complete.json"
    if marker.exists() and score_file.exists():
        assert json.loads(marker.read_text())["n_ctrl"] == n_ctrl
        return
    out = scenario / "scdrs_out"
    out.mkdir(exist_ok=True)
    execute([
        executable, "compute-score",
        "--h5ad_file", str(scenario / "data/scdrs_input_raw.h5ad"),
        "--h5ad_species", "human", "--gs_file", str(scenario / "data/magma_genes.gs"),
        "--gs_species", "human", "--flag_raw_count", "True",
        "--flag_filter_data", "False", "--n_ctrl", str(n_ctrl),
        "--flag_return_ctrl_norm_score", "False", "--out_folder", str(out),
    ], scenario / "scdrs.log")
    assert score_file.exists()
    save_json(marker, {"n_ctrl": n_ctrl, "score": "norm_score", "complete": True})


def run_scdcf(scenario: Path, iterations: int, seed: int):
    out = scenario / "scDCF_output_all"
    marker = out / "benchmark_complete.json"
    if marker.exists():
        previous = json.loads(marker.read_text())
        assert previous["iterations"] == iterations and previous["random_seed"] == seed
        assert previous["source_sha256"] == {
            str(p.relative_to(REPO)): digest(p) for p in (REPO / "scDCF").glob("*.py")
        }, "Cached scDCF results used a different source version; choose a new output directory"
        return
    execute([
        sys.executable, "-m", "scDCF",
        "--csv_file", str(scenario / "data/gene_list_scdcf.csv"),
        "--h5ad_file", str(scenario / "data/scdrs_input_raw.h5ad"),
        "--output_dir", str(out), "--celltype_column", "Cell_Type",
        "--disease_marker", "Condition", "--disease_value", "Disease",
        "--healthy_value", "Healthy", "--rna_count_column", "nCount_RNA",
        "--iterations", str(iterations), "--random_seed", str(seed), "--serial",
    ], scenario / "scdcf.log")
    for ct in ["CellTypeA", "CellTypeB", "CellTypeC"]:
        assert (out / ct / f"{ct}_final_summary.csv").exists()
    save_json(marker, {"iterations": iterations, "random_seed": seed,
                       "expression_input": "raw counts, matching earlier benchmark pipeline",
                       "source_sha256": {str(p.relative_to(REPO)): digest(p)
                                          for p in (REPO / "scDCF").glob("*.py")},
                       "complete": True})


def load_scores(scenario: Path):
    scores = {}
    for method, relative in SCORE_FILES.items():
        file = scenario / relative
        if file.exists():
            scores[method] = pd.read_csv(file).set_index("cell_id").score
    file = scenario / "scdrs_out/trait1.score.gz"
    if file.exists():
        frame = pd.read_csv(file, sep="\t", index_col=0)
        scores["scDRS"] = frame.norm_score
    files = [scenario / "scDCF_output_all" / ct / f"{ct}_final_summary.csv"
             for ct in ["CellTypeA", "CellTypeB", "CellTypeC"]]
    if all(f.exists() for f in files):
        frame = pd.concat([pd.read_csv(f) for f in files]).set_index("cell_id")
        if "combined_minus_log10_p" in frame:
            scores["scDCF"] = frame.combined_minus_log10_p
        else:
            scores["scDCF"] = -np.log10(frame.combined_p_value_fisher.clip(lower=1e-300))
    return scores


def summarize(root: Path, scenarios: list[Path]):
    rows, integrity = [], []
    for scenario in scenarios:
        manifest, meta = validate_scenario(scenario)
        truth = meta.set_index("Cell_ID")
        scores = load_scores(scenario)
        for method in METHODS + ["LegacyFullRank", "VAM_distance"]:
            if method not in scores:
                integrity.append({"scenario": scenario.name, "method": method, "status": "missing"})
                continue
            series = scores[method]
            assert series.index.is_unique
            assert set(series.index) == set(truth.index), f"Incomplete cell matching: {method}"
            aligned = series.reindex(truth.index).astype(float)
            assert np.isfinite(aligned).all()
            integrity.append({"scenario": scenario.name, "method": method, "status": "ok", "n_cells": len(aligned)})
            for compartment, mask in {
                "all_cells": np.ones(len(truth), dtype=bool),
                "CellTypeC": truth.Cell_Type.eq("CellTypeC").to_numpy(),
                "disease_CellTypeC": (truth.Cell_Type.eq("CellTypeC") & truth.Condition.eq("Disease")).to_numpy(),
            }.items():
                y = truth.TruePositive.eq("TP").to_numpy()[mask].astype(int)
                s = aligned.to_numpy()[mask]
                precision, recall, _ = precision_recall_curve(y, s)
                row = {
                    "scenario": scenario.name, "seed": manifest["seed"],
                    "down_fraction": manifest["tp"]["mixed_frac_down"], "method": method,
                    "compartment": compartment, "n_cells": len(y), "n_positive": int(y.sum()),
                    "baseline_precision": y.mean(), "score_direction": "higher_is_positive",
                    "auroc": roc_auc_score(y, s),
                    "auprc_average_precision": average_precision_score(y, s),
                    "auprc_trapezoid": -float((np.trapezoid if hasattr(np, "trapezoid") else np.trapz)(precision, recall)),
                    "mean_score_positive": s[y == 1].mean(), "mean_score_negative": s[y == 0].mean(),
                }
                rows.append(row)
    if not rows:
        raise ValueError("No completed method scores were found; run the scoring stages first")
    frame = pd.DataFrame(rows)
    frame.to_csv(root / "metrics_fixed_direction.csv", index=False)
    table_metrics(frame).to_csv(root / "table_metrics.csv", index=False)
    pd.DataFrame(integrity).to_csv(root / "cell_alignment_audit.csv", index=False)
    primary = frame[(frame.compartment == "all_cells") & frame.method.isin(METHODS)].copy()
    primary.to_csv(root / "primary_results.csv", index=False)
    if not primary.empty:
        group = primary.groupby(["down_fraction", "method"])
        summary = group.agg(n_replicates=("auroc", "count"),
                            auroc_mean=("auroc", "mean"), auroc_sd=("auroc", "std"),
                            auprc_mean=("auprc_average_precision", "mean"),
                            auprc_sd=("auprc_average_precision", "std")).reset_index()
        summary.to_csv(root / "replicate_summary.csv", index=False)
        lines = [r"\begin{tabular}{llrr}", r"\toprule",
                 r"Down-regulated genes & Method & AUROC & AUPRC (AP) \\", r"\midrule"]
        for fraction in sorted(primary.down_fraction.unique()):
            for method in METHODS:
                sub = primary[(primary.down_fraction == fraction) & (primary.method == method)]
                if sub.empty:
                    continue
                if len(sub) == 1:
                    a, p = f"{sub.auroc.iloc[0]:.3f}", f"{sub.auprc_average_precision.iloc[0]:.3f}"
                else:
                    a = f"${sub.auroc.mean():.3f} \\pm {sub.auroc.std():.3f}$"
                    p = f"${sub.auprc_average_precision.mean():.3f} \\pm {sub.auprc_average_precision.std():.3f}$"
                lines.append(f"{fraction * 100:.0f}\\% & {method} & {a} & {p} \\\\")
        lines += [r"\bottomrule", r"\end{tabular}"]
        (root / "results_table.tex").write_text("\n".join(lines) + "\n")
        print(primary[["down_fraction", "method", "auroc", "auprc_average_precision"]].to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=HERE / "output")
    parser.add_argument("--generate", action="store_true")
    parser.add_argument("--seeds", default="46")
    parser.add_argument("--fractions", default="10,20", help="Two increasing down-regulation percentages")
    parser.add_argument("--pair-mode", choices=["nested", "sign_reverse"], default="nested")
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--n-ctrl", type=int, default=1000)
    parser.add_argument("--scdrs-executable", default="scdrs", help="scDRS CLI on PATH, or its full path in a separate environment")
    parser.add_argument("--check-reference", action="store_true", help="Compare table values against the saved seed-46 S7/S8 results")
    parser.add_argument("--stages", default="prepare,reference,scdrs,scdcf,summarize",
                        help="Comma-separated stages; validate checks generated inputs without scoring")
    args = parser.parse_args()
    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    stages = args.stages.split(",")
    assert set(stages) <= {"prepare", "reference", "scdrs", "scdcf", "summarize", "validate"}
    seeds = [int(s) for s in args.seeds.split(",")]
    fractions = [int(f) for f in args.fractions.split(",")]
    if len(fractions) != 2 or not 0 < fractions[0] < fractions[1] < 100:
        parser.error("Provide two increasing fractions between 0 and 100")
    if args.pair_mode == "sign_reverse" and sum(fractions) != 100:
        parser.error("Exact sign reversal requires complementary fractions summing to 100")
    if len(seeds) != len(set(seeds)) or args.iterations < 1 or args.n_ctrl < 1:
        parser.error("Seeds must be unique, and iterations/n-ctrl must be positive")
    if args.generate:
        execute(["Rscript", str(HERE / "generate_scenarios.R"), str(root),
                 args.seeds, args.fractions, args.pair_mode], root / "generation.log")
    scenarios = [root / f"T2_mixed{f}down_clean_dense_geneDelta_seed{seed}"
                 for seed in seeds for f in fractions]
    for scenario in scenarios:
        manifest, _ = validate_scenario(scenario)
        assert manifest.get("pair_mode", "nested") == args.pair_mode
    validate_pairs(root, scenarios)
    save_json(root / "benchmark_configuration.json", {
        "seeds": seeds, "iterations": args.iterations, "scdrs_n_ctrl": args.n_ctrl,
        "fractions_percent": fractions, "pair_mode": args.pair_mode,
        "score_direction": "Fixed: higher is positive; never oriented using labels",
        "primary_auprc": "sklearn.average_precision_score",
        "reference_aucell": "AUCell package with top 5% gene ranks",
        "replicate_interpretation": "Single realization when one seed; SD unavailable in that case",
        "python": sys.executable,
        "python_versions": {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "scipy", "scikit-learn", "anndata"]},
    })
    failures = []
    for scenario in scenarios:
        print(f"\nScenario: {scenario.name}", flush=True)
        for stage in stages:
            if stage == "summarize":
                continue
            try:
                if stage == "prepare": prepare(scenario)
                elif stage == "reference": run_reference(scenario)
                elif stage == "scdrs": run_scdrs(scenario, args.n_ctrl, args.scdrs_executable)
                elif stage == "scdcf": run_scdcf(scenario, args.iterations, seeds[scenarios.index(scenario) // 2])
            except Exception as exc:
                failures.append({"scenario": scenario.name, "stage": stage, "error": str(exc)})
                print(str(exc), file=sys.stderr, flush=True)
    if "summarize" in stages:
        summarize(root, scenarios)
    missing = [{"scenario": scenario.name, "method": method}
               for scenario in scenarios for method in METHODS + ["LegacyFullRank"]
               if method not in load_scores(scenario)]
    save_json(root / "run_status.json", {"complete": not failures and not missing,
                                         "failures": failures, "missing_methods": missing})
    requested = set()
    if "summarize" in stages:
        requested.update(METHODS + ["LegacyFullRank"])
    if "reference" in stages:
        requested.update(SCORE_FILES)
    if "scdrs" in stages:
        requested.add("scDRS")
    if "scdcf" in stages:
        requested.add("scDCF")
    requested_missing = [{"scenario": scenario.name, "method": method}
                         for scenario in scenarios for method in sorted(requested)
                         if method not in load_scores(scenario)]
    if failures or requested_missing:
        raise SystemExit("Requested stages incomplete; see run_status.json and method logs")
    if args.check_reference:
        if "summarize" not in stages:
            parser.error("--check-reference requires the summarize stage")
        check_reference(pd.read_csv(root / "table_metrics.csv"))


if __name__ == "__main__":
    main()
