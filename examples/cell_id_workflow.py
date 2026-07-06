#!/usr/bin/env python3
"""
Example: inspect a final scDCF summary with metadata columns.
"""

from pathlib import Path

import pandas as pd


def main():
    summary_path = Path("results") / "T_cell" / "T_cell_final_summary.csv"
    summary = pd.read_csv(summary_path)

    print(f"Loaded {len(summary)} rows from {summary_path}")
    print("\nColumns:")
    print(summary.columns.tolist())

    significant = summary[summary["scdcf_significant"] == True].copy()
    print(f"\nSignificant cells: {len(significant)}")

    keep_cols = [
        col for col in [
            "cell_id",
            "group",
            "q_cell",
            "combined_p_value_fisher",
            "original_cell_name",
            "sample",
            "batch",
        ]
        if col in significant.columns
    ]

    print("\nTop significant cells:")
    print(significant.sort_values("q_cell").loc[:, keep_cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
