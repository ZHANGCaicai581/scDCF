"""Small independent checks for the paired count-matrix audit."""

import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from audit_matched import audit_count_support
from run_benchmark import check_reference, table_metrics


class CountSupportTests(unittest.TestCase):
    def setUp(self):
        self.plan = pd.DataFrame({"gene_index": [1, 2], "cell_index": [2, 2],
                                  "signed_delta": [0.5, -0.5]})
        self.first = csr_matrix([[3, 8], [2, 2]], dtype=np.int64)
        self.second = csr_matrix([[3, 2], [2, 8]], dtype=np.int64)

    def test_exact_pair_passes(self):
        result = audit_count_support(self.first, self.second, self.plan)
        self.assertEqual(result["changed_gene_cell_pairs"], 2)

    def test_count_ties_are_allowed(self):
        result = audit_count_support(self.first, self.first, self.plan)
        self.assertEqual(result["unchanged_planned_pairs_after_clipping_rounding"], 2)

    def test_change_outside_plan_fails(self):
        altered = csr_matrix([[4, 2], [2, 8]])
        with self.assertRaisesRegex(AssertionError, "outside"):
            audit_count_support(self.first, altered, self.plan)

    def test_wrong_direction_fails(self):
        with self.assertRaisesRegex(AssertionError, "contradict"):
            audit_count_support(self.second, self.first, self.plan)

    def test_duplicate_effects_fail(self):
        duplicate = pd.concat([self.plan, self.plan.iloc[:1]], ignore_index=True)
        with self.assertRaisesRegex(AssertionError, "Duplicated"):
            audit_count_support(self.first, self.second, duplicate)


class TableDefinitionTests(unittest.TestCase):
    def setUp(self):
        self.reference = pd.read_csv(Path(__file__).with_name("reference_metrics.csv"))

    def test_saved_table_values_pass(self):
        check_reference(self.reference)

    def test_missing_method_fails(self):
        with self.assertRaises(ValueError):
            check_reference(self.reference.iloc[1:])

    def test_changed_metric_fails(self):
        wrong = self.reference.copy()
        wrong.loc[0, "average_precision"] += 0.01
        with self.assertRaisesRegex(ValueError, "average_precision"):
            check_reference(wrong)

    def test_reference_aucell_cannot_replace_full_rank(self):
        scores = pd.DataFrame([
            {"seed": 46, "down_fraction": 0.6, "compartment": "all_cells",
             "method": "AUCell", "auroc": 0.927, "auprc_average_precision": 0.406},
            {"seed": 46, "down_fraction": 0.6, "compartment": "all_cells",
             "method": "LegacyFullRank", "auroc": 0.4242, "auprc_average_precision": 0.0135},
        ])
        result = table_metrics(scores)
        self.assertEqual(result.method.tolist(), ["Full-rank enrichment"])
        self.assertEqual(result.source_score.tolist(), ["LegacyFullRank"])
        self.assertEqual(result.auroc.tolist(), [0.4242])


if __name__ == "__main__":
    unittest.main()
