"""Unit tests for BurdenBench core computation functions."""
import unittest
import pandas as pd
import numpy as np
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from run_BurdenBench import (
    compute_burden,
    validate_happy_dataframe,
    format_median_range,
    DataValidationError,
)


class TestPPDerivation(unittest.TestCase):
    """Test percentage-point derivation correctness."""

    def test_recall_gain_pp(self):
        base = pd.DataFrame({
            "Type": ["SNV"],
            "Subtype": ["*"],
            "Subset": ["AllAutosomes"],
            "Filter": ["ALL"],
            "Genotype": ["*"],
            "METRIC.Recall": [0.985],
            "METRIC.Precision": [0.995],
            "QUERY.FP": [2000],
            "QUERY.TP": [395000],
            "TRUTH.TOTAL": [400000],
            "Subset.IS_CONF.Size": [2500000000]
        })
        cons = base.copy()
        cons["METRIC.Recall"] = [0.988]
        cons["QUERY.TP"] = [396200]
        cons["QUERY.FP"] = [2150]

        result, _ = compute_burden(base, cons, "HG001", "bcftools", "AF10")
        self.assertAlmostEqual(result["Recall_gain_pp"].iloc[0], 0.3, places=10)
        self.assertAlmostEqual(result["Recall_gain_pp_check"].iloc[0], 0.3, places=10)

    def test_precision_gain_pp(self):
        base = pd.DataFrame({
            "Type": ["SNV"],
            "Subtype": ["*"],
            "Subset": ["AllAutosomes"],
            "Filter": ["ALL"],
            "Genotype": ["*"],
            "METRIC.Recall": [0.985],
            "METRIC.Precision": [0.995],
            "QUERY.FP": [2000],
            "QUERY.TP": [395000],
            "TRUTH.TOTAL": [400000],
            "Subset.IS_CONF.Size": [2500000000]
        })
        cons = base.copy()
        cons["METRIC.Precision"] = [0.993]
        cons["QUERY.TP"] = [396200]
        cons["QUERY.FP"] = [2150]

        result, _ = compute_burden(base, cons, "HG001", "bcftools", "AF10")
        self.assertAlmostEqual(result["Precision_gain_pp"].iloc[0], -0.2, places=10)


class TestTPFPRatio(unittest.TestCase):
    """Test TP:FP ratio edge cases."""

    def test_normal_ratio(self):
        base = pd.DataFrame({
            "Type": ["SNV"],
            "Subtype": ["*"],
            "Subset": ["AllAutosomes"],
            "Filter": ["ALL"],
            "Genotype": ["*"],
            "METRIC.Recall": [0.985],
            "METRIC.Precision": [0.995],
            "QUERY.FP": [2000],
            "QUERY.TP": [395000],
            "TRUTH.TOTAL": [400000],
            "Subset.IS_CONF.Size": [2500000000]
        })
        cons = base.copy()
        cons["QUERY.TP"] = [396200]  # +1200 TP
        cons["QUERY.FP"] = [2150]    # +150 FP

        result, _ = compute_burden(base, cons, "HG001", "bcftools", "AF10")
        self.assertEqual(result["Additional_TP"].iloc[0], 1200)
        self.assertEqual(result["Additional_FP"].iloc[0], 150)
        self.assertAlmostEqual(result["TP_FP_ratio_raw"].iloc[0], 8.0, places=10)
        self.assertEqual(result["Net_benefit"].iloc[0], 1050)

    def test_pure_gain(self):
        """TP > 0, FP = 0 → inf raw ratio."""
        base = pd.DataFrame({
            "Type": ["SNV"],
            "Subtype": ["*"],
            "Subset": ["AllAutosomes"],
            "Filter": ["ALL"],
            "Genotype": ["*"],
            "METRIC.Recall": [0.985],
            "METRIC.Precision": [0.995],
            "QUERY.FP": [2000],
            "QUERY.TP": [395000],
            "TRUTH.TOTAL": [400000],
            "Subset.IS_CONF.Size": [2500000000]
        })
        cons = base.copy()
        cons["QUERY.TP"] = [396200]  # +1200 TP
        cons["QUERY.FP"] = [2000]    # +0 FP

        result, _ = compute_burden(base, cons, "HG001", "bcftools", "AF10")
        self.assertTrue(np.isinf(result["TP_FP_ratio_raw"].iloc[0]))
        self.assertEqual(result["Net_benefit"].iloc[0], 1200)

    def test_cost_only(self):
        """TP = 0, FP > 0 → 0 raw ratio."""
        base = pd.DataFrame({
            "Type": ["SNV"],
            "Subtype": ["*"],
            "Subset": ["AllAutosomes"],
            "Filter": ["ALL"],
            "Genotype": ["*"],
            "METRIC.Recall": [0.985],
            "METRIC.Precision": [0.995],
            "QUERY.FP": [2000],
            "QUERY.TP": [395000],
            "TRUTH.TOTAL": [400000],
            "Subset.IS_CONF.Size": [2500000000]
        })
        cons = base.copy()
        cons["QUERY.TP"] = [395000]  # +0 TP
        cons["QUERY.FP"] = [2150]    # +150 FP

        result, _ = compute_burden(base, cons, "HG001", "bcftools", "AF10")
        self.assertEqual(result["TP_FP_ratio_raw"].iloc[0], 0.0)
        self.assertEqual(result["Net_benefit"].iloc[0], -150)

    def test_no_change(self):
        """TP = 0, FP = 0 → NaN raw ratio."""
        base = pd.DataFrame({
            "Type": ["SNV"],
            "Subtype": ["*"],
            "Subset": ["AllAutosomes"],
            "Filter": ["ALL"],
            "Genotype": ["*"],
            "METRIC.Recall": [0.985],
            "METRIC.Precision": [0.995],
            "QUERY.FP": [2000],
            "QUERY.TP": [395000],
            "TRUTH.TOTAL": [400000],
            "Subset.IS_CONF.Size": [2500000000]
        })
        cons = base.copy()  # Identical

        result, _ = compute_burden(base, cons, "HG001", "bcftools", "AF10")
        self.assertTrue(np.isnan(result["TP_FP_ratio_raw"].iloc[0]))
        self.assertTrue(np.isnan(result["TP_FP_ratio_adj"].iloc[0]))
        self.assertEqual(result["Net_benefit"].iloc[0], 0)


class TestFormatMedianRange(unittest.TestCase):
    """Test median range formatting."""

    def test_normal(self):
        self.assertEqual(format_median_range(4.75, 2.27, 7.11), "4.75 [2.27-7.11]")

    def test_single_sample(self):
        self.assertEqual(format_median_range(5.0, 5.0, 5.0), "5.0")

    def test_precision_adaptation(self):
        self.assertEqual(format_median_range(150, 100, 200), "150 [100-200]")
        self.assertEqual(format_median_range(15.5, 10.2, 20.8), "15.5 [10.2-20.8]")
        self.assertEqual(format_median_range(1.55, 1.02, 2.08), "1.55 [1.02-2.08]")
        self.assertEqual(format_median_range(0.155, 0.102, 0.208), "0.155 [0.102-0.208]")

    def test_nan(self):
        self.assertEqual(format_median_range(np.nan, 1.0, 2.0), "N/A [1.000-2.000]")


class TestValidation(unittest.TestCase):
    """Test input validation."""

    def test_empty_dataframe(self):
        df = pd.DataFrame()
        with self.assertRaises(DataValidationError):
            validate_happy_dataframe(df, "test")

    def test_missing_columns(self):
        df = pd.DataFrame({"Type": ["SNV"]})
        with self.assertRaises(DataValidationError):
            validate_happy_dataframe(df, "test")


if __name__ == "__main__":
    unittest.main()
