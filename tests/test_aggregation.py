"""Unit tests for BurdenBench aggregation functions."""
import unittest
import pandas as pd
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from run_BurdenBench import (
    aggregate_summary,
    compute_overall_by_caller_region,
    compute_median_of_medians,
    format_median_range,
    AggregationError,
)


class TestAggregateSummary(unittest.TestCase):
    """Test aggregation across samples."""

    def test_basic_aggregation(self):
        df = pd.DataFrame({
            "Sample": ["HG001", "HG001", "HG002", "HG002"],
            "Caller": ["bcftools", "bcftools", "bcftools", "bcftools"],
            "Condition": ["AF10", "AF10", "AF10", "AF10"],
            "Type": ["SNV", "INDEL", "SNV", "INDEL"],
            "Subset": ["AllAutosomes", "AllAutosomes", "AllAutosomes", "AllAutosomes"],
            "Filter": ["ALL", "ALL", "ALL", "ALL"],
            "Genotype": ["*", "*", "*", "*"],
            "Recall_gain_pp": [0.3, 0.5, 0.4, 0.6],
            "Additional_TP": [100, 200, 150, 250],
            "Additional_FP": [50, 100, 75, 125],
            "TP_FP_ratio": [2.0, 2.0, 2.0, 2.0],
            "Net_benefit": [50, 100, 75, 125],
        })

        result = aggregate_summary(df)
        self.assertGreater(len(result), 0)
        # Check that formatted columns exist
        self.assertIn("Recall_gain_pp_formatted", result.columns)
        self.assertIn("Additional_TP_formatted", result.columns)

    def test_empty_groups(self):
        df = pd.DataFrame({"A": [], "Metric": []})
        with self.assertRaises(AggregationError):
            aggregate_summary(df, group_cols=["A"])


class TestOverallAggregation(unittest.TestCase):
    """Test overall (summed) aggregation."""

    def test_caller_region_level(self):
        df = pd.DataFrame({
            "Sample": ["HG001", "HG001", "HG002", "HG002"],
            "Caller": ["bcftools", "bcftools", "bcftools", "bcftools"],
            "Condition": ["AF10", "AF10", "AF10", "AF10"],
            "Type": ["SNV", "INDEL", "SNV", "INDEL"],
            "Subset": ["AllAutosomes", "AllAutosomes", "AllAutosomes", "AllAutosomes"],
            "Filter": ["ALL", "ALL", "ALL", "ALL"],
            "Genotype": ["*", "*", "*", "*"],
            "Additional_TP": [100, 200, 150, 250],
            "Additional_FP": [50, 100, 75, 125],
            "QUERY.TP_base": [395000, 54400, 395000, 54400],
            "QUERY.TP_cons": [396000, 54600, 396150, 54650],
            "QUERY.FP_base": [2000, 8000, 2000, 8000],
            "QUERY.FP_cons": [2050, 8100, 2075, 8125],
            "TRUTH.TOTAL_base": [400000, 80000, 400000, 80000],
            "TRUTH.TOTAL_cons": [400000, 80000, 400000, 80000],
            "Subset.IS_CONF.Size_base": [2500000000, 2500000000, 2500000000, 2500000000],
            "Recall_base_raw": [0.985, 0.680, 0.985, 0.680],
            "Recall_cons_raw": [0.988, 0.690, 0.989, 0.695],
            "Precision_base_raw": [0.995, 0.750, 0.995, 0.750],
            "Precision_cons_raw": [0.994, 0.740, 0.993, 0.735],
        })

        result = compute_overall_by_caller_region(df, levels=["caller-region"])
        self.assertGreater(len(result), 0)
        self.assertIn("Overall_Level", result.columns)


class TestMedianOfMedians(unittest.TestCase):
    """Test median-of-medians sensitivity analysis."""

    def test_basic_mom(self):
        df = pd.DataFrame({
            "Sample": ["HG001", "HG001", "HG002", "HG002"],
            "Caller": ["bcftools", "bcftools", "bcftools", "bcftools"],
            "Condition": ["AF10", "AF10", "AF10", "AF10"],
            "Type": ["SNV", "INDEL", "SNV", "INDEL"],
            "Subset": ["AllAutosomes", "CMRG", "AllAutosomes", "CMRG"],
            "Filter": ["ALL", "ALL", "ALL", "ALL"],
            "Genotype": ["*", "*", "*", "*"],
            "Recall_gain_pp": [0.3, 0.5, 0.4, 0.6],
            "Additional_TP": [100, 10, 150, 15],
            "Additional_FP": [50, 5, 75, 8],
            "TP_FP_ratio": [2.0, 2.0, 2.0, 1.875],
            "Net_benefit": [50, 5, 75, 7],
        })

        result = compute_median_of_medians(df, levels=["caller-condition"])
        self.assertGreater(len(result), 0)
        self.assertIn("MoM_Level", result.columns)


if __name__ == "__main__":
    unittest.main()
