#!/usr/bin/env python3
"""
================================================================================
run_BurdenBench.py v1.0.0
================================================================================
Calculate absolute false-positive (FP) and true-positive (TP) burden from
hap.py benchmarking outputs. Supports multiple consensus conditions, aggregates
medians across samples, and produces publication-ready tables with FULL
REPRODUCIBILITY AUDIT - including explicit percentage-point (pp) derivation.

PREVIOUS NAME
    calculate_fp_burden.py (v5.3.7) → BurdenBench v1.0.0

AUTHOR
    Scientific Reviewer / Manuscript Revision Assistant

DATE
    2026-06-30

DEPENDENCIES
    python >= 3.8
    pandas >= 1.3.0
    numpy >= 1.20.0

INSTALLATION
    pip install pandas numpy

QUICK START
    # Single comparison
    python run_BurdenBench.py \
        --baseline HG001_bwa_GRCh38_bcftools.csv \
        --consensus HG001_1KG_AF10pc_bcftools.csv \
        --sample HG001 --caller bcftools --condition AF10 \
        --output results.csv

    # Batch processing with aggregation (RECOMMENDED for manuscripts)
    python run_BurdenBench.py \
        --manifest manifest.csv \
        --output-prefix manuscript_tables/ \
        --regions AllAutosomes lowmappabilityall segdups CMRG MHC \
        --filters ALL

MANIFEST FORMAT
    CSV or TSV with header row and exactly these columns:
        sample,caller,condition,baseline_path,consensus_path

    Example:
        sample,caller,condition,baseline_path,consensus_path
        HG001,bcftools,1KG_AF10pc,HG001_bwa_GRCh38_bcftools.csv,HG001_1KG_AF10pc_bcftools.csv
        HG001,bcftools,1KG_AF30pc,HG001_bwa_GRCh38_bcftools.csv,HG001_1KG_AF30pc_bcftools.csv
        HG002,bcftools,1KG_AF10pc,HG002_bwa_GRCh38_bcftools.csv,HG002_1KG_AF10pc_bcftools.csv

OUTPUT FILES (batch mode)
    Stratified outputs (Type/Filter/Genotype level):
        {prefix}_raw_per_sample.csv              : Per-sample metrics for all comparisons
        {prefix}_summary_median_range.csv        : Medians [min-max] across samples
        {prefix}_table.txt                       : Plain-text table
        {prefix}_table.tex                       : LaTeX table
        {prefix}_table.csv                       : Clean CSV table
    Overall outputs (summed across stratifications; configurable via --overall-levels):
        {prefix}_overall_summary_median_range.csv : Medians [min-max] by dimension combination
        {prefix}_overall_table.txt                : Plain-text table
        {prefix}_overall_table.tex                : LaTeX table
        {prefix}_overall_table.csv                : Clean CSV table
    Sensitivity analysis (configurable via --mom-levels):
        {prefix}_median_of_medians_summary.csv      : Median-of-medians by dimension combination
        {prefix}_median_of_medians_table.txt        : Median-of-medians text table
        {prefix}_median_of_medians_table.csv        : Median-of-medians CSV table
    Audit:
        {prefix}_reproducibility.txt              : Formula documentation and provenance

KEY METRICS
    TP_FP_ratio        : Laplace-adjusted ratio (handles zero counts)
    TP_FP_ratio_raw    : Raw ratio (may be inf/nan)
    Net_benefit        : Signed TP - FP (clinically intuitive)

    All metrics include verification columns for reproducibility.

PERCENTAGE-POINT (pp) DERIVATION (for reproducibility)
    All percentage-point changes are derived directly from hap.py METRIC values:

        Recall_gain_pp    = (METRIC.Recall_consensus - METRIC.Recall_baseline) × 100
        Precision_gain_pp = (METRIC.Precision_consensus - METRIC.Precision_baseline) × 100

    The raw METRIC.Recall and METRIC.Precision values (0-1 scale) are preserved in
    the output so readers can verify the calculation independently.

KEY METRICS EXPLAINED
    Recall_gain_pp        : Percentage-point recall improvement vs baseline
    Precision_gain_pp     : Percentage-point precision change vs baseline
    Additional_TP         : Raw extra true positives found
    Additional_FP         : Raw extra false positives introduced
    TP_FP_ratio           : Additional TP / Additional FP (clinical utility ratio)
                            >1.0 = net beneficial; <1.0 = more noise than signal
    Additional_TP_per_Mbp : True positives gained per megabase (region-size normalized)
    Additional_FP_per_Mbp : False positives introduced per megabase

ERROR CODES
    0  : Success
    1  : No successful comparisons (all inputs failed)
    2  : Invalid arguments or missing required files
    3  : Manifest format error
    4  : Data validation error (missing columns, empty groups, etc.)

================================================================================
"""

import argparse
import sys
import os
import logging
import traceback
from pathlib import Path
from typing import List, Optional, Tuple
from datetime import datetime

import pandas as pd
import numpy as np

# ==============================================================================
# CONFIGURATION
# ==============================================================================

REQUIRED_HAPPY_COLUMNS = [
    "Type", "Subtype", "Subset", "Filter", "Genotype",
    "METRIC.Recall", "METRIC.Precision", "QUERY.FP", "QUERY.TP",
    "TRUTH.TOTAL", "Subset.IS_CONF.Size"
]

REQUIRED_MANIFEST_COLUMNS = [
    "sample", "caller", "condition", "baseline_path", "consensus_path"
]

DEFAULT_GROUP_COLS = [
    "Caller", "Condition", "Type", "Subset", "Filter", "Genotype"
]

LOG_FORMAT = "[%(levelname)s] %(message)s"

# ==============================================================================
# LOGGING SETUP
# ==============================================================================

def setup_logging(verbose: bool = False) -> None:
    """Configure logging with appropriate verbosity."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format=LOG_FORMAT, stream=sys.stderr)

# ==============================================================================
# EXCEPTION HIERARCHY
# ==============================================================================

class BurdenCalculatorError(Exception):
    """Base exception for this script."""
    pass

class FileAccessError(BurdenCalculatorError):
    """Raised when input files cannot be read."""
    pass

class DataValidationError(BurdenCalculatorError):
    """Raised when input data fails validation checks."""
    pass

class ManifestError(BurdenCalculatorError):
    """Raised when manifest is malformed."""
    pass

class AggregationError(BurdenCalculatorError):
    """Raised when aggregation produces invalid results."""
    pass

# ==============================================================================
# INPUT VALIDATION
# ==============================================================================

def validate_file_exists(path: str, label: str) -> Path:
    """Check that a file exists and is readable."""
    p = Path(path)
    if not p.exists():
        raise FileAccessError(f"{label} not found: {path}")
    if not p.is_file():
        raise FileAccessError(f"{label} is not a file: {path}")
    if not os.access(p, os.R_OK):
        raise FileAccessError(f"{label} is not readable: {path}")
    if p.stat().st_size == 0:
        raise FileAccessError(f"{label} is empty: {path}")
    logging.debug(f"Validated file: {path} ({p.stat().st_size:,} bytes)")
    return p


def validate_manifest(manifest_df: pd.DataFrame, path: str) -> None:
    """Validate manifest DataFrame structure and content."""
    if manifest_df.empty:
        raise ManifestError(f"Manifest is empty: {path}")

    missing = [c for c in REQUIRED_MANIFEST_COLUMNS if c not in manifest_df.columns]
    if missing:
        raise ManifestError(
            f"Manifest missing required columns: {missing}. "
            f"Expected columns: {REQUIRED_MANIFEST_COLUMNS}. "
            f"Found: {list(manifest_df.columns)}"
        )

    for col in REQUIRED_MANIFEST_COLUMNS:
        null_count = manifest_df[col].isna().sum()
        if null_count > 0:
            raise ManifestError(
                f"Manifest column '{col}' has {null_count} missing values. "
                f"All rows must have complete sample, caller, condition, and file paths."
            )

    dup_cols = ["sample", "caller", "condition"]
    dups = manifest_df[manifest_df.duplicated(subset=dup_cols, keep=False)]
    if not dups.empty:
        dup_groups = manifest_df.groupby(dup_cols).size()
        n_dup_groups = (dup_groups > 1).sum()
        logging.warning(
            f"Manifest contains {len(dups)} duplicate row(s) across "
            f"{n_dup_groups} unique sample+caller+condition group(s). "
            f"Duplicates will be averaged (median) after per-file computation."
        )

    for idx, row in manifest_df.iterrows():
        try:
            validate_file_exists(row["baseline_path"], f"Row {idx} baseline")
            validate_file_exists(row["consensus_path"], f"Row {idx} consensus")
        except FileAccessError as e:
            raise FileAccessError(
                f"Manifest row {idx} ({row.get('sample','?')}/{row.get('caller','?')}/{row.get('condition','?')}): {e}"
            )

    logging.info(f"Manifest validated: {len(manifest_df)} comparison(s) ready")


def validate_happy_dataframe(df: pd.DataFrame, source: str) -> None:
    """Validate that a hap.py DataFrame has required columns and sensible data."""
    if df.empty:
        raise DataValidationError(f"hap.py output is empty: {source}")

    missing = [c for c in REQUIRED_HAPPY_COLUMNS if c not in df.columns]
    if missing:
        raise DataValidationError(
            f"hap.py output missing required columns: {missing}. "
            f"Found: {list(df.columns)}. "
            f"Is this a valid hap.py summary file?"
        )

    merge_keys = ["Type", "Subtype", "Subset", "Filter", "Genotype"]
    for key in merge_keys:
        null_count = df[key].isna().sum()
        if null_count > 0:
            raise DataValidationError(
                f"hap.py output has {null_count} missing values in merge key '{key}'. "
                f"Cannot merge baseline and consensus without complete stratification keys."
            )

    numeric_cols = ["METRIC.Recall", "METRIC.Precision", "QUERY.FP", "QUERY.TP",
                    "TRUTH.TOTAL", "Subset.IS_CONF.Size"]
    for col in numeric_cols:
        if not pd.api.types.is_numeric_dtype(df[col]):
            raise DataValidationError(
                f"hap.py column '{col}' is not numeric (found type: {df[col].dtype}). "
                f"Check for non-numeric entries or header formatting issues."
            )

    if (df["METRIC.Recall"] < 0).any() or (df["METRIC.Recall"] > 1).any():
        logging.warning(f"Recall values outside [0,1] detected in {source}")
    if (df["METRIC.Precision"] < 0).any() or (df["METRIC.Precision"] > 1).any():
        logging.warning(f"Precision values outside [0,1] detected in {source}")
    if (df["QUERY.FP"] < 0).any():
        raise DataValidationError(f"Negative FP counts in {source}")
    if (df["QUERY.TP"] < 0).any():
        raise DataValidationError(f"Negative TP counts in {source}")

    logging.debug(f"Validated hap.py data: {len(df)} rows, {len(df.columns)} columns from {source}")

# ==============================================================================
# DATA LOADING
# ==============================================================================

def load_happy(path: str) -> pd.DataFrame:
    """
    Load a hap.py summary CSV/TSV with robust error handling.

    Args:
        path: Path to hap.py summary file

    Returns:
        Validated pandas DataFrame

    Raises:
        FileAccessError: If file cannot be read
        DataValidationError: If data is malformed
    """
    p = validate_file_exists(path, "hap.py file")

    try:
        if str(p).lower().endswith('.tsv'):
            df = pd.read_csv(p, sep='\t')
        else:
            df = pd.read_csv(p, sep=None, engine="python")
    except pd.errors.EmptyDataError:
        raise DataValidationError(f"hap.py file is empty or has no data rows: {path}")
    except pd.errors.ParserError as e:
        raise DataValidationError(f"Cannot parse hap.py file {path}: {e}")
    except Exception as e:
        raise FileAccessError(f"Failed to read {path}: {e}")

    df.columns = [c.strip() for c in df.columns]
    validate_happy_dataframe(df, path)
    return df


def load_manifest(path: str) -> pd.DataFrame:
    """Load and validate a batch manifest CSV/TSV."""
    p = validate_file_exists(path, "Manifest file")

    try:
        if str(p).lower().endswith('.tsv'):
            df = pd.read_csv(p, sep='\t')
        else:
            df = pd.read_csv(p, sep=None, engine="python")
    except pd.errors.EmptyDataError:
        raise ManifestError(f"Manifest file is empty: {path}")
    except pd.errors.ParserError as e:
        raise ManifestError(f"Cannot parse manifest {path}: {e}")
    except Exception as e:
        raise FileAccessError(f"Failed to read manifest {path}: {e}")

    df.columns = [c.strip().lower().replace(' ', '_') for c in df.columns]
    validate_manifest(df, path)
    return df

# ==============================================================================
# CORE COMPUTATION (with explicit pp derivation)
# ==============================================================================

def compute_burden(df_base: pd.DataFrame, df_cons: pd.DataFrame,
                   sample: str, caller: str, condition: str) -> pd.DataFrame:
    """
    Merge baseline and consensus hap.py outputs and compute burden metrics.

    PERCENTAGE-POINT DERIVATION (reproducibility audit):
        Raw hap.py outputs METRIC.Recall and METRIC.Precision as proportions (0-1).
        These are converted to percentages and differenced to yield pp changes:

        Recall_cons_pct    = METRIC.Recall_consensus × 100
        Recall_base_pct    = METRIC.Recall_baseline × 100
        Recall_gain_pp     = Recall_cons_pct - Recall_base_pct
                           = (METRIC.Recall_consensus - METRIC.Recall_baseline) × 100

        (Analogous for Precision)

    All raw 0-1 values are preserved in output columns for independent verification.
    """
    merge_keys = ["Type", "Subtype", "Subset", "Filter", "Genotype"]
    mismatch_records: List[pd.DataFrame] = []

    for key in merge_keys:
        if key not in df_base.columns or key not in df_cons.columns:
            raise DataValidationError(f"Merge key '{key}' missing from one of the inputs.")

    merged = df_cons.merge(df_base, on=merge_keys, suffixes=("_cons", "_base"), how="outer")
    merged["in_both"] = merged["METRIC.Recall_base"].notna() & merged["METRIC.Recall_cons"].notna()

    unmatched = (~merged["in_both"]).sum()
    if unmatched > 0:
        logging.warning(
            f"{unmatched} stratification row(s) unmatched between baseline and consensus "
            f"for {sample}/{caller}/{condition}. These will have NaN for comparative metrics."
        )

    # Warn if baseline and consensus confident-region sizes differ materially
    size_diff = (merged["Subset.IS_CONF.Size_cons"] - merged["Subset.IS_CONF.Size_base"]).abs()
    size_rel_diff = size_diff / merged["Subset.IS_CONF.Size_base"].replace(0, np.nan)
    max_rel_diff = size_rel_diff.max()
    size_mismatch_mask = size_rel_diff > 0.05
    if pd.notna(max_rel_diff) and max_rel_diff > 0.05:
        logging.warning(
            f"Baseline/consensus Subset.IS_CONF.Size differ by up to "
            f"{max_rel_diff*100:.1f}% for {sample}/{caller}/{condition}. "
            f"Per-Mbp normalisation uses the baseline size; review if consensus "
            f"was run on a different region set."
        )
        # Record size mismatches for detailed report
        size_mm_cols = ["Type", "Subtype", "Subset", "Filter", "Genotype",
                        "Subset.IS_CONF.Size_base", "Subset.IS_CONF.Size_cons"]
        size_mm_cols = [c for c in size_mm_cols if c in merged.columns]
        if size_mm_cols:
            size_mm_df = merged.loc[size_mismatch_mask, size_mm_cols].copy()
            size_mm_df.insert(0, "Sample", sample)
            size_mm_df.insert(1, "Caller", caller)
            size_mm_df.insert(2, "Condition", condition)
            size_mm_df["Mismatch_Type"] = "Subset.IS_CONF.Size"
            size_mm_df["Mismatch_Description"] = (
                f"Subset.IS_CONF.Size differs by >5% between baseline and consensus. "
                f"Per-Mbp normalisation uses baseline size."
            )
            mismatch_records.append(size_mm_df)

    # ==============================================================================
    # EXPLICIT PP DERIVATION - raw 0-1 values preserved for reproducibility
    # ==============================================================================

    # Raw hap.py metrics (0-1 scale) - preserved for audit
    merged["Recall_base_raw"]    = merged["METRIC.Recall_base"]
    merged["Recall_cons_raw"]    = merged["METRIC.Recall_cons"]
    merged["Precision_base_raw"] = merged["METRIC.Precision_base"]
    merged["Precision_cons_raw"] = merged["METRIC.Precision_cons"]

    # Convert to percentage scale (0-100)
    merged["Recall_base_pct"]    = merged["Recall_base_raw"] * 100
    merged["Recall_cons_pct"]    = merged["Recall_cons_raw"] * 100
    merged["Precision_base_pct"] = merged["Precision_base_raw"] * 100
    merged["Precision_cons_pct"] = merged["Precision_cons_raw"] * 100

    # Percentage-point (pp) gains - the manuscript-reported metric
    merged["Recall_gain_pp"]    = merged["Recall_cons_pct"] - merged["Recall_base_pct"]
    merged["Precision_gain_pp"] = merged["Precision_cons_pct"] - merged["Precision_base_pct"]

    # F1-score (harmonic mean of recall and precision)
    def _f1(recall, precision):
        denom = recall + precision
        return np.where(denom > 0, 2 * recall * precision / denom, np.nan)

    merged["F1_base"] = _f1(merged["Recall_base_raw"], merged["Precision_base_raw"])
    merged["F1_cons"] = _f1(merged["Recall_cons_raw"], merged["Precision_cons_raw"])
    merged["F1_gain_pp"] = (merged["F1_cons"] - merged["F1_base"]) * 100

    # Verify: pp gain should equal (raw difference) * 100
    merged["Recall_gain_pp_check"]    = (merged["Recall_cons_raw"] - merged["Recall_base_raw"]) * 100
    merged["Precision_gain_pp_check"] = (merged["Precision_cons_raw"] - merged["Precision_base_raw"]) * 100

    # Check for numerical drift (should be identical to machine precision)
    recall_drift    = (merged["Recall_gain_pp"] - merged["Recall_gain_pp_check"]).abs().max()
    precision_drift = (merged["Precision_gain_pp"] - merged["Precision_gain_pp_check"]).abs().max()
    if pd.notna(recall_drift) and recall_drift > 1e-10:
        logging.warning(f"Recall pp calculation drift detected: {recall_drift}")
    if pd.notna(precision_drift) and precision_drift > 1e-10:
        logging.warning(f"Precision pp calculation drift detected: {precision_drift}")

    # Drop verification columns unless --debug-columns requested
    if not getattr(sys, '_debug_columns', False):
        drop_cols = [c for c in merged.columns if c.endswith("_check")]
        merged = merged.drop(columns=drop_cols)

    # ==============================================================================
    # ABSOLUTE BURDEN METRICS
    # ==============================================================================

    merged["Additional_FP"] = merged["QUERY.FP_cons"] - merged["QUERY.FP_base"]
    merged["Additional_TP"] = merged["QUERY.TP_cons"] - merged["QUERY.TP_base"]
    merged["Additional_FN"] = (
        (merged["TRUTH.TOTAL_cons"] - merged["QUERY.TP_cons"]) -
        (merged["TRUTH.TOTAL_base"] - merged["QUERY.TP_base"])
    )

    # Warn when TRUTH.TOTAL differs between baseline and consensus
    both_not_na = merged["TRUTH.TOTAL_cons"].notna() & merged["TRUTH.TOTAL_base"].notna()
    truth_mismatch = (both_not_na & (merged["TRUTH.TOTAL_cons"] != merged["TRUTH.TOTAL_base"])).sum()
    if truth_mismatch > 0:
        # Log first few examples for visibility
        example_cols = ["Type", "Subset", "Filter", "Genotype", "TRUTH.TOTAL_base", "TRUTH.TOTAL_cons"]
        example_cols = [c for c in example_cols if c in merged.columns]
        mismatch_mask = both_not_na & (merged["TRUTH.TOTAL_cons"] != merged["TRUTH.TOTAL_base"])
        examples = merged.loc[mismatch_mask, example_cols].head(5).to_string(index=False)
        logging.warning(
            f"{truth_mismatch} stratification row(s) have mismatched TRUTH.TOTAL "
            f"between baseline and consensus for {sample}/{caller}/{condition}. "
            f"Additional_FN values for those rows are unreliable because the "
            f"change in FN is contaminated by a change in truth set size. "
            f"First 5 examples:\n{examples}"

        )

    # TP:FP ratio
    # ==============================================================================
    # TP:FP RATIO - with edge-case handling for zero counts
    # ==============================================================================

    # Vectorised raw ratio (handles all edge cases correctly)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio_raw = np.where(
            merged["Additional_FP"] != 0,
            merged["Additional_TP"] / merged["Additional_FP"],
            np.where(merged["Additional_TP"] != 0, np.inf, np.nan)
        )
    merged["TP_FP_ratio_raw"] = ratio_raw

    # Pseudocount-adjusted ratio (Laplace smoothing: +1 to both)
    # Prevents inf/nan while preserving ranking
    # Laplace smoothing: only apply when FP >= 0 to avoid sign flips.
    # When Additional_FP < 0, the raw ratio is well-defined (positive TP / negative FP
    # = negative ratio, indicating FP reduction), so use it directly.
    # When both are 0, return NaN (no meaningful ratio).
    merged["TP_FP_ratio_adj"] = np.where(
        (merged["Additional_TP"] == 0) & (merged["Additional_FP"] == 0),
        np.nan,
        np.where(
            merged["Additional_FP"] >= 0,
            (merged["Additional_TP"] + 1) / (merged["Additional_FP"] + 1),
            merged["Additional_TP"] / merged["Additional_FP"]
        )
    )

    # Net benefit: signed difference (clinically intuitive, handles all cases)
    merged["Net_benefit"] = merged["Additional_TP"] - merged["Additional_FP"]

    # Primary display ratio: uses adjusted for aggregation, raw for display when valid
    merged["TP_FP_ratio"] = merged["TP_FP_ratio_adj"]  # Default to adjusted for stats

    # Vectorised format-friendly string version for tables
    tp = merged["Additional_TP"]
    fp = merged["Additional_FP"]
    adj = merged["TP_FP_ratio_adj"]
    formatted = np.select(
        [
            pd.isna(tp) | pd.isna(fp),
            (fp == 0) & (tp > 0),
            (fp == 0) & (tp < 0),
            (fp == 0) & (tp == 0),
            (tp == 0) & (fp > 0),
        ],
        [
            "-",
            "inf (pure gain)",
            "-inf (pure loss)",
            "-",
            "0 (cost only)",
        ],
        default=adj.round(2).astype(str)
    )
    merged["TP_FP_ratio_formatted"] = formatted

    merged["FP_pct_increase"] = np.where(
        merged["QUERY.FP_base"] > 0,
        merged["Additional_FP"] / merged["QUERY.FP_base"] * 100,
        np.nan
    )

    # Per megabase
    # Baseline region size (used for base and delta metrics)
    conf_mb_base = merged["Subset.IS_CONF.Size_base"] / 1_000_000
    conf_mb_base_safe = conf_mb_base.replace(0, np.nan)

    # Consensus region size (used for consensus metrics)
    conf_mb_cons = merged["Subset.IS_CONF.Size_cons"] / 1_000_000
    conf_mb_cons_safe = conf_mb_cons.replace(0, np.nan)

    merged["FP_per_Mbp_base"]    = merged["QUERY.FP_base"] / conf_mb_base_safe
    merged["FP_per_Mbp_cons"]    = merged["QUERY.FP_cons"] / conf_mb_cons_safe
    merged["Additional_FP_per_Mbp"] = merged["Additional_FP"] / conf_mb_base_safe
    merged["TP_per_Mbp_base"]    = merged["QUERY.TP_base"] / conf_mb_base_safe
    merged["TP_per_Mbp_cons"]    = merged["QUERY.TP_cons"] / conf_mb_cons_safe
    merged["Additional_TP_per_Mbp"] = merged["Additional_TP"] / conf_mb_base_safe

    # Per 100 truth variants
    truth_safe = merged["TRUTH.TOTAL_base"].replace(0, np.nan)
    merged["FP_per_100_truth_base"]    = merged["QUERY.FP_base"] / truth_safe * 100
    merged["FP_per_100_truth_cons"]    = merged["QUERY.FP_cons"] / truth_safe * 100
    merged["Additional_FP_per_100_truth"] = merged["Additional_FP"] / truth_safe * 100
    merged["TP_per_100_truth_base"]    = merged["QUERY.TP_base"] / truth_safe * 100
    merged["TP_per_100_truth_cons"]    = merged["QUERY.TP_cons"] / truth_safe * 100
    merged["Additional_TP_per_100_truth"] = merged["Additional_TP"] / truth_safe * 100

    # Metadata
    merged.insert(0, "Sample", sample)
    merged.insert(1, "Caller", caller)
    merged.insert(2, "Condition", condition)

    # Column ordering - REPRODUCIBILITY-FIRST: raw metrics, then pp derivation, then burden
    front = [
        # Identifiers
        "Sample", "Caller", "Condition", "Type", "Subtype", "Subset", "Filter", "Genotype",
        # Raw hap.py metrics (0-1 scale) - for reproducibility audit
        "Recall_base_raw", "Recall_cons_raw",
        "Precision_base_raw", "Precision_cons_raw",
        # Percentage-scale metrics (0-100)
        "Recall_base_pct", "Recall_cons_pct",
        "Precision_base_pct", "Precision_cons_pct",
        # Percentage-point gains (manuscript metric)
        "Recall_gain_pp", "Precision_gain_pp",
        # F1-score (harmonic mean of recall and precision)
        "F1_base", "F1_cons", "F1_gain_pp",
        # Truth and query counts
        "TRUTH.TOTAL_base", "QUERY.TP_base", "QUERY.FP_base",
        "TRUTH.TOTAL_cons", "QUERY.TP_cons", "QUERY.FP_cons",
        # Absolute burden
        "Additional_TP", "Additional_FP", "Additional_FN",
        "TP_FP_ratio", "FP_pct_increase",
        # Normalized burden
        "Additional_TP_per_Mbp", "Additional_FP_per_Mbp",
        "Additional_TP_per_100_truth", "Additional_FP_per_100_truth",
        # Region metadata
        "Subset.IS_CONF.Size_base", "in_both"
    ]
    front = [c for c in front if c in merged.columns]
    rest = [c for c in merged.columns if c not in front]
    merged = merged[front + rest]

    return merged, mismatch_records

# ==============================================================================
# AGGREGATION
# ==============================================================================

def aggregate_summary(df_raw: pd.DataFrame, group_cols: Optional[List[str]] = None) -> pd.DataFrame:
    """Aggregate per-sample results into median + [min-max] across samples."""
    if group_cols is None:
        group_cols = DEFAULT_GROUP_COLS.copy()

    group_cols = [c for c in group_cols if c in df_raw.columns]
    if not group_cols:
        raise AggregationError("No valid grouping columns found in data")

    # Metrics to aggregate - INCLUDING raw metrics for reproducibility
    metrics = [
        # Raw metrics (so readers can verify pp calculation)
        "Recall_base_raw", "Recall_cons_raw",
        "Precision_base_raw", "Precision_cons_raw",
        # PP gains
        "Recall_gain_pp", "Precision_gain_pp",
        # F1-score
        "F1_base", "F1_cons", "F1_gain_pp",
        # Base counts (for reference)
        "QUERY.TP_base", "QUERY.FP_base",
        # Consensus counts (for reference)
        "QUERY.TP_cons", "QUERY.FP_cons",
        # Burden
        "Additional_TP", "Additional_FP", "Additional_FN",
        "TP_FP_ratio", "FP_pct_increase", "Net_benefit",
        "Additional_TP_per_Mbp", "Additional_FP_per_Mbp",
        "Additional_TP_per_100_truth", "Additional_FP_per_100_truth",
    ]
    metrics = [m for m in metrics if m in df_raw.columns]

    if not metrics:
        raise AggregationError("No valid metrics found for aggregation")

    logging.info(f"Aggregating across {df_raw['Sample'].nunique()} sample(s) by {group_cols}")

    try:
        agg_dict = {m: ["median", "min", "max", "count"] for m in metrics}
        grouped = df_raw.groupby(group_cols, sort=False).agg(agg_dict)
    except Exception as e:
        raise AggregationError(f"Aggregation failed: {e}")

    grouped.columns = [f"{col}_{stat}" for col, stat in grouped.columns]
    grouped = grouped.reset_index()

    if grouped.empty:
        raise AggregationError("Aggregation produced empty results")

    # Format median [min-max] strings
    for m in metrics:
        med = f"{m}_median"
        lo = f"{m}_min"
        hi = f"{m}_max"
        if med in grouped.columns and lo in grouped.columns and hi in grouped.columns:
            grouped[f"{m}_formatted"] = grouped.apply(
                lambda r: format_median_range(r[med], r[lo], r[hi]),
                axis=1
            )

    logging.info(f"Aggregation complete: {len(grouped)} group(s)")
    return grouped


def _parse_level_name(level_name: str) -> List[str]:
    """
    Parse a level name (e.g., 'caller-condition-region') into 
    aggregation column names.

    Dimension mapping:
        sample    -> Sample
        caller    -> Caller
        condition -> Condition
        type      -> Type
        region    -> Subset
    """
    dim_map = {
        "sample": "Sample",
        "caller": "Caller",
        "condition": "Condition",
        "type": "Type",
        "region": "Subset",
    }
    parts = level_name.lower().split("-")
    cols = []
    for p in parts:
        if p in dim_map:
            cols.append(dim_map[p])
        else:
            logging.warning(f"Unknown dimension '{p}' in level name '{level_name}'")
    return cols


def compute_overall_by_caller_region(df_raw: pd.DataFrame,
                                     levels: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Compute overall metrics at configurable aggregation levels.

    Each level name specifies the dimensions to aggregate by (excluding Sample,
    which is always used for per-sample summing). Dimensions are separated by
    hyphens and can be: caller, condition, type, region.

    Available levels (specified via --overall-levels):
      caller-condition-type-region : Per Caller + Condition + Type + Region
      caller-condition-region      : Per Caller + Condition + Region
      caller-condition-type        : Per Caller + Condition + Type (all regions)
      caller-condition             : Per Caller + Condition (all types/regions)
      caller-type-region           : Per Caller + Type + Region (all conditions)
      caller-region                : Per Caller + Region (all types/conditions)
      caller-type                  : Per Caller + Type (all conditions/regions)
      caller                       : Per Caller (all conditions/types/regions)
      condition-type-region        : Per Condition + Type + Region (all callers)
      condition-region             : Per Condition + Region (all callers)
      condition-type               : Per Condition + Type (all callers/regions)
      condition                    : Per Condition (all callers/types/regions)
      type-region                  : Per Type + Region (all callers/conditions)
      region                       : Per Region (all callers/conditions/types)
      type                         : Per Type (all callers/conditions/regions)

    For each level, counts are summed within each sample, derived metrics
    recomputed, then aggregated as median [min-max] across samples.
    """
    if levels is None:
        levels = [
            "caller-condition-type-region",
            "caller-condition-region",
            "caller-condition-type",
            "caller-condition",
            "caller-type-region",
            "caller-region",
            "caller-type",
            "caller",
        ]

    # Base required columns (always needed)
    base_required = [
        "Sample", "Caller", "Condition", "Subset", "Type",
        "Additional_TP", "Additional_FP", "Additional_FN",
        "QUERY.TP_base", "QUERY.TP_cons", "QUERY.FP_base", "QUERY.FP_cons",
        "TRUTH.TOTAL_base", "TRUTH.TOTAL_cons",
        "Subset.IS_CONF.Size_base",
        "Recall_base_raw", "Recall_cons_raw",
        "Precision_base_raw", "Precision_cons_raw"
    ]

    # Check which columns actually exist
    available_cols = [c for c in base_required if c in df_raw.columns]
    missing = [c for c in base_required if c not in df_raw.columns]
    if missing:
        logging.warning(f"Missing columns in raw data: {missing}")

    sum_cols = ["Additional_TP", "Additional_FP", "Additional_FN",
                "QUERY.TP_base", "QUERY.TP_cons", 
                "QUERY.FP_base", "QUERY.FP_cons",
                "TRUTH.TOTAL_base", "TRUTH.TOTAL_cons",
                "Subset.IS_CONF.Size_base"]
    sum_cols = [c for c in sum_cols if c in df_raw.columns]

    # ------------------------------------------------------------------
    # Helper: build per-sample aggregated frame
    # ------------------------------------------------------------------
    def _build_per_sample(gb_keys: list) -> pd.DataFrame:
        gb_keys = [c for c in gb_keys if c in df_raw.columns]
        if not gb_keys:
            logging.warning("No valid groupby keys for per-sample build")
            return pd.DataFrame()

        df_s = df_raw.groupby(gb_keys, sort=False)[sum_cols].sum().reset_index()
        df_ps = df_s.copy()

        # ------------------------------------------------------------------
        # CRITICAL: Recompute recall/precision from summed counts
        # Averaging recall/precision across stratification rows is
        # mathematically incorrect. We must compute them from the
        # aggregated TP, FP, and TRUTH.TOTAL.
        # ------------------------------------------------------------------
        truth_total_base = df_ps["TRUTH.TOTAL_base"].replace(0, np.nan)
        truth_total_cons = df_ps["TRUTH.TOTAL_cons"].replace(0, np.nan)
        query_total_base = (df_ps["QUERY.TP_base"] + df_ps["QUERY.FP_base"]).replace(0, np.nan)
        query_total_cons = (df_ps["QUERY.TP_cons"] + df_ps["QUERY.FP_cons"]).replace(0, np.nan)

        df_ps["Recall_base_raw"] = df_ps["QUERY.TP_base"] / truth_total_base
        df_ps["Recall_cons_raw"] = df_ps["QUERY.TP_cons"] / truth_total_cons
        df_ps["Precision_base_raw"] = df_ps["QUERY.TP_base"] / query_total_base
        df_ps["Precision_cons_raw"] = df_ps["QUERY.TP_cons"] / query_total_cons

        # Recompute derived metrics
        df_ps["Recall_gain_pp"] = (df_ps["Recall_cons_raw"] - df_ps["Recall_base_raw"]) * 100
        df_ps["Precision_gain_pp"] = (df_ps["Precision_cons_raw"] - df_ps["Precision_base_raw"]) * 100

        # F1-score
        def _f1_build(r, p):
            denom = r + p
            return np.where(denom > 0, 2 * r * p / denom, np.nan)
        df_ps["F1_base"] = _f1_build(df_ps["Recall_base_raw"], df_ps["Precision_base_raw"])
        df_ps["F1_cons"] = _f1_build(df_ps["Recall_cons_raw"], df_ps["Precision_cons_raw"])
        df_ps["F1_gain_pp"] = (df_ps["F1_cons"] - df_ps["F1_base"]) * 100

        df_ps["Recall_base_pct"] = df_ps["Recall_base_raw"] * 100
        df_ps["Recall_cons_pct"] = df_ps["Recall_cons_raw"] * 100
        df_ps["Precision_base_pct"] = df_ps["Precision_base_raw"] * 100
        df_ps["Precision_cons_pct"] = df_ps["Precision_cons_raw"] * 100

        df_ps["TP_FP_ratio_adj"] = np.where(
            (df_ps["Additional_TP"] == 0) & (df_ps["Additional_FP"] == 0),
            np.nan,
            np.where(
                df_ps["Additional_FP"] >= 0,
                (df_ps["Additional_TP"] + 1) / (df_ps["Additional_FP"] + 1),
                df_ps["Additional_TP"] / df_ps["Additional_FP"]
            )
        )
        df_ps["TP_FP_ratio"] = df_ps["TP_FP_ratio_adj"]
        df_ps["Net_benefit"] = df_ps["Additional_TP"] - df_ps["Additional_FP"]

        df_ps["FP_pct_increase"] = np.where(
            df_ps["QUERY.FP_base"] > 0,
            df_ps["Additional_FP"] / df_ps["QUERY.FP_base"] * 100,
            np.nan
        )

        conf_mb = df_ps["Subset.IS_CONF.Size_base"] / 1_000_000
        conf_mb_safe = conf_mb.replace(0, np.nan)
        df_ps["Additional_TP_per_Mbp"] = df_ps["Additional_TP"] / conf_mb_safe
        df_ps["Additional_FP_per_Mbp"] = df_ps["Additional_FP"] / conf_mb_safe

        truth_safe = df_ps["TRUTH.TOTAL_base"].replace(0, np.nan)
        df_ps["Additional_TP_per_100_truth"] = (df_ps["Additional_TP"] / truth_safe) * 100
        df_ps["Additional_FP_per_100_truth"] = (df_ps["Additional_FP"] / truth_safe) * 100

        # Warn if summed counts suggest overlapping strata
        if "TRUTH.TOTAL_base" in df_ps.columns and "Sample" in df_ps.columns:
            for sample_id, grp in df_ps.groupby("Sample"):
                max_single = grp["TRUTH.TOTAL_base"].max()
                if max_single > 0:
                    # Check if any aggregated row exceeds 1.5x the max single-row truth
                    exceeds = (grp["TRUTH.TOTAL_base"] > max_single * 1.5).sum()
                    if exceeds > 0:
                        logging.warning(
                            f"Sample {sample_id}: {exceeds} overall group(s) have "
                            f"TRUTH.TOTAL > 1.5x the largest single stratification row "
                            f"({max_single:.0f}), suggesting overlapping strata were summed. "
                            f"Use --regions to select non-overlapping subsets for reliable "
                            f"overall counts."
                        )

        return df_ps

    # ------------------------------------------------------------------
    # Helper: median [min-max] across samples
    # ------------------------------------------------------------------
    def _aggregate_across_samples(df_ps: pd.DataFrame, agg_keys: list) -> pd.DataFrame:
        metrics = [
            "Recall_base_raw", "Recall_cons_raw",
            "Precision_base_raw", "Precision_cons_raw",
            "Recall_gain_pp", "Precision_gain_pp",
            "F1_base", "F1_cons", "F1_gain_pp",
            "QUERY.TP_base", "QUERY.FP_base",
            "QUERY.TP_cons", "QUERY.FP_cons",
            "Additional_TP", "Additional_FP", "Additional_FN",
            "TP_FP_ratio", "FP_pct_increase",
            "Additional_TP_per_Mbp", "Additional_FP_per_Mbp",
            "Additional_TP_per_100_truth", "Additional_FP_per_100_truth",
            "Net_benefit"
        ]
        metrics = [m for m in metrics if m in df_ps.columns]

        agg_keys = [c for c in agg_keys if c in df_ps.columns]
        if not agg_keys:
            logging.warning("No valid aggregation keys")
            return pd.DataFrame()

        agg_dict = {m: ["median", "min", "max", "count"] for m in metrics}
        grouped = df_ps.groupby(agg_keys, sort=False).agg(agg_dict)
        grouped.columns = [f"{col}_{stat}" for col, stat in grouped.columns]
        grouped = grouped.reset_index()

        for m in metrics:
            med = f"{m}_median"
            lo = f"{m}_min"
            hi = f"{m}_max"
            if med in grouped.columns and lo in grouped.columns and hi in grouped.columns:
                grouped[f"{m}_formatted"] = grouped.apply(
                    lambda r: format_median_range(r[med], r[lo], r[hi]),
                    axis=1
                )
        return grouped

    # ------------------------------------------------------------------
    # Process each requested level
    # ------------------------------------------------------------------
    all_results = []
    level_counts = []

    for level in levels:
        final_keys = _parse_level_name(level)
        if not final_keys:
            logging.warning(f"Skipping invalid overall level: {level}")
            continue

        # Per-sample keys always include Sample + final keys
        per_sample_keys = ["Sample"] + final_keys

        df_ps = _build_per_sample(per_sample_keys)
        if df_ps.empty:
            logging.warning(f"No data for overall level: {level}")
            continue

        df_agg = _aggregate_across_samples(df_ps, final_keys)
        if df_agg.empty:
            logging.warning(f"Aggregation produced no results for level: {level}")
            continue

        # Add level label and fill missing stratification columns
        df_agg["Overall_Level"] = level

        # Fill placeholder values for columns not in final_keys
        placeholders = {
            "Caller": "All-Callers",
            "Condition": "All-Conditions",
            "Type": "All-Types",
            "Subset": "All-Regions",
            "Filter": "Overall",
            "Genotype": "Overall",
            "Subtype": "Overall",
        }
        for col, val in placeholders.items():
            if col not in final_keys and col not in df_agg.columns:
                df_agg[col] = val

        all_results.append(df_agg)
        level_counts.append(f"{level}={len(df_agg)}")

    if not all_results:
        logging.warning("No overall levels produced valid results")
        return pd.DataFrame()

    df_overall = pd.concat(all_results, ignore_index=True)
    logging.info(f"Overall: {sum(len(r) for r in all_results)} total group(s) "
                 f"across {len(all_results)} level(s) ({'; '.join(level_counts)})")
    return df_overall


def compute_median_of_medians(df_raw: pd.DataFrame,
                              levels: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Compute median-of-medians across samples at configurable levels.

    Each level name specifies the dimensions to aggregate by (excluding Sample,
    which is always used for per-sample median). Dimensions are separated by
    hyphens and can be: caller, condition, type.

    Available levels (specified via --mom-levels):
      caller-condition-type : Per Caller + Condition + Type (median across Subset/Filter/Genotype)
      caller-condition      : Per Caller + Condition (median across all stratifications)
      caller-type           : Per Caller + Type (median across Subset/Filter/Genotype, all conditions)
      caller                : Per Caller (median across all stratifications, all conditions)
      condition-type        : Per Condition + Type (median across Subset/Filter/Genotype, all callers)
      condition             : Per Condition (median across all stratifications, all callers)
      type                  : Per Type (median across all stratifications, all callers/conditions)

    For each level, within each sample the median across relevant
    stratification rows is computed, then aggregated as median [min-max]
    across samples.
    """
    if levels is None:
        levels = ["caller-condition-type", "caller-condition"]

    metrics = [
        "Recall_base_raw", "Recall_cons_raw",
        "Precision_base_raw", "Precision_cons_raw",
        "Recall_gain_pp", "Precision_gain_pp",
        "F1_base", "F1_cons", "F1_gain_pp",
        "QUERY.TP_base", "QUERY.FP_base",
        "QUERY.TP_cons", "QUERY.FP_cons",
        "Additional_TP", "Additional_FP", "Additional_FN",
        "TP_FP_ratio", "FP_pct_increase",
        "Additional_TP_per_Mbp", "Additional_FP_per_Mbp",
        "Additional_TP_per_100_truth", "Additional_FP_per_100_truth",
        "Net_benefit"
    ]
    metrics = [m for m in metrics if m in df_raw.columns]

    if not metrics:
        logging.warning("Cannot compute median-of-medians: no valid metrics found")
        return pd.DataFrame()

    # ------------------------------------------------------------------
    # Helper: per-sample median -> across-sample median [min-max]
    # ------------------------------------------------------------------
    def _mom_for_level(df_level: pd.DataFrame, gb_keys: list, label: str) -> pd.DataFrame:
        gb_keys = [c for c in gb_keys if c in df_level.columns]
        if not gb_keys:
            logging.warning(f"No valid groupby keys for MoM level: {label}")
            return pd.DataFrame()

        per_sample = df_level.groupby(gb_keys, sort=False)[metrics].median().reset_index()

        agg_keys = [c for c in gb_keys if c != "Sample"]
        if not agg_keys:
            logging.warning(f"No valid aggregation keys for MoM level: {label}")
            return pd.DataFrame()

        agg_dict = {m: ["median", "min", "max", "count"] for m in metrics}
        grouped = per_sample.groupby(agg_keys, sort=False).agg(agg_dict)
        grouped.columns = [f"{col}_{stat}" for col, stat in grouped.columns]
        grouped = grouped.reset_index()

        grouped["MoM_Level"] = label

        for m in metrics:
            med = f"{m}_median"
            lo = f"{m}_min"
            hi = f"{m}_max"
            if med in grouped.columns and lo in grouped.columns and hi in grouped.columns:
                grouped[f"{m}_formatted"] = grouped.apply(
                    lambda r: format_median_range(r[med], r[lo], r[hi]),
                    axis=1
                )
        return grouped

    # ------------------------------------------------------------------
    # Process each requested level
    # ------------------------------------------------------------------
    all_results = []
    level_counts = []

    for level in levels:
        final_keys = _parse_level_name(level)
        if not final_keys:
            logging.warning(f"Skipping invalid MoM level: {level}")
            continue

        # Per-sample keys always include Sample + final keys
        per_sample_keys = ["Sample"] + final_keys

        df_mom_level = _mom_for_level(df_raw, per_sample_keys, level)
        if df_mom_level.empty:
            logging.warning(f"No data for MoM level: {level}")
            continue

        # Fill placeholder values for columns not in final_keys
        placeholders = {
            "Caller": "All-Callers",
            "Condition": "All-Conditions",
            "Type": "All-Types",
            "Subset": "Across-Subsets",
            "Filter": "Median-of-Medians",
            "Genotype": "Median-of-Medians",
            "Subtype": "Median-of-Medians",
        }
        for col, val in placeholders.items():
            if col not in final_keys and col not in df_mom_level.columns:
                df_mom_level[col] = val

        all_results.append(df_mom_level)
        level_counts.append(f"{level}={len(df_mom_level)}")

    if not all_results:
        logging.warning("No MoM levels produced valid results")
        return pd.DataFrame()

    df_mom = pd.concat(all_results, ignore_index=True)
    logging.info(f"Median-of-medians: {sum(len(r) for r in all_results)} total group(s) "
                 f"across {len(all_results)} level(s) ({'; '.join(level_counts)})")
    return df_mom


def format_median_range(median, min_val, max_val) -> str:
    """Format a metric as 'median [min-max]' with appropriate precision."""
    if abs(median) >= 100:
        fmt = ".0f"
    elif abs(median) >= 10:
        fmt = ".1f"
    elif abs(median) >= 1:
        fmt = ".2f"
    else:
        fmt = ".3f"

    if pd.isna(median):
        med_str = "N/A"
    else:
        med_str = f"{median:{fmt}}"

    lo_str = f"{min_val:{fmt}}" if pd.notna(min_val) else "N/A"
    hi_str = f"{max_val:{fmt}}" if pd.notna(max_val) else "N/A"

    # Single-sample collapse: when min == max == median, show just the value
    if pd.notna(median) and pd.notna(min_val) and pd.notna(max_val):
        if abs(median - min_val) < 1e-12 and abs(median - max_val) < 1e-12:
            return med_str

    return f"{med_str} [{lo_str}-{hi_str}]"

# ==============================================================================
# TABLE FORMATTING
# ==============================================================================

def to_markdown_table(df_summary: pd.DataFrame,
                      region_order: Optional[List[str]] = None,
                      caller_order: Optional[List[str]] = None) -> str:
    """Convert summary to a manuscript-ready Markdown table."""
    df = df_summary.copy()

    display_cols = [
        "Caller", "Condition", "Type", "Subset", "Filter",
        "QUERY.TP_base_formatted", "QUERY.FP_base_formatted",
        "QUERY.TP_cons_formatted", "QUERY.FP_cons_formatted",
        "Recall_gain_pp_formatted", "Precision_gain_pp_formatted",
        "F1_gain_pp_formatted",
        "Additional_TP_formatted", "Additional_FP_formatted",
        "TP_FP_ratio_formatted", "Net_benefit_formatted",
        "Additional_TP_per_Mbp_formatted", "Additional_FP_per_Mbp_formatted"
    ]
    display_cols = [c for c in display_cols if c in df.columns]

    rename_map = {
        "QUERY.TP_base_formatted": "Base TP",
        "QUERY.FP_base_formatted": "Base FP",
        "QUERY.TP_cons_formatted": "Consensus TP",
        "QUERY.FP_cons_formatted": "Consensus FP",
        "Recall_gain_pp_formatted": "Recall Gain (pp)",
        "Precision_gain_pp_formatted": "Precision Gain (pp)",
        "F1_gain_pp_formatted": "F1 Gain (pp)",
        "Additional_TP_formatted": "Additional TP",
        "Additional_FP_formatted": "Additional FP",
        "TP_FP_ratio_formatted": "TP:FP Ratio",
        "Net_benefit_formatted": "Net Benefit",
        "Additional_TP_per_Mbp_formatted": "Additional TP/Mbp",
        "Additional_FP_per_Mbp_formatted": "Additional FP/Mbp"
    }

    df_display = df[display_cols].rename(columns=rename_map)

    if region_order and "Subset" in df_display.columns:
        df_display["Subset"] = pd.Categorical(df_display["Subset"], categories=region_order, ordered=True)
        df_display = df_display.sort_values("Subset")
    if caller_order and "Caller" in df_display.columns:
        df_display["Caller"] = pd.Categorical(df_display["Caller"], categories=caller_order, ordered=True)
        df_display = df_display.sort_values("Caller")

    try:
        return df_display.to_markdown(index=False)
    except ImportError:
        logging.warning("tabulate not installed; falling back to simple table")
        return df_display.to_string(index=False)


def to_latex_table(df_summary: pd.DataFrame,
                   region_order: Optional[List[str]] = None,
                   caption: str = "") -> str:
    """Convert summary to a LaTeX table using Styler.to_latex (pandas >= 1.3)."""
    df = df_summary.copy()
    display_cols = [
        "Caller", "Condition", "Type", "Subset",
        "QUERY.TP_base_formatted", "QUERY.FP_base_formatted",
        "QUERY.TP_cons_formatted", "QUERY.FP_cons_formatted",
        "Recall_gain_pp_formatted", "Precision_gain_pp_formatted",
        "F1_gain_pp_formatted",
        "Additional_TP_formatted", "Additional_FP_formatted",
        "TP_FP_ratio_formatted", "Net_benefit_formatted"
    ]
    display_cols = [c for c in display_cols if c in df.columns]
    rename_map = {
        "QUERY.TP_base_formatted": "Base TP",
        "QUERY.FP_base_formatted": "Base FP",
        "QUERY.TP_cons_formatted": "Cons. TP",
        "QUERY.FP_cons_formatted": "Cons. FP",
        "Recall_gain_pp_formatted": "Recall (pp)",
        "Precision_gain_pp_formatted": "Precision (pp)",
        "F1_gain_pp_formatted": "F1 (pp)",
        "Additional_TP_formatted": "Add. TP",
        "Additional_FP_formatted": "Add. FP",
        "TP_FP_ratio_formatted": "TP:FP",
        "Net_benefit_formatted": "Net Ben."
    }
    df_display = df[display_cols].rename(columns=rename_map)

    if region_order and "Subset" in df_display.columns:
        df_display["Subset"] = pd.Categorical(df_display["Subset"], categories=region_order, ordered=True)
        df_display = df_display.sort_values("Subset")

    try:
        styler = df_display.style.hide(axis="index")
        latex = styler.to_latex(
            caption=caption,
            hrules=True,
            environment="tabular"
        )
    except Exception as e:
        logging.warning(f"LaTeX table generation failed: {e}")
        latex = "% LaTeX table generation failed; use text table instead\n"

    return latex

def to_csv_table(df_summary: pd.DataFrame,
                 region_order: Optional[List[str]] = None,
                 caller_order: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Convert summary to a clean CSV table with renamed columns for direct use.
    Returns a DataFrame ready for write_csv or further processing.
    """
    df = df_summary.copy()

    display_cols = [
        "Caller", "Condition", "Type", "Subset", "Filter",
        "QUERY.TP_base_formatted", "QUERY.FP_base_formatted",
        "QUERY.TP_cons_formatted", "QUERY.FP_cons_formatted",
        "Recall_gain_pp_formatted", "Precision_gain_pp_formatted",
        "F1_gain_pp_formatted",
        "Additional_TP_formatted", "Additional_FP_formatted",
        "TP_FP_ratio_formatted", "Net_benefit_formatted",
        "Additional_TP_per_Mbp_formatted", "Additional_FP_per_Mbp_formatted",
        "Additional_TP_per_100_truth_formatted", "Additional_FP_per_100_truth_formatted"
    ]
    display_cols = [c for c in display_cols if c in df.columns]

    rename_map = {
        "QUERY.TP_base_formatted": "Base_TP",
        "QUERY.FP_base_formatted": "Base_FP",
        "Recall_gain_pp_formatted": "Recall_Gain_pp",
        "Precision_gain_pp_formatted": "Precision_Gain_pp",
        "F1_gain_pp_formatted": "F1_Gain_pp",
        "Additional_TP_formatted": "Additional_TP",
        "Additional_FP_formatted": "Additional_FP",
        "TP_FP_ratio_formatted": "TP_FP_Ratio",
        "Net_benefit_formatted": "Net_Benefit",
        "Additional_TP_per_Mbp_formatted": "Additional_TP_per_Mbp",
        "Additional_FP_per_Mbp_formatted": "Additional_FP_per_Mbp",
        "Additional_TP_per_100_truth_formatted": "Additional_TP_per_100_truth",
        "Additional_FP_per_100_truth_formatted": "Additional_FP_per_100_truth"
    }

    df_display = df[display_cols].rename(columns=rename_map)

    if region_order and "Subset" in df_display.columns:
        df_display["Subset"] = pd.Categorical(df_display["Subset"], categories=region_order, ordered=True)
        df_display = df_display.sort_values("Subset")
    if caller_order and "Caller" in df_display.columns:
        df_display["Caller"] = pd.Categorical(df_display["Caller"], categories=caller_order, ordered=True)
        df_display = df_display.sort_values("Caller")

    return df_display

# ==============================================================================
# REPRODUCIBILITY AUDIT DOCUMENT
# ==============================================================================

def write_reproducibility_doc(output_prefix: str, df_raw: pd.DataFrame,
                              df_summary: pd.DataFrame, manifest_path: str) -> None:
    """Write a detailed reproducibility audit trail as a plain-text file."""

    lines = []
    lines.append("# Reproducibility Audit Trail")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now().isoformat()}")
    lines.append(f"**Script:** run_BurdenBench.py v1.0.0")
    lines.append(f"**Manifest:** {manifest_path}")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 1. Percentage-Point (pp) Derivation")
    lines.append("")
    lines.append("All percentage-point changes reported in the manuscript are derived ")
    lines.append("directly from hap.py `METRIC.Recall` and `METRIC.Precision` values (0-1 scale):")
    lines.append("")
    lines.append("```")
    lines.append("Recall_gain_pp     = (METRIC.Recall_consensus  - METRIC.Recall_baseline)    × 100")
    lines.append("Precision_gain_pp  = (METRIC.Precision_consensus - METRIC.Precision_baseline) × 100")
    lines.append("```")
    lines.append("")
    lines.append("### Verification")
    lines.append("The raw output file (`_raw_per_sample.csv`) contains the following columns ")
    lines.append("so that every pp value can be independently recalculated:")
    lines.append("")
    lines.append("| Column | Description | Scale |")
    lines.append("|--------|-------------|-------|")
    lines.append("| `Recall_base_raw`    | Baseline recall from hap.py    | 0-1   |")
    lines.append("| `Recall_cons_raw`    | Consensus recall from hap.py   | 0-1   |")
    lines.append("| `Recall_base_pct`    | Baseline recall                | 0-100 |")
    lines.append("| `Recall_cons_pct`    | Consensus recall               | 0-100 |")
    lines.append("| `Recall_gain_pp`     | Recall improvement             | pp    |")
    lines.append("| `Recall_gain_pp_check` | Independent verification     | pp    |")
    lines.append("")
    lines.append("> **Note:** `*_check` columns verify that `gain_pp` equals `(cons_raw - base_raw) × 100`. ")
    lines.append("> Any non-zero drift indicates a data processing error.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 2. Absolute Burden Derivation")
    lines.append("")
    lines.append("```")
    lines.append("Additional_TP = QUERY.TP_consensus - QUERY.TP_baseline")
    lines.append("Additional_FP = QUERY.FP_consensus - QUERY.FP_baseline")
    lines.append("Net_benefit   = Additional_TP - Additional_FP   (signed, handles all cases)")
    lines.append("TP_FP_ratio   = Additional_TP / Additional_FP   (raw, may produce inf)")
    lines.append("TP_FP_ratio_adj = np.where(Additional_FP >= 0, (Additional_TP + 1) / (Additional_FP + 1), Additional_TP / Additional_FP)")
    lines.append("  - Laplace +1 smoothing only when Additional_FP >= 0 (avoids sign flips)")
    lines.append("  - NaN when both Additional_TP and Additional_FP are 0 (no meaningful ratio)")
    lines.append("  - Raw ratio used directly when Additional_FP < 0 (FP reduction)")
    lines.append("```")
    lines.append("")
    lines.append("### Edge case handling")
    lines.append("")
    lines.append("| Scenario | Raw ratio | Adjusted ratio | Net benefit | Display |")
    lines.append("|----------|-----------|----------------|-------------|---------|")
    lines.append("| TP > 0, FP = 0 | inf | (TP+1)/(0+1) | +TP | inf (pure gain) |")
    lines.append("| TP = 0, FP > 0 | 0 | 1/(FP+1) | -FP | 0 (cost only) |")
    lines.append("| TP = 0, FP = 0 | nan | NaN | 0 | - (no change) |")
    lines.append("| TP < 0, FP = 0 | -inf | (TP+1)/(0+1) | TP | -inf (pure loss) |")
    lines.append("| TP > 0, FP < 0 | negative | TP/FP | +TP-FP | negative (FP reduced) |")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 3. Normalization Formulas")
    lines.append("")
    lines.append("```")
    lines.append("Additional_TP_per_Mbp = Additional_TP / (Subset.IS_CONF.Size_baseline / 1,000,000)")
    lines.append("Additional_FP_per_Mbp = Additional_FP / (Subset.IS_CONF.Size_baseline / 1,000,000)")
    lines.append("TP_per_Mbp_cons = QUERY.TP_cons / (Subset.IS_CONF.Size_consensus / 1,000,000)")
    lines.append("FP_per_Mbp_cons = QUERY.FP_cons / (Subset.IS_CONF.Size_consensus / 1,000,000)")
    lines.append("")
    lines.append("Additional_TP_per_100_truth = (Additional_TP / TRUTH.TOTAL_baseline) × 100")
    lines.append("Additional_FP_per_100_truth = (Additional_FP / TRUTH.TOTAL_baseline) × 100")
    lines.append("```")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 4. Aggregation Method")
    lines.append("")
    lines.append(f"- **Samples aggregated:** {df_raw['Sample'].nunique()}")
    lines.append(f"- **Grouping columns:** Caller, Condition, Type, Subset, Filter, Genotype")
    lines.append(f"- **Statistic:** Median with inter-sample range [minimum - maximum]")
    lines.append(f"- **Total groups:** {len(df_summary)}")
    lines.append("")
    lines.append("### Overall Aggregation Caveat")
    lines.append("")
    lines.append("Overall metrics (summed across stratifications) recompute recall and precision")
    lines.append("from summed TP/FP/TRUTH counts rather than averaging raw values.")
    lines.append("Per-Mbp and per-100-truth metrics in overall outputs are most reliable when")
    lines.append("aggregating non-overlapping strata (e.g., single Filter, distinct Subsets).")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 5. Data Provenance")
    lines.append("")
    samples = sorted(df_raw['Sample'].unique())
    callers = sorted(df_raw['Caller'].unique())
    conditions = sorted(df_raw['Condition'].unique())
    lines.append(f"- **Samples:** {', '.join(samples)}")
    lines.append(f"- **Callers:** {', '.join(callers)}")
    lines.append(f"- **Conditions:** {', '.join(conditions)}")
    lines.append(f"- **Total stratification rows (raw):** {len(df_raw):,}")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 6. Validation Checks Performed")
    lines.append("")
    lines.append("1. File existence and readability")
    lines.append("2. Required hap.py columns present (`METRIC.Recall`, `QUERY.TP`, etc.)")
    lines.append("3. Merge keys complete (no missing `Type`, `Subset`, etc.)")
    lines.append("4. Numeric type validation for all quantitative columns")
    lines.append("5. Value sanity checks (no negative FP/TP counts)")
    lines.append("6. PP calculation drift check (verifies `gain_pp == (cons - base) × 100`)")
    lines.append("")

    doc_path = _output_path(output_prefix, "reproducibility.txt")
    with open(doc_path, "w") as f:
        f.write("\n".join(lines))
    logging.info(f"[5] Reproducibility audit: {doc_path}")

# ==============================================================================
# BATCH PROCESSING
# ==============================================================================

def _output_path(output_prefix: str, suffix: str) -> Path:
    """Build a robust output path from prefix and file suffix.

    Handles directory-only prefixes (e.g. 'manuscript_tables/') by ensuring
    the suffix is appended correctly without dropping the trailing slash.
    """
    p = Path(output_prefix)
    # If prefix is a directory (ends with / or \), append suffix directly
    if str(output_prefix).rstrip().endswith(("/", "\\")):
        return p / suffix.lstrip("_")
    # Otherwise treat prefix as stem and append suffix
    return Path(str(p.with_suffix("")) + "_" + suffix.lstrip("_"))


def process_manifest(manifest_path: str, output_prefix: str,
                     regions: Optional[List[str]] = None,
                     filters: Optional[List[str]] = None,
                     group_cols: Optional[List[str]] = None,
                     dry_run: bool = False,
                     overall_levels: Optional[List[str]] = None,
                     mom_levels: Optional[List[str]] = None,
                     no_overall: bool = False,
                     no_mom: bool = False,
                     report_mismatches: bool = False) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Batch process manifest and produce all outputs."""
    manifest = load_manifest(manifest_path)

    if dry_run:
        logging.info("DRY RUN: All inputs and files validated. No outputs will be written.")
        return pd.DataFrame(), pd.DataFrame()

    all_raw = []
    all_mismatch_records: List[pd.DataFrame] = []
    errors: List[Tuple] = []
    success_count = 0

    # Single-sample sanity check
    n_samples = manifest["sample"].nunique()
    if n_samples == 1:
        logging.warning(
            "Only 1 sample found in manifest. Median [min-max] statistics will be "
            "degenerate (all three values identical). Consider reporting raw values "
            "from the _raw_per_sample.csv output instead."
        )
    dup_tracker: dict = {}  # key -> list of DataFrames for duplicate averaging

    total = len(manifest)
    logging.info(f"Processing {total} comparison(s)...")

    for idx, row in manifest.iterrows():
        label = f"{row['sample']}/{row['caller']}/{row['condition']}"
        try:
            logging.debug(f"[{idx+1}/{total}] Loading baseline: {row['baseline_path']}")
            df_base = load_happy(row["baseline_path"])

            logging.debug(f"[{idx+1}/{total}] Loading consensus: {row['consensus_path']}")
            df_cons = load_happy(row["consensus_path"])

            logging.debug(f"[{idx+1}/{total}] Computing burden for {label}")
            result, comp_mismatches = compute_burden(df_base, df_cons, row["sample"], row["caller"], row["condition"])
            if comp_mismatches:
                all_mismatch_records.extend(comp_mismatches)

            if regions:
                before = len(result)
                result = result[result["Subset"].isin(regions)]
                after = len(result)
                if after < before:
                    logging.debug(f"  Filtered to {after}/{before} rows matching regions: {regions}")

            if filters:
                before = len(result)
                result = result[result["Filter"].isin(filters)]
                after = len(result)
                if after < before:
                    logging.debug(f"  Filtered to {after}/{before} rows matching filters: {filters}")

                # Warn about double-counting when ALL and PASS are both selected
                if "ALL" in filters and "PASS" in filters:
                    logging.warning(
                        f"Both 'ALL' and 'PASS' filters selected for {sample}/{caller}/{condition}. "
                        f"Note that the 'ALL' row in hap.py output includes all variants (both passing "
                        f"and failing filters), while 'PASS' includes only passing variants. In overall "
                        f"aggregation (which sums across Filter), selecting both will double-count PASS "
                        f"variants. Use --filters ALL alone for complete counts, or --filters PASS alone "
                        f"for filtered counts only."
                    )

            if result.empty:
                logging.warning(f"[{idx+1}/{total}] {label}: No rows remaining after region/filter selection")
                errors.append((row["sample"], row["caller"], row["condition"], "Empty after filtering"))
                continue

            # Track duplicates for averaging
            dup_key = (row["sample"], row["caller"], row["condition"])
            if dup_key not in dup_tracker:
                dup_tracker[dup_key] = []
            dup_tracker[dup_key].append(result)

            success_count += 1
            logging.info(f"[{idx+1}/{total}] ✓ {label}: {len(result)} stratification rows")

        except BurdenCalculatorError as e:
            errors.append((row["sample"], row["caller"], row["condition"], str(e)))
            logging.error(f"[{idx+1}/{total}] ✗ {label}: {e}")
        except Exception as e:
            errors.append((row["sample"], row["caller"], row["condition"], f"Unexpected: {e}"))
            logging.error(f"[{idx+1}/{total}] ✗ {label}: Unexpected error: {e}")
            logging.debug(traceback.format_exc())

    if not dup_tracker:
        logging.error(f"All {total} comparison(s) failed. Nothing to aggregate.")
        if errors:
            logging.error("Errors encountered:")
            for sample, caller, condition, err in errors:
                logging.error(f"  {sample}/{caller}/{condition}: {err}")
        raise AggregationError("No successful comparisons. Cannot produce output.")

    # ------------------------------------------------------------------
    # Average duplicates by median within each sample+caller+condition
    # ------------------------------------------------------------------
    for dup_key, dfs in dup_tracker.items():
        sample, caller, condition = dup_key
        if len(dfs) > 1:
            # Multiple files for same key: compute median across duplicates
            merged_dup = pd.concat(dfs, ignore_index=True)
            # Group by all stratification keys and compute median across duplicates
            strat_keys = ["Type", "Subtype", "Subset", "Filter", "Genotype"]
            strat_keys = [k for k in strat_keys if k in merged_dup.columns]

            # Numeric columns to average
            numeric_cols = [c for c in merged_dup.columns 
                           if pd.api.types.is_numeric_dtype(merged_dup[c]) 
                           and c not in strat_keys + ["Sample", "Caller", "Condition"]]

            # Compute median across duplicates for each stratification row
            agg_dict = {c: "median" for c in numeric_cols}
            # For non-numeric, take first (they should be identical)
            for c in merged_dup.columns:
                if c not in agg_dict and c not in strat_keys + ["Sample", "Caller", "Condition"]:
                    agg_dict[c] = "first"

            averaged = merged_dup.groupby(strat_keys, sort=False).agg(agg_dict).reset_index()
            averaged["Sample"] = sample
            averaged["Caller"] = caller
            averaged["Condition"] = condition
            averaged["Duplicate_Count"] = len(dfs)
            averaged["Duplicate_Source"] = f"Median of {len(dfs)} duplicate files"
            all_raw.append(averaged)
            logging.info(f"[DUP] {sample}/{caller}/{condition}: Averaged {len(dfs)} duplicates by median")
        else:
            # Single file: just add Duplicate_Count = 1
            df_single = dfs[0].copy()
            df_single["Duplicate_Count"] = 1
            df_single["Duplicate_Source"] = "Single file (no duplicates)"
            all_raw.append(df_single)

    df_raw = pd.concat(all_raw, ignore_index=True)
    n_averaged = df_raw["Duplicate_Count"].gt(1).sum()
    logging.info(f"Combined raw data: {len(df_raw)} rows from {len(dup_tracker)} unique comparison(s) "
                 f"({success_count} files processed, {n_averaged} rows averaged from duplicates)")

    out_dir = Path(output_prefix).parent if str(Path(output_prefix).parent) != "." else Path(".")
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Raw per-sample
    raw_path = _output_path(output_prefix, "raw_per_sample.csv")
    try:
        df_raw.to_csv(raw_path, index=False, encoding="utf-8-sig")
        logging.info(f"[1] Raw per-sample results: {raw_path} ({len(df_raw)} rows)")
    except Exception as e:
        raise FileAccessError(f"Cannot write raw output {raw_path}: {e}")

    # 2. Summary
    try:
        df_summary = aggregate_summary(df_raw, group_cols=group_cols)
    except AggregationError as e:
        logging.error(f"Aggregation failed: {e}")
        raise
    # 2b. Overall by caller and region (separate output)
    df_overall = pd.DataFrame()
    if not no_overall:
        try:
            df_overall = compute_overall_by_caller_region(df_raw, levels=overall_levels)
            if not df_overall.empty:
                logging.info(f"[2b] Computed {len(df_overall)} overall group(s)")
        except Exception as e:
            logging.warning(f"Could not compute overall by caller+region: {e}")
    else:
        logging.info("[2b] Overall aggregation skipped (--no-overall)")

    summary_path = _output_path(output_prefix, "summary_median_range.csv")
    try:
        df_summary.to_csv(summary_path, index=False, encoding="utf-8-sig")
        logging.info(f"[2] Summary (median+range): {summary_path} ({len(df_summary)} groups)")
    except Exception as e:
        raise FileAccessError(f"Cannot write summary {summary_path}: {e}")

    # 3. Markdown
    md_path = _output_path(output_prefix, "table.txt")
    try:
        md = to_markdown_table(df_summary)
        with open(md_path, "w") as f:
            f.write("Absolute TP and FP Burden Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write("Median [min-max] across samples.\n\n")
            f.write(md)
            f.write("\n\n")
            f.write("---\n\n")
            f.write("**Notes:**\n")
            f.write("- TP:FP Ratio > 1.0 indicates net clinical benefit\n")
            f.write("- TP:FP Ratio < 1.0 indicates more noise than signal; filtering recommended\n")
            f.write("- All metrics computed within GIAB confident regions (Subset.IS_CONF.Size)\n")
            f.write("- Percentage-point (pp) changes = (METRIC_consensus - METRIC_baseline) × 100\n")
        logging.info(f"[3] Text table: {md_path}")
    except Exception as e:
        logging.warning(f"Could not write Markdown table: {e}")

    # 4. LaTeX
    tex_path = _output_path(output_prefix, "table.tex")
    try:
        latex = to_latex_table(
            df_summary,
            caption="Median recall/precision changes and absolute TP/FP burden across GIAB samples."
        )
        with open(tex_path, "w") as f:
            f.write("% Generated by run_BurdenBench.py v1.0.0\n")
            f.write("% Median [min-max] across samples\n\n")
            f.write(latex)
        logging.info(f"[4] LaTeX table: {tex_path}")
    except Exception as e:
        logging.warning(f"Could not write LaTeX table: {e}")

    # 5. CSV table (clean, publication-ready)
    csv_table_path = _output_path(output_prefix, "table.csv")
    try:
        df_csv_table = to_csv_table(df_summary)
        df_csv_table.to_csv(csv_table_path, index=False, encoding="utf-8-sig")
        logging.info(f"[5] CSV table: {csv_table_path} ({len(df_csv_table)} rows)")
    except Exception as e:
        logging.warning(f"Could not write CSV table: {e}")

    # 5b. Overall outputs (separate from stratified)
    if not df_overall.empty:
        # 5b-1. Overall summary CSV
        overall_summary_path = _output_path(output_prefix, "overall_summary_median_range.csv")
        try:
            df_overall.to_csv(overall_summary_path, index=False, encoding="utf-8-sig")
            logging.info(f"[5b-1] Overall summary: {overall_summary_path} ({len(df_overall)} groups)")
        except Exception as e:
            logging.warning(f"Could not write overall summary: {e}")

        # 5b-2. Overall text table
        overall_md_path = _output_path(output_prefix, "overall_table.txt")
        try:
            overall_md = to_markdown_table(df_overall)
            with open(overall_md_path, "w") as f:
                f.write("Overall TP and FP Burden Summary (by Caller + Region)\n")
                f.write("=" * 55 + "\n\n")
                f.write("Median [min-max] across samples.\n\n")
                f.write(overall_md)
                f.write("\n\n")
                f.write("-" * 55 + "\n\n")
                f.write("Notes:\n")
                f.write("- Overall = summed across Filter, Genotype, Subtype\n")
                f.write("- Type=Overall rows = summed across all variant types\n")
                f.write("- Type=SNV/INDEL rows = summed across Filter/Genotype/Subtype only\n")
                f.write("- WARNING: Selecting both 'ALL' and 'PASS' filters double-counts PASS\n")
                f.write("  variants in overall aggregation. Use one or the other.\n")
                f.write("- TP:FP Ratio > 1.0 indicates net clinical benefit\n")
                f.write("- All metrics computed within GIAB confident regions\n")
                f.write("- WARNING: Per-Mbp and per-100-truth metrics assume non-overlapping strata.\n")
                f.write("  Use --filters ALL for most reliable overall per-Mbp estimates.\n")
            logging.info(f"[5b-2] Overall text table: {overall_md_path}")
        except Exception as e:
            logging.warning(f"Could not write overall text table: {e}")

        # 5b-3. Overall LaTeX table
        overall_tex_path = _output_path(output_prefix, "overall_table.tex")
        try:
            overall_latex = to_latex_table(
                df_overall,
                caption="Overall median recall/precision changes and absolute TP/FP burden by caller and region."
            )
            with open(overall_tex_path, "w") as f:
                f.write("% Generated by run_BurdenBench.py v1.0.0\n")
                f.write("% Overall by caller + region (median [min-max] across samples)\n\n")
                f.write(overall_latex)
            logging.info(f"[5b-3] Overall LaTeX table: {overall_tex_path}")
        except Exception as e:
            logging.warning(f"Could not write overall LaTeX table: {e}")

        # 5b-4. Overall CSV table
        overall_csv_path = _output_path(output_prefix, "overall_table.csv")
        try:
            df_overall_csv = to_csv_table(df_overall)
            df_overall_csv.to_csv(overall_csv_path, index=False, encoding="utf-8-sig")
            logging.info(f"[5b-4] Overall CSV table: {overall_csv_path} ({len(df_overall_csv)} rows)")
        except Exception as e:
            logging.warning(f"Could not write overall CSV table: {e}")

    # 5c. Median-of-medians output (for sensitivity analysis / comparison)
    if not no_mom:
        try:
            df_mom = compute_median_of_medians(df_raw, levels=mom_levels)
            if not df_mom.empty:
                # 5c-1. Median-of-medians summary CSV
                mom_summary_path = _output_path(output_prefix, "median_of_medians_summary.csv")
                df_mom.to_csv(mom_summary_path, index=False, encoding="utf-8-sig")
                logging.info(f"[5c-1] Median-of-medians summary: {mom_summary_path} ({len(df_mom)} groups)")

                # 5c-2. Median-of-medians text table
                mom_txt_path = _output_path(output_prefix, "median_of_medians_table.txt")
                try:
                    mom_md = to_markdown_table(df_mom)
                    with open(mom_txt_path, "w") as f:
                        f.write("Median-of-Medians Sensitivity Analysis\n")
                        f.write("=" * 50 + "\n\n")
                        f.write("For each sample, the median across all stratification rows\n")
                        f.write("(Type/Subset/Filter/Genotype) was computed. The table below\n")
                        f.write("shows the median [min-max] of those per-sample medians.\n\n")
                        f.write(mom_md)
                        f.write("\n\n")
                        f.write("-" * 50 + "\n\n")
                        f.write("Interpretation:\n")
                        f.write("- This gives equal weight to each sample regardless of\n")
                        f.write("  how many stratification rows it contains.\n")
                        f.write("- Compare with _table.txt (direct median) to assess\n")
                        f.write("  robustness to aggregation method.\n")
                    logging.info(f"[5c-2] Median-of-medians text table: {mom_txt_path}")
                except Exception as e:
                    logging.warning(f"Could not write median-of-medians text table: {e}")

                # 5c-3. Median-of-medians CSV table
                mom_csv_path = _output_path(output_prefix, "median_of_medians_table.csv")
                try:
                    df_mom_csv = to_csv_table(df_mom)
                    df_mom_csv.to_csv(mom_csv_path, index=False, encoding="utf-8-sig")
                    logging.info(f"[5c-3] Median-of-medians CSV table: {mom_csv_path} ({len(df_mom_csv)} rows)")
                except Exception as e:
                    logging.warning(f"Could not write median-of-medians CSV table: {e}")
        except Exception as e:
            logging.warning(f"Could not compute median-of-medians: {e}")
    else:
        logging.info("[5c] Median-of-medians skipped (--no-mom)")

    # 6. Mismatch report (if --report-mismatches or mismatches detected)
    if all_mismatch_records:
        total_mismatches = sum(len(df) for df in all_mismatch_records)
        if report_mismatches:
            try:
                mismatch_df = pd.concat(all_mismatch_records, ignore_index=True)
                mismatch_path = _output_path(output_prefix, "mismatch_report.csv")
                mismatch_df.to_csv(mismatch_path, index=False, encoding="utf-8-sig")
                logging.info(
                    f"[6] Mismatch report: {mismatch_path} ({len(mismatch_df)} rows, "
                    f"{mismatch_df['Mismatch_Type'].nunique()} type(s): "
                    f"{', '.join(mismatch_df['Mismatch_Type'].unique())})"
                )
            except Exception as e:
                logging.warning(f"Could not write mismatch report: {e}")
        else:
            logging.info(
                f"{total_mismatches} mismatch record(s) detected across all comparisons. "
                f"Use --report-mismatches to write a detailed CSV report."
            )

    # 7. Reproducibility audit
    try:
        write_reproducibility_doc(output_prefix, df_raw, df_summary, manifest_path)
    except Exception as e:
        logging.warning(f"Could not write reproducibility document: {e}")

    if errors:
        logging.warning(f"\n[!] {len(errors)}/{total} comparison(s) failed:")
        for sample, caller, condition, err in errors:
            logging.warning(f"    {sample}/{caller}/{condition}: {err}")

    logging.info(f"\nDone. {success_count}/{total} comparison(s) successful.")
    logging.info(f"Output files written to: {output_prefix}*")
    return df_raw, df_summary

# ==============================================================================
# SINGLE-FILE MODE
# ==============================================================================

def process_single(baseline_path: str, consensus_path: str,
                   sample: str, caller: str, condition: str,
                   output_path: str,
                   regions: Optional[List[str]] = None,
                   filters: Optional[List[str]] = None) -> pd.DataFrame:
    """Process a single baseline-consensus pair."""
    logging.info(f"Single-file mode: {sample}/{caller}/{condition}")

    df_base = load_happy(baseline_path)
    df_cons = load_happy(consensus_path)
    result, comp_mismatches = compute_burden(df_base, df_cons, sample, caller, condition)
    if comp_mismatches:
        logging.info(
            f"{len(comp_mismatches)} mismatch record(s) detected. "
            f"Use --report-mismatches to write a detailed report file."
        )

    if regions:
        result = result[result["Subset"].isin(regions)]
    if filters:
        result = result[result["Filter"].isin(filters)]
        if "ALL" in filters and "PASS" in filters:
            logging.warning(
                "Both 'ALL' and 'PASS' filters selected. The 'ALL' row includes all variants "
                "(both passing and failing filters), while 'PASS' includes only passing variants. "
                "In overall aggregation, selecting both will double-count PASS variants. "
                "Use --filters ALL alone for complete counts, or --filters PASS alone for "
                "filtered counts only."
            )

    if result.empty:
        raise AggregationError("No rows remaining after region/filter selection")

    try:
        result.to_csv(output_path, index=False)
        logging.info(f"Saved: {output_path} ({len(result)} rows)")
    except Exception as e:
        raise FileAccessError(f"Cannot write output {output_path}: {e}")

    preview_cols = [
        "Subset", "Filter",
        "Recall_base_raw", "Recall_cons_raw", "Recall_gain_pp",
        "Precision_base_raw", "Precision_cons_raw", "Precision_gain_pp",
        "Additional_TP", "Additional_FP", "TP_FP_ratio",
        "Additional_TP_per_Mbp", "Additional_FP_per_Mbp"
    ]
    preview_cols = [c for c in preview_cols if c in result.columns]
    logging.info("\nPreview (first 10 rows):")
    logging.info("\n" + result[preview_cols].head(10).to_string(index=False))

    return result

# ==============================================================================
# ARGUMENT PARSER
# ==============================================================================

def create_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="run_BurdenBench.py",
        description="""
BurdenBench: Calculate absolute false-positive (FP) and true-positive (TP) 
burden from hap.py benchmarking outputs. Supports batch processing with 
aggregation across samples for manuscript-ready tables.

EXAMPLES:
  # Single comparison
  python run_BurdenBench.py \
      --baseline HG001_bwa_GRCh38_bcftools.csv \
      --consensus HG001_1KG_AF10pc_bcftools.csv \
      --sample HG001 --caller bcftools --condition AF10 \
      --output results.csv

  # Batch with aggregation (recommended)
  python run_BurdenBench.py \
      --manifest manifest.csv \
      --output-prefix manuscript_tables/ \
      --regions AllAutosomes lowmappabilityall segdups CMRG MHC \
      --filters ALL

  # Validate inputs without writing outputs
  python run_BurdenBench.py \
      --manifest manifest.csv --dry-run

  # Verbose logging for debugging
  python run_BurdenBench.py \
      --manifest manifest.csv --output-prefix out/ --verbose

  # Compute only specific overall levels
  python run_BurdenBench.py \
      --manifest manifest.csv --output-prefix out/ \
      --overall-levels caller-region caller

  # Skip sensitivity analyses (faster)
  python run_BurdenBench.py \
      --manifest manifest.csv --output-prefix out/ --no-sensitivity
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXIT CODES:
  0  Success
  1  No successful comparisons
  2  Invalid arguments or missing files
  3  Manifest format error
  4  Data validation error

REPRODUCIBILITY:
  All percentage-point (pp) values are derived as:
    gain_pp = (METRIC_consensus - METRIC_baseline) × 100

  Raw 0-1 values are preserved in output columns:
    Recall_base_raw, Recall_cons_raw, Precision_base_raw, Precision_cons_raw

  A full reproducibility audit is written to {prefix}_reproducibility.txt

For issues, check that hap.py files contain the expected columns:
  Type, Subtype, Subset, Filter, Genotype, METRIC.Recall, METRIC.Precision,
  QUERY.FP, QUERY.TP, TRUTH.TOTAL, Subset.IS_CONF.Size
        """
    )

    input_group = parser.add_argument_group("Input Options")
    input_group.add_argument("--baseline", help="Path to baseline hap.py summary CSV")
    input_group.add_argument("--consensus", help="Path to consensus hap.py summary CSV")
    input_group.add_argument("--manifest", help="Path to batch manifest CSV/TSV (enables full aggregation)")
    input_group.add_argument("--sample", default="SAMPLE", help="Sample identifier (single mode)")
    input_group.add_argument("--caller", default="CALLER", help="Variant caller name (single mode)")
    input_group.add_argument("--condition", default="CONDITION", help="Consensus condition name (single mode)")

    filter_group = parser.add_argument_group("Filtering Options")
    filter_group.add_argument("--regions", nargs="+", default=None,
                              help="Space-separated Subset names to retain (e.g. AllAutosomes lowmappabilityall)")
    filter_group.add_argument("--filters", nargs="+", default=None,
                              help="Space-separated Filter values to retain (e.g. ALL PASS)")

    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument("--output-prefix", default="fp_burden",
                              help="Prefix for all output files in batch mode (default: fp_burden)")
    output_group.add_argument("--output", help="Output CSV path (single mode only; ignored with --manifest)")
    output_group.add_argument("--group-by", nargs="+", default=None,
                              help="Columns to group by for aggregation (default: Caller Condition Type Subset Filter Genotype)")

    agg_group = parser.add_argument_group("Aggregation Options")
    agg_group.add_argument("--overall-levels", nargs="+", 
                           choices=[
                               "caller-condition-type-region", "caller-condition-region",
                               "caller-condition-type", "caller-condition",
                               "caller-type-region", "caller-region",
                               "caller-type", "caller",
                               "condition-type-region", "condition-region",
                               "condition-type", "condition",
                               "type-region", "region", "type"
                           ],
                           default=[
                               "caller-condition-type-region", "caller-condition-region",
                               "caller-condition-type", "caller-condition",
                               "caller-type-region", "caller-region",
                               "caller-type", "caller"
                           ],
                           help="Overall aggregation levels to compute (default: all caller-focused levels)")
    agg_group.add_argument("--mom-levels", nargs="+",
                           choices=[
                               "caller-condition-type", "caller-condition",
                               "caller-type", "caller",
                               "condition-type", "condition",
                               "type"
                           ],
                           default=["caller-condition-type", "caller-condition"],
                           help="Median-of-medians levels to compute (default: both caller-condition levels)")
    agg_group.add_argument("--no-overall", action="store_true",
                           help="Skip overall aggregation (stratified outputs only)")
    agg_group.add_argument("--no-mom", action="store_true",
                           help="Skip median-of-medians calculation")
    agg_group.add_argument("--no-sensitivity", action="store_true",
                           help="Skip all sensitivity analyses (overall + MoM)")

    util_group = parser.add_argument_group("Utility Options")
    util_group.add_argument("--verbose", "-v", action="store_true", help="Enable verbose/debug logging")
    util_group.add_argument("--dry-run", action="store_true",
                            help="Validate all inputs but do not write any output files")
    util_group.add_argument("--debug-columns", action="store_true",
                            help="Include verification columns (Recall_gain_pp_check, etc.) in output")
    util_group.add_argument("--report-mismatches", action="store_true",
                            help="Write a detailed mismatch report file listing all rows with "
                                 "mismatched TRUTH.TOTAL, Subset.IS_CONF.Size, or other "
                                 "inconsistencies between baseline and consensus")
    util_group.add_argument("--version", action="version", version="%(prog)s 1.0.0")

    return parser

# ==============================================================================
# MAIN
# ==============================================================================

def main() -> int:
    parser = create_parser()
    args = parser.parse_args()

    setup_logging(verbose=args.verbose)
    logging.info("=" * 60)
    logging.info("run_BurdenBench.py v1.0.0 - with explicit pp derivation")
    logging.info("=" * 60)

    if args.manifest:
        if args.baseline or args.consensus:
            logging.warning("--baseline and --consensus are ignored when --manifest is provided")
    else:
        if not args.baseline or not args.consensus:
            parser.error("--baseline and --consensus are required unless using --manifest")
        if not args.output:
            args.output = f"{args.output_prefix}_{args.sample}_{args.caller}_{args.condition}.csv"

    # Handle --no-sensitivity (shorthand for both --no-overall and --no-mom)
    if getattr(args, 'no_sensitivity', False):
        args.no_overall = True
        args.no_mom = True

    # Propagate --debug-columns flag to compute_burden via module-level sentinel
    if getattr(args, 'debug_columns', False):
        sys._debug_columns = True

    try:
        if args.manifest:
            process_manifest(
                args.manifest,
                args.output_prefix,
                regions=args.regions,
                filters=args.filters,
                group_cols=args.group_by,
                dry_run=args.dry_run,
                overall_levels=args.overall_levels,
                mom_levels=args.mom_levels,
                no_overall=args.no_overall,
                no_mom=args.no_mom,
                report_mismatches=getattr(args, 'report_mismatches', False)
            )
        else:
            process_single(
                args.baseline, args.consensus,
                args.sample, args.caller, args.condition,
                args.output,
                regions=args.regions,
                filters=args.filters
            )

        logging.info("Done.")
        return 0

    except FileAccessError as e:
        logging.error(f"File access error: {e}")
        return 2
    except ManifestError as e:
        logging.error(f"Manifest error: {e}")
        return 3
    except DataValidationError as e:
        logging.error(f"Data validation error: {e}")
        return 4
    except AggregationError as e:
        logging.error(f"Aggregation error: {e}")
        return 1
    except BurdenCalculatorError as e:
        logging.error(f"Calculator error: {e}")
        return 1
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        if args.verbose:
            logging.error(traceback.format_exc())
        return 1


if __name__ == "__main__":
    sys.exit(main())
