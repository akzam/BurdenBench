# Usage Manual — calculate_fp_burden.py

**Version:** 5.1  
**Purpose:** Calculate absolute false-positive (FP) and true-positive (TP) burden from hap.py benchmarking outputs, with explicit percentage-point (pp) derivation, Laplace-adjusted ratios, and full reproducibility audit.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Installation](#installation)
3. [Input Requirements](#input-requirements)
4. [Running the Script](#running-the-script)
5. [Output Files](#output-files)
6. [Interpreting Results](#interpreting-results)
7. [Troubleshooting](#troubleshooting)
8. [Examples](#examples)

---

## Quick Start

### Single comparison (one sample, one caller, one condition)

```bash
python calculate_fp_burden.py \
    --baseline HG001_bwa_GRCh38_bcftools.csv \
    --consensus HG001_1KG_AF10pc_bcftools.csv \
    --sample HG001 \
    --caller bcftools \
    --condition 1KG_AF10pc \
    --output results.csv
```

### Batch processing with aggregation (recommended for manuscripts)

```bash
python calculate_fp_burden.py \
    --manifest manifest.csv \
    --output-prefix manuscript_tables/ \
    --regions AllAutosomes lowmappabilityall segdups CMRG MHC \
    --filters ALL
```

---

## Installation

### Requirements

- Python >= 3.8
- pandas >= 1.3.0
- numpy >= 1.20.0

### Install dependencies

```bash
pip install pandas numpy
```

### Download

```bash
git clone https://github.com/YOUR_USERNAME/IUPAC-consensus-reference.git
cd IUPAC-consensus-reference/scripts
```

---

## Input Requirements

### hap.py summary files

The script consumes **hap.py summary CSV/TSV files** (the standard output of `hap.py` benchmarking). Each file must contain these columns:

| Column | Description | Required |
|--------|-------------|----------|
| `Type` | Variant type (SNV, INDEL, etc.) | ✅ |
| `Subtype` | Variant subtype | ✅ |
| `Subset` | Stratification region name | ✅ |
| `Filter` | Filter status (ALL, PASS, etc.) | ✅ |
| `Genotype` | Genotype category | ✅ |
| `METRIC.Recall` | Recall proportion (0–1) | ✅ |
| `METRIC.Precision` | Precision proportion (0–1) | ✅ |
| `QUERY.FP` | False positive count | ✅ |
| `QUERY.TP` | True positive count | ✅ |
| `TRUTH.TOTAL` | Total truth variants | ✅ |
| `Subset.IS_CONF.Size` | Confident region size (bp) | ✅ |

> **Tip:** These columns are standard in hap.py v0.3.15+ outputs. If your files are missing any, re-run hap.py with the `--write-vcf` and `--stratification` flags.

### Manifest file (for batch mode)

A CSV or TSV with **exactly** these columns:

```csv
sample,caller,condition,baseline_path,consensus_path
HG001,bcftools,1KG_AF10pc,HG001_bwa_GRCh38_bcftools.csv,HG001_1KG_AF10pc_bcftools.csv
HG001,bcftools,1KG_AF30pc,HG001_bwa_GRCh38_bcftools.csv,HG001_1KG_AF30pc_bcftools.csv
HG001,freebayes,1KG_AF10pc,HG001_bwa_GRCh38_freebayes.csv,HG001_1KG_AF10pc_freebayes.csv
HG002,bcftools,1KG_AF10pc,HG002_bwa_GRCh38_bcftools.csv,HG002_1KG_AF10pc_bcftools.csv
HG005,bcftools,1KG_AF10pc,HG005_bwa_GRCh38_bcftools.csv,HG005_1KG_AF10pc_bcftools.csv
```

| Column | Description |
|--------|-------------|
| `sample` | Sample identifier (e.g., HG001) |
| `caller` | Variant caller name (e.g., bcftools, gatk, freebayes) |
| `condition` | Consensus/reference condition (e.g., 1KG_AF10pc, novo_GRCh38) |
| `baseline_path` | Path to baseline hap.py CSV |
| `consensus_path` | Path to consensus hap.py CSV |

> **Note:** The `condition` column can be any label. It appears in output tables to distinguish AF thresholds, aligners, or population-specific references.

---

## Running the Script

### Command-line arguments

```
usage: calculate_fp_burden.py [-h] [--baseline BASELINE] [--consensus CONSENSUS]
                                [--manifest MANIFEST] [--sample SAMPLE] [--caller CALLER]
                                [--condition CONDITION] [--regions REGIONS [REGIONS ...]]
                                [--filters FILTERS [FILTERS ...]] [--output-prefix OUTPUT_PREFIX]
                                [--output OUTPUT] [--group-by GROUP_BY [GROUP_BY ...]]
                                [--verbose] [--dry-run] [--version]
```

### Input Options

| Argument | Required | Description |
|----------|----------|-------------|
| `--baseline` | Yes* | Path to baseline hap.py CSV |
| `--consensus` | Yes* | Path to consensus hap.py CSV |
| `--manifest` | Yes* | Path to batch manifest CSV/TSV |
| `--sample` | No | Sample ID (default: `SAMPLE`) |
| `--caller` | No | Caller name (default: `CALLER`) |
| `--condition` | No | Condition label (default: `CONDITION`) |

*Either `--manifest` OR both `--baseline` and `--consensus` are required.

### Filtering Options

| Argument | Description |
|----------|-------------|
| `--regions` | Space-separated Subset names to retain (e.g., `AllAutosomes lowmappabilityall`) |
| `--filters` | Space-separated Filter values to retain (e.g., `ALL PASS`) |

### Output Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--output-prefix` | `fp_burden` | Prefix for all output files in batch mode |
| `--output` | Auto | Output CSV path (single mode only) |
| `--group-by` | See below | Columns for aggregation grouping |

Default grouping: `Caller`, `Condition`, `Type`, `Subset`, `Filter`, `Genotype`

### Utility Options

| Argument | Description |
|----------|-------------|
| `--verbose`, `-v` | Enable debug-level logging |
| `--dry-run` | Validate all inputs without writing outputs |
| `--version` | Show version and exit |

---

## Output Files

### Batch mode produces 6 files:

| # | File | Description | Use case |
|---|------|-------------|----------|
| 1 | `{prefix}_raw_per_sample.csv` | Every metric for every comparison | Supplementary data, custom analysis |
| 2 | `{prefix}_summary_median_range.csv` | Median [min–max] across samples | Statistical reporting, further aggregation |
| 3 | `{prefix}_table.md` | Markdown table | Paste into Word, GitHub, or supplements |
| 4 | `{prefix}_table.tex` | LaTeX table | Overleaf, journal submission |
| 5 | `{prefix}_table.csv` | Clean CSV table | Direct import to Excel/Word/R |
| 6 | `{prefix}_reproducibility.md` | Full formula audit | Methods section, reviewer response |

### File details

#### 1. Raw per-sample CSV

Every metric for every sample-caller-condition combination. Contains 40+ columns including:
- Raw hap.py metrics (0–1 scale): `Recall_base_raw`, `Recall_cons_raw`
- Percentage-scale metrics (0–100): `Recall_base_pct`, `Recall_cons_pct`
- PP gains: `Recall_gain_pp`, `Precision_gain_pp`
- Verification columns: `Recall_gain_pp_check`, `Precision_gain_pp_check`
- Absolute counts: `Additional_TP`, `Additional_FP`, `Additional_FN`
- **Ratios (v5.1):** `TP_FP_ratio` (adjusted), `TP_FP_ratio_raw`, `TP_FP_ratio_adj`
- **Net benefit (v5.1):** `Net_benefit` — signed TP − FP
- Normalized: `Additional_TP_per_Mbp`, `Additional_FP_per_Mbp`, `Additional_TP_per_100_truth`
- Metadata: `Subset.IS_CONF.Size_base`, `in_both`

#### 2. Summary CSV

One row per aggregation group with statistics:
- `{metric}_median` — median across samples
- `{metric}_min` — minimum across samples
- `{metric}_max` — maximum across samples
- `{metric}_count` — number of samples
- `{metric}_formatted` — human-readable string (e.g., `4.75 [2.27–7.11]`)

#### 3. Markdown table (`_table.md`)

GitHub-flavored Markdown with renamed columns:
- `Recall Gain (pp)`, `Precision Gain (pp)`
- `Additional TP`, `Additional FP`
- `TP:FP Ratio` (adjusted, with edge-case labels)
- `Net Benefit` (signed, clinically intuitive)
- `Additional TP/Mbp`, `Additional FP/Mbp`

Includes interpretation notes at the bottom.

#### 4. LaTeX table (`_table.tex`)

Standard `tabular` environment with:
- `escape=False` (preserves `–` en-dash in ranges)
- Caption included
- Compatible with `longtable` if needed

#### 5. CSV table (`_table.csv`)

Clean CSV with machine-friendly column names (no spaces, no special characters):

| Column | Example value |
|--------|---------------|
| `Caller` | bcftools |
| `Condition` | 1KG_AF10pc |
| `Type` | INDEL |
| `Subset` | lowmappabilityall |
| `Filter` | ALL |
| `Recall_Gain_pp` | 4.75 [2.27–7.11] |
| `Precision_Gain_pp` | -3.38 [-3.39–-3.20] |
| `Additional_TP` | 166 [79–249] |
| `Additional_FP` | 215 [173–251] |
| `TP_FP_Ratio` | 0.77 [0.46–0.99] or "∞ (pure gain)" |
| `Net_Benefit` | -49 [-172–+74] |
| `Additional_TP_per_Mbp` | 3.32 [1.58–4.98] |
| `Additional_FP_per_Mbp` | 4.30 [3.46–5.02] |
| `Additional_TP_per_100_truth` | 1.11 [0.53–1.66] |
| `Additional_FP_per_100_truth` | 1.43 [1.15–1.67] |

**Use cases:**
- **Excel:** Direct import, no column renaming needed
- **R:** `read.csv()` → `ggplot()` immediately
- **Word:** Mail merge or table import
- **Python:** `pd.read_csv()` → analysis pipeline

#### 6. Reproducibility audit (`_reproducibility.md`)

Self-documenting formula reference including:
- PP derivation formula with verification procedure
- Absolute burden derivation
- **Edge-case handling (v5.1):** How zero counts are managed
- Normalization formulas
- Aggregation method
- Data provenance (samples, callers, conditions, row counts)
- Validation checks performed

---

## Interpreting Results

### TP:FP Ratio — clinical utility metric

The **adjusted ratio** uses Laplace smoothing (+1 to both TP and FP) to handle zero counts gracefully:

| Adjusted Ratio | Meaning | Clinical implication |
|----------------|---------|----------------------|
| **> 1.0** | More signal than noise | Net beneficial — adopt consensus |
| **~1.0** | Equal trade-off | Marginal — depends on filtering cost |
| **< 1.0** | More noise than signal | Questionable without post-filtering |

**Raw ratio** (in `_raw_per_sample.csv`) may show `inf` when FP = 0. This is **good** — it means pure gain at zero cost.

### Net Benefit — the most intuitive metric

$$
	ext{Net Benefit} = 	ext{Additional TP} - 	ext{Additional FP}$$

| Net Benefit | Interpretation |
|-------------|----------------|
| **Positive** | Consensus recovers more true variants than false positives introduced |
| **Zero** | Break-even |
| **Negative** | More false positives than true variants gained |

**Example:**

| Region | Additional TP | Additional FP | TP:FP Ratio | Net Benefit |
|--------|-------------|---------------|-------------|-------------|
| CMRG (INDEL) | +12 | +14 | 0.87 | **−2** |
| Low-mappability (INDEL) | +166 | +215 | 0.78 | **−49** |
| AllAutosomes (SNP) | +2628 | +4446 | 0.61 | **−1818** |

> **CMRG:** Near break-even (−2), but every extra TP is potentially pathogenic. **Consider with filtering.**  
> **Low-mappability:** Net negative (−49). **Needs filtering before clinical use.**  
> **AllAutosomes:** Strongly net negative. **Only viable with aggressive post-filtering.**

### Edge-case handling for ratios

When FP = 0 (or TP = 0), the raw ratio is mathematically undefined. The script handles this as follows:

| Scenario | Raw ratio | Adjusted ratio | Net benefit | Display |
|----------|-----------|----------------|-------------|---------|
| TP > 0, FP = 0 | `inf` | (TP+1)/(0+1) | +TP | **∞ (pure gain)** |
| TP = 0, FP > 0 | 0 | 1/(FP+1) | −FP | **0 (cost only)** |
| TP = 0, FP = 0 | `nan` | 1.0 | 0 | **0 (no change)** |
| TP < 0, FP = 0 | `-inf` | (TP+1)/(0+1) | TP | **−∞ (pure loss)** |

> **Note:** The adjusted ratio (used in summary statistics) prevents `inf`/`nan` while preserving ranking. The raw ratio is preserved in per-sample output for transparency.

### Percentage-point verification

All pp values in the output can be independently verified from raw columns:

```
Recall_gain_pp = (Recall_cons_raw - Recall_base_raw) × 100
```

The `*_check` columns confirm this to machine precision. Any discrepancy triggers a warning.

---

## Troubleshooting

### Error: "hap.py file not found"

- Check that paths in manifest are correct relative to working directory
- Use absolute paths if unsure: `/home/user/data/HG001_baseline.csv`

### Error: "Missing required columns"

- Ensure hap.py was run with `--write-vcf` and stratification BED files
- Check that the summary file is not the VCF itself — look for `.csv` or `.summary.csv`

### Error: "No rows remaining after region/filter selection"

- Verify that `--regions` names exactly match `Subset` values in hap.py output
- Common names: `AllAutosomes`, `AllHomopolymers_ge7bp_imperfectge11bp_slop5`, `lowmappabilityall`, `segdups`, `CMRG`, `MHC`
- Use `--verbose` to see available regions before filtering

### Warning: "Recall values outside [0,1]"

- Usually harmless if slightly above 1.0 due to floating-point rounding
- If severely outside range, check for data corruption in hap.py output

### All comparisons fail

- Run with `--dry-run` first to validate manifest
- Run with `--verbose` to see exactly which row fails
- Check that baseline and consensus files have matching `Type`, `Subset`, `Filter` values

### TP:FP ratio shows "∞ (pure gain)" or "0 (cost only)"

- This is **expected** when one sample has zero FP or zero TP
- Check `_raw_per_sample.csv` to see which sample caused it
- The adjusted ratio and net benefit columns provide interpretable alternatives

---

## Examples

### Example 1: Validate manifest without running

```bash
python calculate_fp_burden.py \
    --manifest manifest.csv \
    --dry-run
```

**Expected output:**
```
[INFO] Manifest validated: 18 comparison(s) ready
[INFO] DRY RUN: All inputs validated. No outputs will be written.
```

### Example 2: Process only SNVs in difficult regions

```bash
python calculate_fp_burden.py \
    --manifest manifest.csv \
    --output-prefix snv_difficult/ \
    --regions lowmappabilityall segdups CMRG MHC \
    --filters ALL
```

### Example 3: Compare AF thresholds side-by-side

```bash
python calculate_fp_burden.py \
    --manifest manifest.csv \
    --output-prefix af_comparison/ \
    --regions AllAutosomes \
    --filters ALL \
    --group-by Caller Condition Type Subset Filter
```

### Example 4: Import CSV table into R for plotting

```r
# After running the script
df <- read.csv("manuscript_tables_table.csv")

# Plot Net Benefit by region
library(ggplot2)
ggplot(df, aes(x = Subset, y = Net_Benefit, fill = Type)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~Caller) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Net Benefit of IUPAC Consensus Reference",
       subtitle = "Positive = net gain, Negative = net cost",
       y = "Net Benefit (TP - FP)", x = "Genomic Region")
```

### Example 5: Import CSV table into Excel

1. Open Excel → Data → From Text/CSV
2. Select `manuscript_tables_table.csv`
3. Excel auto-detects comma delimiter
4. Columns are ready to use (no spaces or special characters)

### Example 6: Debug a failing comparison

```bash
python calculate_fp_burden.py \
    --baseline HG001_baseline.csv \
    --consensus HG001_consensus.csv \
    --sample HG001 --caller bcftools --condition test \
    --output debug.csv \
    --verbose
```

---

## Citation

If you use this script in published work, please cite:

> Saidin A, Ricos MG, Dibbens LM. *IUPAC consensus references improve variant detection in clinically challenging genomic regions.* [Journal TBD]. 2026.

And reference the GitHub repository:

```
https://github.com/akzam/IUPAC-consensus-reference
```

---

## License

MIT License — see LICENSE file in repository.

---

## Contact

For issues or feature requests, open a GitHub Issue at:
https://github.com/akzam/IUPAC-consensus-reference/issues
