# BurdenBench

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Calculate absolute false-positive (FP) and true-positive (TP) burden from hap.py benchmarking outputs.**

BurdenBench supports multiple consensus conditions, aggregates medians across samples, and produces publication-ready tables with a full reproducibility audit — including explicit percentage-point (pp) derivation.

> **Previous name:** `calculate_fp_burden.py` (v5.3.7)

---

## Quick Start

### Single comparison
```bash
python run_BurdenBench.py \
    --baseline HG001_bwa_GRCh38_bcftools.csv \
    --consensus HG001_1KG_AF10pc_bcftools.csv \
    --sample HG001 --caller bcftools --condition AF10 \
    --output results.csv
```

### Batch processing with aggregation (recommended for manuscripts)
```bash
python run_BurdenBench.py \
    --manifest manifest.csv \
    --output-prefix manuscript_tables/ \
    --regions AllAutosomes lowmappabilityall segdups CMRG MHC \
    --filters ALL
```

---

## Installation

```bash
pip install pandas numpy
```

**Requirements:**
- Python >= 3.8
- pandas >= 1.3.0
- numpy >= 1.20.0

---

## Documentation

- [Usage Manual](docs/USAGE_MANUAL.md) — Command-line reference, troubleshooting, and examples
- [Technical Manual](docs/TECHNICAL_MANUAL.md) — Architecture, algorithms, API reference, and testing

---

## Key Features

- **Explicit percentage-point (pp) derivation** with raw 0–1 value preservation
- **Laplace-adjusted TP:FP ratio** (handles zero counts gracefully)
- **Net Benefit metric** (`Additional TP − Additional FP`) for clinical interpretation
- **Batch aggregation** with median [min–max] across samples
- **Multiple output formats**: raw CSV, summary CSV, Markdown, LaTeX, and clean CSV tables
- **Full reproducibility audit** auto-generated for every run
- **Overall aggregation** (summed across stratifications) at configurable levels
- **Median-of-medians sensitivity analysis**

---

## Output Files

| # | File | Description |
|---|------|-------------|
| 1 | `{prefix}_raw_per_sample.csv` | Every metric for every comparison |
| 2 | `{prefix}_summary_median_range.csv` | Median [min–max] across samples |
| 3 | `{prefix}_table.txt` | Markdown table |
| 4 | `{prefix}_table.tex` | LaTeX table |
| 5 | `{prefix}_table.csv` | Clean CSV table |
| 6 | `{prefix}_reproducibility.txt` | Full formula audit |
| 7 | `{prefix}_overall_summary_median_range.csv` | Overall aggregated metrics |
| 8 | `{prefix}_overall_table.*` | Overall tables (txt/tex/csv) |
| 9 | `{prefix}_median_of_medians_summary.csv` | Sensitivity analysis |
| 10 | `{prefix}_median_of_medians_table.*` | MoM tables (txt/csv) |

---

## Citation

If you use BurdenBench in published work, please cite:

> Saidin A, Ricos MG, Dibbens LM. *IUPAC consensus references improve variant detection in clinically challenging genomic regions.* [In preparation]. 2026.

---

## License

MIT License — see [LICENSE](LICENSE).

---

## Contact

For issues or feature requests, please open a GitHub Issue.
