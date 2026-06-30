# BurdenBench Examples

This directory contains example input files for testing BurdenBench.

## Files

- `manifest_example.csv` — Example batch manifest with 6 comparisons
- `example_happy_baseline.csv` — Example baseline hap.py output (3 regions × 2 types)
- `example_happy_consensus.csv` — Example consensus hap.py output (3 regions × 2 types)

## Quick test

```bash
# Single comparison
python run_BurdenBench.py \
    --baseline examples/example_happy_baseline.csv \
    --consensus examples/example_happy_consensus.csv \
    --sample HG001 --caller bcftools --condition AF10 \
    --output examples/single_output.csv

# Batch processing
python run_BurdenBench.py \
    --manifest examples/manifest_example.csv \
    --output-prefix examples/output/ \
    --regions AllAutosomes lowmappabilityall CMRG \
    --filters ALL
```
