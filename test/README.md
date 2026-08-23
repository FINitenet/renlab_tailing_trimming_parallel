# Minimal smoke test

This directory contains a synthetic, non-biological example that exercises the
miR3End profile stage without downloading a reference genome or running an
external aligner. It verifies that a precomputed SAM file can be converted into
the expected profile, summary, 5GMC and sequence-logo inputs.

From the repository root, after creating the Conda environment, run:

```bash
conda env create -f environment.yml
conda activate mir3end-0.1.0
python test/run_smoke_test.py
```

The test uses a temporary output directory and leaves the repository unchanged.
The synthetic files are only for software validation and must not be used for
biological interpretation.
