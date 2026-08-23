# Changelog

All notable changes to miR3End are documented here.

## [0.1.0] - 2026-08-23

### Added

- First public release candidate of the parallel and resumable miRNA 3′-end
  trimming/tailing workflow.
- Versioned Conda environment and executable `--version` output.
- Synthetic profile-stage test data, smoke-test runner, and GitHub Actions CI.
- MIT license, citation metadata, Zenodo metadata, and release checklist.
- Reference-resource documentation independent of the original server layout.

### Changed

- Removed machine-specific default paths for indexes, annotations, ShortStack,
  featureCounts, and Rscript.
- External tools are now resolved from `PATH`; reference paths can be supplied
  by command-line options or `MIR3END_*` environment variables.
