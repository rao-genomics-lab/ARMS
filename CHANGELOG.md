# Changelog

All notable changes to the ARMS (Adaptive Resolution Multiscale Spatial DNA sequencing) package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- **Critical bug in pseudobulk nucleotide aggregation** (2025-12-15)
  - Fixed incorrect nucleotide counting in `create_pseudobulk_from_bins()` and `create_bulk_from_bins()`
  - **Root cause**: When aggregating multiple samples into pseudobulk profiles, the `records` field (read counts) was correctly summed across all samples, but the `nucleotides` field was only copied from the first sample instead of being summed
  - **Impact**: This mismatch caused ASCAT.sc normalization calculations to fail with `Error in if (any(nonround < 0)) { : missing value where TRUE/FALSE needed` due to invalid `records/nucleotides` ratios producing NA values
  - **Solution**: Both `records` and `nucleotides` are now properly summed across all samples during pseudobulk aggregation
  - **Files modified**: `R/pseudobulk.R` (lines 93-116 and lines 243-266)
  - This fix ensures correct normalization when running `run_ascat_on_pseudobulk_tracks()` and `run_ascat_on_bulk_track()`

## [Initial Release]

### Added

- Initial implementation of ARMS package
- ASCAT.sc integration for single-cell copy number calling
- Hierarchical clustering with dynamic tree cutting
- Pseudobulk analysis from bin counts (fast mode)
- BAM merging workflow for pseudobulk creation
- Spatial tile matching and visualization
- Publication-ready plotting functions (heatmaps, UMAP, spatial maps)
- Comprehensive vignette with example workflow
