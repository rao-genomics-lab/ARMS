# Changelog

All notable changes to the ARMS (Adaptive Resolution Multiscale Spatial DNA sequencing) package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- **Sensible clustering defaults matching the typical workflow** (2026-07-05)
  - Function defaults now reproduce the canonical analysis (`test_workflow_2.Rmd`) without needing to override the same arguments on every call. `weighted_hierarchical_clustering()` is now the recommended default entry point in the README and vignette.
  - `weighted_hierarchical_clustering()` default changes: `discretize` `FALSE` → `TRUE`, `loss_threshold` `-0.2` → `-0.1`, `use_size_weighting` `TRUE` → `FALSE`, `min_cluster_size` `10` → `3`, `silhouette_threshold` `0.05` → `0.001`, `pvclust_threshold` `0.95` → `0.3`.
  - `hierarchical_clustering()` default changes: `min_cluster_size` `10` → `3`, `silhouette_threshold` `0.05` → `0.001`, `pvclust_threshold` `0.95` → `0.3`.
  - `split_clusters_manually()` weighted-path defaults aligned: `discretize` `FALSE` → `TRUE`, `loss_threshold` `-0.2` → `-0.1`, `use_size_weighting` `TRUE` → `FALSE`.
  - **Behavior change**: `weighted_hierarchical_clustering()` now discretizes logR by default (handles tumor-purity differences), and size weighting is off by default. Because `use_size_weighting` now defaults to `FALSE`, a bare `weighted_hierarchical_clustering(logr_matrix)` call no longer requires `grid` and no longer errors.
  - **Experimental features de-emphasized**: correlation-based `refine_clusters_by_correlation()` and the size-weighting machinery are now labeled experimental in their docs and removed from the default README/vignette path (both remain exported and callable).
  - **Verified**: a bare `weighted_hierarchical_clustering(logr_matrix)` call produces cluster assignments identical to the previous explicit-argument call, and reproduces the reference clustering partition on the canonical dataset (per-sample co-assignment agreement 1.0).
  - **Known limitation (pre-existing)**: the recursive split step evaluates a gap statistic and pvclust bootstraps without a fixed seed, so the pipeline is not bit-for-bit reproducible unless the caller sets `set.seed()` (the split decisions are silhouette-driven in practice, so the final partition is stable). These bootstraps run whenever a cluster exceeds `max_module_size` and can be slow at fine bin resolution.
  - **Files modified**: `R/clustering.R`, `R/analysis.R`, `README.md`, `vignettes/example_workflow.Rmd`, `CLAUDE.md`

### Added

- **Optional sample-name row labels for LogR heatmap** (2026-07-02)
  - `plot_logr_heatmap()` gains two new parameters: `show_sample_names` (default `FALSE`) to display sample/tile ID row labels, and `row_names_fontsize` (default `6`) to control their size
  - Row labels render on the right side of the heatmap and do not disturb the existing left-side cluster color bar, ploidy annotation, or per-cluster dendrogram
  - Default behavior is unchanged for existing callers (`show_sample_names = FALSE`)
  - **Files modified**: `R/visualization.R`

- **Manual cluster splitting** (2026-02-20)
  - New exported function `split_clusters_manually()` in `R/clustering.R` is the inverse of `merge_clusters_manually()`
  - Takes specified cluster(s), applies hierarchical clustering to subdivide them, and integrates new subclusters back into the original clustering object
  - **Method selection**:
    - `method = "hierarchical"`: Standard hierarchical clustering with Ward's D2 linkage
    - `method = "weighted"`: Weighted clustering with optional discretization and size-based weighting
  - **Clustering parameters**:
    - `k`: Number of subclusters (NULL = auto-detect via dynamic tree cut)
    - `use_dynamic`: Use dynamic tree cut for auto k detection
    - `min_cluster_size`: Minimum subcluster size (default: 3)
  - **Weighted method parameters**: `discretize`, `gain_threshold`, `loss_threshold`, `use_size_weighting`, `min_bins_full_weight`, `weight_function`
  - **Output control**:
    - `renumber`: Renumber clusters sequentially after splitting (default: TRUE)
  - **Edge case handling**:
    - Warns and skips clusters with < 2 samples
    - Warns if clustering produces only 1 subcluster
    - Automatically reduces k if greater than sample count
  - Returns updated clustering object with `$split_info` documenting which clusters were split
  - **Files modified**: `R/clustering.R`, `NAMESPACE`

- **Weighted hierarchical clustering** (2026-02-19)
  - New exported function `weighted_hierarchical_clustering()` in `R/clustering.R` addresses two common challenges in copy number clustering:
    1. **Purity differences**: Optional discretization converts continuous logR to gain/neutral/loss categories, so samples with identical CNA profiles but different tumor purity cluster together
    2. **Small aberrations**: Optional size-based weighting reduces influence of small altered regions (1-2 bins) that may be sequencing artifacts
  - **Discretization parameters**:
    - `discretize`: Enable logR discretization (default: FALSE)
    - `gain_threshold` / `loss_threshold`: Thresholds for gain/loss classification (default: ±0.2)
    - `discretize_mode`: "ternary" (+1/0/-1) or "binary" (1/0)
  - **Size weighting parameters**:
    - `use_size_weighting`: Enable size-based weighting (default: TRUE)
    - `min_bins_full_weight`: Minimum contiguous bins for full weight (default: 5)
    - `weight_function`: "linear", "sigmoid", or "step"
    - `weight_scope`: "per_sample" or "consensus"
  - Returns extended output including `weight_matrix`, `discretized_matrix`, and `parameters`
  - Four new internal helper functions: `discretize_logr()`, `identify_altered_regions()`, `compute_size_weights()`, `weighted_euclidean_distance()`
  - **Files modified**: `R/clustering.R`

- **Raw counts import for ASCAT.sc pipeline** (2026-02-04)
  - New exported function `create_res_from_counts()` in `R/counts-to-ascat.R` enables users with pre-binned read count data to run the full ASCAT.sc copy number calling pipeline without BAM files
  - Accepts a tab-separated file with bins as rows (`chr:start-end` format, e.g. `1:1-500000`) and samples as columns
  - Computes GC content from a BSgenome object for GC correction
  - Builds ASCAT.sc-compatible track structures, runs segmentation, purity/ploidy fitting, and profile generation
  - Writes standard ASCAT.sc profile files and saves the result `.Rda` object for downstream ARMS analysis
  - Two internal helpers: `parse_counts_file()` (parses bin coordinates, auto-detects binsize) and `compute_gc_for_bins()` (efficient GC computation via `Biostrings::Views()`)
  - Supports configurable purity/ploidy grids, segmentation alpha, and parallel processing via `MC.CORES`
  - **Files added**: `R/counts-to-ascat.R`
  - **Dependencies added**: `BSgenome` (Imports); `BSgenome.Hsapiens.UCSC.hg19`, `BSgenome.Hsapiens.UCSC.hg38` (Suggests)

- **hg19 genomic grid** (2026-02-04)
  - New exported function `create_genomic_grid_hg19()` in `R/data-loading.R` provides a fixed grid based on hg19 chromosome lengths (NCBI-style naming without "chr" prefix)
  - Mirrors the existing `create_genomic_grid_hg38()` structure

### Changed

- **Auto-detect chromosome naming in `create_genomic_grid()`** (2026-02-04)
  - `create_genomic_grid()` now detects whether profile data uses "chr1" or "1" style chromosome names instead of hardcoding `paste0("chr", c(1:22, "X"))`
  - Enables seamless grid creation for both hg19 (NCBI-style) and hg38 (UCSC-style) profiles
  - **Files modified**: `R/data-loading.R`

- **Manual cluster merging** (2026-02-02)
  - New exported function `merge_clusters_manually()` in `R/clustering.R` allows users to manually merge specified clusters
  - Accepts a single vector (e.g. `c(1, 3, 5)` — merge 3 and 5 into 1) or a list of vectors for multiple independent merges
  - `renumber` parameter (default `TRUE`) renumbers clusters sequentially after merging; set to `FALSE` to preserve original target IDs
  - Validates that all cluster IDs exist, no ID appears in multiple merge groups, and warns on single-element groups
  - Returns the same clustering object structure with updated `$clusters`, compatible with all downstream functions
  - **Files modified**: `R/clustering.R`, `NAMESPACE`

- **Row-mean NA imputation for LogR matrix** (2026-01-29)
  - New exported function `impute_missing()` in `R/utils.R` provides two imputation strategies:
    - `"row_mean"` (new default): Two-pass approach — first replaces each NA with the mean of that sample's non-NA values, then fills any remaining NAs (from entirely-NA rows) with the overall matrix mean. This preserves each sample's baseline LogR level rather than forcing missing positions to the diploid baseline.
    - `"scalar"`: Replaces all NAs with a fixed value (previous behaviour, default 0).
  - `load_ascat_profiles()` and `align_profiles_to_grid()` now accept an `imputation_method` parameter (default `"row_mean"`)
  - The previous behaviour is available via `imputation_method = "scalar"`; the existing `logr_missing_value` parameter continues to control the fill value in scalar mode
  - Informational message now reports the number of NAs imputed and the method used
  - **Files modified**: `R/utils.R` (new `impute_missing()`), `R/data-loading.R` (both loading functions), `NAMESPACE`
  - **Rationale**: Scalar 0 imputation assumes missing regions are diploid, which creates sharp discontinuities for samples with widespread CN alterations. Row-mean imputation preserves each sample's overall ploidy level, matching the strategy used by the ASCAT workflow scripts. See `ascat_workflow_vs_ARMS_differences.md` section 3.4 for detailed comparison.

- **Manual cluster assignment from GeoJSON groups** (2025-12-15)
  - New function `create_clusters_from_geojson_groups()` allows manual definition of spatial groups using multiple GeoJSON files
  - Each GeoJSON file represents a distinct cluster/group for pseudobulk analysis
  - Automatically matches simplified tile names (e.g., "plate4_11E") to full sample names using existing normalization infrastructure
  - Returns cluster assignments in same format as `hierarchical_clustering()`, enabling seamless integration with pseudobulk workflow
  - **Files modified**: `R/clustering.R` (lines 328-475)
  - **Use case**: Manually select tumor regions, normal regions, or other spatial groups for targeted pseudobulk copy number profiling

### Changed

- **Enhanced ASCAT pseudobulk/bulk functions with custom purity/ploidy grids** (2025-12-16)
  - `run_ascat_on_pseudobulk_tracks()` and `run_ascat_on_bulk_track()` now accept optional `purs` and `ploidies` parameters
  - When specified, these custom grids are used for ASCAT fitting instead of the default grids from the original ASCAT run
  - If not specified (NULL), functions use the original behavior (extract grids from ascat_rda_file)
  - Useful for refining purity/ploidy search ranges based on initial results or known sample characteristics
  - **Files modified**: `R/ascat-core.R` (function signatures and fitting logic)

- **Improved spatial cluster visualization** (2025-12-15)
  - `plot_spatial_clusters()` now fills actual GeoJSON tile polygons with cluster colors by default (more accurate spatial representation)
  - Added `use_polygons` parameter: set to `FALSE` to use the old behavior (colored square points)
  - Added `flipped` parameter (default `TRUE`): reverses the y-axis so the origin is top-left, matching typical image coordinates; set to `FALSE` for standard mathematical orientation (origin bottom-left)
  - **Files modified**: `R/visualization.R` (lines 460-568)

### Fixed

- **Chromosome label rendering in heatmap annotations** (2026-07-02)
  - `create_chr_annotation()` (internal helper used by `plot_logr_heatmap()`) now uses `ComplexHeatmap::anno_text()` instead of `anno_mark()` to render chromosome number labels at column midpoints
  - **Root cause**: `anno_mark()` draws labels with connector lines designed for sparse annotations and could misplace or omit labels when applied densely across chromosome midpoints on the column axis
  - **Solution**: Build a per-column label vector (empty string except at each chromosome's midpoint) and render with `anno_text()`, which places labels directly without connector-line logic
  - **Files modified**: `R/visualization.R`

- **Critical bug in spatial matching** (2025-12-15)
  - Fixed sample name normalization failures in `match_spatial_tiles_to_samples()` and `match_samples_to_tiles()`
  - **Root cause**: Missing `@importFrom stringr str_match str_split` directives caused internal functions to silently fail
  - **Impact**: All sample names were returned unchanged instead of being normalized, resulting in 0% match rate
  - **Solution**: Added proper stringr imports to exported functions
  - **Files modified**: `R/spatial-matching.R` (lines 130, 263)

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
