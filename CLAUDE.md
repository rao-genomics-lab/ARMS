# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

ARMS (Adaptive Resolution Multiscale Spatial DNA sequencing) is an R package for detecting copy number alterations in spatial genomics data. It wraps ASCAT.sc for single-cell copy number calling, performs hierarchical clustering to identify clonal populations, and generates publication-ready visualizations.

## Development Commands

```r
# Load package for development (use this during iterative development)
devtools::load_all()

# Generate documentation from roxygen2 comments (run after changing @param/@export etc.)
devtools::document()

# Check package (R CMD check equivalent)
devtools::check()

# Build and install
devtools::build()
devtools::install()
```

**Note:** No test suite exists yet. The `tests/` directory has not been created despite `testthat` being in Suggests. When tests are added, run with `devtools::test()` or `testthat::test_file("tests/testthat/test-*.R")`.

## Architecture

### Five-Stage Pipeline

1. **ASCAT.sc Integration** (`R/ascat-core.R`) - Wraps ASCAT.sc's `run_sc_sequencing()` for copy number calling. `run_ascat_on_pseudobulk_tracks()` and `run_ascat_on_bulk_track()` accept optional custom `purs`/`ploidies` grids for refined fitting.

2. **Data Loading** (`R/data-loading.R`) - `load_ascat_profiles()` reads ASCAT.sc profile files (`*.ASCAT.scprofile.txt`) into matrices aligned to a common genomic grid (default 1Mb). Returns `copy_number_matrix`, `logr_matrix`, `grid`, `sample_names`. Supports `imputation_method` parameter (`"row_mean"` default, or `"scalar"` for legacy behavior).

3. **Clustering** (`R/clustering.R`) - `hierarchical_clustering()` uses Ward's D2 with `dynamicTreeCut::cutreeDynamic()`. `refine_clusters_by_correlation()` merges similar clusters. `create_clusters_from_geojson_groups()` enables manual cluster assignment from GeoJSON spatial selections.

4. **Pseudobulk** (`R/pseudobulk.R`) - `create_pseudobulk_from_bins()` (fast, preferred) and `merge_bams_by_cluster()` (BAM-level, slower). Aggregates cells within clusters for improved signal-to-noise.

5. **Visualization** (`R/visualization.R`, largest file ~974 lines) - `plot_logr_heatmap()` (ComplexHeatmap), `plot_umap()`, `plot_spatial_clusters()` (GeoJSON polygon fills by default), and various QC plots.

### Supporting Modules

- **Counts Import** (`R/counts-to-ascat.R`) - `create_res_from_counts()` enables running the full ASCAT.sc pipeline from pre-binned read count files (no BAM files needed). Computes GC content from BSgenome, builds track structures, and runs segmentation + purity/ploidy fitting.
- **Spatial Matching** (`R/spatial-matching.R`) - Links sample IDs to spatial coordinates via complex filename parsing. `normalize_identifier()` (internal) handles multiple barcode schemes (16-barcode and 24-barcode layouts). Parses patterns like `plate2-1A-2H_lcmad_007` to `P2_A1`.
- **Utilities** (`R/utils.R`) - `impute_missing()` (row-mean or scalar NA imputation), `compute_umap()`, `calculate_morans_i()` (spatial autocorrelation), `plot_lisa()`.
- **Analysis** (`R/analysis.R`) - Contains `compare_sc_vs_pseudobulk()` for comparing single-cell and pseudobulk profiles.

### Key Data Structures

- **ASCAT Result Object** (`res`): Loaded from `.Rda` files. Contains `allTracks`, `chr`, `sex`, `binsize`, `lGCT`, `lSe`, `segmentation_alpha`. Each track has `lCTS.tumour` (bin counts) and `profile` (copy number).
- **Matrices**: `logr_matrix` and `copy_number_matrix` are samples x genomic bins. `grid` provides genomic coordinates (chr, start, end, pos) for columns.
- **Clustering Output**: `clusters` (integer vector), `dendrogram` (hclust object), `quality_metrics` (silhouette scores, gap statistics).

## Important Details

### Genomic Resolution
- Default analysis resolution: 1Mb. ASCAT.sc binsize typically 30kb.
- Profiles aligned via `create_genomic_grid()` and `align_profiles_to_grid()`.
- Chromosomes: "chr1"-"chr22" + "chrX".

### Sample Naming
- BAM → profile: `sample.bam` → `sample.bam.ASCAT.scprofile.txt`
- Refitted: `sample.bam_refitted.ASCAT.scprofile.txt`
- Pseudobulk: `cluster_1_pseudobulk.ASCAT.scprofile.txt`
- Multiple naming variations handled by `normalize_identifier()` in `R/spatial-matching.R`.

### Unused Shiny Dependencies
The DESCRIPTION imports many Shiny packages (`shiny`, `shinydashboard`, `shinyFiles`, etc.) but no Shiny app currently exists in the codebase. These may be planned for a future interactive interface.

## Code Conventions

- All exported functions have roxygen2 docs with `@param`, `@return`, `@export`, `@examples`
- `@importFrom` for all external package functions (no `library()` calls inside functions)
- Use `\dontrun{}` for examples requiring external data
- Helper functions prefixed with verbs: `create_`, `calculate_`, `plot_`
- Internal functions use `@keywords internal`
- dplyr pipes (`%>%`) for data transformations
- `readr::read_tsv()` for file I/O
- `stop()` for critical errors, `warning()` for non-fatal issues, `message()` for progress
