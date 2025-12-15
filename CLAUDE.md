# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

ARMS (Adaptive Resolution Multiscale Spatial DNA sequencing) is an R/Bioconductor package for detecting copy number alterations in spatial genomics data. It wraps ASCAT.sc for single-cell copy number calling, performs hierarchical clustering to identify clonal populations, and generates publication-ready visualizations.

## Development Commands

### Building and Checking
```r
# Build package
devtools::build()

# Check package (R CMD check)
devtools::check()

# Install package locally
devtools::install()

# Load package for development
devtools::load_all()
```

### Documentation
```r
# Generate documentation from roxygen2 comments
devtools::document()

# Build vignettes
devtools::build_vignettes()
```

### Testing
```r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-clustering.R")
```

## Architecture

### Core Workflow Pipeline

The package implements a five-stage analysis pipeline:

1. **ASCAT.sc Integration** (`R/ascat-core.R`)
   - Wraps ASCAT.sc's `run_sc_sequencing()` for initial single-cell copy number calling
   - `run_ascat_on_pseudobulk_tracks()`: Reruns ASCAT.sc on aggregated pseudobulk profiles
   - `run_ascat_on_bulk_track()`: Processes bulk/reference samples

2. **Data Loading** (`R/data-loading.R`)
   - `load_ascat_profiles()`: Primary entry point - reads ASCAT.sc profile files (*.ASCAT.scprofile.txt) into matrices
   - Creates aligned genomic grids at specified resolution (default: 1Mb)
   - Handles both refitted and non-refitted profiles
   - Returns: `copy_number_matrix`, `logr_matrix`, `grid`, `sample_names`

3. **Clustering** (`R/clustering.R`)
   - `hierarchical_clustering()`: Ward's D2 linkage with dynamic tree cutting
   - `refine_clusters_by_correlation()`: Post-processing to merge similar clones
   - Quality metrics: silhouette scores, gap statistics, pvclust AU p-values
   - Uses `dynamicTreeCut::cutreeDynamic()` for automatic cluster detection

4. **Pseudobulk Analysis** (`R/pseudobulk.R`)
   - Two approaches:
     - `create_pseudobulk_from_bins()`: Fast aggregation from bin counts (preferred)
     - `merge_bams_by_cluster()`: BAM-level merging (slower, for resequencing)
   - Improves signal-to-noise ratio by aggregating cells within clusters
   - Feeds back into ASCAT.sc for refined copy number calling

5. **Visualization** (`R/visualization.R`)
   - `plot_logr_heatmap()`: Primary visualization using ComplexHeatmap
   - `plot_umap()`: Dimensionality reduction with uwot
   - `plot_spatial_clusters()`: GeoJSON-based spatial mapping
   - `plot_cluster_correlation()`: Inter-cluster similarity
   - Quality plots: ploidy distribution, chromosomal instability, logR variance

### Spatial Integration

**Spatial Matching** (`R/spatial-matching.R`)
- `match_spatial_tiles_to_samples()`: Links sample IDs to spatial coordinates
- Complex filename parsing for plate/well identifiers:
  - Handles multiple barcode schemes (16-barcode and 24-barcode layouts)
  - Parses patterns like "plate2-1A-2H_lcmad_007" → "P2_A1"
- `load_geojson_tiles()`: Reads spatial geometry data
- Enables spatial autocorrelation analysis (Moran's I, LISA)

### Key Data Structures

**ASCAT Result Object** (`res`)
- Loaded from `.Rda` files produced by ASCAT.sc
- Contains: `allTracks`, `chr`, `sex`, `binsize`, `lGCT`, `lSe`, `segmentation_alpha`
- Each track includes: `lCTS.tumour` (bin counts), `profile` (copy number)

**Matrices**
- `logr_matrix`: Samples × genomic bins, log2 ratios
- `copy_number_matrix`: Samples × genomic bins, integer copy numbers
- `grid`: Genomic coordinates (chr, start, end, pos) for matrix columns

**Clustering Output**
- `clusters`: Vector of cluster assignments (integers)
- `dendrogram`: hclust object for visualization
- `quality_metrics`: Silhouette scores, gap statistics

## Important Implementation Details

### Genomic Resolution
- Default resolution: 1Mb (1,000,000 bp)
- ASCAT.sc binsize typically: 30kb (30,000 bp)
- Profiles are aligned to a common grid using `create_genomic_grid()` and `align_profiles_to_grid()`

### Chromosome Handling
- Works with autosomes (1-22) and sex chromosome X
- Chromosome names: "chr1", "chr2", ..., "chr22", "chrX"
- Sex-specific analysis supported through ASCAT.sc's `sex` parameter

### Sample Naming Conventions
- BAM files → ASCAT profiles: `sample.bam` → `sample.bam.ASCAT.scprofile.txt`
- Refitted profiles: `sample.bam_refitted.ASCAT.scprofile.txt`
- Pseudobulk: `cluster_1_pseudobulk.ASCAT.scprofile.txt`
- The package handles multiple naming variations through `normalize_identifier()` (`R/utils.R`)

### Performance Considerations
- Pseudobulk from bins is 10-100x faster than BAM merging
- Use `min_cluster_size` parameter to filter small clusters
- Clustering quality metrics use reduced bootstrap samples for speed (default: B=10)
- Parallel processing via `parallel::mclapply()` in ASCAT.sc wrapper

## Typical Workflow Usage

```r
# 1. Run ASCAT.sc (external, or via wrapper)
result <- ASCAT.sc::run_sc_sequencing(tumour_bams, ...)

# 2. Load profiles into matrices
profiles <- load_ascat_profiles("output_dir/", resolution = 1e6)

# 3. Cluster samples
clustering <- hierarchical_clustering(profiles$logr_matrix, min_cluster_size = 3)

# 4. Visualize
plot_logr_heatmap(profiles$logr_matrix, clustering$clusters, profiles$grid, "heatmap.pdf")
plot_umap(profiles$logr_matrix, clustering$clusters, "umap.pdf")

# 5. Optional: Pseudobulk refinement
pb_tracks <- create_pseudobulk_from_bins("result.Rda", clustering$clusters,
                                         profiles$sample_names, "pb_output/")
pb_result <- run_ascat_on_pseudobulk_tracks(pb_tracks, "result.Rda", "pb_output/")

# 6. Optional: Spatial analysis (if GeoJSON available)
matches <- match_spatial_tiles_to_samples(profiles$sample_names, "tiles.geojson")
plot_spatial_clusters(matches, clustering$clusters, "spatial_map.pdf")
```

## Dependencies

### Critical External Packages
- **ASCAT.sc**: VanLoo-lab implementation (GitHub-only, specified in Remotes)
- **ComplexHeatmap**: Bioconductor, for advanced heatmap generation
- **uwot**: UMAP dimensionality reduction
- **dynamicTreeCut**: Automatic dendrogram cutting
- **sf/geojsonio**: Spatial data handling

### Package Structure
- `R/`: All source code (8 files organized by function)
- `man/`: Auto-generated Roxygen2 documentation
- `vignettes/`: Example workflow Rmd showing complete analysis
- `inst/extdata/`: Example data (ASCAT profiles, GeoJSON tiles)

## Code Conventions

### Function Organization
- Exported functions: User-facing API
- Internal functions: `@keywords internal`, no export
- Helper functions: Prefixed with verb (e.g., `create_`, `calculate_`, `plot_`)

### Roxygen Documentation
- All exported functions have `@param`, `@return`, `@export`, `@examples`
- `@importFrom` for all external package functions (no library() calls in functions)
- Use `\dontrun{}` for examples requiring external data

### Data Handling
- Use `readr::read_tsv()` for ASCAT profile loading (handles large files efficiently)
- dplyr pipes (`%>%`) for data transformations
- Preserve sample order through clustering and visualization

### Error Handling
- Use `stop()` with informative messages for critical errors
- `warning()` for non-fatal issues (e.g., grid mismatch)
- `message()` for progress updates in long-running functions

## Testing Strategy

Tests should cover:
- Profile loading with different file patterns
- Clustering with varying parameters
- Matrix dimension consistency through pipeline
- Sample name matching and normalization
- Edge cases: single cluster, small sample sizes
