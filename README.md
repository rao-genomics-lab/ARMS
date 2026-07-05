# ARMS

**A**daptive **R**esolution **M**ultiscale **S**patial DNA sequencing — an R package for detecting copy number alterations (CNAs) in spatial genomics data.

## Overview

ARMS wraps [ASCAT.sc](https://github.com/VanLoo-lab/ASCAT.sc) for single-cell copy number calling, performs hierarchical clustering to identify clonal populations, and generates publication-ready visualizations (LogR heatmaps, UMAP, spatial maps). It also supports pseudobulk analysis for improved signal-to-noise ratio.

## Installation

Requires R >= 4.3.0. This repository is **private**, so installing directly requires either cloning with an SSH key that has access, or a GitHub personal access token (PAT).

**1. Install Bioconductor dependencies**

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "IRanges", "GenomicRanges",
                        "DNAcopy", "BSgenome", "Biostrings"))
```

**2. Install ASCAT.sc**

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("VanLoo-lab/ASCAT.sc")
```

**3. Install ARMS**

Recommended — clone and install locally:

```bash
git clone git@github.com:rao-genomics-lab/ARMS.git
cd ARMS
```

```r
devtools::install(".")
```

Alternative — install directly from GitHub using a PAT with access to this repo:

```r
devtools::install_github("rao-genomics-lab/ARMS", auth_token = "<your_github_pat>")
```

## Quick Start

A minimal end-to-end workflow. See [`vignettes/example_workflow.Rmd`](vignettes/example_workflow.Rmd) for the full walkthrough (QC plots, spatial mapping, pseudobulk, clone refinement, bulk simulation).

**1. Run ASCAT.sc on BAM files** to generate per-sample copy number profiles:

```r
ascat_result <- ASCAT.sc::run_sc_sequencing(
  tumour_bams = bam_files,
  allchr = paste0("chr", c(1:22, "X")),
  sex = rep("male", length(bam_files)),
  binsize = 30000,
  outdir = "ascat_output/",
  projectname = "myproject",
  build = "hg38"
)
```

**2. Load the profiles** into R matrices:

```r
library(ARMS)
profiles <- load_ascat_profiles(
  profile_dir = "ascat_output/",
  pattern = "*refitted.ASCAT.scprofile.txt",
  resolution = 30000
)
```

**3. Cluster samples** into clonal populations:

```r
clustering <- weighted_hierarchical_clustering(
  logr_matrix = profiles$logr_matrix
)
```

`weighted_hierarchical_clustering()` is the recommended default: it discretizes logR
values (gain/neutral/loss) before clustering so samples with the same copy-number
profile but different tumor purity group together. Use `hierarchical_clustering()` if
you want to cluster on continuous logR instead.

**4. Visualize** the results:

```r
plot_logr_heatmap(
  logr_matrix = profiles$logr_matrix,
  grid = profiles$grid,
  clusters = clustering$clusters,
  output_file = "heatmap.pdf"
)

umap_coords <- compute_umap(profiles$logr_matrix)
plot_umap(umap_coords, clustering$clusters, "umap.pdf")
```

**5. Optional next steps** — spatial mapping and pseudobulk aggregation are supported, and clusters can be refined manually with `split_clusters_manually()` and `merge_clusters_manually()`; see the vignette for details.

## Documentation

- Full worked example: [`vignettes/example_workflow.Rmd`](vignettes/example_workflow.Rmd)
- More detailed guides and troubleshooting: [project wiki](https://github.com/rao-genomics-lab/ARMS/wiki) (in progress)

## License

GPL-3
