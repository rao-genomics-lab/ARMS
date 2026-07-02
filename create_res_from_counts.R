#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  # Keep startup noise down for batch runs
  NULL
}))

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript create_res_from_counts.R <counts_file> <output_dir> --genome <BSgenomePkg> [options]\n\n",
    "Required:\n",
    "  counts_file   Tab-separated counts file (first column: chr:start-end)\n",
    "  output_dir    Output directory for ASCAT.sc profiles and result .Rda\n",
    "  --genome      BSgenome package name (e.g., BSgenome.Hsapiens.UCSC.hg19)\n\n",
    "Options:\n",
    "  --build        Genome build label (default: hg19)\n",
    "  --sex          Sex for all samples: male|female (default: male)\n",
    "  --binsize      Bin size in bp (default: auto-detect from counts)\n",
    "  --projectname  Project name used in outputs (default: counts)\n",
    "  --cores        Number of cores (default: 1)\n",
    "  --seg-alpha    CBS segmentation alpha (default: 0.01)\n",
    "  --purs         Purity grid, comma-separated (default: 0.1..1 by 0.01)\n",
    "  --ploidies     Ploidy grid, comma-separated (default: 1.5..5.5 by 0.01)\n\n",
    "Example:\n",
    "  Rscript create_res_from_counts.R counts.tsv ascat_out \\\n",
    "    --genome BSgenome.Hsapiens.UCSC.hg19 --build hg19 --sex male --cores 8\n\n",
    sep = ""
  )
}

if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  usage()
  quit(status = 1)
}

counts_file <- NULL
output_dir <- NULL
genome <- NULL
build <- "hg19"
sex <- "male"
binsize <- NA_real_
projectname <- "counts"
cores <- 1
seg_alpha <- 0.01
purs <- seq(0.1, 1, 0.01)
ploidies <- seq(1.5, 5.5, 0.01)

idx <- 1
while (idx <= length(args)) {
  arg <- args[[idx]]
  if (is.null(counts_file) && !startsWith(arg, "--")) {
    counts_file <- arg
    idx <- idx + 1
    next
  }
  if (is.null(output_dir) && !startsWith(arg, "--")) {
    output_dir <- arg
    idx <- idx + 1
    next
  }

  if (arg == "--genome") {
    genome <- args[[idx + 1]]
    idx <- idx + 2
  } else if (arg == "--build") {
    build <- args[[idx + 1]]
    idx <- idx + 2
  } else if (arg == "--sex") {
    sex <- args[[idx + 1]]
    idx <- idx + 2
  } else if (arg == "--binsize") {
    binsize <- as.numeric(args[[idx + 1]])
    idx <- idx + 2
  } else if (arg == "--projectname") {
    projectname <- args[[idx + 1]]
    idx <- idx + 2
  } else if (arg == "--cores") {
    cores <- as.integer(args[[idx + 1]])
    idx <- idx + 2
  } else if (arg == "--seg-alpha") {
    seg_alpha <- as.numeric(args[[idx + 1]])
    idx <- idx + 2
  } else if (arg == "--purs") {
    purs <- as.numeric(strsplit(args[[idx + 1]], ",", fixed = TRUE)[[1]])
    idx <- idx + 2
  } else if (arg == "--ploidies") {
    ploidies <- as.numeric(strsplit(args[[idx + 1]], ",", fixed = TRUE)[[1]])
    idx <- idx + 2
  } else {
    stop("Unknown argument: ", arg)
  }
}

if (is.null(counts_file) || is.null(output_dir) || is.null(genome)) {
  usage()
  stop("counts_file, output_dir, and --genome are required")
}

if (!file.exists(counts_file)) {
  stop("Counts file not found: ", counts_file)
}

if (!is.finite(cores) || cores < 1) {
  stop("--cores must be a positive integer")
}

if (!is.na(binsize) && (!is.finite(binsize) || binsize <= 0)) {
  stop("--binsize must be a positive number if set")
}

if (!sex %in% c("male", "female")) {
  stop("--sex must be 'male' or 'female'")
}

# Resolve project root from script location when possible
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1])))
  }
  NA_character_
}

script_path <- get_script_path()
project_root <- if (!is.na(script_path)) dirname(script_path) else getwd()

local_counts_path <- file.path(project_root, "R", "counts-to-ascat.R")
if (file.exists(local_counts_path)) {
  # Prefer local source to pick up recent fixes
  source(local_counts_path)
} else if (requireNamespace("ARMS", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ARMS))
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(project_root, quiet = TRUE)
} else {
  stop("ARMS not available and local source not found at: ", local_counts_path)
}

message("Running create_res_from_counts...")
res <- create_res_from_counts(
  counts_file = counts_file,
  output_dir = output_dir,
  genome = genome,
  build = build,
  sex = sex,
  binsize = if (is.na(binsize)) NULL else binsize,
  purs = purs,
  ploidies = ploidies,
  segmentation_alpha = seg_alpha,
  projectname = projectname,
  MC.CORES = cores
)

res_path <- file.path(output_dir, paste0(projectname, "_res.Rda"))
save(res, file = res_path)
message("Saved result object: ", res_path)
