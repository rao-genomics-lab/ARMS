#!/usr/bin/env Rscript

suppressWarnings(suppressMessages({
  # Keep startup noise down for batch runs
  NULL
}))

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript heatmap_from_counts.R <counts_file> [--output <pdf>] [--pseudocount <n>] [--width <in>] [--height <in>] [--logr-range <n>]\n\n",
    "Arguments:\n",
    "  counts_file   Tab-separated counts file (first column: chr:start-end)\n\n",
    "Options:\n",
    "  --output       Output PDF path (default: <counts_file>_heatmap.pdf)\n",
    "  --pseudocount  Pseudocount added before log2 (default: 1)\n",
    "  --width        Heatmap width in inches (default: 20)\n",
    "  --height       Heatmap height in inches (default: auto by sample count)\n",
    "  --logr-range   Color range for LogR (default: 99th percentile of |logR|)\n\n",
    sep = ""
  )
}

if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  usage()
  quit(status = 1)
}

# Simple flag parser
counts_file <- NULL
output_file <- NULL
pseudocount <- 1
width <- 20
height <- NA_real_
logr_range <- NA_real_

idx <- 1
while (idx <= length(args)) {
  arg <- args[[idx]]
  if (is.null(counts_file) && !startsWith(arg, "--")) {
    counts_file <- arg
    idx <- idx + 1
    next
  }

  if (arg == "--output") {
    output_file <- args[[idx + 1]]
    idx <- idx + 2
  } else if (arg == "--pseudocount") {
    pseudocount <- as.numeric(args[[idx + 1]])
    idx <- idx + 2
  } else if (arg == "--width") {
    width <- as.numeric(args[[idx + 1]])
    idx <- idx + 2
  } else if (arg == "--height") {
    height <- as.numeric(args[[idx + 1]])
    idx <- idx + 2
  } else if (arg == "--logr-range") {
    logr_range <- as.numeric(args[[idx + 1]])
    idx <- idx + 2
  } else {
    stop("Unknown argument: ", arg)
  }
}

if (is.null(counts_file)) {
  usage()
  stop("counts_file is required")
}

if (!file.exists(counts_file)) {
  stop("Counts file not found: ", counts_file)
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

# Load ARMS functions
parse_counts <- NULL
plot_heatmap <- NULL

if (requireNamespace("ARMS", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ARMS))
  # parse_counts_file is internal, so access via namespace
  parse_counts <- getFromNamespace("parse_counts_file", "ARMS")
  plot_heatmap <- getFromNamespace("plot_logr_heatmap", "ARMS")
} else if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(project_root, quiet = TRUE)
  parse_counts <- get("parse_counts_file", envir = .GlobalEnv)
  plot_heatmap <- get("plot_logr_heatmap", envir = .GlobalEnv)
} else {
  source(file.path(project_root, "R", "counts-to-ascat.R"))
  source(file.path(project_root, "R", "visualization.R"))
  parse_counts <- get("parse_counts_file", envir = .GlobalEnv)
  plot_heatmap <- get("plot_logr_heatmap", envir = .GlobalEnv)
}

if (is.null(parse_counts) || is.null(plot_heatmap)) {
  stop("Failed to load ARMS functions. Ensure ARMS is installed or run from the repo root.")
}

message("Parsing counts file...")
parsed <- parse_counts(counts_file)

counts <- parsed$counts
if (!is.finite(pseudocount) || pseudocount < 0) {
  stop("pseudocount must be a non-negative number")
}

# Build grid for chromosome annotations (one row per bin)
grid <- data.frame(
  chr = parsed$bins$chr,
  start = parsed$bins$start,
  end = parsed$bins$end,
  pos = (parsed$bins$start + parsed$bins$end) / 2,
  stringsAsFactors = FALSE
)

# Convert counts to LogR-like matrix: log2(counts + pseudocount), centered per sample
log_counts <- log2(counts + pseudocount)
logr_matrix <- t(log_counts)
rownames(logr_matrix) <- parsed$sample_names

row_medians <- apply(logr_matrix, 1, median, na.rm = TRUE)
logr_matrix <- logr_matrix - row_medians

if (is.na(logr_range)) {
  logr_range <- as.numeric(stats::quantile(abs(logr_matrix), 0.99, na.rm = TRUE))
  if (!is.finite(logr_range) || logr_range <= 0) {
    logr_range <- 2
  }
}

logr_range <- 2 # manually setting logr range

n_samples <- nrow(logr_matrix)
if (is.na(height)) {
  height <- max(6, min(0.2 * n_samples, 60))
}

if (is.null(output_file)) {
  output_file <- paste0(tools::file_path_sans_ext(basename(counts_file)), "_heatmap.pdf")
}

clusters <- rep(1, n_samples)

message(sprintf("Generating heatmap for %d samples and %d bins", n_samples, ncol(logr_matrix)))
plot_heatmap(
  logr_matrix = logr_matrix,
  grid = grid,
  clusters = clusters,
  output_file = output_file,
  logr_range = logr_range,
  width = width,
  height = height
)

message("Done. Output written to: ", output_file)
