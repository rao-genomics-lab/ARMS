
#' Parse a tab-separated counts file into bins and count matrix
#'
#' Reads a file where rows are genomic bins (format \code{chr:start-end}) and
#' columns are samples.  Returns parsed bin coordinates, a count matrix, and
#' derived metadata.
#'
#' @param counts_file Path to the tab-separated counts file.
#' @return A list with elements \code{bins} (data frame with chr, start, end),
#'   \code{counts} (numeric matrix, bins x samples), \code{sample_names},
#'   \code{allchr} (unique chromosomes in genome order), and \code{binsize}.
#' @importFrom readr read_tsv
#' @keywords internal
parse_counts_file <- function(counts_file) {
  message("Reading counts file: ", counts_file)
  raw <- readr::read_tsv(counts_file, show_col_types = FALSE)

  # First column is "bin" (chr:start-end)
  bin_col <- raw[[1]]
  sample_names <- colnames(raw)[-1]
  counts_raw <- as.matrix(raw[, -1])
  storage.mode(counts_raw) <- "numeric"

  # Parse bin coordinates
  m <- regmatches(bin_col, regexec("^([^:]+):(\\d+)-(\\d+)$", bin_col))
  parsed <- do.call(rbind, lapply(m, function(x) {
    if (length(x) != 4) stop("Failed to parse bin: ", x[1])
    x[2:4]
  }))

  bins <- data.frame(
    chr = parsed[, 1],
    start = as.numeric(parsed[, 2]),
    end = as.numeric(parsed[, 3]),
    stringsAsFactors = FALSE
  )

  # Sort by chromosome (numeric order, then X, Y) and start position
  chr_order <- c(as.character(1:22), "X", "Y")
  bins$chr_rank <- match(bins$chr, chr_order)
  if (any(is.na(bins$chr_rank))) {
    # Try stripping "chr" prefix
    bins$chr_rank <- match(sub("^chr", "", bins$chr), chr_order)
  }
  sort_idx <- order(bins$chr_rank, bins$start)
  bins <- bins[sort_idx, ]
  bins$chr_rank <- NULL
  counts_raw <- counts_raw[sort_idx, , drop = FALSE]
  rownames(bins) <- NULL

  # Detect binsize from most common bin width
  widths <- bins$end - bins$start
  binsize <- as.numeric(names(sort(table(widths), decreasing = TRUE))[1])
  message(sprintf("Detected binsize: %d bp", binsize))

  # Get unique chromosomes in genome order
  allchr <- unique(bins$chr)
  message(sprintf("Parsed %d bins across %d chromosomes for %d samples",
                  nrow(bins), length(allchr), length(sample_names)))

  list(
    bins = bins,
    counts = counts_raw,
    sample_names = sample_names,
    allchr = allchr,
    binsize = binsize
  )
}


#' Compute GC content for genomic bins using a BSgenome object
#'
#' @param lSe Named list of bin coordinates per chromosome, each with
#'   \code{$starts} and \code{$ends} vectors.
#' @param genome A BSgenome object.
#' @param allchr Character vector of chromosome names (NCBI-style, e.g. "1").
#' @param chr_prefix Prefix to prepend when accessing BSgenome sequences
#'   (default \code{"chr"}).
#' @return Named list of numeric vectors (GC fractions) per chromosome.
#' @importFrom Biostrings letterFrequency Views
#' @importFrom BSgenome getSeq
#' @keywords internal
compute_gc_for_bins <- function(lSe, genome, allchr, chr_prefix = "chr") {
  message("Computing GC content for bins...")
  lGCT <- list()

  for (chr in allchr) {
    genome_chr_name <- paste0(chr_prefix, chr)
    chr_seq <- genome[[genome_chr_name]]
    chr_len <- length(chr_seq)

    starts <- lSe[[chr]]$starts
    ends <- pmin(lSe[[chr]]$ends, chr_len)

    v <- Biostrings::Views(chr_seq, start = starts, end = ends)
    gc_counts <- Biostrings::letterFrequency(v, letters = c("G", "C"))
    total <- Biostrings::letterFrequency(v, letters = c("A", "C", "G", "T"))
    gc_frac <- rowSums(gc_counts) / rowSums(total)
    gc_frac[!is.finite(gc_frac)] <- 0.5  # fallback for any zero-length bins

    lGCT[[chr]] <- gc_frac
  }

  message(sprintf("Computed GC content for %d chromosomes", length(lGCT)))
  lGCT
}


#' Create ASCAT.sc result object from a raw counts file
#'
#' Takes a tab-separated file of pre-binned read counts (bins as rows, samples
#' as columns) and runs the full ASCAT.sc copy number calling pipeline.
#' This enables users who have already binned their data (e.g. at 500kb
#' resolution) to use ARMS without BAM files.
#'
#' @param counts_file Path to tab-separated counts file.  First column
#'   \code{"bin"} with format \code{chr:start-end} (e.g. \code{1:1-500000}),
#'   remaining columns are sample counts.
#' @param output_dir Directory for output profile files and the saved result
#'   object.
#' @param genome A \code{BSgenome} object or the name of a BSgenome package
#'   (e.g. \code{"BSgenome.Hsapiens.UCSC.hg19"}).
#' @param build Genome build identifier stored in the result object
#'   (default \code{"hg19"}).
#' @param sex Sex of all samples: \code{"male"} or \code{"female"}
#'   (default \code{"male"}).
#' @param binsize Bin size in bp.  If \code{NULL} (default), auto-detected from
#'   the counts file.
#' @param purs Purity grid for ASCAT fitting
#'   (default \code{seq(0.1, 1, 0.01)}).
#' @param ploidies Ploidy grid for ASCAT fitting
#'   (default \code{seq(1.5, 5.5, 0.01)}).
#' @param maxtumourpsi Maximum tumour psi for fitting (default 5).
#' @param segmentation_alpha Alpha for CBS segmentation (default 0.01).
#' @param projectname Project name used when writing output files
#'   (default \code{"counts"}).
#' @param MC.CORES Number of cores for parallel processing (default 1).
#' @return An ASCAT.sc result object (list) with profiles, solutions, and
#'   tracks.  Profile files are written to \code{output_dir}.
#' @importFrom ASCAT.sc getTrackForAll searchGrid getProfile fitProfile predictRefit_all printResults_all
#' @importFrom parallel mclapply
#' @importFrom DNAcopy getbdry
#' @importFrom readr read_tsv
#' @importFrom Biostrings letterFrequency Views
#' @importFrom BSgenome getSeq
#' @export
#' @examples
#' \dontrun{
#' res <- create_res_from_counts(
#'   counts_file = "binned_counts.tsv",
#'   output_dir = "ascat_output",
#'   genome = "BSgenome.Hsapiens.UCSC.hg19",
#'   build = "hg19",
#'   sex = "male",
#'   projectname = "my_project"
#' )
#' # Load resulting profiles into ARMS
#' profiles <- load_ascat_profiles("ascat_output/")
#' }
create_res_from_counts <- function(counts_file,
                                   output_dir,
                                   genome,
                                   build = "hg19",
                                   sex = "male",
                                   binsize = NULL,
                                   purs = seq(0.1, 1, 0.01),
                                   ploidies = seq(1.5, 5.5, 0.01),
                                   maxtumourpsi = 5,
                                   segmentation_alpha = 0.01,
                                   projectname = "counts",
                                   MC.CORES = 1) {

  message("=== ARMS: Creating ASCAT.sc result from raw counts ===")

  library(ASCAT.sc)
  library(parallel)
  library(DNAcopy)

  # --- Load genome ---
  if (is.character(genome)) {
    message(sprintf("Loading genome package: %s", genome))
    if (!requireNamespace(genome, quietly = TRUE)) {
      stop("Genome package '", genome, "' is not installed. ",
           "Install with: BiocManager::install('", genome, "')")
    }
    genome <- get(genome, envir = asNamespace(genome))
  }

  # --- Parse counts file ---
  parsed <- parse_counts_file(counts_file)
  allchr <- parsed$allchr
  if (is.null(binsize)) {
    binsize <- parsed$binsize
  }
  message(sprintf("Using binsize: %d bp", binsize))

  # --- Build lSe (bin coordinates per chromosome) ---
  lSe <- list()
  for (chr in allchr) {
    chr_bins <- parsed$bins[parsed$bins$chr == chr, ]
    lSe[[chr]] <- list(
      starts = chr_bins$start,
      ends = chr_bins$end
    )
  }

  # --- Compute GC content ---
  # Determine chr_prefix: check if genome uses "chr" prefix
  genome_seqnames <- if (requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    GenomeInfoDb::seqnames(genome)
  } else {
    names(genome)
  }
  if (any(grepl("^chr", genome_seqnames))) {
    chr_prefix <- "chr"
  } else {
    chr_prefix <- ""
  }
  lGCT <- compute_gc_for_bins(lSe, genome, allchr, chr_prefix = chr_prefix)

  # --- Build allTracks ---
  message("Building track structures for all samples...")
  n_samples <- length(parsed$sample_names)
  allTracks <- list()

  for (s in seq_len(n_samples)) {
    lCTS <- list()
    for (chr in allchr) {
      chr_bins <- parsed$bins[parsed$bins$chr == chr, ]
      chr_idx <- which(parsed$bins$chr == chr)
      records <- parsed$counts[chr_idx, s]

      lCTS[[chr]] <- list(
        space = rep(chr, nrow(chr_bins)),
        start = chr_bins$start,
        end = chr_bins$end,
        width = chr_bins$end - chr_bins$start,
        file = parsed$sample_names[s],
        records = records,
        nucleotides = records  # not used in core analysis
      )
    }
    allTracks[[s]] <- list(
      lCTS.tumour = lCTS,
      nlCTS.tumour = lCTS  # already at target resolution
    )
  }
  names(allTracks) <- parsed$sample_names
  message(sprintf("Built tracks for %d samples", n_samples))

  # --- Get SBDRY for segmentation ---
  data("SBDRYs_precomputed", package = "ASCAT.sc")
  if (!as.character(segmentation_alpha) %in% names(SBDRYs)) {
    nperms <- 10000
    max.ones <- floor(nperms * segmentation_alpha) + 1
    SBDRY <- DNAcopy::getbdry(eta = 0.05, nperm = nperms, max.ones = max.ones)
  } else {
    SBDRY <- SBDRYs[[as.character(segmentation_alpha)]]
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Process tracks (GC correction, smoothing, segmentation) ---
  message("Processing tracks (GC correction + segmentation)...")
  allTracks.processed <- mclapply(seq_len(n_samples), function(x) {
    message(sprintf("  Processing %s (%d/%d)...", names(allTracks)[x], x, n_samples))

    lCT_raw <- allTracks[[x]]$nlCTS.tumour

    # Convert to data frame format expected by getTrackForAll
    lCT_formatted <- lapply(allchr, function(chr) {
      chr_data <- lCT_raw[[chr]]
      target_length <- length(lGCT[[chr]])

      df <- data.frame(
        space = chr_data$space,
        start = chr_data$start,
        end = chr_data$end,
        width = chr_data$width,
        file = chr_data$file,
        records = chr_data$records,
        nucleotides = chr_data$nucleotides,
        stringsAsFactors = FALSE
      )

      # Trim to match lGCT length if needed
      if (nrow(df) > target_length) {
        df <- df[1:target_length, ]
      }

      df
    })
    names(lCT_formatted) <- allchr

    ASCAT.sc::getTrackForAll(
      bamfile = NULL,
      window = NULL,
      lCT = lCT_formatted,
      lSe = lSe,
      lGCT = lGCT,
      lNormals = NULL,
      allchr = allchr,
      sdNormalise = 0,
      SBDRY = SBDRY,
      svinput = NULL,
      segmentation_alpha = segmentation_alpha
    )
  }, mc.cores = MC.CORES)
  names(allTracks.processed) <- names(allTracks)

  # --- Fit purity/ploidy ---
  message("Fitting purity/ploidy...")
  message(sprintf("  Purity grid: %d values (%.2f to %.2f)", length(purs), min(purs), max(purs)))
  message(sprintf("  Ploidy grid: %d values (%.2f to %.2f)", length(ploidies), min(ploidies), max(ploidies)))

  allSols <- mclapply(seq_len(n_samples), function(x) {
    message(sprintf("  Fitting %s (%d/%d)...", names(allTracks.processed)[x], x, n_samples))
    sol <- try(ASCAT.sc::searchGrid(
      allTracks.processed[[x]],
      purs = purs,
      ploidies = ploidies,
      maxTumourPhi = maxtumourpsi,
      ismale = (sex == "male"),
      isPON = FALSE
    ), silent = FALSE)
    sol
  }, mc.cores = MC.CORES)
  names(allSols) <- names(allTracks.processed)

  # --- Generate profiles ---
  message("Generating profiles...")
  allProfiles <- mclapply(seq_len(n_samples), function(x) {
    message(sprintf("  Profile for %s (%d/%d)...", names(allTracks.processed)[x], x, n_samples))
    try(ASCAT.sc::getProfile(
      ASCAT.sc::fitProfile(
        allTracks.processed[[x]],
        purity = allSols[[x]]$purity,
        ploidy = allSols[[x]]$ploidy,
        ismale = (sex == "male")
      ),
      CHRS = allchr
    ), silent = FALSE)
  }, mc.cores = MC.CORES)
  names(allProfiles) <- names(allTracks.processed)

  # --- Compile result object ---
  message("Compiling result object...")
  res <- list(
    allTracks.processed = allTracks.processed,
    allTracks = allTracks,
    allSolutions = allSols,
    allProfiles = allProfiles,
    chr = allchr,
    chrstring_bam = allchr,  # same as allchr for counts-based input
    purs = purs,
    ploidies = ploidies,
    maxtumourpsi = maxtumourpsi,
    build = build,
    binsize = binsize,
    sex = rep(sex, n_samples),
    segmentation_alpha = segmentation_alpha,
    multipcf = FALSE,
    lSe = lSe,
    lGCT = lGCT,
    isPON = FALSE,
    mode = "counts"
  )

  # --- Predict refit ---
  message("Predicting refit...")
  res <- ASCAT.sc::predictRefit_all(res)

  # --- Write profile files ---
  message("Writing output files...")
  res <- ASCAT.sc::printResults_all(res, svinput = NULL, lSVinput = NULL,
                                     outdir = output_dir, projectname = projectname)

  # --- Save result object ---
  rda_file <- file.path(output_dir, paste0(projectname, "_ascat_result.Rda"))
  save(res, file = rda_file)
  message(sprintf("Saved result object to: %s", rda_file))

  # --- Summary ---
  n_success <- sum(!sapply(allProfiles, inherits, "try-error"))
  message(sprintf("\n=== ASCAT.sc analysis complete ==="))
  message(sprintf("  Samples processed: %d", n_samples))
  message(sprintf("  Profiles generated: %d", n_success))
  message(sprintf("  Output directory: %s", output_dir))

  return(res)
}
