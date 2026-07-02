

#' Run full ASCAT.sc on pseudobulk tracks (mirrors run_sc_sequencing workflow)
#'
#' @param pseudobulk_tracks List of pseudobulk track data
#' @param ascat_rda_file Path to original ASCAT result
#' @param output_dir Output directory
#' @param purs Purity grid for ASCAT fitting (default: NULL, uses values from ascat_rda_file)
#' @param ploidies Ploidy grid for ASCAT fitting (default: NULL, uses values from ascat_rda_file)
#' @return ASCAT result object with profiles
#' @importFrom ASCAT.sc getTrackForAll searchGrid getProfile fitProfile predictRefit_all printResults_all
#' @importFrom parallel mclapply
#' @importFrom DNAcopy getbdry
#' @export
#' @examples
#' \dontrun{
#' # Use default purity/ploidy grids from original ASCAT run
#' pb_res <- run_ascat_on_pseudobulk_tracks(pb_tracks, "ascat.Rda", "output")
#'
#' # Specify custom purity/ploidy grids
#' pb_res <- run_ascat_on_pseudobulk_tracks(
#'   pb_tracks, "ascat.Rda", "output",
#'   purs = seq(0.1, 1, 0.01),
#'   ploidies = seq(1.5, 5, 0.01)
#' )
#' }
run_ascat_on_pseudobulk_tracks <- function(pseudobulk_tracks, ascat_rda_file, output_dir,
                                            purs = NULL, ploidies = NULL) {

  message("Running full ASCAT.sc on pseudobulk tracks...")

  library(ASCAT.sc)
  library(parallel)
  library(DNAcopy)

  # Load original ASCAT result to get parameters
  load(ascat_rda_file)  # Loads 'res'

  message(sprintf("Processing %d pseudobulk samples...", length(pseudobulk_tracks)))

  # Extract parameters from original run
  allchr <- res$chr
  chrstring_bam <- res$chrstring_bam
  sex <- res$sex[1]  # Use first sex value
  segmentation_alpha <- res$segmentation_alpha
  multipcf <- res$multipcf
  binsize <- res$binsize
  isPON <- res$isPON

  # Get reference data - at same resolution, use lSe and lGCT directly
  # (nlSe/nlGCT are only needed when normalizing to different binsize)
  lGCT <- res$lGCT
  lSe <- res$lSe

  # Get SBDRY for segmentation
  data("SBDRYs_precomputed", package = "ASCAT.sc")
  if (!as.character(segmentation_alpha) %in% names(SBDRYs)) {
    nperms <- 10000
    max.ones <- floor(nperms * segmentation_alpha) + 1
    SBDRY <- DNAcopy::getbdry(eta = 0.05, nperm = nperms, max.ones = max.ones)
  } else {
    SBDRY <- SBDRYs[[as.character(segmentation_alpha)]]
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Build allTracks structure matching run_sc_sequencing format
  allTracks <- list()
  for (i in seq_along(pseudobulk_tracks)) {
    cluster_name <- names(pseudobulk_tracks)[i]
    track <- pseudobulk_tracks[[cluster_name]]

    # Already at target binsize, so nlCTS.tumour = lCTS.tumour
    allTracks[[i]] <- list(
      lCTS.tumour = track$lCTS.tumour,
      nlCTS.tumour = track$lCTS.tumour  # Same since already at target binsize
    )
  }
  names(allTracks) <- names(pseudobulk_tracks)

  # Process tracks using getTrackForAll (same as run_sc_sequencing)
  message("Smoothing and segmenting tracks...")
  allTracks.processed <- mclapply(1:length(allTracks), function(x) {
    message(sprintf("  Processing %s...", names(allTracks)[x]))

    # Get the raw track data - list with one element per chromosome
    lCT_raw <- allTracks[[x]]$nlCTS.tumour

    # Convert list of lists to proper data frame format expected by getTrackForAll
    # Also trim to match lGCT length (lCTS may have 1 extra bin per chr)
    lCT_formatted <- lapply(seq_along(lCT_raw), function(i) {
      chr_data <- lCT_raw[[i]]
      chr_name <- allchr[i]
      target_length <- length(lGCT[[chr_name]])

      # Create data frame
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

    getTrackForAll(
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
  }, mc.cores = 1)
  names(allTracks.processed) <- names(allTracks)

  # Fit purity/ploidy
  message("Fitting purity/ploidy...")
  # Use provided purs/ploidies if specified, otherwise extract from original ASCAT result
  if (is.null(purs)) {
    # res$purs and res$ploidies are lists (one per sample), extract the grid from first element
    purs <- if (is.list(res$purs)) res$purs[[1]] else res$purs
    message(sprintf("Using purity grid from original ASCAT run: %d values (%.2f to %.2f)",
                    length(purs), min(purs), max(purs)))
  } else {
    message(sprintf("Using custom purity grid: %d values (%.2f to %.2f)",
                    length(purs), min(purs), max(purs)))
  }

  if (is.null(ploidies)) {
    ploidies <- if (is.list(res$ploidies)) res$ploidies[[1]] else res$ploidies
    message(sprintf("Using ploidy grid from original ASCAT run: %d values (%.2f to %.2f)",
                    length(ploidies), min(ploidies), max(ploidies)))
  } else {
    message(sprintf("Using custom ploidy grid: %d values (%.2f to %.2f)",
                    length(ploidies), min(ploidies), max(ploidies)))
  }

  maxtumourpsi <- res$maxtumourpsi

  allSols <- mclapply(1:length(allTracks.processed), function(x) {
    message(sprintf("  Fitting %s...", names(allTracks.processed)[x]))
    sol <- try(searchGrid(
      allTracks.processed[[x]],
      purs = purs,
      ploidies = ploidies,
      maxTumourPhi = maxtumourpsi,
      ismale = (sex == "male"),
      isPON = isPON
    ), silent = FALSE)
  }, mc.cores = 1)
  names(allSols) <- names(allTracks.processed)

  # Generate profiles
  message("Generating profiles...")
  allProfiles <- mclapply(1:length(allTracks.processed), function(x) {
    message(sprintf("  Profile for %s...", names(allTracks.processed)[x]))
    try(getProfile(
      fitProfile(
        allTracks.processed[[x]],
        purity = allSols[[x]]$purity,
        ploidy = allSols[[x]]$ploidy,
        ismale = (sex == "male")
      ),
      CHRS = allchr
    ), silent = FALSE)
  }, mc.cores = 1)
  names(allProfiles) <- names(allTracks.processed)

  # Create result object
  pb_res <- list(
    allTracks.processed = allTracks.processed,
    allTracks = allTracks,
    allSolutions = allSols,
    allProfiles = allProfiles,
    chr = allchr,
    chrstring_bam = chrstring_bam,
    purs = purs,
    ploidies = ploidies,
    maxtumourpsi = maxtumourpsi,
    build = res$build,
    binsize = binsize,
    sex = rep(sex, length(allTracks)),
    segmentation_alpha = segmentation_alpha,
    multipcf = multipcf,
    lSe = lSe,
    lGCT = lGCT,
    isPON = isPON,
    mode = "pseudobulk"
  )

  # Predict refit
  message("Predicting refit...")
  pb_res <- predictRefit_all(pb_res)

  # Print results (writes profile files)
  message("Writing output files...")
  pb_res <- printResults_all(pb_res, svinput = NULL, lSVinput = NULL,
                             outdir = output_dir, projectname = "pseudobulk")

  message(sprintf("\nFull ASCAT.sc analysis complete!"))
  message(sprintf("Created %d profiles in: %s", length(allProfiles), output_dir))

  return(pb_res)
}





#' Run full ASCAT.sc workflow on aggregated bulk track
#'
#' @param bulk_track Bulk track data (list with lCTS.tumour)
#' @param ascat_rda_file Path to original ASCAT result for parameters
#' @param output_dir Output directory for results
#' @param purs Purity grid for ASCAT fitting (default: NULL, uses values from ascat_rda_file)
#' @param ploidies Ploidy grid for ASCAT fitting (default: NULL, uses values from ascat_rda_file)
#' @return ASCAT result object with bulk profile
#' @importFrom ASCAT.sc getTrackForAll searchGrid getProfile fitProfile predictRefit_all printResults_all
#' @importFrom parallel mclapply
#' @importFrom DNAcopy getbdry
#' @export
#' @examples
#' \dontrun{
#' # Use default purity/ploidy grids
#' bulk_res <- run_ascat_on_bulk_track(bulk_track, "ascat.Rda", "output/")
#'
#' # Specify custom grids
#' bulk_res <- run_ascat_on_bulk_track(
#'   bulk_track, "ascat.Rda", "output/",
#'   purs = seq(0.1, 1, 0.01),
#'   ploidies = seq(1.5, 5, 0.01)
#' )
#' }
run_ascat_on_bulk_track <- function(bulk_track, ascat_rda_file, output_dir,
                                     purs = NULL, ploidies = NULL) {

  message("Running full ASCAT.sc on bulk track...")

  library(ASCAT.sc)
  library(parallel)
  library(DNAcopy)

  # Load original ASCAT result to get parameters
  load(ascat_rda_file)  # Loads 'res'

  # Extract parameters from original run
  allchr <- res$chr
  chrstring_bam <- res$chrstring_bam
  sex <- res$sex[1]
  segmentation_alpha <- res$segmentation_alpha
  multipcf <- res$multipcf
  binsize <- res$binsize
  isPON <- res$isPON

  # Get reference data
  lGCT <- res$lGCT
  lSe <- res$lSe

  # Get SBDRY for segmentation
  data("SBDRYs_precomputed", package = "ASCAT.sc")
  if (!as.character(segmentation_alpha) %in% names(SBDRYs)) {
    nperms <- 10000
    max.ones <- floor(nperms * segmentation_alpha) + 1
    SBDRY <- DNAcopy::getbdry(eta = 0.05, nperm = nperms, max.ones = max.ones)
  } else {
    SBDRY <- SBDRYs[[as.character(segmentation_alpha)]]
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Build allTracks structure
  allTracks <- list()
  allTracks[[1]] <- list(
    lCTS.tumour = bulk_track$lCTS.tumour,
    nlCTS.tumour = bulk_track$lCTS.tumour
  )
  names(allTracks) <- "bulk"

  # Process track using getTrackForAll
  message("Smoothing and segmenting bulk track...")

  # Get the raw track data
  lCT_raw <- allTracks[[1]]$nlCTS.tumour

  # Convert list of lists to proper data frame format
  lCT_formatted <- lapply(seq_along(lCT_raw), function(i) {
    chr_data <- lCT_raw[[i]]
    chr_name <- allchr[i]
    target_length <- length(lGCT[[chr_name]])

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

  allTracks.processed <- list()
  allTracks.processed[[1]] <- getTrackForAll(
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
  names(allTracks.processed) <- "bulk"

  # Fit purity/ploidy
  message("Fitting purity/ploidy...")
  # Use provided purs/ploidies if specified, otherwise extract from original ASCAT result
  if (is.null(purs)) {
    purs <- if (is.list(res$purs)) res$purs[[1]] else res$purs
    message(sprintf("Using purity grid from original ASCAT run: %d values (%.2f to %.2f)",
                    length(purs), min(purs), max(purs)))
  } else {
    message(sprintf("Using custom purity grid: %d values (%.2f to %.2f)",
                    length(purs), min(purs), max(purs)))
  }

  if (is.null(ploidies)) {
    ploidies <- if (is.list(res$ploidies)) res$ploidies[[1]] else res$ploidies
    message(sprintf("Using ploidy grid from original ASCAT run: %d values (%.2f to %.2f)",
                    length(ploidies), min(ploidies), max(ploidies)))
  } else {
    message(sprintf("Using custom ploidy grid: %d values (%.2f to %.2f)",
                    length(ploidies), min(ploidies), max(ploidies)))
  }

  maxtumourpsi <- res$maxtumourpsi

  allSols <- list()
  allSols[[1]] <- try(searchGrid(
    allTracks.processed[[1]],
    purs = purs,
    ploidies = ploidies,
    maxTumourPhi = maxtumourpsi,
    ismale = (sex == "male"),
    isPON = isPON
  ), silent = FALSE)
  names(allSols) <- "bulk"

  # Generate profile
  message("Generating bulk profile...")
  allProfiles <- list()
  allProfiles[[1]] <- try(getProfile(
    fitProfile(
      allTracks.processed[[1]],
      purity = allSols[[1]]$purity,
      ploidy = allSols[[1]]$ploidy,
      ismale = (sex == "male")
    ),
    CHRS = allchr
  ), silent = FALSE)
  names(allProfiles) <- "bulk"

  # Create result object
  bulk_res <- list(
    allTracks.processed = allTracks.processed,
    allTracks = allTracks,
    allSolutions = allSols,
    allProfiles = allProfiles,
    chr = allchr,
    chrstring_bam = chrstring_bam,
    purs = purs,
    ploidies = ploidies,
    maxtumourpsi = maxtumourpsi,
    build = res$build,
    binsize = binsize,
    sex = sex,
    segmentation_alpha = segmentation_alpha,
    multipcf = multipcf,
    lSe = lSe,
    lGCT = lGCT,
    isPON = isPON,
    mode = "bulk"
  )

  # Predict refit
  message("Predicting refit...")
  bulk_res <- predictRefit_all(bulk_res)

  # Print results (writes profile files)
  message("Writing output files...")
  bulk_res <- printResults_all(bulk_res, svinput = NULL, lSVinput = NULL,
                               outdir = output_dir, projectname = "bulk")

  message(sprintf("\nBulk ASCAT.sc analysis complete!"))
  message(sprintf("Bulk profile saved in: %s", output_dir))

  return(bulk_res)
}
