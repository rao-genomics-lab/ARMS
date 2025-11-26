#' Create pseudobulk profiles from bin counts (fast, no BAM merging)
#'
#' @param ascat_rda_file Path to ASCAT.sc result Rda file
#' @param clusters Cluster assignments
#' @param sample_names Sample names matching cluster assignments
#' @param output_dir Output directory for pseudobulk profiles
#' @param min_samples Minimum samples per cluster (default: 3)
#' @return List of pseudobulk track data
#' @export
#' @examples
#' \dontrun{
#' pb_tracks <- create_pseudobulk_from_bins("ascat.Rda", clusters, sample_names, "output")
#' }
create_pseudobulk_from_bins <- function(ascat_rda_file, clusters, sample_names,
                                        output_dir, min_samples = 3) {

  message("Creating pseudobulk profiles from bin counts (fast mode)...")

  # Load ASCAT result
  load(ascat_rda_file)  # Loads 'res' object

  if (!exists("res")) {
    stop("ASCAT result object 'res' not found in ", ascat_rda_file)
  }

  # Get track structure
  if (!"allTracks" %in% names(res)) {
    stop("allTracks not found in ASCAT result")
  }

  n_tracks <- length(res$allTracks)
  message(sprintf("Found %d tracks (samples) in ASCAT result", n_tracks))

  # Extract sample names from tracks
  track_sample_names <- sapply(seq_along(res$allTracks), function(i) {
    if (length(res$allTracks[[i]]$lCTS.tumour) > 0) {
      return(res$allTracks[[i]]$lCTS.tumour[[1]]$file[1])
    }
    return(NA)
  })

  # Clean track sample names
  track_sample_names <- gsub("\\.bam$", "", track_sample_names)

  # Clean input sample names
  sample_names_clean <- gsub("\\.bam_ascat_profiles$", "", sample_names)
  sample_names_clean <- gsub("_ascat_profiles$", "", sample_names_clean)
  sample_names_clean <- gsub("\\.bam$", "", sample_names_clean)

  message(sprintf("First 5 track samples: %s", paste(head(track_sample_names), collapse = ", ")))
  message(sprintf("First 5 input samples (original): %s", paste(head(sample_names), collapse = ", ")))
  message(sprintf("First 5 input samples (cleaned): %s", paste(head(sample_names_clean), collapse = ", ")))

  # Match tracks to sample names
  sample_to_track <- match(sample_names_clean, track_sample_names)

  if (all(is.na(sample_to_track))) {
    stop("Could not match any sample names to tracks")
  }

  matched_count <- sum(!is.na(sample_to_track))
  message(sprintf("Matched %d/%d samples to tracks", matched_count, length(sample_names)))

  # Create pseudobulk by summing bin counts per cluster
  cluster_ids <- sort(unique(clusters))
  pseudobulk_tracks <- list()

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (cluster_id in cluster_ids) {
    cluster_mask <- clusters == cluster_id
    cluster_sample_names <- sample_names[cluster_mask]
    cluster_track_ids <- sample_to_track[cluster_mask]
    cluster_track_ids <- cluster_track_ids[!is.na(cluster_track_ids)]

    if (length(cluster_track_ids) < min_samples) {
      message(sprintf("Skipping cluster %d (only %d samples)",
                      cluster_id, length(cluster_track_ids)))
      next
    }

    message(sprintf("Processing cluster %d with %d samples...", cluster_id, length(cluster_track_ids)))

    # Sum bin counts across all chromosomes
    pseudobulk_track <- list()
    pseudobulk_track$lCTS.tumour <- list()

    # Initialize with first track structure
    first_track <- res$allTracks[[cluster_track_ids[1]]]

    for (chr_idx in seq_along(first_track$lCTS.tumour)) {
      # Start with first sample's bins
      chr_data <- first_track$lCTS.tumour[[chr_idx]]
      summed_records <- chr_data$records

      # Add counts from remaining samples
      for (track_id in cluster_track_ids[-1]) {
        track <- res$allTracks[[track_id]]
        summed_records <- summed_records + track$lCTS.tumour[[chr_idx]]$records
      }

      # Create summed track for this chromosome
      pseudobulk_track$lCTS.tumour[[chr_idx]] <- list(
        space = chr_data$space,
        start = chr_data$start,
        end = chr_data$end,
        width = chr_data$width,
        file = sprintf("cluster_%d", cluster_id),
        records = summed_records,
        nucleotides = chr_data$nucleotides
      )
    }

    # Store pseudobulk track
    pseudobulk_tracks[[sprintf("cluster_%d", cluster_id)]] <- pseudobulk_track
    message(sprintf("  Created pseudobulk for cluster %d", cluster_id))
  }

  message(sprintf("Created %d pseudobulk profiles", length(pseudobulk_tracks)))

  # Save pseudobulk tracks
  pseudobulk_file <- file.path(output_dir, "pseudobulk_tracks.Rda")
  save(pseudobulk_tracks, file = pseudobulk_file)
  message(sprintf("Saved pseudobulk tracks to: %s", pseudobulk_file))

  return(pseudobulk_tracks)
}

#' Merge BAM files by cluster
#'
#' @param bam_dir Directory containing BAM files
#' @param sample_names Sample names
#' @param clusters Cluster assignments
#' @param output_dir Output directory for merged BAMs
#' @param min_samples Minimum samples per cluster
#' @return Character vector of merged BAM file paths
#' @export
#' @examples
#' \dontrun{
#' merged_bams <- merge_bams_by_cluster("bams/", sample_names, clusters, "merged/")
#' }
merge_bams_by_cluster <- function(bam_dir, sample_names, clusters,
                                  output_dir, min_samples = 3) {

  message("Merging BAM files by cluster...")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cluster_ids <- sort(unique(clusters))
  merged_bams <- character()

  for (cluster_id in cluster_ids) {
    cluster_samples <- sample_names[clusters == cluster_id]

    if (length(cluster_samples) < min_samples) {
      message(sprintf("Skipping cluster %d (only %d samples)",
                      cluster_id, length(cluster_samples)))
      next
    }

    # Find BAM files - strip project suffix
    # Sample names may be like: "plate13-11A-12H_lcmad_000_dupmarked.bam_p15_analysis"
    # BAM files are like: "plate13-11A-12H_lcmad_000_dupmarked.bam"
    bam_names <- sub("(\\.bam).*", "\\1", cluster_samples)  # Keep everything up to and including .bam

    cluster_bams <- file.path(bam_dir, bam_names)
    cluster_bams <- cluster_bams[file.exists(cluster_bams)]

    if (length(cluster_bams) == 0) {
      message(sprintf("No BAM files found for cluster %d", cluster_id))
      next
    }

    # Merge using samtools
    output_bam <- file.path(output_dir, sprintf("cluster_%d.bam", cluster_id))

    cmd <- sprintf("samtools merge -f %s %s", output_bam,
                   paste(cluster_bams, collapse = " "))

    message(sprintf("Merging %d BAMs for cluster %d...",
                    length(cluster_bams), cluster_id))

    ret <- system(cmd, wait = TRUE)
    if (ret != 0) {
      stop(sprintf("samtools merge failed for cluster %d (exit code: %d)", cluster_id, ret))
    }
    message(sprintf("  ✓ Cluster %d merged", cluster_id))

    # Index
    ret <- system(sprintf("samtools index %s", output_bam), wait = TRUE)
    if (ret != 0) {
      stop(sprintf("samtools index failed for cluster %d (exit code: %d)", cluster_id, ret))
    }
    message(sprintf("  ✓ Cluster %d indexed", cluster_id))

    merged_bams <- c(merged_bams, output_bam)
  }

  message(sprintf("Created %d merged BAM files", length(merged_bams)))

  return(merged_bams)
}




#' Create bulk profile by aggregating all single-cell bin counts
#'
#' @param ascat_rda_file Path to ASCAT result Rda file
#' @param output_dir Output directory
#' @return Bulk track data (list with lCTS.tumour)
#' @export
#' @examples
#' \dontrun{
#' bulk_track <- create_bulk_from_bins("ascat.Rda", "output/")
#' }
create_bulk_from_bins <- function(ascat_rda_file, output_dir) {

  message("Creating bulk profile from all sample bin counts...")

  # Load ASCAT result
  load(ascat_rda_file)  # Loads 'res' object

  if (!exists("res")) {
    stop("ASCAT result object 'res' not found in ", ascat_rda_file)
  }

  n_tracks <- length(res$allTracks)
  message(sprintf("Merging %d samples into bulk profile", n_tracks))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Initialize with first track structure
  first_track <- res$allTracks[[1]]
  bulk_track <- list()
  bulk_track$lCTS.tumour <- list()

  for (chr_idx in seq_along(first_track$lCTS.tumour)) {
    # Start with first sample's bins
    chr_data <- first_track$lCTS.tumour[[chr_idx]]
    summed_records <- chr_data$records

    # Add counts from all remaining samples
    for (track_id in 2:n_tracks) {
      track <- res$allTracks[[track_id]]
      summed_records <- summed_records + track$lCTS.tumour[[chr_idx]]$records
    }

    # Create summed track for this chromosome
    bulk_track$lCTS.tumour[[chr_idx]] <- list(
      space = chr_data$space,
      start = chr_data$start,
      end = chr_data$end,
      width = chr_data$width,
      file = "bulk",
      records = summed_records,
      nucleotides = chr_data$nucleotides
    )
  }

  # Save bulk track
  bulk_file <- file.path(output_dir, "bulk_track.Rda")
  save(bulk_track, file = bulk_file)
  message(sprintf("Saved bulk track to: %s", bulk_file))

  return(bulk_track)
}


