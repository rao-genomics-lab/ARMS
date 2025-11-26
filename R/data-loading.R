#' Load ASCAT.sc profile files and create matrices
#'
#' @param profile_dir Directory containing profile files
#' @param pattern File pattern (default: "*refitted.ASCAT.scprofile.txt")
#' @param resolution Genomic resolution in bp
#' @param logr_missing_value Value for missing logR data
#' @return List with copy_number_matrix, logr_matrix, grid, sample_names
#' @importFrom readr read_tsv
#' @importFrom dplyr %>%
#' @importFrom stringr str_replace
#' @export
#' @examples
#' \dontrun{
#' profiles <- load_ascat_profiles("path/to/profiles", resolution = 1e6)
#' }
load_ascat_profiles <- function(profile_dir,
                                pattern = "*refitted.ASCAT.scprofile.txt",
                                resolution = 1000000,
                                logr_missing_value = 0) {

  message("Loading ASCAT.sc profile files...")

  # Find profile files
  profile_files <- list.files(profile_dir, pattern = glob2rx(pattern), full.names = TRUE)

  if (length(profile_files) == 0) {
    # Try non-refitted profiles
    pattern <- "*ASCAT.scprofile.txt"
    profile_files <- list.files(profile_dir, pattern = glob2rx(pattern), full.names = TRUE)
  }

  if (length(profile_files) == 0) {
    stop("No profile files found in ", profile_dir)
  }

  # Sort files numerically by cluster ID if they are cluster files
  file_basenames <- basename(profile_files)
  cluster_ids <- as.numeric(gsub(".*cluster_(\\d+).*", "\\1", file_basenames))
  if (!all(is.na(cluster_ids))) {
    profile_files <- profile_files[order(cluster_ids)]
  }

  message(sprintf("Found %d profile files", length(profile_files)))

  # Read first file to get chromosomes
  first_profile <- read_tsv(profile_files[1], show_col_types = FALSE)
  chromosomes <- paste0("chr", c(1:22, "X"))

  # Create genomic grid
  grid <- create_genomic_grid(first_profile, resolution)
  message(sprintf("Created genomic grid with %d positions", nrow(grid)))

  # Initialize matrices
  n_samples <- length(profile_files)
  n_positions <- nrow(grid)

  copy_number_matrix <- matrix(NA, nrow = n_samples, ncol = n_positions)
  logr_matrix <- matrix(NA, nrow = n_samples, ncol = n_positions)

  sample_names <- basename(profile_files) %>%
    str_replace("_refitted.ASCAT.scprofile.txt", "") %>%
    str_replace(".ASCAT.scprofile.txt", "")

  rownames(copy_number_matrix) <- sample_names
  rownames(logr_matrix) <- sample_names

  # Load each profile - iterate through segments (LEGACY APPROACH for speed!)
  # Key insight: ~150 segments per profile vs ~3000 grid positions
  # So iterating through segments and finding overlapping grid positions is much faster
  for (i in seq_along(profile_files)) {
    if (i %% 50 == 0) message(sprintf("Processing file %d/%d", i, n_samples))

    profile <- read_tsv(profile_files[i], show_col_types = FALSE)

    if (nrow(profile) == 0) {
      warning(sprintf("Empty profile file: %s", profile_files[i]))
      next
    }

    # Handle different column naming conventions
    if ("chromosome" %in% colnames(profile)) {
      # New format: chromosome, start, end, total_copy_number, logr
      # ITERATE THROUGH SEGMENTS (not grid positions) - this is the key!
      for (j in 1:nrow(profile)) {
        seg <- profile[j, ]

        # Find all grid positions that overlap this segment
        overlaps_idx <- which(
          grid$chr == seg$chromosome &
            grid$position >= seg$start &
            grid$position <= seg$end
        )

        if (length(overlaps_idx) > 0) {
          if ("total_copy_number" %in% colnames(profile)) {
            copy_number_matrix[i, overlaps_idx] <- seg$total_copy_number
          }
          if ("logr" %in% colnames(profile)) {
            logr_matrix[i, overlaps_idx] <- seg$logr
          }
        }
      }
    } else {
      # Old format: chr, startpos, endpos, nMajor, nMinor, LogR
      # ITERATE THROUGH SEGMENTS (not grid positions)
      for (j in 1:nrow(profile)) {
        seg <- profile[j, ]

        # Find all grid positions that overlap this segment
        overlaps_idx <- which(
          grid$chr == seg$chr &
            grid$position >= seg$startpos &
            grid$position <= seg$endpos
        )

        if (length(overlaps_idx) > 0) {
          copy_number_matrix[i, overlaps_idx] <- seg$nMajor + seg$nMinor
          if ("LogR" %in% colnames(profile)) {
            logr_matrix[i, overlaps_idx] <- seg$LogR
          }
        }
      }
    }
  }

  # Replace NA with missing value for logR
  logr_matrix[is.na(logr_matrix)] <- logr_missing_value

  message(sprintf("Final data dimensions: %d samples x %d genomic positions",
                  nrow(copy_number_matrix), ncol(copy_number_matrix)))
  message(sprintf("LogR data range: %.3f to %.3f",
                  min(logr_matrix, na.rm = TRUE), max(logr_matrix, na.rm = TRUE)))

  return(list(
    copy_number_matrix = copy_number_matrix,
    logr_matrix = logr_matrix,
    grid = grid,
    sample_names = sample_names
  ))
}

#' Align profiles to an existing genomic grid (using IRanges for speed)
#'
#' @param profile_dir Directory containing profile files
#' @param pattern File pattern to match
#' @param grid Existing genomic grid to align to
#' @param logr_missing_value Value for missing LogR data
#' @return List with copy_number_matrix and logr_matrix
#' @importFrom readr read_tsv
#' @importFrom dplyr %>% filter
#' @importFrom stringr str_replace
#' @importFrom IRanges IRanges findOverlaps
#' @export
#' @examples
#' \dontrun{
#' aligned <- align_profiles_to_grid("path/to/profiles", "*.txt", grid)
#' }
align_profiles_to_grid <- function(profile_dir, pattern, grid, logr_missing_value = 0) {

  message("Loading profiles and aligning to grid...")

  # Find profile files
  profile_files <- list.files(profile_dir, pattern = glob2rx(pattern), full.names = TRUE)

  if (length(profile_files) == 0) {
    pattern <- "*ASCAT.scprofile.txt"
    profile_files <- list.files(profile_dir, pattern = glob2rx(pattern), full.names = TRUE)
  }

  if (length(profile_files) == 0) {
    stop("No profile files found in ", profile_dir)
  }

  # Sort files numerically by cluster ID (cluster_1, cluster_2, ..., cluster_10, cluster_11)
  file_basenames <- basename(profile_files)
  cluster_ids <- as.numeric(gsub(".*cluster_(\\d+).*", "\\1", file_basenames))
  if (!all(is.na(cluster_ids))) {
    profile_files <- profile_files[order(cluster_ids)]
  }

  message(sprintf("Found %d profile files", length(profile_files)))

  # Initialize matrices
  n_samples <- length(profile_files)
  n_positions <- nrow(grid)

  copy_number_matrix <- matrix(NA, nrow = n_samples, ncol = n_positions)
  logr_matrix <- matrix(NA, nrow = n_samples, ncol = n_positions)

  sample_names <- basename(profile_files) %>%
    str_replace("_refitted.ASCAT.scprofile.txt", "") %>%
    str_replace(".ASCAT.scprofile.txt", "")

  rownames(copy_number_matrix) <- sample_names
  rownames(logr_matrix) <- sample_names

  # Load each profile and align to grid using vectorized operations
  for (i in seq_along(profile_files)) {
    message(sprintf("Processing file %d/%d", i, n_samples))

    profile <- read_tsv(profile_files[i], show_col_types = FALSE)

    if (nrow(profile) == 0) {
      warning(sprintf("Empty profile file: %s", profile_files[i]))
      next
    }

    # Standardize column names
    if ("chromosome" %in% colnames(profile)) {
      profile$chr_std <- profile$chromosome
      profile$start_std <- profile$start
      profile$end_std <- profile$end
      profile$cn_std <- if ("total_copy_number" %in% colnames(profile)) profile$total_copy_number else NA
      profile$logr_std <- if ("logr" %in% colnames(profile)) profile$logr else NA
    } else {
      profile$chr_std <- profile$chr
      profile$start_std <- profile$startpos
      profile$end_std <- profile$endpos
      profile$cn_std <- profile$nMajor + profile$nMinor
      profile$logr_std <- if ("LogR" %in% colnames(profile)) profile$LogR else NA
    }

    # Process each chromosome separately for efficiency
    for (chr in unique(grid$chr)) {
      grid_chr <- grid[grid$chr == chr, ]
      profile_chr <- profile[profile$chr_std == chr, ]

      if (nrow(profile_chr) == 0) next

      # Use IRanges for fast overlap detection
      grid_ranges <- IRanges::IRanges(start = grid_chr$position, width = 1)
      profile_ranges <- IRanges::IRanges(start = profile_chr$start_std, end = profile_chr$end_std)

      # Find overlaps: which profile bin does each grid position fall into?
      overlaps <- IRanges::findOverlaps(grid_ranges, profile_ranges, select = "first")

      # Assign values to grid positions
      valid_overlaps <- !is.na(overlaps)
      if (any(valid_overlaps)) {
        grid_indices <- which(grid$chr == chr)[valid_overlaps]
        profile_indices <- overlaps[valid_overlaps]

        copy_number_matrix[i, grid_indices] <- profile_chr$cn_std[profile_indices]
        logr_matrix[i, grid_indices] <- profile_chr$logr_std[profile_indices]
      }
    }
  }

  # Fill NAs with missing value
  logr_matrix[is.na(logr_matrix)] <- logr_missing_value

  message(sprintf("Aligned %d samples to grid with %d positions", n_samples, n_positions))

  return(list(
    copy_number_matrix = copy_number_matrix,
    logr_matrix = logr_matrix,
    sample_names = sample_names
  ))
}

#' Create genomic grid at specified resolution
#'
#' @param profile First profile for chromosome info
#' @param resolution Bin size in bp
#' @return Data frame with chr, position
#' @importFrom dplyr %>% filter bind_rows
#' @export
#' @examples
#' \dontrun{
#' profile <- read_tsv("sample.ASCAT.scprofile.txt")
#' grid <- create_genomic_grid(profile, resolution = 1e6)
#' }
create_genomic_grid <- function(profile, resolution) {
  chromosomes <- paste0("chr", c(1:22, "X"))

  grid_list <- lapply(chromosomes, function(chr) {
    # Handle both column naming conventions
    if ("chromosome" %in% colnames(profile)) {
      chr_data <- profile %>% filter(chromosome == !!chr)
      chr_start <- min(chr_data$start, na.rm = TRUE)
      chr_end <- max(chr_data$end, na.rm = TRUE)
    } else {
      chr_data <- profile %>% filter(chr == !!chr)
      chr_start <- min(chr_data$startpos, na.rm = TRUE)
      chr_end <- max(chr_data$endpos, na.rm = TRUE)
    }

    if (nrow(chr_data) == 0 || !is.finite(chr_start) || !is.finite(chr_end)) {
      return(NULL)
    }

    positions <- seq(chr_start, chr_end, by = resolution)

    data.frame(
      chr = chr,
      position = positions,
      stringsAsFactors = FALSE
    )
  })

  grid <- bind_rows(grid_list)
  return(grid)
}
