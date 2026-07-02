#' Load ASCAT.sc profile files and create matrices
#'
#' @param profile_dir Directory containing profile files
#' @param pattern File pattern (default: "*refitted.ASCAT.scprofile.txt")
#' @param resolution Genomic resolution in bp
#' @param logr_missing_value Value for missing logR data (used when
#'   \code{imputation_method = "scalar"}).
#' @param ascat_compat Logical. If TRUE, mimic ASCAT workflow behavior:
#'   recursive file discovery with refitted preference, hg38-length grid,
#'   sample filtering (<10% coverage removed), and scalar imputation
#'   (CN=2, LogR=logr_missing_value).
#' @param imputation_method How to fill NA values in the LogR matrix.
#'   \code{"row_mean"} (default) uses a two-pass approach: first replaces each
#'   NA with the row (sample) mean, then fills any remaining NAs with the
#'   overall matrix mean. \code{"scalar"} replaces all NAs with
#'   \code{logr_missing_value} (the previous default behaviour).
#'   See \code{\link{impute_missing}} for details.
#' @return List with copy_number_matrix, logr_matrix, grid, sample_names
#' @importFrom readr read_tsv
#' @importFrom dplyr %>%
#' @importFrom stringr str_replace
#' @export
#' @examples
#' \dontrun{
#' profiles <- load_ascat_profiles("path/to/profiles", resolution = 1e6)
#' # Use legacy scalar imputation:
#' profiles <- load_ascat_profiles("path/to/profiles", imputation_method = "scalar")
#' }
load_ascat_profiles <- function(profile_dir,
                                pattern = "*refitted.ASCAT.scprofile.txt",
                                resolution = 1000000,
                                logr_missing_value = 0,
                                ascat_compat = FALSE,
                                imputation_method = c("row_mean", "scalar")) {

  imputation_method <- match.arg(imputation_method)

  message("Loading ASCAT.sc profile files...")

  # Find profile files
  if (ascat_compat) {
    if (dir.exists(profile_dir)) {
      profile_files <- list.files(profile_dir, pattern = glob2rx(pattern),
                                  full.names = TRUE, recursive = TRUE)
    } else {
      profile_files <- list.files(dirname(profile_dir), pattern = glob2rx(pattern),
                                  full.names = TRUE, recursive = TRUE)
    }

    # Prefer refitted profiles if any exist
    refitted_files <- profile_files[grepl("refitted", profile_files)]
    if (length(refitted_files) > 0) {
      profile_files <- refitted_files
    }
  } else {
    profile_files <- list.files(profile_dir, pattern = glob2rx(pattern), full.names = TRUE)

    if (length(profile_files) == 0) {
      # Try non-refitted profiles
      pattern <- "*ASCAT.scprofile.txt"
      profile_files <- list.files(profile_dir, pattern = glob2rx(pattern), full.names = TRUE)
    }
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

  # Create genomic grid
  if (ascat_compat) {
    grid <- create_genomic_grid_hg38(resolution)
  } else {
    first_profile <- read_tsv(profile_files[1], show_col_types = FALSE)
    grid <- create_genomic_grid(first_profile, resolution)
  }
  message(sprintf("Created genomic grid with %d positions", nrow(grid)))

  # Initialize matrices
  n_samples <- length(profile_files)
  n_positions <- nrow(grid)

  copy_number_matrix <- matrix(NA, nrow = n_samples, ncol = n_positions)
  logr_matrix <- matrix(NA, nrow = n_samples, ncol = n_positions)

  if (ascat_compat) {
    base_filenames <- sub("\\.ASCAT\\.scprofile\\.txt$", "", basename(profile_files))
    sample_names <- simplify_plate_well(base_filenames)
  } else {
    sample_names <- basename(profile_files) %>%
      str_replace("_refitted.ASCAT.scprofile.txt", "") %>%
      str_replace(".ASCAT.scprofile.txt", "")
  }

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
        if (ascat_compat) {
          overlaps_idx <- grid$genomic_pos[
            grid$chromosome == seg$chromosome &
              grid$end >= seg$start &
              grid$start <= seg$end
          ]
        } else {
          overlaps_idx <- which(
            grid$chr == seg$chromosome &
              grid$position >= seg$start &
              grid$position <= seg$end
          )
        }

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

  if (ascat_compat) {
    # Remove samples with too much missing data (ASCAT workflow behavior)
    valid_samples <- rowSums(!is.na(copy_number_matrix)) > (ncol(copy_number_matrix) * 0.1)
    copy_number_matrix <- copy_number_matrix[valid_samples, , drop = FALSE]
    logr_matrix <- logr_matrix[valid_samples, , drop = FALSE]
    sample_names <- sample_names[valid_samples]

    rownames(copy_number_matrix) <- sample_names
    rownames(logr_matrix) <- sample_names

    if (nrow(copy_number_matrix) < 10) {
      stop("Not enough valid samples for clustering analysis")
    }

    # Scalar imputation (ASCAT workflow behavior)
    copy_number_matrix[is.na(copy_number_matrix)] <- 2
    logr_matrix[is.na(logr_matrix)] <- logr_missing_value
  } else {
    # Impute missing LogR values
    n_na <- sum(is.na(logr_matrix))
    if (n_na > 0) {
      message(sprintf("Imputing %d missing LogR values using method: %s", n_na, imputation_method))
      logr_matrix <- impute_missing(logr_matrix, method = imputation_method,
                                    scalar_value = logr_missing_value)
    }
  }

  message(sprintf("Final data dimensions: %d samples x %d genomic positions",
                  nrow(copy_number_matrix), ncol(copy_number_matrix)))
  message(sprintf("LogR data range: %.3f to %.3f",
                  min(logr_matrix, na.rm = TRUE), max(logr_matrix, na.rm = TRUE)))

  if (ascat_compat) {
    colnames(copy_number_matrix) <- paste0(grid$chromosome, ":", grid$start)
    colnames(logr_matrix) <- paste0(grid$chromosome, ":", grid$start)
  }

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
#' @param logr_missing_value Value for missing LogR data (used when
#'   \code{imputation_method = "scalar"}).
#' @param imputation_method How to fill NA values in the LogR matrix.
#'   \code{"row_mean"} (default) or \code{"scalar"}.
#'   See \code{\link{impute_missing}} for details.
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
align_profiles_to_grid <- function(profile_dir, pattern, grid,
                                   logr_missing_value = 0,
                                   imputation_method = c("row_mean", "scalar")) {

  imputation_method <- match.arg(imputation_method)

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

  # Impute missing LogR values
  n_na <- sum(is.na(logr_matrix))
  if (n_na > 0) {
    message(sprintf("Imputing %d missing LogR values using method: %s", n_na, imputation_method))
    logr_matrix <- impute_missing(logr_matrix, method = imputation_method,
                                  scalar_value = logr_missing_value)
  }

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
  # Auto-detect chromosome naming convention from profile data
  if ("chromosome" %in% colnames(profile)) {
    profile_chrs <- unique(profile$chromosome)
  } else if ("chr" %in% colnames(profile)) {
    profile_chrs <- unique(profile$chr)
  } else {
    stop("Profile must have a 'chromosome' or 'chr' column")
  }
  has_chr_prefix <- any(grepl("^chr", profile_chrs))
  if (has_chr_prefix) {
    chromosomes <- paste0("chr", c(1:22, "X"))
  } else {
    chromosomes <- c(as.character(1:22), "X")
  }

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

#' Create fixed hg38 genomic grid at specified resolution
#'
#' Internal helper for ASCAT workflow compatibility.
#'
#' @param resolution Bin size in bp
#' @return Data frame with chromosome, start, end, genomic_pos
#' @keywords internal
create_genomic_grid_hg38 <- function(resolution) {
  chr_lengths <- c(
    chr1 = 248956422, chr2 = 242193529, chr3 = 198295559, chr4 = 190214555,
    chr5 = 181538259, chr6 = 170805979, chr7 = 159345973, chr8 = 145138636,
    chr9 = 138394717, chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
    chr13 = 114364328, chr14 = 107043718, chr15 = 101991189, chr16 = 90338345,
    chr17 = 83257441, chr18 = 80373285, chr19 = 58617616, chr20 = 64444167,
    chr21 = 46709983, chr22 = 50818468, chrX = 156040895
  )

  grid <- data.frame()
  cumulative_pos <- 0

  for (chr in names(chr_lengths)) {
    chr_grid <- data.frame(
      chromosome = chr,
      start = seq(0, chr_lengths[chr], by = resolution),
      stringsAsFactors = FALSE
    )
    chr_grid$end <- pmin(chr_grid$start + resolution - 1, chr_lengths[chr])
    chr_grid$genomic_pos <- cumulative_pos + ceiling(chr_grid$start / resolution) + 1
    chr_grid <- chr_grid[chr_grid$start < chr_lengths[chr], , drop = FALSE]

    grid <- rbind(grid, chr_grid)
    cumulative_pos <- max(chr_grid$genomic_pos)
  }

  rownames(grid) <- NULL
  return(grid)
}

#' Create fixed hg19 genomic grid at specified resolution
#'
#' Internal helper that provides a fixed grid based on hg19 chromosome lengths.
#' Uses NCBI-style chromosome names (no "chr" prefix) to match ASCAT.sc hg19
#' convention.
#'
#' @param resolution Bin size in bp
#' @return Data frame with chromosome, start, end, genomic_pos
#' @export
create_genomic_grid_hg19 <- function(resolution) {
  chr_lengths <- c(
    "1" = 249250621, "2" = 243199373, "3" = 198022430, "4" = 191154276,
    "5" = 180915260, "6" = 171115067, "7" = 159138663, "8" = 146364022,
    "9" = 141213431, "10" = 135534747, "11" = 135006516, "12" = 133851895,
    "13" = 115169878, "14" = 107349540, "15" = 102531392, "16" = 90354753,
    "17" = 81195210, "18" = 78077248, "19" = 59128983, "20" = 63025520,
    "21" = 48129895, "22" = 51304566, "X" = 155270560
  )

  grid <- data.frame()
  cumulative_pos <- 0

  for (chr in names(chr_lengths)) {
    chr_grid <- data.frame(
      chromosome = chr,
      start = seq(0, chr_lengths[chr], by = resolution),
      stringsAsFactors = FALSE
    )
    chr_grid$end <- pmin(chr_grid$start + resolution - 1, chr_lengths[chr])
    chr_grid$genomic_pos <- cumulative_pos + ceiling(chr_grid$start / resolution) + 1
    chr_grid <- chr_grid[chr_grid$start < chr_lengths[chr], , drop = FALSE]

    grid <- rbind(grid, chr_grid)
    cumulative_pos <- max(chr_grid$genomic_pos)
  }

  rownames(grid) <- NULL
  return(grid)
}
