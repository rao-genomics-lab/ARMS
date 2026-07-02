

#' Parse plate and well identifiers from complex filename patterns
#'
#' Internal helper function to parse complex filename patterns.
#'
#' @param files Character vector of filenames
#' @return Character vector of simplified plate_well identifiers
#' @keywords internal
simplify_plate_well <- function(files) {
  row_letters <- LETTERS[1:8]  # A-H

  parse_one <- function(fn) {
    # Extract plate number
    plate <- str_match(fn, "plate\\s*([0-9]+)")[,2]
    if (is.na(plate)) {
      plate <- str_match(fn, "(?:^|[\\-_])P([0-9]+)(?=[\\-_])")[,2]
    }
    if (is.na(plate)) {
      plate_match <- str_match(fn, "_P([0-9]+)_")
      if (!is.na(plate_match[1])) {
        plate <- plate_match[,2]
      }
    }
    if (is.na(plate)) {
      plate_match <- str_match(fn, "P([0-9]+)")
      if (!is.na(plate_match[1])) {
        plate <- plate_match[,2]
      }
    }
    if (is.na(plate)) return(fn)

    # Detect barcode scheme
    m_rows <- str_match(fn, "([A-H])1-([A-H])12")
    m_cols <- str_match(fn, "(\\d{1,2})A-(\\d{1,2})H")

    if (!is.na(m_cols[1])) {
      # 16 barcodes, 0-based
      idx_m <- str_match(fn, "lcmad_0*([0-9]{1,3})")
      if (is.na(idx_m[1])) return(fn)
      idx0 <- as.integer(idx_m[,2]) %% 16
      col_start <- as.integer(m_cols[,2])
      col_end   <- as.integer(m_cols[,3])
      col <- if (idx0 < 8) col_start else col_end
      row <- row_letters[(idx0 %% 8) + 1]
    } else if (!is.na(m_rows[1])) {
      # 24 barcodes, 0-based
      idx_m <- str_match(fn, "ILMN_0*([0-9]{1,3})(?!\\d)")
      if (is.na(idx_m[1])) {
        idx_m <- str_match(fn, "ILMN_[A-Za-z0-9]+_0*([0-9]{1,3})")
      }
      if (is.na(idx_m[1])) return(fn)
      idx0 <- as.integer(idx_m[,2]) %% 24
      row1 <- m_rows[,2]
      row2 <- m_rows[,3]
      row <- if (idx0 < 12) row1 else row2
      col <- (idx0 %% 12) + 1
    } else {
      # 96 barcodes, 1-based
      idx_m <- str_match(fn, "ILMN_[A-Za-z0-9]+_0*([0-9]{1,3})")
      if (is.na(idx_m[1])) {
        idx_m <- str_match(fn, ".*?([0-9]{1,3})(?!.*[0-9])")
      }
      if (is.na(idx_m[1])) return(fn)
      idx1 <- ((as.integer(idx_m[,2]) - 1) %% 96) + 1
      row <- row_letters[((idx1 - 1) %/% 12) + 1]
      col <- ((idx1 - 1) %% 12) + 1
    }

    sprintf("plate%s_%d%s", plate, col, row)
  }

  sapply(files, function(f) {
    tryCatch(parse_one(f), error = function(e) f)
  }, USE.NAMES = FALSE)
}

#' Normalize identifiers to standard plate_column-row format
#'
#' Internal helper function to normalize identifiers.
#'
#' @param identifiers Character vector of identifiers to normalize
#' @return Character vector of normalized identifiers in format plate{N}_{column}{row}
#' @keywords internal
normalize_identifier <- function(identifiers) {
  sapply(identifiers, function(x) {
    # Handle simple GeoJSON formats
    if (grepl("^\\d+_\\d*[A-H]+\\d*$", x)) {
      parts <- str_split(x, "_")[[1]]
      plate_num <- parts[1]
      well <- parts[2]
      well_match <- str_match(well, "^(\\d*)([A-H]+)(\\d*)$")
      if (!is.na(well_match[1])) {
        col_prefix <- well_match[,2]
        row_letters <- well_match[,3]
        col_suffix <- well_match[,4]
        if (col_prefix == "" && col_suffix != "") {
          standardized_well <- paste0(col_suffix, row_letters)
        } else if (col_prefix != "" && col_suffix == "") {
          standardized_well <- well
        } else {
          standardized_well <- well
        }
        return(sprintf("plate%s_%s", plate_num, standardized_well))
      }
    }

    # Check if already in correct format
    if (grepl("^plate\\d+_\\d+[A-H]$", x)) {
      return(x)
    }

    # Use simplify_plate_well for complex filenames
    simplified <- simplify_plate_well(x)
    if (simplified != x && grepl("^plate\\d+_\\d+[A-H]$", simplified)) {
      return(simplified)
    }

    # Fallback: return original
    return(x)
  }, USE.NAMES = FALSE)
}

#' Match spatial tiles to samples using robust matching strategies
#'
#' @param sample_ids Character vector of sample identifiers
#' @param tile_names Character vector of tile names from GeoJSON
#' @param verbose Logical, whether to print matching information
#' @return List with matches (data.frame), unmatched_samples, unmatched_tiles, stats
#' @importFrom stringr str_match str_split
#' @export
#' @examples
#' \dontrun{
#' sample_ids <- c("plate1_1A", "plate1_2B")
#' tile_names <- c("1_1A", "1_2B")
#' matches <- match_spatial_tiles_to_samples(sample_ids, tile_names)
#' }
match_spatial_tiles_to_samples <- function(sample_ids, tile_names, verbose = TRUE) {
  if (verbose) {
    message(sprintf("Matching %d sample IDs to %d tile names...",
                    length(sample_ids), length(tile_names)))
  }

  # Normalize both sets
  normalized_samples <- normalize_identifier(sample_ids)
  normalized_tiles <- normalize_identifier(tile_names)

  if (verbose) {
    message("Normalized sample identifiers (first 5):")
    message(paste(head(paste(sample_ids, "->", normalized_samples), 5), collapse = "\n"))
    message("Normalized tile identifiers (first 5):")
    message(paste(head(paste(tile_names, "->", normalized_tiles), 5), collapse = "\n"))
  }

  # Create lookup tables
  sample_lookup <- data.frame(
    original = sample_ids,
    normalized = normalized_samples,
    stringsAsFactors = FALSE
  )

  tile_lookup <- data.frame(
    original = tile_names,
    normalized = normalized_tiles,
    stringsAsFactors = FALSE
  )

  # Direct exact match
  direct_matches <- merge(sample_lookup, tile_lookup,
                          by = "normalized", suffixes = c("_sample", "_tile"))

  all_matches <- data.frame()

  if (nrow(direct_matches) > 0) {
    if (verbose) message(sprintf("Direct matches found: %d", nrow(direct_matches)))

    # Rename columns
    if ("original.x" %in% names(direct_matches)) names(direct_matches)[names(direct_matches) == "original.x"] <- "original_sample"
    if ("original.y" %in% names(direct_matches)) names(direct_matches)[names(direct_matches) == "original.y"] <- "original_tile"

    if ("original_sample" %in% names(direct_matches) && "original_tile" %in% names(direct_matches)) {
      direct_matches_formatted <- data.frame(
        normalized = direct_matches$normalized,
        original_sample = direct_matches$original_sample,
        original_tile = direct_matches$original_tile,
        match_type = "direct",
        stringsAsFactors = FALSE
      )
      all_matches <- rbind(all_matches, direct_matches_formatted)
    }
  }

  # Calculate final stats
  matched_samples <- if (nrow(all_matches) > 0) all_matches$original_sample else character(0)
  matched_tiles <- if (nrow(all_matches) > 0) all_matches$original_tile else character(0)
  final_unmatched_samples <- sample_ids[!sample_ids %in% matched_samples]
  final_unmatched_tiles <- tile_names[!tile_names %in% matched_tiles]

  stats <- list(
    total_samples = length(sample_ids),
    total_tiles = length(tile_names),
    matches_found = nrow(all_matches),
    direct_matches = sum(all_matches$match_type == "direct", na.rm = TRUE),
    unmatched_samples = length(final_unmatched_samples),
    unmatched_tiles = length(final_unmatched_tiles),
    match_rate = if (min(length(sample_ids), length(tile_names)) > 0) round(nrow(all_matches) / min(length(sample_ids), length(tile_names)) * 100, 1) else 0
  )

  if (verbose) {
    message("=== Matching Statistics ===")
    message(sprintf("Total matches: %d/%d (%.1f%%)", stats$matches_found, min(stats$total_samples, stats$total_tiles), stats$match_rate))
    message(sprintf("Direct matches: %d", stats$direct_matches))
    message(sprintf("Unmatched samples: %d", stats$unmatched_samples))
    message(sprintf("Unmatched tiles: %d", stats$unmatched_tiles))
  }

  return(list(
    matches = all_matches,
    unmatched_samples = final_unmatched_samples,
    unmatched_tiles = final_unmatched_tiles,
    stats = stats
  ))
}


#' Load and process GeoJSON tile data with spatial coordinates
#'
#' @param geojson_file Path to GeoJSON file
#' @return sf object with tile geometries and centroids
#' @importFrom sf st_read st_set_crs st_centroid st_coordinates
#' @export
#' @examples
#' \dontrun{
#' tiles <- load_geojson_tiles("path/to/tiles.geojson")
#' }
load_geojson_tiles <- function(geojson_file) {

  message("Loading GeoJSON tile data...")

  tiles <- st_read(geojson_file, quiet = TRUE)

  # Set CRS to NA to avoid coordinate issues
  tiles <- st_set_crs(tiles, NA)

  # Extract centroids if meanx/meany don't exist
  if (!"meanx" %in% colnames(tiles) || !"meany" %in% colnames(tiles)) {
    centroids <- st_centroid(tiles$geometry)
    coords <- st_coordinates(centroids)
    tiles$meanx <- coords[, "X"]
    tiles$meany <- coords[, "Y"]
  }

  message(sprintf("Loaded %d tiles", nrow(tiles)))

  return(tiles)
}

#' Match sample names to spatial tiles using robust identifier matching
#'
#' @param sample_names Sample names from profiles
#' @param tiles sf object with tile data
#' @return Data frame with sample, tile, x, y coordinates
#' @importFrom stringr str_match str_split
#' @export
#' @examples
#' \dontrun{
#' tiles <- load_geojson_tiles("tiles.geojson")
#' matches <- match_samples_to_tiles(sample_names, tiles)
#' }
match_samples_to_tiles <- function(sample_names, tiles) {

  # Extract tile names from GeoJSON
  if ("name" %in% colnames(tiles)) {
    tile_names <- tiles$name
  } else if ("id" %in% colnames(tiles)) {
    tile_names <- tiles$id
  } else {
    stop("GeoJSON must have 'name' or 'id' column")
  }

  # Use legacy matching function from utils.R
  # Function returns a list with: matches, unmatched_samples, unmatched_tiles, stats
  matching_result_list <- match_spatial_tiles_to_samples(sample_names, tile_names, verbose = FALSE)
  matching_results <- matching_result_list$matches

  # Create output dataframe in the same order as input sample_names
  matches <- data.frame(
    sample = sample_names,
    tile = NA_character_,
    x = NA_real_,
    y = NA_real_,
    stringsAsFactors = FALSE
  )

  # Map coordinates using the matching results
  # matching_results has columns: normalized, original_sample, original_tile, match_type
  for (i in seq_along(sample_names)) {
    # Find this sample in matching_results using which()
    match_idx <- which(matching_results$original_sample == sample_names[i])

    if (length(match_idx) > 0) {
      tile_name <- matching_results$original_tile[match_idx[1]]

      # Find this tile in the tiles dataframe using which()
      tile_idx <- which(tiles$name == tile_name)

      if (length(tile_idx) > 0) {
        matches$tile[i] <- tile_name
        matches$x[i] <- tiles$meanx[tile_idx[1]]
        matches$y[i] <- tiles$meany[tile_idx[1]]
      }
    }
  }

  message(sprintf("Matched %d/%d samples to spatial coordinates",
                  sum(!is.na(matches$x)), length(sample_names)))

  return(matches)
}
