

#' Calculate silhouette score for clustering
#'
#' @param data Data matrix
#' @param clusters Cluster assignments
#' @return Mean silhouette score
#' @importFrom cluster silhouette
#' @keywords internal
calculate_silhouette_score <- function(data, clusters) {
  if (length(unique(clusters)) < 2) return(0)

  tryCatch({
    dist_matrix <- dist(data)
    silhouette_result <- silhouette(clusters, dist_matrix)
    return(mean(silhouette_result[, 3]))
  }, error = function(e) {
    return(0)
  })
}

#' Calculate gap statistic for clustering (optimized)
#'
#' @param data Data matrix
#' @param clusters Cluster assignments
#' @param B Number of bootstrap samples (reduced for speed)
#' @return Gap statistic value
#' @importFrom cluster clusGap
#' @keywords internal
calculate_gap_statistic <- function(data, clusters, B = 10) {
  tryCatch({
    gap_result <- clusGap(data, FUN = function(x, k) list(cluster = clusters),
                          K.max = length(unique(clusters)), B = B)
    return(max(gap_result$Tab[, "gap"]))
  }, error = function(e) {
    return(0)
  })
}

#' Calculate pvclust AU p-value for a split (optimized)
#'
#' @param data Data matrix
#' @param split_indices Indices for the split
#' @param nboot Number of bootstrap replicates (reduced for speed)
#' @param fast_mode If TRUE, use simplified approximation
#' @return AU p-value
#' @importFrom pvclust pvclust
#' @keywords internal
calculate_pvclust_pvalue <- function(data, split_indices, nboot = 10, fast_mode = FALSE) {
  if (length(split_indices) < 3 || (nrow(data) - length(split_indices)) < 3) {
    return(0)  # Not enough samples for meaningful test
  }

  if (fast_mode) {
    # Use simple approximation instead of full pvclust
    dist_matrix <- as.matrix(dist(data))
    within_cluster_dist <- mean(dist_matrix[split_indices, split_indices], na.rm = TRUE)
    between_cluster_dist <- mean(dist_matrix[split_indices, -split_indices], na.rm = TRUE)

    # Simple separation ratio as proxy for significance
    separation_ratio <- between_cluster_dist / (within_cluster_dist + 1e-10)
    return(min(separation_ratio / 2, 1.0))  # Normalize to [0,1]
  }

  tryCatch({
    # Highly optimized pvclust with minimal bootstrap iterations
    pv_result <- pvclust(t(data), method.hclust = "ward.D2",
                         method.dist = "euclidean", nboot = nboot, quiet = TRUE)

    au_values <- pv_result$edges[, "au"]
    return(max(au_values, na.rm = TRUE))

  }, error = function(e) {
    return(0)
  })
}

#' Merge small clusters into nearest large clusters
#'
#' @param clusters Vector of cluster assignments
#' @param data Data matrix
#' @param min_size Minimum cluster size threshold
#' @return Vector of merged cluster assignments
#' @keywords internal
merge_small_clusters <- function(clusters, data, min_size) {
  cluster_sizes <- table(clusters)
  small_clusters <- which(cluster_sizes < min_size)

  if (length(small_clusters) == 0) {
    return(clusters)
  }

  # Calculate centroids
  unique_clusters <- sort(unique(clusters))
  centroids <- matrix(0, nrow = length(unique_clusters), ncol = ncol(data))

  for (i in seq_along(unique_clusters)) {
    cluster_id <- unique_clusters[i]
    cluster_indices <- which(clusters == cluster_id)
    if (length(cluster_indices) > 1) {
      centroids[i, ] <- colMeans(data[cluster_indices, , drop = FALSE])
    } else {
      centroids[i, ] <- data[cluster_indices, ]
    }
  }

  # Merge small clusters with nearest large clusters
  merged_clusters <- clusters

  for (small_cluster_id in names(small_clusters)) {
    small_idx <- which(unique_clusters == as.numeric(small_cluster_id))

    # Find distances to all other clusters
    distances <- numeric(length(unique_clusters))
    for (j in seq_along(unique_clusters)) {
      if (j != small_idx) {
        distances[j] <- dist(rbind(centroids[small_idx, ], centroids[j, ]))
      } else {
        distances[j] <- Inf
      }
    }

    # Find nearest cluster that is large enough
    large_clusters <- which(cluster_sizes >= min_size)
    valid_targets <- intersect(seq_along(unique_clusters),
                               which(unique_clusters %in% names(large_clusters)))

    if (length(valid_targets) > 0) {
      nearest_large <- valid_targets[which.min(distances[valid_targets])]
      target_cluster <- unique_clusters[nearest_large]

      # Merge small cluster into target cluster
      merged_clusters[merged_clusters == as.numeric(small_cluster_id)] <- target_cluster

      message(sprintf("Merged small cluster %s (%d samples) into cluster %d",
                      small_cluster_id, cluster_sizes[small_cluster_id], target_cluster))
    }
  }

  return(merged_clusters)
}

#' Recursive dynamic hierarchical clustering with quality-based splitting
#'
#' Internal recursive function used by hierarchical_clustering.
#'
#' @param data Data matrix (samples x features)
#' @param min_cluster_size Minimum cluster size threshold
#' @param max_module_size Maximum module size for recursive splitting
#' @param silhouette_threshold Minimum silhouette improvement threshold
#' @param gap_threshold Minimum gap statistic threshold
#' @param pvclust_threshold Minimum AU p-value threshold
#' @param max_recursion_depth Maximum recursion depth
#' @param hierarchical_k Fallback k for fixed cut
#' @param depth Current recursion depth (internal)
#' @param dist_matrix Pre-computed distance matrix for performance
#' @return Vector of cluster assignments
#' @importFrom dynamicTreeCut cutreeDynamic
#' @keywords internal
recursive_dynamic_hac <- function(data,
                                  min_cluster_size = 10,
                                  max_module_size = 100,
                                  silhouette_threshold = 0.05,
                                  gap_threshold = 0.01,
                                  pvclust_threshold = 0.95,
                                  max_recursion_depth = 4,
                                  hierarchical_k = 20,
                                  depth = 0,
                                  dist_matrix = NULL) {

  n_samples <- nrow(data)

  # Early termination
  if (n_samples < min_cluster_size * 2 || depth >= max_recursion_depth) {
    return(rep(1, n_samples))
  }

  message(sprintf("Recursive HAC depth %d: Processing %d samples", depth, n_samples))

  # Compute or reuse distance matrix
  if (is.null(dist_matrix)) {
    dist_matrix <- dist(data, method = "euclidean")
  }

  # Perform hierarchical clustering
  hclust_result <- hclust(dist_matrix, method = "ward.D2")

  # Apply dynamic tree cut
  tryCatch({
    dynamic_clusters <- cutreeDynamic(
      dendro = hclust_result,
      distM = as.matrix(dist_matrix),
      minClusterSize = min_cluster_size,
      method = "hybrid",
      deepSplit = 2,
      verbose = 0
    )

    # If dynamic tree cut fails or produces only one cluster, fall back to fixed cut
    if (length(unique(dynamic_clusters)) <= 1) {
      k_clusters <- min(hierarchical_k, floor(n_samples / min_cluster_size))
      k_clusters <- max(k_clusters, 2)
      dynamic_clusters <- cutree(hclust_result, k = k_clusters)
    }

  }, error = function(e) {
    message(sprintf("Dynamic tree cut failed at depth %d, using fixed cut", depth))
    k_clusters <- min(hierarchical_k, floor(n_samples / min_cluster_size))
    k_clusters <- max(k_clusters, 2)
    dynamic_clusters <- cutree(hclust_result, k = k_clusters)
  })

  # Merge small clusters
  dynamic_clusters <- merge_small_clusters(dynamic_clusters, data, min_cluster_size)

  # Get cluster sizes
  unique_clusters <- unique(dynamic_clusters)
  cluster_sizes <- table(dynamic_clusters)

  message(sprintf("Dynamic cut produced %d clusters with sizes: %s",
                  length(unique_clusters), paste(cluster_sizes, collapse = ", ")))

  # Check for recursive splitting
  large_clusters <- which(cluster_sizes > max_module_size)

  if (length(large_clusters) == 0 || depth >= max_recursion_depth) {
    return(dynamic_clusters)
  }

  # Process each large cluster recursively
  final_clusters <- dynamic_clusters
  next_cluster_id <- max(dynamic_clusters) + 1

  for (cluster_id in names(large_clusters)) {
    cluster_indices <- which(dynamic_clusters == as.numeric(cluster_id))
    cluster_data <- data[cluster_indices, , drop = FALSE]

    # Skip if too small to split
    if (nrow(cluster_data) < min_cluster_size * 3) {
      next
    }

    message(sprintf("Evaluating cluster %s (%d samples) for recursive splitting",
                    cluster_id, length(cluster_indices)))

    # Calculate current quality
    current_silhouette <- calculate_silhouette_score(cluster_data, rep(1, nrow(cluster_data)))

    # Subset distance matrix
    sub_dist_matrix <- as.dist(as.matrix(dist_matrix)[cluster_indices, cluster_indices])

    # Get subclusters
    sub_clusters <- recursive_dynamic_hac(
      data = cluster_data,
      min_cluster_size = min_cluster_size,
      max_module_size = max_module_size,
      silhouette_threshold = silhouette_threshold,
      gap_threshold = gap_threshold,
      pvclust_threshold = pvclust_threshold,
      max_recursion_depth = max_recursion_depth,
      hierarchical_k = hierarchical_k,
      depth = depth + 1,
      dist_matrix = sub_dist_matrix
    )

    # Only proceed if we get more than one subcluster
    if (length(unique(sub_clusters)) > 1) {

      # Calculate quality metrics for the split
      new_silhouette <- calculate_silhouette_score(cluster_data, sub_clusters)

      # Use fast mode for deeper levels
      use_fast_mode <- depth > 1

      gap_statistic <- if (use_fast_mode) {
        # Simple approximation
        old_wcss <- sum(apply(cluster_data, 2, var))
        new_wcss <- sum(by(cluster_data, sub_clusters, function(x) sum(apply(x, 2, var))))
        old_wcss - new_wcss
      } else {
        calculate_gap_statistic(cluster_data, sub_clusters, B = 10)
      }

      # Calculate pvclust p-value
      au_pvalue <- calculate_pvclust_pvalue(cluster_data,
                                            which(sub_clusters == 1),
                                            nboot = ifelse(use_fast_mode, 5, 10),
                                            fast_mode = use_fast_mode)

      # Apply stopping criteria
      silhouette_improved <- (new_silhouette - current_silhouette) >= silhouette_threshold
      gap_significant <- gap_statistic >= gap_threshold
      au_significant <- au_pvalue >= pvclust_threshold

      message(sprintf("Split evaluation - Silhouette: %.3f -> %.3f (improved: %s), Gap: %.3f (significant: %s), AU p-value: %.3f (significant: %s)",
                      current_silhouette, new_silhouette, silhouette_improved,
                      gap_statistic, gap_significant, au_pvalue, au_significant))

      # Accept split if criteria are met
      accept_split <- silhouette_improved || (gap_significant && au_significant)

      if (accept_split) {
        message(sprintf("Accepting split for cluster %s", cluster_id))

        # Update cluster assignments
        for (i in seq_along(cluster_indices)) {
          original_idx <- cluster_indices[i]
          sub_cluster <- sub_clusters[i]

          if (sub_cluster == 1) {
            final_clusters[original_idx] <- as.numeric(cluster_id)
          } else {
            final_clusters[original_idx] <- next_cluster_id + sub_cluster - 2
          }
        }

        next_cluster_id <- next_cluster_id + length(unique(sub_clusters)) - 1

      } else {
        message(sprintf("No meaningful subclusters found for cluster %s", cluster_id))
      }
    }
  }

  return(final_clusters)
}

#' Create cluster assignments from manually defined GeoJSON groups
#'
#' This function allows you to manually define spatial groups using multiple GeoJSON files,
#' where each file represents a distinct cluster/group. The function assigns cluster IDs
#' based on which GeoJSON file contains each sample, creating cluster assignments compatible
#' with the pseudobulk workflow.
#'
#' @param geojson_files Character vector of paths to GeoJSON files (each file = one cluster)
#' @param sample_names Character vector of full sample names from your dataset
#' @param unassigned_value Cluster ID for samples not in any GeoJSON (default: 0 = excluded from pseudobulk)
#' @param verbose Logical, whether to print matching statistics (default: TRUE)
#' @return Named integer vector of cluster assignments (names = sample names, values = cluster IDs)
#' @importFrom sf st_read
#' @export
#' @examples
#' \dontrun{
#' # Define manual groups via GeoJSON files
#' geojson_groups <- c(
#'   "tumor_region1.geojson",
#'   "tumor_region2.geojson",
#'   "normal_region.geojson"
#' )
#'
#' # Create cluster assignments
#' manual_clusters <- create_clusters_from_geojson_groups(
#'   geojson_files = geojson_groups,
#'   sample_names = profiles$sample_names
#' )
#'
#' # Use with pseudobulk workflow
#' pseudobulk <- create_pseudobulk_from_bins(
#'   "result.Rda",
#'   clusters = manual_clusters,
#'   sample_names = profiles$sample_names,
#'   output_dir = "pseudobulk_output"
#' )
#' }
create_clusters_from_geojson_groups <- function(geojson_files,
                                                 sample_names,
                                                 unassigned_value = 0,
                                                 verbose = TRUE) {

  if (length(geojson_files) == 0) {
    stop("geojson_files must contain at least one GeoJSON file path")
  }

  if (length(sample_names) == 0) {
    stop("sample_names must contain at least one sample name")
  }

  # Initialize cluster vector with unassigned value
  clusters <- rep(unassigned_value, length(sample_names))
  names(clusters) <- sample_names

  # Normalize sample names once for efficiency
  normalized_samples <- normalize_identifier(sample_names)

  if (verbose) {
    message(sprintf("Processing %d GeoJSON files for %d samples...",
                    length(geojson_files), length(sample_names)))
  }

  # Process each GeoJSON file
  total_assigned <- 0

  for (i in seq_along(geojson_files)) {
    geojson_file <- geojson_files[i]
    cluster_id <- i

    if (!file.exists(geojson_file)) {
      warning(sprintf("GeoJSON file not found: %s (skipping)", geojson_file))
      next
    }

    if (verbose) {
      message(sprintf("\nProcessing file %d/%d: %s", i, length(geojson_files),
                      basename(geojson_file)))
    }

    # Load GeoJSON
    tryCatch({
      tiles <- st_read(geojson_file, quiet = !verbose)

      # Extract tile names
      if ("name" %in% colnames(tiles)) {
        tile_names <- tiles$name
      } else if ("id" %in% colnames(tiles)) {
        tile_names <- tiles$id
      } else {
        warning(sprintf("GeoJSON file %s has no 'name' or 'id' column (skipping)",
                        basename(geojson_file)))
        next
      }

      # Normalize tile names
      normalized_tiles <- normalize_identifier(tile_names)

      # Match normalized tiles to normalized samples
      matches <- match(normalized_tiles, normalized_samples)
      valid_matches <- !is.na(matches)

      if (sum(valid_matches) == 0) {
        warning(sprintf("No matches found for GeoJSON file %s", basename(geojson_file)))
        next
      }

      # Assign cluster ID to matched samples
      matched_sample_indices <- matches[valid_matches]
      clusters[matched_sample_indices] <- cluster_id

      if (verbose) {
        message(sprintf("  Cluster %d: Matched %d/%d tiles to samples",
                        cluster_id, sum(valid_matches), length(tile_names)))
      }

      total_assigned <- total_assigned + sum(valid_matches)

    }, error = function(e) {
      warning(sprintf("Error processing GeoJSON file %s: %s",
                      basename(geojson_file), e$message))
    })
  }

  # Summary statistics
  if (verbose) {
    message("\n=== Cluster Assignment Summary ===")
    cluster_counts <- table(clusters)

    if (unassigned_value %in% names(cluster_counts)) {
      n_unassigned <- cluster_counts[as.character(unassigned_value)]
      message(sprintf("Unassigned samples: %d (%.1f%%)",
                      n_unassigned, 100 * n_unassigned / length(sample_names)))
    }

    assigned_clusters <- setdiff(unique(clusters), unassigned_value)
    if (length(assigned_clusters) > 0) {
      message(sprintf("Assigned samples: %d (%.1f%%)",
                      total_assigned, 100 * total_assigned / length(sample_names)))
      message(sprintf("Number of clusters: %d", length(assigned_clusters)))
      message("Cluster sizes:")
      for (cid in sort(assigned_clusters)) {
        message(sprintf("  Cluster %d: %d samples", cid, sum(clusters == cid)))
      }
    }
  }

  return(clusters)
}

#' Manually merge specified clusters
#'
#' Merge two or more clusters into one by reassigning cluster IDs. The first
#' element of each merge group becomes the target ID; all other IDs in the group
#' are reassigned to that target.
#'
#' @param clustering_obj List as returned by \code{hierarchical_clustering()},
#'   containing at minimum a \code{clusters} element (named integer vector).
#' @param clusters_to_merge Either a single integer vector (e.g. \code{c(1, 3, 5)})
#'   where the first element is the target, or a list of such vectors for
#'   multiple independent merges.
#' @param renumber Logical. If \code{TRUE} (default), renumber clusters
#'   sequentially (1, 2, 3, ...) after merging. If \code{FALSE}, preserve the
#'   original target IDs (may leave gaps).
#' @return The same clustering object with \code{$clusters} updated.
#' @export
#' @examples
#' \dontrun{
#' # Merge clusters 3 and 5 into cluster 1
#' result <- merge_clusters_manually(clustering_obj, c(1, 3, 5))
#'
#' # Multiple independent merges
#' result <- merge_clusters_manually(clustering_obj,
#'   list(c(1, 3, 5), c(2, 4)))
#'
#' # Keep original IDs (no renumbering)
#' result <- merge_clusters_manually(clustering_obj, c(1, 3), renumber = FALSE)
#' }
merge_clusters_manually <- function(clustering_obj, clusters_to_merge, renumber = TRUE) {

  # Validate clustering_obj

  if (!is.list(clustering_obj) || is.null(clustering_obj$clusters)) {
    stop("clustering_obj must be a list with a 'clusters' element")
  }

  clusters <- clustering_obj$clusters

  # Normalise input: wrap a plain vector in a list

  if (!is.list(clusters_to_merge)) {
    clusters_to_merge <- list(clusters_to_merge)
  }

  existing_ids <- unique(clusters)

  # --- Validation ---
  # Check all IDs exist
  all_ids <- unlist(clusters_to_merge)
  missing <- setdiff(all_ids, existing_ids)
  if (length(missing) > 0) {
    stop(sprintf("Cluster IDs not found in data: %s",
                 paste(missing, collapse = ", ")))
  }

  # Check no ID appears in multiple merge groups
  if (length(all_ids) != length(unique(all_ids))) {
    stop("A cluster ID must not appear in more than one merge group")
  }

  # Warn about single-element groups
  for (i in seq_along(clusters_to_merge)) {
    if (length(clusters_to_merge[[i]]) < 2) {
      warning(sprintf("Merge group %d has only one element (%s) — nothing to merge",
                      i, clusters_to_merge[[i]]))
    }
  }

  # --- Perform merges ---
  for (group in clusters_to_merge) {
    if (length(group) < 2) next
    target <- group[1]
    sources <- group[-1]
    clusters[clusters %in% sources] <- target
    message(sprintf("Merged cluster(s) %s into cluster %d",
                    paste(sources, collapse = ", "), target))
  }

  # --- Renumber ---
  if (renumber) {
    old_ids <- sort(unique(clusters))
    new_ids <- seq_along(old_ids)
    mapping <- setNames(new_ids, old_ids)
    clusters <- unname(mapping[as.character(clusters)])
    names(clusters) <- names(clustering_obj$clusters)
    if (!identical(as.integer(old_ids), new_ids)) {
      message(sprintf("Renumbered clusters: %s",
                      paste(sprintf("%d->%d", old_ids, new_ids), collapse = ", ")))
    }
  }

  message(sprintf("Result: %d clusters with sizes: %s",
                  length(unique(clusters)),
                  paste(table(clusters), collapse = ", ")))

  clustering_obj$clusters <- clusters
  return(clustering_obj)
}

#' Manually split specified clusters into subclusters
#'
#' The inverse of \code{merge_clusters_manually()}. Takes specified cluster(s),
#' applies hierarchical clustering to subdivide them, and integrates the new
#' subclusters back into the original clustering object.
#'
#' @param clustering_obj List as returned by \code{hierarchical_clustering()} or
#'   \code{weighted_hierarchical_clustering()}, containing at minimum a
#'   \code{clusters} element (named integer vector).
#' @param clusters_to_split Integer vector of cluster IDs to split (e.g., \code{c(2, 5)}).
#' @param logr_matrix The original logR matrix used for clustering (samples x genomic bins).
#' @param method Clustering method: \code{"hierarchical"} (standard) or \code{"weighted"}
#'   (uses size weighting and optional discretization). Default: \code{"hierarchical"}.
#' @param grid Genomic grid data frame. Required if \code{method = "weighted"} with
#'   \code{use_size_weighting = TRUE}.
#' @param k Number of subclusters per split cluster. If \code{NULL} (default), uses
#'   dynamic tree cut to determine automatically.
#' @param use_dynamic Logical. Use dynamic tree cut for determining subclusters.
#'   Default: \code{TRUE}.
#' @param min_cluster_size Minimum subcluster size. Default: 3 (smaller than usual
#'   since splits often involve smaller groups).
#' @param discretize Logical. If \code{TRUE} and \code{method = "weighted"}, discretize
#'   logR values before clustering. Default: \code{FALSE}.
#' @param gain_threshold LogR value above which is classified as gain. Default: 0.2.
#' @param loss_threshold LogR value below which is classified as loss. Default: -0.2.
#' @param use_size_weighting Logical. If \code{TRUE} and \code{method = "weighted"},
#'   weight bins by contiguous region size. Default: \code{TRUE}.
#' @param min_bins_full_weight Minimum contiguous bins for full weight. Default: 5.
#' @param weight_function Weight function: \code{"linear"}, \code{"sigmoid"}, or
#'   \code{"step"}. Default: \code{"linear"}.
#' @param renumber Logical. If \code{TRUE} (default), renumber all clusters
#'   sequentially (1, 2, 3, ...) after splitting. If \code{FALSE}, original cluster
#'   IDs are preserved where possible and new subclusters get IDs starting from
#'   max(existing) + 1.
#'
#' @return The same clustering object with:
#' \describe{
#'   \item{clusters}{Updated cluster assignments}
#'   \item{split_info}{List documenting which clusters were split and resulting subclusters}
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' # Initial clustering
#' result <- hierarchical_clustering(logr_matrix, k = 5)
#'
#' # Split cluster 2 into 3 subclusters using hierarchical clustering
#' result_split <- split_clusters_manually(
#'   clustering_obj = result,
#'   clusters_to_split = 2,
#'   logr_matrix = logr_matrix,
#'   method = "hierarchical",
#'   k = 3
#' )
#'
#' # Split multiple clusters using weighted clustering
#' result_split <- split_clusters_manually(
#'   clustering_obj = result,
#'   clusters_to_split = c(2, 4),
#'   logr_matrix = logr_matrix,
#'   method = "weighted",
#'   grid = profiles$grid,
#'   discretize = TRUE,
#'   k = NULL  # auto-detect via dynamic tree cut
#' )
#'
#' # Split without renumbering (keeps original IDs where possible)
#' result_split <- split_clusters_manually(
#'   clustering_obj = result,
#'   clusters_to_split = 2,
#'   logr_matrix = logr_matrix,
#'   method = "hierarchical",
#'   k = 2,
#'   renumber = FALSE
#' )
#' }
split_clusters_manually <- function(clustering_obj,
                                     clusters_to_split,
                                     logr_matrix,
                                     method = "hierarchical",
                                     grid = NULL,
                                     k = NULL,
                                     use_dynamic = TRUE,
                                     min_cluster_size = 3,
                                     discretize = FALSE,
                                     gain_threshold = 0.2,
                                     loss_threshold = -0.2,
                                     use_size_weighting = TRUE,
                                     min_bins_full_weight = 5,
                                     weight_function = "linear",
                                     renumber = TRUE) {

  # --- Input validation ---

  # Validate clustering_obj

  if (!is.list(clustering_obj) || is.null(clustering_obj$clusters)) {
    stop("clustering_obj must be a list with a 'clusters' element")
  }

  clusters <- clustering_obj$clusters

  # Validate clusters_to_split
  if (length(clusters_to_split) == 0) {
    message("No clusters specified to split. Returning unchanged clustering_obj.")
    return(clustering_obj)
  }

  existing_ids <- unique(clusters)
  missing <- setdiff(clusters_to_split, existing_ids)
  if (length(missing) > 0) {
    stop(sprintf("Cluster IDs not found in data: %s",
                 paste(missing, collapse = ", ")))
  }

  # Validate logr_matrix
  if (!is.matrix(logr_matrix) && !is.data.frame(logr_matrix)) {
    stop("logr_matrix must be a matrix or data frame")
  }
  logr_matrix <- as.matrix(logr_matrix)

  # Check sample name alignment
  cluster_samples <- names(clusters)
  matrix_samples <- rownames(logr_matrix)

  if (is.null(cluster_samples) || is.null(matrix_samples)) {
    stop("Both clustering_obj$clusters and logr_matrix must have named samples")
  }

  if (!all(cluster_samples %in% matrix_samples)) {
    missing_samples <- setdiff(cluster_samples, matrix_samples)
    stop(sprintf("Samples in clustering_obj not found in logr_matrix: %s",
                 paste(head(missing_samples, 5), collapse = ", ")))
  }

  # Validate method
  if (!method %in% c("hierarchical", "weighted")) {
    stop("method must be 'hierarchical' or 'weighted'")
  }

  # Validate grid requirement for weighted clustering with size weighting
  if (method == "weighted" && use_size_weighting && is.null(grid)) {
    stop("grid is required when method = 'weighted' with use_size_weighting = TRUE")
  }

  # --- Perform splits ---

  split_info <- list()
  next_cluster_id <- max(clusters) + 1

  for (cluster_id in clusters_to_split) {
    # Get samples in this cluster
    sample_indices <- which(clusters == cluster_id)
    sample_names_in_cluster <- names(clusters)[sample_indices]
    n_samples <- length(sample_names_in_cluster)

    message(sprintf("Splitting cluster %d (%d samples) using %s clustering...",
                    cluster_id, n_samples, method))

    # Check if cluster has enough samples to split
    if (n_samples < 2) {
      warning(sprintf("Cluster %d has only %d sample(s) - cannot split. Skipping.",
                      cluster_id, n_samples))
      next
    }

    # Subset logr_matrix to samples in this cluster
    subset_matrix <- logr_matrix[sample_names_in_cluster, , drop = FALSE]

    # Determine k for this split
    effective_k <- k
    if (!is.null(k) && k >= n_samples) {
      warning(sprintf("k=%d >= number of samples (%d) in cluster %d. Reducing k to %d.",
                      k, n_samples, cluster_id, n_samples - 1))
      effective_k <- n_samples - 1
    }

    # Apply clustering
    if (method == "hierarchical") {
      sub_result <- tryCatch({
        hierarchical_clustering(
          logr_matrix = subset_matrix,
          method = "ward.D2",
          use_dynamic = use_dynamic && is.null(effective_k),
          k = effective_k,
          min_cluster_size = min_cluster_size,
          max_module_size = ceiling(n_samples / 2),
          max_recursion_depth = 2
        )
      }, error = function(e) {
        warning(sprintf("Clustering failed for cluster %d: %s. Skipping.",
                        cluster_id, e$message))
        return(NULL)
      })
    } else {
      # method == "weighted"
      sub_result <- tryCatch({
        weighted_hierarchical_clustering(
          logr_matrix = subset_matrix,
          grid = grid,
          discretize = discretize,
          gain_threshold = gain_threshold,
          loss_threshold = loss_threshold,
          use_size_weighting = use_size_weighting,
          min_bins_full_weight = min_bins_full_weight,
          weight_function = weight_function,
          method = "ward.D2",
          use_dynamic = use_dynamic && is.null(effective_k),
          k = effective_k,
          min_cluster_size = min_cluster_size,
          max_module_size = ceiling(n_samples / 2),
          max_recursion_depth = 2
        )
      }, error = function(e) {
        warning(sprintf("Weighted clustering failed for cluster %d: %s. Skipping.",
                        cluster_id, e$message))
        return(NULL)
      })
    }

    if (is.null(sub_result)) {
      next
    }

    sub_clusters <- sub_result$clusters
    n_subclusters <- length(unique(sub_clusters))

    # Check if split produced more than one subcluster
    if (n_subclusters <= 1) {
      warning(sprintf("Cluster %d could not be split (produced only 1 subcluster). Keeping original.",
                      cluster_id))
      next
    }

    # Record split info
    subcluster_sizes <- table(sub_clusters)
    message(sprintf("  Cluster %d split into %d subclusters with sizes: %s",
                    cluster_id, n_subclusters,
                    paste(subcluster_sizes, collapse = ", ")))

    # Integrate subclusters back into main clustering
    # First subcluster keeps the original cluster ID
    # Additional subclusters get new IDs
    unique_sub_ids <- sort(unique(sub_clusters))
    new_ids <- c(cluster_id, seq(from = next_cluster_id,
                                   length.out = n_subclusters - 1))
    id_mapping <- setNames(new_ids, as.character(unique_sub_ids))

    # Update cluster assignments
    for (sample_name in sample_names_in_cluster) {
      old_sub_id <- sub_clusters[sample_name]
      new_id <- id_mapping[as.character(old_sub_id)]
      clusters[sample_name] <- new_id
    }

    # Update next available cluster ID
    next_cluster_id <- max(clusters) + 1

    # Store split info
    split_info[[as.character(cluster_id)]] <- list(
      original_cluster = cluster_id,
      n_samples = n_samples,
      n_subclusters = n_subclusters,
      subcluster_ids = unname(new_ids),
      subcluster_sizes = as.vector(subcluster_sizes),
      method = method
    )
  }

  # --- Renumber if requested ---
  if (renumber && length(split_info) > 0) {
    old_ids <- sort(unique(clusters))
    new_ids <- seq_along(old_ids)
    mapping <- setNames(new_ids, old_ids)
    old_clusters <- clusters
    clusters <- unname(mapping[as.character(clusters)])
    names(clusters) <- names(old_clusters)

    if (!identical(as.integer(old_ids), new_ids)) {
      message(sprintf("Renumbered clusters: %s",
                      paste(sprintf("%d->%d", old_ids, new_ids), collapse = ", ")))
    }
  }

  # --- Final summary ---
  message(sprintf("Result: %d clusters with sizes: %s",
                  length(unique(clusters)),
                  paste(table(clusters), collapse = ", ")))

  # --- Update and return clustering object ---
  clustering_obj$clusters <- clusters
  clustering_obj$split_info <- split_info

  return(clustering_obj)
}

#' Discretize logR values to categorical states
#'
#' Converts continuous logR values to discrete gain/neutral/loss categories.
#' This helps handle purity differences by treating samples with the same
#' CNA profile but different purity as identical.
#'
#' @param logr_matrix LogR matrix (samples x genomic bins)
#' @param gain_threshold LogR value above which is classified as gain (default: 0.2)
#' @param loss_threshold LogR value below which is classified as loss (default: -0.2)
#' @param mode Discretization mode: "ternary" (+1/0/-1) or "binary" (1/0 for altered/neutral)
#' @return Matrix of discretized values (same dimensions as input)
#' @keywords internal
discretize_logr <- function(logr_matrix,
                            gain_threshold = 0.2,
                            loss_threshold = -0.2,
                            mode = "ternary") {

  if (!mode %in% c("ternary", "binary")) {
    stop("mode must be 'ternary' or 'binary'")
  }

  result <- matrix(0, nrow = nrow(logr_matrix), ncol = ncol(logr_matrix))
  rownames(result) <- rownames(logr_matrix)
  colnames(result) <- colnames(logr_matrix)

  if (mode == "ternary") {
    result[logr_matrix > gain_threshold] <- 1
    result[logr_matrix < loss_threshold] <- -1
  } else {
    # Binary mode: any alteration = 1
    result[logr_matrix > gain_threshold | logr_matrix < loss_threshold] <- 1
  }

  return(result)
}

#' Identify contiguous altered regions per sample
#'
#' For each sample, identifies runs of contiguous altered bins (non-zero values
#' in discretized data). Returns region boundaries respecting chromosome boundaries.
#'
#' @param data_matrix Discretized or original logR matrix (samples x genomic bins)
#' @param grid Genomic grid data frame with chr, start, end columns
#' @param altered_threshold Threshold for considering a bin altered (for non-discretized data)
#' @return List with one element per sample, each containing a data frame of regions
#'   with columns: start_bin, end_bin, size, chromosome
#' @keywords internal
identify_altered_regions <- function(data_matrix, grid, altered_threshold = 0.2) {

  n_samples <- nrow(data_matrix)
  n_bins <- ncol(data_matrix)

  # Get chromosome boundaries
  chr_boundaries <- which(diff(as.numeric(factor(grid$chr))) != 0)
  chr_starts <- c(1, chr_boundaries + 1)
  chr_ends <- c(chr_boundaries, n_bins)

  regions_list <- vector("list", n_samples)
  names(regions_list) <- rownames(data_matrix)

  for (i in seq_len(n_samples)) {
    sample_data <- data_matrix[i, ]
    sample_regions <- data.frame(
      start_bin = integer(0),
      end_bin = integer(0),
      size = integer(0),
      chromosome = character(0),
      stringsAsFactors = FALSE
    )

    # Process each chromosome separately
    for (chr_idx in seq_along(chr_starts)) {
      start <- chr_starts[chr_idx]
      end <- chr_ends[chr_idx]
      chr_name <- grid$chr[start]

      # Identify altered bins in this chromosome
      chr_data <- sample_data[start:end]
      is_altered <- abs(chr_data) > altered_threshold

      if (!any(is_altered)) next

      # Find runs of altered bins using rle
      runs <- rle(is_altered)
      run_ends <- cumsum(runs$lengths)
      run_starts <- c(1, run_ends[-length(run_ends)] + 1)

      # Extract altered regions
      altered_runs <- which(runs$values)
      for (run_idx in altered_runs) {
        region_start <- start + run_starts[run_idx] - 1
        region_end <- start + run_ends[run_idx] - 1
        region_size <- region_end - region_start + 1

        sample_regions <- rbind(sample_regions, data.frame(
          start_bin = region_start,
          end_bin = region_end,
          size = region_size,
          chromosome = chr_name,
          stringsAsFactors = FALSE
        ))
      }
    }

    regions_list[[i]] <- sample_regions
  }

  return(regions_list)
}

#' Compute size-based weights for genomic bins
#'
#' Assigns weights to each bin based on the size of the contiguous altered region
#' it belongs to. Smaller regions get reduced weight to minimize influence of
#' potential sequencing artifacts or missed calls.
#'
#' @param data_matrix Data matrix (samples x genomic bins), discretized or original
#' @param grid Genomic grid data frame with chr, start, end columns
#' @param min_bins_full_weight Minimum number of contiguous bins for full weight (default: 5)
#' @param weight_function Weight function type: "linear", "sigmoid", or "step"
#' @param scope Weight scope: "per_sample" computes weights for each sample independently,
#'   "consensus" uses the union of altered regions across samples
#' @param altered_threshold Threshold for considering a bin altered (for non-discretized data)
#' @return Matrix of weights (samples x genomic bins), values in [0, 1]
#' @keywords internal
compute_size_weights <- function(data_matrix,
                                  grid,
                                  min_bins_full_weight = 5,
                                  weight_function = "linear",
                                  scope = "per_sample",
                                  altered_threshold = 0.2) {

  n_samples <- nrow(data_matrix)
  n_bins <- ncol(data_matrix)

  # Initialize weight matrix with 1s (full weight for neutral regions)
  weight_matrix <- matrix(1, nrow = n_samples, ncol = n_bins)
  rownames(weight_matrix) <- rownames(data_matrix)
  colnames(weight_matrix) <- colnames(data_matrix)

  # Weight calculation function
  calc_weight <- function(size, min_bins, func_type) {
    if (func_type == "linear") {
      return(min(size / min_bins, 1.0))
    } else if (func_type == "sigmoid") {
      # Smooth S-curve centered at min_bins/2
      x <- (size - min_bins / 2) / (min_bins / 4)
      return(1 / (1 + exp(-x)))
    } else if (func_type == "step") {
      return(ifelse(size >= min_bins, 1.0, 0.0))
    } else {
      return(1.0)
    }
  }

  if (scope == "per_sample") {
    # Identify regions for each sample
    regions_list <- identify_altered_regions(data_matrix, grid, altered_threshold)

    for (i in seq_len(n_samples)) {
      regions <- regions_list[[i]]
      if (nrow(regions) == 0) next

      for (j in seq_len(nrow(regions))) {
        region <- regions[j, ]
        weight <- calc_weight(region$size, min_bins_full_weight, weight_function)
        weight_matrix[i, region$start_bin:region$end_bin] <- weight
      }
    }

  } else if (scope == "consensus") {
    # Compute consensus altered regions across all samples
    # A bin is "consensus altered" if altered in any sample
    consensus_altered <- apply(abs(data_matrix) > altered_threshold, 2, any)

    # Create a pseudo-sample with consensus alteration pattern
    consensus_matrix <- matrix(as.numeric(consensus_altered), nrow = 1)
    colnames(consensus_matrix) <- colnames(data_matrix)
    rownames(consensus_matrix) <- "consensus"

    consensus_regions <- identify_altered_regions(consensus_matrix, grid, 0.5)[[1]]

    # Apply consensus weights to all samples
    if (nrow(consensus_regions) > 0) {
      for (j in seq_len(nrow(consensus_regions))) {
        region <- consensus_regions[j, ]
        weight <- calc_weight(region$size, min_bins_full_weight, weight_function)
        # Apply weight only to bins that are altered in each sample
        for (i in seq_len(n_samples)) {
          bins_in_region <- region$start_bin:region$end_bin
          altered_in_sample <- abs(data_matrix[i, bins_in_region]) > altered_threshold
          weight_matrix[i, bins_in_region[altered_in_sample]] <- weight
        }
      }
    }
  }

  return(weight_matrix)
}

#' Compute weighted Euclidean distance matrix
#'
#' Computes pairwise weighted Euclidean distances between samples.
#' For each pair, the weight at each bin is the geometric mean of
#' the individual sample weights.
#'
#' @param data_matrix Data matrix (samples x genomic bins)
#' @param weight_matrix Weight matrix (samples x genomic bins), values in [0, 1]
#' @return dist object containing weighted distances
#' @keywords internal
weighted_euclidean_distance <- function(data_matrix, weight_matrix) {

  n_samples <- nrow(data_matrix)
  n_bins <- ncol(data_matrix)

  # Initialize distance matrix
  dist_mat <- matrix(0, nrow = n_samples, ncol = n_samples)
  rownames(dist_mat) <- rownames(data_matrix)
  colnames(dist_mat) <- rownames(data_matrix)

  for (i in seq_len(n_samples - 1)) {
    for (j in (i + 1):n_samples) {
      # Compute combined weights (geometric mean)
      combined_weights <- sqrt(weight_matrix[i, ] * weight_matrix[j, ])

      # Compute weighted squared differences
      diff_sq <- (data_matrix[i, ] - data_matrix[j, ])^2
      weighted_diff_sq <- combined_weights * diff_sq

      # Normalized weighted Euclidean distance
      total_weight <- sum(combined_weights)
      if (total_weight > 0) {
        dist_val <- sqrt(sum(weighted_diff_sq) / total_weight)
      } else {
        dist_val <- 0
      }

      dist_mat[i, j] <- dist_val
      dist_mat[j, i] <- dist_val
    }
  }

  return(as.dist(dist_mat))
}

#' Perform weighted hierarchical clustering on LogR data
#'
#' A clustering function that addresses two common challenges in copy number analysis:
#' \enumerate{
#'   \item \strong{Purity differences}: Samples with identical CNA profiles but different
#'     tumor purity show different logR magnitudes. Discretization converts continuous
#'     values to gain/neutral/loss categories, making purity-different but CNA-identical
#'     samples cluster together.
#'   \item \strong{Small aberrations}: Small altered regions (1-2 bins) may represent
#'     sequencing artifacts or missed calls. Size weighting reduces the influence of
#'     these small regions on clustering.
#' }
#'
#' @param logr_matrix LogR matrix (samples x genomic bins). Required.
#' @param grid Genomic grid data frame from \code{load_ascat_profiles()}.
#'   Required if \code{use_size_weighting = TRUE}. Contains chr, start, end columns
#'   for identifying chromosome boundaries.
#' @param discretize Logical. If TRUE, discretize logR values to gain/neutral/loss
#'   categories before clustering. Default: FALSE.
#' @param gain_threshold LogR value above which is classified as gain. Default: 0.2.
#' @param loss_threshold LogR value below which is classified as loss. Default: -0.2.
#' @param discretize_mode Discretization mode: "ternary" (+1/0/-1 for gain/neutral/loss)
#'   or "binary" (1/0 for altered/neutral). Default: "ternary".
#' @param use_size_weighting Logical. If TRUE, weight bins by the size of their
#'   contiguous altered region. Default: TRUE.
#' @param min_bins_full_weight Minimum number of contiguous altered bins for full weight.
#'   Regions smaller than this get reduced weight. Default: 5.
#' @param weight_function Weight function for region sizes: "linear" (weight = size/min_bins,
#'   capped at 1), "sigmoid" (smooth S-curve transition), or "step" (0 below threshold, 1 above).
#'   Default: "linear".
#' @param weight_scope How to compute region sizes: "per_sample" (independently per sample)
#'   or "consensus" (using union of altered regions across samples). Default: "per_sample".
#' @param method Linkage method for hierarchical clustering. Default: "ward.D2".
#' @param use_dynamic Logical. If TRUE, use dynamic tree cut. Default: TRUE.
#' @param k Number of clusters for fixed cut (if \code{use_dynamic = FALSE}). Default: 20.
#' @param min_cluster_size Minimum cluster size for dynamic tree cut. Default: 10.
#' @param max_module_size Maximum module size before recursive splitting. Default: 100.
#' @param silhouette_threshold Minimum silhouette improvement for accepting splits. Default: 0.05.
#' @param gap_threshold Minimum gap statistic for accepting splits. Default: 0.01
#' @param pvclust_threshold Minimum AU p-value for accepting splits. Default: 0.95.
#' @param max_recursion_depth Maximum depth for recursive splitting. Default: 4.
#'
#' @return A list containing:
#' \describe{
#'   \item{clusters}{Integer vector of cluster assignments (named by sample)}
#'   \item{hclust_obj}{The hclust object (dendrogram)}
#'   \item{dist_matrix}{The (weighted) distance matrix used for clustering}
#'   \item{weight_matrix}{Matrix of size weights (samples x bins), or NULL if not used}
#'   \item{discretized_matrix}{Discretized data matrix, or NULL if not used}
#'   \item{parameters}{List of all parameter values used}
#' }
#'
#' @importFrom cluster silhouette
#' @importFrom dynamicTreeCut cutreeDynamic
#' @export
#' @examples
#' \dontrun{
#' profiles <- load_ascat_profiles("path/to/profiles", resolution = 1e6)
#'
#' # Basic weighted clustering with size weighting
#' result <- weighted_hierarchical_clustering(
#'   logr_matrix = profiles$logr_matrix,
#'   grid = profiles$grid,
#'   use_size_weighting = TRUE,
#'   min_bins_full_weight = 5
#' )
#'
#' # Clustering with discretization (handles purity differences)
#' result_discrete <- weighted_hierarchical_clustering(
#'   logr_matrix = profiles$logr_matrix,
#'   grid = profiles$grid,
#'   discretize = TRUE,
#'   gain_threshold = 0.2,
#'   loss_threshold = -0.2
#' )
#'
#' # Combined approach: discretization + size weighting
#' result_both <- weighted_hierarchical_clustering(
#'   logr_matrix = profiles$logr_matrix,
#'   grid = profiles$grid,
#'   discretize = TRUE,
#'   use_size_weighting = TRUE
#' )
#'
#' # Access results
#' clusters <- result$clusters
#' plot(result$hclust_obj)
#' }
weighted_hierarchical_clustering <- function(logr_matrix,
                                              grid = NULL,
                                              discretize = FALSE,
                                              gain_threshold = 0.2,
                                              loss_threshold = -0.2,
                                              discretize_mode = "ternary",
                                              use_size_weighting = TRUE,
                                              min_bins_full_weight = 5,
                                              weight_function = "linear",
                                              weight_scope = "per_sample",
                                              method = "ward.D2",
                                              use_dynamic = TRUE,
                                              k = 20,
                                              min_cluster_size = 10,
                                              max_module_size = 100,
                                              silhouette_threshold = 0.05,
                                              gap_threshold = 0.01,
                                              pvclust_threshold = 0.95,
                                              max_recursion_depth = 4) {

  # Input validation
  if (!is.matrix(logr_matrix) && !is.data.frame(logr_matrix)) {
    stop("logr_matrix must be a matrix or data frame")
  }
  logr_matrix <- as.matrix(logr_matrix)

  if (nrow(logr_matrix) < 2) {
    stop("logr_matrix must have at least 2 samples (rows)")
  }

  if (use_size_weighting && is.null(grid)) {
    stop("grid is required when use_size_weighting = TRUE")
  }

  if (!is.null(grid) && ncol(logr_matrix) != nrow(grid)) {
    stop("Number of columns in logr_matrix must match number of rows in grid")
  }

  if (!weight_function %in% c("linear", "sigmoid", "step")) {
    stop("weight_function must be 'linear', 'sigmoid', or 'step'")
  }

  if (!weight_scope %in% c("per_sample", "consensus")) {
    stop("weight_scope must be 'per_sample' or 'consensus'")
  }

  if (!discretize_mode %in% c("ternary", "binary")) {
    stop("discretize_mode must be 'ternary' or 'binary'")
  }

  message("Performing weighted hierarchical clustering...")


  # Store parameters
  parameters <- list(
    discretize = discretize,
    gain_threshold = gain_threshold,
    loss_threshold = loss_threshold,
    discretize_mode = discretize_mode,
    use_size_weighting = use_size_weighting,
    min_bins_full_weight = min_bins_full_weight,
    weight_function = weight_function,
    weight_scope = weight_scope,
    method = method,
    use_dynamic = use_dynamic,
    k = k,
    min_cluster_size = min_cluster_size,
    max_module_size = max_module_size,
    silhouette_threshold = silhouette_threshold,
    gap_threshold = gap_threshold,
    pvclust_threshold = pvclust_threshold,
    max_recursion_depth = max_recursion_depth
  )

  # Step 1: Discretization (optional)
  discretized_matrix <- NULL
  working_matrix <- logr_matrix

  if (discretize) {
    message(sprintf("Discretizing logR values (mode: %s, gain > %.2f, loss < %.2f)...",
                    discretize_mode, gain_threshold, loss_threshold))
    discretized_matrix <- discretize_logr(
      logr_matrix,
      gain_threshold = gain_threshold,
      loss_threshold = loss_threshold,
      mode = discretize_mode
    )
    working_matrix <- discretized_matrix

    # Report discretization summary
    n_gains <- sum(discretized_matrix > 0)
    n_losses <- sum(discretized_matrix < 0)
    n_neutral <- sum(discretized_matrix == 0)
    total <- length(discretized_matrix)
    message(sprintf("  Gains: %d (%.1f%%), Losses: %d (%.1f%%), Neutral: %d (%.1f%%)",
                    n_gains, 100 * n_gains / total,
                    n_losses, 100 * n_losses / total,
                    n_neutral, 100 * n_neutral / total))
  }

  # Step 2: Size-based weighting (optional)
  weight_matrix <- NULL

  if (use_size_weighting) {
    message(sprintf("Computing size-based weights (min_bins: %d, function: %s, scope: %s)...",
                    min_bins_full_weight, weight_function, weight_scope))

    # Use appropriate threshold based on whether data is discretized
    altered_threshold <- if (discretize) 0.5 else max(abs(gain_threshold), abs(loss_threshold))

    weight_matrix <- compute_size_weights(
      data_matrix = working_matrix,
      grid = grid,
      min_bins_full_weight = min_bins_full_weight,
      weight_function = weight_function,
      scope = weight_scope,
      altered_threshold = altered_threshold
    )

    # Report weight summary
    n_full_weight <- sum(weight_matrix == 1)
    n_reduced_weight <- sum(weight_matrix > 0 & weight_matrix < 1)
    n_zero_weight <- sum(weight_matrix == 0)
    total <- length(weight_matrix)
    message(sprintf("  Full weight: %d (%.1f%%), Reduced: %d (%.1f%%), Zero: %d (%.1f%%)",
                    n_full_weight, 100 * n_full_weight / total,
                    n_reduced_weight, 100 * n_reduced_weight / total,
                    n_zero_weight, 100 * n_zero_weight / total))
  }

  # Step 3: Compute distance matrix
  if (use_size_weighting && !is.null(weight_matrix)) {
    message("Computing weighted Euclidean distance matrix...")
    dist_matrix <- weighted_euclidean_distance(working_matrix, weight_matrix)
  } else {
    message("Computing standard Euclidean distance matrix...")
    dist_matrix <- dist(working_matrix, method = "euclidean")
  }

  # Step 4: Hierarchical clustering
  message(sprintf("Performing hierarchical clustering (method: %s)...", method))
  hc <- hclust(dist_matrix, method = method)

  # Step 5: Cut tree
  if (use_dynamic) {
    message("Using dynamic tree cut with recursive splitting...")
    clusters <- recursive_dynamic_hac(
      data = working_matrix,
      min_cluster_size = min_cluster_size,
      max_module_size = max_module_size,
      silhouette_threshold = silhouette_threshold,
      gap_threshold = gap_threshold,
      pvclust_threshold = pvclust_threshold,
      max_recursion_depth = max_recursion_depth,
      hierarchical_k = ifelse(is.null(k), 20, k),
      depth = 0,
      dist_matrix = dist_matrix
    )
  } else if (!is.null(k)) {
    message(sprintf("Using fixed tree cut with k=%d", k))
    clusters <- cutree(hc, k = k)
  } else {
    # Use silhouette to find optimal k
    message("Finding optimal k using silhouette analysis...")
    sil_scores <- sapply(2:min(10, nrow(logr_matrix) - 1), function(k_val) {
      clust <- cutree(hc, k = k_val)
      sil <- silhouette(clust, dist_matrix)
      mean(sil[, 3])
    })

    optimal_k <- which.max(sil_scores) + 1
    message(sprintf("Optimal k by silhouette: %d", optimal_k))
    clusters <- cutree(hc, k = optimal_k)
  }

  # Name clusters by sample names
  names(clusters) <- rownames(logr_matrix)

  message(sprintf("Final result: %d clusters", length(unique(clusters))))
  message(sprintf("Cluster sizes: %s", paste(table(clusters), collapse = ", ")))

  return(list(
    clusters = clusters,
    hclust_obj = hc,
    dist_matrix = dist_matrix,
    weight_matrix = weight_matrix,
    discretized_matrix = discretized_matrix,
    parameters = parameters
  ))
}

#' Perform hierarchical clustering on LogR data with dynamic tree cut
#'
#' @param logr_matrix LogR matrix (samples x positions)
#' @param method Linkage method (default: "ward.D2")
#' @param use_dynamic Use dynamic tree cut (default: TRUE)
#' @param k Number of clusters for fixed cut (if use_dynamic = FALSE)
#' @param min_cluster_size Minimum cluster size
#' @param max_module_size Maximum module size for recursive splitting
#' @param silhouette_threshold Minimum silhouette improvement threshold
#' @param gap_threshold Minimum gap statistic threshold
#' @param pvclust_threshold Minimum AU p-value threshold
#' @param max_recursion_depth Maximum recursion depth
#' @param ascat_compat Logical. If TRUE, default to ASCAT workflow settings
#'   when arguments are omitted (e.g., hierarchical_k = 9).
#' @return List with clusters, hclust_obj, and dist_matrix
#' @importFrom cluster silhouette
#' @importFrom dynamicTreeCut cutreeDynamic
#' @export
#' @examples
#' \dontrun{
#' result <- hierarchical_clustering(logr_matrix, min_cluster_size = 10)
#' clusters <- result$clusters
#' }
hierarchical_clustering <- function(logr_matrix,
                                    method = "ward.D2",
                                    use_dynamic = TRUE,
                                    k = 20,
                                    min_cluster_size = 10,
                                    max_module_size = 100,
                                    silhouette_threshold = 0.05,
                                    gap_threshold = 0.01,
                                    pvclust_threshold = 0.95,
                                    max_recursion_depth = 4,
                                    ascat_compat = FALSE) {

  if (ascat_compat) {
    if (missing(method)) method <- "ward.D2"
    if (missing(use_dynamic)) use_dynamic <- TRUE
    if (missing(k)) k <- 9
  }

  message("Performing hierarchical clustering...")

  # Calculate distance matrix
  dist_matrix <- dist(logr_matrix, method = "euclidean")

  # Perform hierarchical clustering
  hc <- hclust(dist_matrix, method = method)

  # Cut tree
  if (use_dynamic) {
    message("Using dynamic tree cut with recursive splitting...")
    clusters <- recursive_dynamic_hac(
      data = logr_matrix,
      min_cluster_size = min_cluster_size,
      max_module_size = max_module_size,
      silhouette_threshold = silhouette_threshold,
      gap_threshold = gap_threshold,
      pvclust_threshold = pvclust_threshold,
      max_recursion_depth = max_recursion_depth,
      hierarchical_k = ifelse(is.null(k), 20, k),
      depth = 0,
      dist_matrix = dist_matrix
    )
  } else if (!is.null(k)) {
    message(sprintf("Using fixed tree cut with k=%d", k))
    clusters <- cutree(hc, k = k)
  } else {
    # Use silhouette to find optimal k
    message("Finding optimal k using silhouette analysis...")
    sil_scores <- sapply(2:min(10, nrow(logr_matrix) - 1), function(k) {
      clust <- cutree(hc, k = k)
      sil <- silhouette(clust, dist_matrix)
      mean(sil[, 3])
    })

    optimal_k <- which.max(sil_scores) + 1
    message(sprintf("Optimal k by silhouette: %d", optimal_k))
    clusters <- cutree(hc, k = optimal_k)
  }

  message(sprintf("Final result: %d clusters", length(unique(clusters))))

  return(list(clusters = clusters, hclust_obj = hc, dist_matrix = dist_matrix))
}
