

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
                                    max_recursion_depth = 4) {

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
