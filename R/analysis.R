
#' Merge similar clusters based on correlation
#'
#' @param logr_matrix LogR matrix
#' @param clusters Initial cluster assignments
#' @param cor_threshold Correlation threshold for merging (default: 0.9)
#' @return Refined cluster assignments
#' @importFrom stats cor
#' @export
#' @examples
#' \dontrun{
#' refined <- refine_clusters_by_correlation(logr_matrix, clusters, cor_threshold = 0.9)
#' }
refine_clusters_by_correlation <- function(logr_matrix, clusters,
                                           cor_threshold = 0.9) {

  message("Refining clusters by correlation...")

  # Calculate cluster centroids
  cluster_ids <- sort(unique(clusters))
  n_clusters <- length(cluster_ids)

  centroids <- matrix(NA, nrow = n_clusters, ncol = ncol(logr_matrix))
  rownames(centroids) <- paste0("Cluster_", cluster_ids)

  for (i in seq_along(cluster_ids)) {
    cid <- cluster_ids[i]
    cluster_samples <- logr_matrix[clusters == cid, , drop = FALSE]
    centroids[i, ] <- colMeans(cluster_samples, na.rm = TRUE)
  }

  # Calculate pairwise correlations between centroids
  cor_matrix <- cor(t(centroids), use = "pairwise.complete.obs")

  # Iterative merging: merge highest correlated pairs one at a time
  # and recalculate centroids after each merge
  refined_clusters <- clusters
  merged <- TRUE

  while (merged) {
    merged <- FALSE

    # Get current cluster IDs
    current_cluster_ids <- sort(unique(refined_clusters))
    n_current <- length(current_cluster_ids)

    if (n_current <= 1) break

    # Recalculate centroids based on current clustering
    current_centroids <- matrix(NA, nrow = n_current, ncol = ncol(logr_matrix))
    rownames(current_centroids) <- paste0("Cluster_", current_cluster_ids)

    for (i in seq_along(current_cluster_ids)) {
      cid <- current_cluster_ids[i]
      cluster_samples <- logr_matrix[refined_clusters == cid, , drop = FALSE]
      current_centroids[i, ] <- colMeans(cluster_samples, na.rm = TRUE)
    }

    # Recalculate correlations
    current_cor_matrix <- cor(t(current_centroids), use = "pairwise.complete.obs")

    # Find highest correlation pair above threshold
    # Set diagonal and lower triangle to -1 to avoid duplicates
    search_matrix <- current_cor_matrix
    diag(search_matrix) <- -1
    search_matrix[lower.tri(search_matrix)] <- -1

    max_cor <- max(search_matrix, na.rm = TRUE)

    if (max_cor > cor_threshold) {
      # Find the indices of max correlation
      max_idx <- which(search_matrix == max_cor, arr.ind = TRUE)[1, ]
      cluster1 <- current_cluster_ids[max_idx[1]]
      cluster2 <- current_cluster_ids[max_idx[2]]

      message(sprintf("Merging cluster %d and %d (correlation: %.3f)",
                      cluster1, cluster2, max_cor))

      # Merge cluster2 into cluster1
      refined_clusters[refined_clusters == cluster2] <- cluster1
      merged <- TRUE
    }
  }

  # Renumber clusters sequentially
  unique_clusters <- sort(unique(refined_clusters))
  cluster_mapping <- setNames(1:length(unique_clusters), unique_clusters)
  refined_clusters <- cluster_mapping[as.character(refined_clusters)]

  message(sprintf("Refined from %d to %d clusters",
                  n_clusters, length(unique(refined_clusters))))

  return(refined_clusters)
}
