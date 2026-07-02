
#' Impute missing values in a matrix
#'
#' Two methods are available:
#' \describe{
#'   \item{"row_mean"}{Two-pass imputation. First, each NA is replaced with the
#'     mean of the non-NA values in its row (i.e., that sample's average). Then,
#'     any remaining NAs (from rows that were entirely NA) are replaced with the
#'     overall matrix mean. This preserves each sample's baseline ploidy/LogR level.}
#'   \item{"scalar"}{All NAs are replaced with a single fixed value
#'     (default 0, i.e., the diploid LogR baseline).}
#' }
#'
#' @param mat Numeric matrix (samples x positions) potentially containing NAs.
#' @param method Character: \code{"row_mean"} (default) or \code{"scalar"}.
#' @param scalar_value Numeric value used when \code{method = "scalar"}.
#'   Ignored when \code{method = "row_mean"}. Default 0.
#' @return Matrix with NAs replaced according to the chosen method.
#' @export
#' @examples
#' \dontrun{
#' mat <- matrix(c(1, NA, 3, NA, 5, 6), nrow = 2)
#' impute_missing(mat, method = "row_mean")
#' impute_missing(mat, method = "scalar", scalar_value = 0)
#' }
impute_missing <- function(mat, method = c("row_mean", "scalar"), scalar_value = 0) {
  method <- match.arg(method)

  if (!any(is.na(mat))) {
    return(mat)
  }

  if (method == "scalar") {
    mat[is.na(mat)] <- scalar_value
  } else {
    # Pass 1: replace NAs with per-row (per-sample) mean
    for (i in seq_len(nrow(mat))) {
      row_na <- is.na(mat[i, ])
      if (any(row_na) && !all(row_na)) {
        mat[i, row_na] <- mean(mat[i, !row_na], na.rm = TRUE)
      }
    }
    # Pass 2: any remaining NAs (entirely-NA rows) get the overall mean
    if (any(is.na(mat))) {
      overall_mean <- mean(mat, na.rm = TRUE)
      if (is.finite(overall_mean)) {
        mat[is.na(mat)] <- overall_mean
      } else {
        mat[is.na(mat)] <- 0
      }
    }
  }

  mat
}

#' Compute UMAP embedding
#'
#' @param logr_matrix LogR matrix
#' @param n_neighbors Number of neighbors
#' @param min_dist Minimum distance
#' @return Matrix with UMAP coordinates
#' @importFrom umap umap umap.defaults
#' @export
#' @examples
#' \dontrun{
#' umap_coords <- compute_umap(logr_matrix, n_neighbors = 15)
#' }
compute_umap <- function(logr_matrix, n_neighbors = 15, min_dist = 0.1) {

  message("Computing UMAP embedding...")

  umap_config <- umap.defaults
  umap_config$n_neighbors <- n_neighbors
  umap_config$min_dist <- min_dist

  umap_result <- umap(logr_matrix, config = umap_config)

  return(umap_result$layout)
}

#' Calculate Local Indicators of Spatial Association (LISA) for a cluster
#'
#' Internal function used by calculate_morans_i.
#'
#' @param spatial_data Data frame with sample_name, x, y, cluster columns
#' @param cluster_id Cluster ID to analyze
#' @param W_std Row-standardized spatial weights matrix
#' @param n_permutations Number of permutations for significance testing
#' @return Data frame with Local Moran's I, p-values, and quadrant classification
#' @keywords internal
calculate_local_morans_i <- function(spatial_data, cluster_id, W_std, n_permutations = 999) {

  # Create binary indicator for cluster membership
  z <- ifelse(spatial_data$cluster == cluster_id, 1, 0)
  n <- length(z)
  z_mean <- mean(z)
  z_dev <- z - z_mean

  if (sum(z_dev^2) == 0) return(NULL)

  # Variance
  m2 <- sum(z_dev^2) / n

  # Calculate Local Moran's I for each location i
  local_i <- numeric(n)
  for (i in 1:n) {
    # Local Moran's I_i = (z_i - mean(z)) / m2 * sum_j(w_ij * (z_j - mean(z)))
    local_i[i] <- (z_dev[i] / m2) * sum(W_std[i, ] * z_dev)
  }

  # Permutation test for each location
  p_values <- numeric(n)
  for (i in 1:n) {
    # Permute neighbors of location i only
    I_perm <- replicate(n_permutations, {
      z_perm <- z
      neighbors <- which(W_std[i, ] > 0)
      if (length(neighbors) > 0) {
        z_perm[neighbors] <- sample(z[neighbors])
      }
      z_dev_perm <- z_perm - mean(z_perm)
      m2_perm <- sum(z_dev_perm^2) / n
      (z_dev_perm[i] / m2_perm) * sum(W_std[i, ] * z_dev_perm)
    })

    p_values[i] <- sum(abs(I_perm) >= abs(local_i[i])) / n_permutations
  }

  # Classify LISA quadrants
  # HH: High-High (cluster member surrounded by cluster members)
  # HL: High-Low (cluster member surrounded by non-members)
  # LH: Low-High (non-member surrounded by cluster members)
  # LL: Low-Low (non-member surrounded by non-members)

  spatial_lag <- as.numeric(W_std %*% z)
  quadrant <- rep("NS", n)  # Not Significant

  for (i in 1:n) {
    if (p_values[i] < 0.05) {
      if (z[i] > z_mean && spatial_lag[i] > mean(spatial_lag)) {
        quadrant[i] <- "HH"  # Hotspot
      } else if (z[i] > z_mean && spatial_lag[i] < mean(spatial_lag)) {
        quadrant[i] <- "HL"  # Spatial outlier (high)
      } else if (z[i] < z_mean && spatial_lag[i] > mean(spatial_lag)) {
        quadrant[i] <- "LH"  # Spatial outlier (low)
      } else {
        quadrant[i] <- "LL"  # Coldspot
      }
    }
  }

  data.frame(
    sample_name = spatial_data$sample_name,
    x = spatial_data$x,
    y = spatial_data$y,
    cluster = spatial_data$cluster,
    is_member = z,
    local_i = local_i,
    p_value = p_values,
    spatial_lag = spatial_lag,
    quadrant = quadrant,
    stringsAsFactors = FALSE
  )
}

#' Calculate global and local Moran's I for spatial autocorrelation
#'
#' @param spatial_data Data frame with x, y, cluster columns
#' @param max_distance Maximum distance for neighbor definition
#' @param n_permutations Number of permutations for significance testing
#' @param calculate_local Whether to calculate Local Moran's I (LISA)
#' @return List with global results (data frame) and local results (list)
#' @importFrom dplyr %>% filter bind_rows
#' @export
#' @examples
#' \dontrun{
#' morans <- calculate_morans_i(spatial_data, max_distance = 1000)
#' }
calculate_morans_i <- function(spatial_data, max_distance = NULL, n_permutations = 100, calculate_local = TRUE) {

  # Filter to samples with spatial coordinates
  spatial_data <- spatial_data %>%
    filter(!is.na(x), !is.na(y), !is.na(cluster))

  if (nrow(spatial_data) < 3) {
    message("Not enough samples with spatial coordinates for Moran's I")
    return(NULL)
  }

  # Convert x, y to numeric
  spatial_data$x <- as.numeric(spatial_data$x)
  spatial_data$y <- as.numeric(spatial_data$y)

  # Calculate distance matrix
  coords <- as.matrix(spatial_data[, c("x", "y")])
  dist_matrix <- as.matrix(dist(coords))

  # Create spatial weights matrix (inverse distance)
  # W[i,j] = 1/d[i,j] if d < max_distance, 0 otherwise
  if (is.null(max_distance)) {
    max_distance <- max(dist_matrix) * 0.5  # Use half of max distance as default
  }

  W <- 1 / (dist_matrix + diag(nrow(dist_matrix)) * 1e10)  # Avoid division by zero on diagonal
  W[dist_matrix > max_distance] <- 0
  diag(W) <- 0

  # Row-standardize weights
  row_sums <- rowSums(W)
  W_std <- W / ifelse(row_sums > 0, row_sums, 1)

  # Helper function to calculate Moran's I
  calc_morans_i <- function(z, W_std) {
    n <- length(z)
    z_mean <- mean(z)
    z_dev <- z - z_mean
    if (sum(z_dev^2) == 0) return(NA)

    # Calculate numerator: sum of weighted cross-products
    numerator <- 0
    for (i in 1:n) {
      for (j in 1:n) {
        numerator <- numerator + W_std[i, j] * z_dev[i] * z_dev[j]
      }
    }

    # Calculate denominator: variance * n
    denominator <- sum(z_dev^2)

    # Moran's I = (n / S0) * (numerator / denominator)
    # With row-standardized weights, S0 = n
    I <- (n / n) * (numerator / denominator)

    return(I)
  }

  # Calculate Moran's I for each cluster
  cluster_ids <- sort(unique(spatial_data$cluster))

  results <- lapply(cluster_ids, function(cid) {
    # Create binary indicator: 1 if cluster matches, 0 otherwise
    z <- ifelse(spatial_data$cluster == cid, 1, 0)
    n_samples <- sum(z)

    if (n_samples < 2) {
      return(data.frame(
        cluster = cid,
        morans_i = NA,
        random_mean = NA,
        random_sd = NA,
        p_value = NA,
        n_samples = n_samples,
        interpretation = "Too few samples"
      ))
    }

    # Calculate observed Moran's I
    I_obs <- calc_morans_i(z, W_std)

    if (is.na(I_obs)) {
      return(data.frame(
        cluster = cid,
        morans_i = NA,
        random_mean = NA,
        random_sd = NA,
        p_value = NA,
        n_samples = n_samples,
        interpretation = "Insufficient data"
      ))
    }

    # Permutation test: shuffle cluster assignments
    I_random <- replicate(n_permutations, {
      z_perm <- sample(z)
      calc_morans_i(z_perm, W_std)
    })
    I_random <- I_random[!is.na(I_random)]

    if (length(I_random) == 0) {
      return(data.frame(
        cluster = cid,
        morans_i = I_obs,
        random_mean = NA,
        random_sd = NA,
        p_value = NA,
        n_samples = n_samples,
        interpretation = "Permutation failed"
      ))
    }

    # Calculate empirical p-value
    p_value <- sum(abs(I_random) >= abs(I_obs)) / length(I_random)

    # Interpretation
    random_mean <- mean(I_random)
    random_sd <- sd(I_random)

    interpretation <- if (p_value < 0.001) {
      if (I_obs > random_mean) "Highly clustered***" else "Highly dispersed***"
    } else if (p_value < 0.01) {
      if (I_obs > random_mean) "Clustered**" else "Dispersed**"
    } else if (p_value < 0.05) {
      if (I_obs > random_mean) "Clustered*" else "Dispersed*"
    } else {
      "Random"
    }

    data.frame(
      cluster = cid,
      morans_i = I_obs,
      random_mean = random_mean,
      random_sd = random_sd,
      p_value = p_value,
      n_samples = n_samples,
      interpretation = interpretation
    )
  })

  results_df <- bind_rows(results)

  # Calculate Local Moran's I (LISA) if requested
  local_results <- NULL
  if (calculate_local && !is.null(spatial_data$sample_name)) {
    message("Calculating Local Moran's I (LISA) for each cluster...")
    local_results <- lapply(cluster_ids, function(cid) {
      n_samples <- sum(spatial_data$cluster == cid)
      if (n_samples < 2) return(NULL)

      lisa <- calculate_local_morans_i(
        spatial_data = spatial_data,
        cluster_id = cid,
        W_std = W_std,
        n_permutations = min(n_permutations, 999)  # LISA needs more permutations
      )

      if (!is.null(lisa)) {
        lisa$cluster_id <- cid
      }
      return(lisa)
    })
    names(local_results) <- paste0("cluster_", cluster_ids)
    local_results <- local_results[!sapply(local_results, is.null)]
  }

  # Return both global and local results
  return(list(
    global = results_df,
    local = local_results
  ))
}

#' Plot Local Indicators of Spatial Association (LISA) results
#'
#' @param lisa_results List of data frames with Local Moran's I results
#' @param tiles GeoJSON tiles for spatial context
#' @param output_file Output PDF path
#' @return Summary table of LISA results
#' @importFrom ggplot2 ggplot aes geom_sf geom_point geom_hline geom_vline scale_color_manual scale_size_continuous facet_wrap theme_minimal labs theme element_rect annotate
#' @importFrom dplyr bind_rows %>% filter group_by mutate ungroup summarise arrange
#' @importFrom tidyr %>%
#' @importFrom grDevices pdf dev.off
#' @importFrom grid grid.newpage grid.text gpar
#' @importFrom gridExtra grid.table
#' @export
#' @examples
#' \dontrun{
#' summary <- plot_lisa(lisa_results, tiles, "lisa.pdf")
#' }
plot_lisa <- function(lisa_results, tiles, output_file) {

  message("Generating LISA (Local Moran's I) plots...")

  if (is.null(lisa_results) || length(lisa_results) == 0) {
    message("No LISA results to plot")
    return(NULL)
  }

  # Combine all LISA results
  lisa_combined <- bind_rows(lisa_results, .id = "cluster_name")

  # Filter to significant locations only
  lisa_sig <- lisa_combined %>%
    filter(p_value < 0.05)

  if (nrow(lisa_sig) == 0) {
    message("No significant local spatial associations found")
    return(NULL)
  }

  # LISA quadrant colors
  quadrant_colors <- c(
    "HH" = "#d7191c",  # Red - Hotspot (cluster surrounded by cluster)
    "HL" = "#fdae61",  # Orange - Outlier (cluster surrounded by non-cluster)
    "LH" = "#abd9e9",  # Light blue - Outlier (non-cluster surrounded by cluster)
    "LL" = "#2c7bb6",  # Blue - Coldspot (non-cluster surrounded by non-cluster)
    "NS" = "#d9d9d9"   # Gray - Not significant
  )

  # Plot 1: LISA quadrants on spatial map
  p1 <- ggplot() +
    geom_sf(data = tiles, fill = "gray95", color = "gray80", size = 0.2) +
    geom_point(data = lisa_sig, aes(x = x, y = y, color = quadrant, size = abs(local_i)),
               alpha = 0.8) +
    scale_color_manual(
      values = quadrant_colors,
      name = "LISA Category",
      labels = c(
        "HH" = "Hotspot (HH)",
        "HL" = "Outlier (HL)",
        "LH" = "Outlier (LH)",
        "LL" = "Coldspot (LL)"
      )
    ) +
    scale_size_continuous(name = "|Local I|", range = c(1, 5)) +
    facet_wrap(~ cluster_name, ncol = 3) +
    theme_minimal() +
    labs(
      title = "Local Indicators of Spatial Association (LISA)",
      subtitle = "Significant locations only (p < 0.05)",
      x = "X coordinate",
      y = "Y coordinate"
    ) +
    theme(
      legend.position = "bottom",
      panel.border = element_rect(color = "gray80", fill = NA)
    )

  # Plot 2: Moran scatterplot (standardized value vs spatial lag)
  # Calculate standardized values for visualization
  lisa_combined <- lisa_combined %>%
    group_by(cluster_name) %>%
    mutate(
      z_standardized = (is_member - mean(is_member)) / sd(is_member),
      spatial_lag_standardized = (spatial_lag - mean(spatial_lag)) / sd(spatial_lag)
    ) %>%
    ungroup()

  p2 <- ggplot(lisa_combined, aes(x = z_standardized, y = spatial_lag_standardized, color = quadrant)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = abs(local_i)), alpha = 0.6) +
    scale_color_manual(values = quadrant_colors, name = "LISA Category") +
    scale_size_continuous(name = "|Local I|", range = c(1, 4)) +
    facet_wrap(~ cluster_name, ncol = 3) +
    theme_minimal() +
    labs(
      title = "Moran Scatterplot",
      subtitle = "Standardized value vs. spatial lag (neighbor average)",
      x = expression("Standardized value ("*z[i]*")"),
      y = expression("Spatial lag ("*sum(w[ij]*z[j], j)*")")
    ) +
    theme(
      legend.position = "bottom",
      panel.border = element_rect(color = "gray80", fill = NA)
    ) +
    annotate("text", x = -Inf, y = Inf, label = "LH\n(Outlier)", hjust = -0.1, vjust = 1.5, color = "gray50", size = 3) +
    annotate("text", x = Inf, y = Inf, label = "HH\n(Hotspot)", hjust = 1.1, vjust = 1.5, color = "gray50", size = 3) +
    annotate("text", x = -Inf, y = -Inf, label = "LL\n(Coldspot)", hjust = -0.1, vjust = -0.5, color = "gray50", size = 3) +
    annotate("text", x = Inf, y = -Inf, label = "HL\n(Outlier)", hjust = 1.1, vjust = -0.5, color = "gray50", size = 3)

  # Plot 3: Local Moran's I values bar chart (top samples only)
  # Simplify sample names - extract just the well identifier
  lisa_sig <- lisa_sig %>%
    mutate(
      simple_name = sub(".*_(\\d+[A-H])$", "\\1", sample_name),  # Extract well (e.g., "1A", "12H")
      abs_local_i = abs(local_i)
    )

  # Keep only top 20% of samples per cluster by |Local I|
  lisa_top <- lisa_sig %>%
    group_by(cluster_name) %>%
    arrange(desc(abs_local_i)) %>%
    slice_head(prop = 0.20) %>%  # Top 20%
    ungroup()

  if (nrow(lisa_top) == 0) {
    message("No samples to display in bar chart after filtering")
    p3 <- ggplot() + theme_void() +
      ggtitle("Local Moran's I Values", subtitle = "No significant samples after filtering")
  } else {
    p3 <- ggplot(lisa_top, aes(x = reorder(simple_name, local_i), y = local_i, fill = quadrant)) +
      geom_col() +
      scale_fill_manual(values = quadrant_colors, name = "LISA Category") +
      coord_flip() +
      facet_wrap(~ cluster_name, scales = "free_y", ncol = 2) +
      theme_minimal() +
      labs(
        title = "Local Moran's I Values (Top 20% per cluster)",
        subtitle = paste0("Samples with strongest spatial associations (n=", nrow(lisa_top), ")"),
        x = "Well",
        y = "Local Moran's I"
      ) +
      theme(
        legend.position = "bottom",
        axis.text.y = element_text(size = 8)
      )
  }

  # Save plots
  pdf(output_file, width = 14, height = 10)
  print(p1)
  print(p2)
  print(p3)

  # Summary table
  summary_table <- lisa_sig %>%
    group_by(cluster_name, quadrant) %>%
    summarise(
      n = n(),
      mean_local_i = mean(local_i),
      .groups = "drop"
    ) %>%
    arrange(cluster_name, quadrant)

  grid.newpage()
  grid.text("LISA Summary Table", y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  grid.table(summary_table, rows = NULL)

  dev.off()

  message(sprintf("Saved LISA plots: %s", output_file))
  return(summary_table)
}
