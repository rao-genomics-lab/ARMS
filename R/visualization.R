
#' Create consistent colorblind-friendly cluster color palette using RColorBrewer
#'
#' @param clusters Vector of cluster assignments
#' @param max_clusters Maximum number of clusters to support
#' @return Named vector of colors for each cluster ID
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @examples
#' clusters <- c(1, 1, 2, 2, 3, 3)
#' colors <- create_cluster_colors(clusters)
create_cluster_colors <- function(clusters, max_clusters = 24) {
  cluster_ids <- sort(unique(clusters))

  # Create a fixed palette for cluster numbers 1-24 (or max_clusters)
  # This ensures cluster 1 always gets the same color, cluster 2 always gets the same color, etc.
  full_palette <- c(
    RColorBrewer::brewer.pal(8, "Set1"),
    RColorBrewer::brewer.pal(8, "Set2"),
    RColorBrewer::brewer.pal(8, "Dark2")
  )[1:max_clusters]

  # Map cluster IDs to colors from the fixed palette
  # Use cluster ID as index (cluster 1 = color 1, cluster 2 = color 2, etc.)
  colors <- full_palette[cluster_ids]
  names(colors) <- as.character(cluster_ids)

  return(colors)
}

#' Create chromosome annotation for heatmaps with alternating colors
#'
#' @param grid Genomic grid with chr and position columns
#' @param matrix Data matrix for column dimension matching
#' @return HeatmapAnnotation object with chromosome labels
#' @importFrom ComplexHeatmap HeatmapAnnotation anno_mark
#' @importFrom grid gpar unit
#' @keywords internal
create_chr_annotation <- function(grid, matrix) {

  # Ensure grid matches matrix dimensions
  if (nrow(grid) != ncol(matrix)) {
    # Truncate or subset grid to match matrix
    if (nrow(grid) > ncol(matrix)) {
      grid <- grid[1:ncol(matrix), ]
    } else {
      warning("Grid has fewer positions than matrix columns")
      return(NULL)
    }
  }

  # Create a vector of chromosome labels for each position
  chr_vector <- grid$chr

  # Create alternating colors for chromosomes
  unique_chrs <- unique(chr_vector)
  chr_colors <- rep(c("grey90", "grey70"), length.out = length(unique_chrs))
  names(chr_colors) <- unique_chrs

  # Find chromosome boundaries and midpoints for labels
  chr_changes <- which(diff(as.numeric(factor(chr_vector, levels = unique_chrs))) != 0)
  chr_starts <- c(1, chr_changes + 1)
  chr_ends <- c(chr_changes, length(chr_vector))
  chr_midpoints <- (chr_starts + chr_ends) / 2

  # Create chromosome labels (strip "chr" prefix for cleaner display)
  chr_labels <- gsub("chr", "", unique_chrs)

  # Create simple block annotation with labels
  chr_annotation <- HeatmapAnnotation(
    Chromosome = chr_vector,
    col = list(Chromosome = chr_colors),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    simple_anno_size = unit(3, "mm"),
    annotation_label = list(Chromosome = ""),
    # Add chromosome number labels
    foo = anno_mark(
      at = chr_midpoints,
      labels = chr_labels,
      which = "column",
      side = "top",
      labels_gp = gpar(fontsize = 8),
      link_height = unit(2, "mm"),
      padding = unit(0.5, "mm")
    )
  )

  return(chr_annotation)
}

#' Plot LogR heatmap with hierarchical clustering
#'
#' @param logr_matrix LogR matrix
#' @param grid Genomic grid
#' @param clusters Cluster assignments
#' @param output_file Output PDF path
#' @param copy_number_matrix Copy number matrix (optional, for ploidy annotation)
#' @param logr_range Color range for LogR
#' @param width Width in inches
#' @param height Height in inches
#' @importFrom ComplexHeatmap Heatmap rowAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit
#' @importFrom grDevices pdf dev.off
#' @export
#' @examples
#' \dontrun{
#' plot_logr_heatmap(logr_matrix, grid, clusters, "heatmap.pdf")
#' }
plot_logr_heatmap <- function(logr_matrix, grid, clusters, output_file,
                              copy_number_matrix = NULL,
                              logr_range = 2.0, width = 20, height = 8) {

  message("Generating LogR heatmap...")
  message(sprintf("Matrix dimensions: %d rows x %d cols", nrow(logr_matrix), ncol(logr_matrix)))
  message(sprintf("Grid dimensions: %d positions", nrow(grid)))

  # Create consistent colorblind-friendly cluster colors
  cluster_colors <- create_cluster_colors(clusters)

  # Order rows by cluster to group samples together (numerically 1, 2, 3, ...)
  # Convert to numeric to ensure proper ordering
  cluster_order <- order(as.numeric(clusters), rownames(logr_matrix))
  logr_matrix_ordered <- logr_matrix[cluster_order, , drop = FALSE]
  clusters_ordered <- clusters[cluster_order]

  # Calculate ploidy (mean copy number) if copy_number_matrix provided
  ploidy_ordered <- NULL
  if (!is.null(copy_number_matrix)) {
    ploidy <- rowMeans(copy_number_matrix, na.rm = TRUE)
    ploidy_ordered <- ploidy[cluster_order]
  }

  # Create annotations
  if (!is.null(ploidy_ordered)) {
    cluster_annotation <- rowAnnotation(
      Cluster = as.character(clusters_ordered),
      Ploidy = ploidy_ordered,
      col = list(
        Cluster = cluster_colors,
        Ploidy = colorRamp2(c(1.5, 2, 2.5, 3, 4), c("blue", "lightblue", "white", "pink", "red"))
      ),
      show_annotation_name = TRUE,
      annotation_name_side = "top",
      annotation_legend_param = list(
        Cluster = list(title = "Cluster"),
        Ploidy = list(title = "Ploidy")
      )
    )
  } else {
    cluster_annotation <- rowAnnotation(
      Cluster = as.character(clusters_ordered),
      col = list(Cluster = cluster_colors),
      show_annotation_name = TRUE,
      annotation_name_side = "top",
      annotation_legend_param = list(
        Cluster = list(title = "Cluster")
      )
    )
  }

  # Create chromosome annotation (or skip if dimensions don't match)
  chr_annotation <- NULL
  if (nrow(grid) == ncol(logr_matrix_ordered)) {
    chr_annotation <- create_chr_annotation(grid, logr_matrix_ordered)
  } else {
    message(sprintf("Warning: Grid/matrix dimension mismatch (%d vs %d), skipping chromosome annotation",
                    nrow(grid), ncol(logr_matrix_ordered)))
  }

  # Create row split to group by cluster with dendrogram within each cluster
  cluster_levels <- sort(unique(as.numeric(clusters_ordered)))
  row_split <- factor(clusters_ordered, levels = cluster_levels)

  # Create heatmap with clusters grouped and dendrogram within each
  ht <- Heatmap(
    logr_matrix_ordered,
    name = "LogR",
    col = colorRamp2(c(-logr_range, 0, logr_range), c("blue", "white", "red")),
    use_raster = TRUE,
    raster_quality = 2,
    raster_device = "png",
    left_annotation = cluster_annotation,
    top_annotation = chr_annotation,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    row_split = row_split,
    row_title = "C%s",
    row_title_side = "left",
    row_title_gp = gpar(fontsize = 10),
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title = "Genomic Coordinates",
    column_title_gp = gpar(fontsize = 12),
    column_title_side = "bottom",
    heatmap_width = unit(width, "inches"),
    heatmap_height = unit(height, "inches"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 10)
    )
  )

  # Save - match legacy exactly
  pdf(output_file, width = width + 2, height = height + 2)
  draw(ht)
  dev.off()

  message(sprintf("Saved heatmap: %s", output_file))
}

#' Plot interleaved single-cell and pseudobulk comparison heatmap
#'
#' @param sc_logr_matrix Single-cell LogR matrix
#' @param pb_logr_matrix Pseudobulk LogR matrix
#' @param grid Genomic grid
#' @param sc_clusters Single-cell cluster assignments
#' @param output_file Output PDF path
#' @param sc_cn_matrix Single-cell copy number matrix (for ploidy)
#' @param pb_cn_matrix Pseudobulk copy number matrix (for ploidy)
#' @param logr_range Color range for LogR
#' @param width Width in inches
#' @param height Height in inches
#' @importFrom ComplexHeatmap Heatmap rowAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom grDevices pdf dev.off
#' @export
#' @examples
#' \dontrun{
#' plot_sc_pb_comparison_heatmap(sc_logr, pb_logr, grid, clusters, "comparison.pdf")
#' }
plot_sc_pb_comparison_heatmap <- function(sc_logr_matrix, pb_logr_matrix, grid,
                                          sc_clusters, output_file,
                                          sc_cn_matrix = NULL, pb_cn_matrix = NULL,
                                          logr_range = 2.0, width = 20, height = 8) {

  message("Generating single-cell vs pseudobulk comparison heatmap...")

  # Create consistent colorblind-friendly cluster colors
  cluster_colors <- create_cluster_colors(sc_clusters)

  # Get unique clusters in numeric order (1, 2, 3, ..., 9)
  cluster_ids <- sort(unique(as.numeric(sc_clusters)))

  # Extract cluster IDs from pseudobulk rownames and sort matrix by cluster ID
  pb_rownames <- rownames(pb_logr_matrix)
  pb_cluster_ids <- as.numeric(gsub("cluster_(\\d+).*", "\\1", pb_rownames))

  # Sort pseudobulk matrix rows by cluster ID to ensure correct order
  pb_order <- order(pb_cluster_ids)
  pb_logr_matrix <- pb_logr_matrix[pb_order, , drop = FALSE]
  if (!is.null(pb_cn_matrix)) {
    pb_cn_matrix <- pb_cn_matrix[pb_order, , drop = FALSE]
  }
  pb_cluster_ids <- pb_cluster_ids[pb_order]

  # Create named mapping from cluster ID to PB row index (now in sorted order)
  pb_cluster_map <- setNames(1:nrow(pb_logr_matrix), pb_cluster_ids)

  # Calculate ploidy for SC and PB if matrices provided
  sc_ploidy <- NULL
  pb_ploidy <- NULL
  if (!is.null(sc_cn_matrix)) {
    sc_ploidy <- rowMeans(sc_cn_matrix, na.rm = TRUE)
  }
  if (!is.null(pb_cn_matrix)) {
    pb_ploidy <- rowMeans(pb_cn_matrix, na.rm = TRUE)
  }

  # Build interleaved matrix: for each cluster, add SC samples then PB sample
  interleaved_rows <- list()
  interleaved_labels <- c()
  interleaved_types <- c()
  interleaved_clusters <- c()
  interleaved_ploidy <- c()

  for (cid in cluster_ids) {
    # Get single-cell samples for this cluster
    sc_samples_idx <- which(sc_clusters == cid)
    sc_subset <- sc_logr_matrix[sc_samples_idx, , drop = FALSE]

    # Add SC samples
    for (i in 1:nrow(sc_subset)) {
      interleaved_rows[[length(interleaved_rows) + 1]] <- sc_subset[i, ]
      interleaved_labels <- c(interleaved_labels, paste0("SC_", cid, "_", i))
      interleaved_types <- c(interleaved_types, "Single-cell")
      interleaved_clusters <- c(interleaved_clusters, cid)

      # Add ploidy if available
      if (!is.null(sc_ploidy)) {
        interleaved_ploidy <- c(interleaved_ploidy, sc_ploidy[sc_samples_idx[i]])
      }
    }

    # Add pseudobulk sample - find matching PB row by cluster ID
    pb_row_idx <- pb_cluster_map[as.character(cid)]
    if (!is.na(pb_row_idx)) {
      n_sc_in_cluster <- nrow(sc_subset)
      n_pb_copies <- max(2, round(n_sc_in_cluster * 0.25))

      for (rep in 1:n_pb_copies) {
        interleaved_rows[[length(interleaved_rows) + 1]] <- pb_logr_matrix[pb_row_idx, ]
        interleaved_labels <- c(interleaved_labels, paste0("PB_", cid, "_copy", rep))
        interleaved_types <- c(interleaved_types, "Pseudobulk")
        interleaved_clusters <- c(interleaved_clusters, cid)

        # Add ploidy if available
        if (!is.null(pb_ploidy)) {
          interleaved_ploidy <- c(interleaved_ploidy, pb_ploidy[pb_row_idx])
        }
      }
    } else {
      message(sprintf("Warning: No pseudobulk profile found for cluster %d", cid))
    }
  }

  # Combine into matrix
  combined_matrix <- do.call(rbind, interleaved_rows)
  rownames(combined_matrix) <- interleaved_labels

  message(sprintf("Combined matrix dimensions: %d rows x %d cols",
                  nrow(combined_matrix), ncol(combined_matrix)))

  # Create annotations
  type_colors <- c("Single-cell" = "lightblue", "Pseudobulk" = "orange")

  if (length(interleaved_ploidy) > 0) {
    ha <- rowAnnotation(
      Cluster = as.character(interleaved_clusters),
      Type = interleaved_types,
      Ploidy = interleaved_ploidy,
      col = list(
        Cluster = cluster_colors,
        Type = type_colors,
        Ploidy = colorRamp2(c(1.5, 2, 2.5, 3, 4), c("blue", "lightblue", "white", "pink", "red"))
      ),
      show_annotation_name = TRUE,
      annotation_name_side = "top",
      annotation_legend_param = list(
        Cluster = list(title = "Cluster"),
        Type = list(title = "Type"),
        Ploidy = list(title = "Ploidy")
      )
    )
  } else {
    ha <- rowAnnotation(
      Cluster = as.character(interleaved_clusters),
      Type = interleaved_types,
      col = list(
        Cluster = cluster_colors,
        Type = type_colors
      ),
      show_annotation_name = TRUE,
      annotation_name_side = "top",
      annotation_legend_param = list(
        Cluster = list(title = "Cluster"),
        Type = list(title = "Type")
      )
    )
  }

  # Create chromosome annotation
  chr_annotation <- NULL
  if (nrow(grid) == ncol(combined_matrix)) {
    chr_annotation <- create_chr_annotation(grid, combined_matrix)
  } else {
    message(sprintf("Warning: Grid/matrix dimension mismatch (%d vs %d), skipping chromosome annotation",
                    nrow(grid), ncol(combined_matrix)))
  }

  # Create row split to group by cluster
  cluster_levels <- sort(unique(as.numeric(interleaved_clusters)))
  row_split <- factor(interleaved_clusters, levels = cluster_levels)

  # Create heatmap with clusters grouped (SC/PB pairs together within each cluster)
  ht <- Heatmap(
    combined_matrix,
    name = "LogR",
    col = colorRamp2(c(-logr_range, 0, logr_range), c("blue", "white", "red")),
    use_raster = TRUE,
    raster_quality = 2,
    raster_device = "png",
    left_annotation = ha,
    top_annotation = chr_annotation,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = row_split,
    row_title = "Cluster %s",
    row_title_side = "left",
    row_title_gp = gpar(fontsize = 10),
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_title = "Genomic Coordinates (SC + PB per cluster)",
    column_title_gp = gpar(fontsize = 12),
    column_title_side = "bottom",
    heatmap_width = unit(width, "inches"),
    heatmap_height = unit(height, "inches"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 10)
    )
  )

  # Save - match legacy
  pdf(output_file, width = width + 2, height = height + 2)
  draw(ht)
  dev.off()

  message(sprintf("Saved comparison heatmap: %s", output_file))
}


#' Plot UMAP with cluster colors
#'
#' @param umap_coords UMAP coordinates matrix
#' @param clusters Cluster assignments
#' @param output_file Output PDF path
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_minimal theme labs element_blank ggsave
#' @export
#' @examples
#' \dontrun{
#' plot_umap(umap_coords, clusters, "umap.pdf")
#' }
plot_umap <- function(umap_coords, clusters, output_file) {

  message("Generating UMAP plot...")

  umap_df <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    Cluster = factor(clusters, levels = sort(unique(clusters)))
  )

  # Use consistent colorblind-friendly palette
  cluster_colors <- create_cluster_colors(clusters)

  p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = cluster_colors) +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "UMAP Projection",
      x = "UMAP 1",
      y = "UMAP 2"
    )

  ggsave(output_file, plot = p, width = 10, height = 8)

  message(sprintf("Saved UMAP plot: %s", output_file))
}




#' Plot spatial distribution of clusters on tissue coordinates
#'
#' @param spatial_data Data frame with x, y, cluster columns
#' @param tiles sf object with tile geometries for context
#' @param output_file Output PDF path
#' @importFrom ggplot2 ggplot aes geom_polygon geom_point scale_fill_manual coord_fixed theme_minimal theme labs ggsave element_blank
#' @importFrom sf st_coordinates
#' @importFrom dplyr %>% filter
#' @export
#' @examples
#' \dontrun{
#' plot_spatial_clusters(spatial_data, tiles, "spatial.pdf")
#' }
plot_spatial_clusters <- function(spatial_data, tiles, output_file) {

  message("Generating spatial cluster plot...")

  # Filter to matched samples
  spatial_data <- spatial_data %>% filter(!is.na(x), !is.na(y))

  message(sprintf("Matched %d samples with %d clusters to spatial coordinates",
                  nrow(spatial_data),
                  length(unique(spatial_data$cluster))))

  # Ensure x and y are numeric (they may have been converted to character)
  spatial_data$x <- as.numeric(spatial_data$x)
  spatial_data$y <- as.numeric(spatial_data$y)

  # Convert cluster to character for simpler handling
  spatial_data$cluster_char <- as.character(spatial_data$cluster)

  # Create explicit color palette for clusters
  # Use consistent colorblind-friendly palette
  cluster_colors <- create_cluster_colors(spatial_data$cluster)

  message(sprintf("Using %d colors for %d clusters: %s",
                  length(cluster_colors), length(cluster_colors),
                  paste(names(cluster_colors), collapse = ", ")))

  # Convert sf tiles to regular polygons for plotting
  tiles_coords <- st_coordinates(tiles)
  tiles_df <- data.frame(
    x = tiles_coords[, "X"],
    y = tiles_coords[, "Y"],
    group = tiles_coords[, "L2"]
  )

  # Create plot with filled squares instead of points
  p <- ggplot(spatial_data, aes(x = x, y = y)) +
    geom_polygon(data = tiles_df,
                 aes(x = x, y = y, group = group),
                 fill = NA, color = "gray80", linewidth = 0.3,
                 inherit.aes = FALSE) +
    geom_point(aes(fill = cluster_char), size = 4, shape = 22, color = "black", stroke = 0.3) +
    scale_fill_manual(values = cluster_colors, name = "Cluster") +
    coord_fixed() +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.position = "right"
    ) +
    labs(title = "Spatial Distribution of Clusters")

  ggsave(output_file, plot = p, width = 10, height = 8)

  message(sprintf("Saved spatial plot: %s", output_file))
}





#' Plot total ploidy distribution
#'
#' @param copy_number_matrix Copy number matrix
#' @param output_file Output PDF path
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline theme_minimal labs theme element_text ggsave
#' @export
#' @examples
#' \dontrun{
#' plot_ploidy_distribution(cn_matrix, "ploidy.pdf")
#' }
plot_ploidy_distribution <- function(copy_number_matrix, output_file) {

  message("Generating ploidy distribution plot...")

  # Calculate total ploidy per sample
  total_ploidy <- rowMeans(copy_number_matrix, na.rm = TRUE)

  ploidy_df <- data.frame(
    sample = rownames(copy_number_matrix),
    total_ploidy = total_ploidy
  )

  p <- ggplot(ploidy_df, aes(x = total_ploidy)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
    geom_vline(xintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
    theme_minimal() +
    labs(
      title = "Total Ploidy Distribution",
      x = "Mean Total Copy Number",
      y = "Number of Samples"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  ggsave(output_file, plot = p, width = 10, height = 6)

  message(sprintf("Saved ploidy plot: %s", output_file))
}

#' Plot LogR variance per sample
#'
#' @param logr_matrix LogR matrix
#' @param output_file Output PDF path
#' @importFrom ggplot2 ggplot aes geom_col theme_minimal labs theme element_blank element_text ggsave
#' @export
#' @examples
#' \dontrun{
#' plot_logr_variance(logr_matrix, "variance.pdf")
#' }
plot_logr_variance <- function(logr_matrix, output_file) {

  message("Generating LogR variance plot...")

  # Calculate variance per sample
  logr_var <- apply(logr_matrix, 1, var, na.rm = TRUE)

  var_df <- data.frame(
    sample = rownames(logr_matrix),
    logr_variance = logr_var
  )

  p <- ggplot(var_df, aes(x = reorder(sample, logr_variance), y = logr_variance)) +
    geom_col(fill = "darkorange", alpha = 0.7) +
    theme_minimal() +
    labs(
      title = "LogR Variance per Sample",
      x = "Sample",
      y = "LogR Variance"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  ggsave(output_file, plot = p, width = 12, height = 6)

  message(sprintf("Saved variance plot: %s", output_file))
}

#' Plot chromosomal instability index
#'
#' @param copy_number_matrix Copy number matrix
#' @param grid Genomic grid
#' @param output_file Output PDF path
#' @importFrom ggplot2 ggplot aes geom_histogram theme_minimal labs theme element_text ggsave
#' @export
#' @examples
#' \dontrun{
#' plot_chromosomal_instability(cn_matrix, grid, "instability.pdf")
#' }
plot_chromosomal_instability <- function(copy_number_matrix, grid, output_file) {

  message("Generating chromosomal instability plot...")

  # Calculate fraction of genome with altered copy number (not 2)
  ci_index <- apply(copy_number_matrix, 1, function(x) {
    sum(x != 2, na.rm = TRUE) / sum(!is.na(x))
  })

  ci_df <- data.frame(
    sample = rownames(copy_number_matrix),
    ci_index = ci_index
  )

  p <- ggplot(ci_df, aes(x = ci_index)) +
    geom_histogram(bins = 30, fill = "purple", color = "white", alpha = 0.7) +
    theme_minimal() +
    labs(
      title = "Chromosomal Instability Index",
      x = "Fraction of Genome Altered",
      y = "Number of Samples"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  ggsave(output_file, plot = p, width = 10, height = 6)

  message(sprintf("Saved instability plot: %s", output_file))
}




#' Compare single-cell vs pseudobulk quality metrics
#'
#' @param sc_profile_dir Single-cell profile directory
#' @param pb_profile_dir Pseudobulk profile directory
#' @param output_file Output PDF path
#' @return Data frame with comparison metrics
#' @importFrom readr read_tsv
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_minimal labs theme ggsave
#' @importFrom gridExtra grid.arrange
#' @importFrom grid textGrob gpar
#' @export
#' @examples
#' \dontrun{
#' metrics <- compare_sc_vs_pseudobulk("sc_profiles/", "pb_profiles/", "comparison.pdf")
#' }
compare_sc_vs_pseudobulk <- function(sc_profile_dir, pb_profile_dir, output_file) {

  message("Comparing single-cell vs pseudobulk quality...")

  # Load single-cell profiles (sample a subset for speed)
  sc_files <- list.files(sc_profile_dir, pattern = "*ASCAT.scprofile.txt",
                         full.names = TRUE)
  sc_sample <- sample(sc_files, min(50, length(sc_files)))

  # Calculate metrics for single-cell
  sc_metrics <- lapply(sc_sample, function(f) {
    prof <- read_tsv(f, show_col_types = FALSE)
    # Handle both LogR and logr column names
    logr_col <- if ("LogR" %in% colnames(prof)) prof$LogR else if ("logr" %in% colnames(prof)) prof$logr else NULL
    if (is.null(logr_col)) return(NULL)

    # Calculate segment lengths if position columns exist
    if ("start" %in% colnames(prof) && "end" %in% colnames(prof)) {
      seg_lengths <- prof$end - prof$start
      mean_seg_length <- mean(seg_lengths, na.rm = TRUE)
    } else {
      mean_seg_length <- NA
    }

    data.frame(
      type = "Single-cell",
      logr_variance = var(logr_col, na.rm = TRUE),
      logr_mad = mad(logr_col, na.rm = TRUE),  # Median absolute deviation
      logr_range = diff(range(logr_col, na.rm = TRUE)),
      n_segments = nrow(prof),
      mean_segment_length = mean_seg_length
    )
  })
  sc_metrics <- bind_rows(sc_metrics[!sapply(sc_metrics, is.null)])

  # Load pseudobulk profiles
  pb_files <- list.files(pb_profile_dir, pattern = "*ASCAT.scprofile.txt",
                         full.names = TRUE)

  if (length(pb_files) == 0) {
    message("No pseudobulk profiles found")
    return(NULL)
  }

  # Calculate metrics for pseudobulk
  pb_metrics <- lapply(pb_files, function(f) {
    prof <- read_tsv(f, show_col_types = FALSE)
    # Handle both LogR and logr column names
    logr_col <- if ("LogR" %in% colnames(prof)) prof$LogR else if ("logr" %in% colnames(prof)) prof$logr else NULL
    if (is.null(logr_col)) return(NULL)

    # Calculate segment lengths if position columns exist
    if ("start" %in% colnames(prof) && "end" %in% colnames(prof)) {
      seg_lengths <- prof$end - prof$start
      mean_seg_length <- mean(seg_lengths, na.rm = TRUE)
    } else {
      mean_seg_length <- NA
    }

    data.frame(
      type = "Pseudobulk",
      logr_variance = var(logr_col, na.rm = TRUE),
      logr_mad = mad(logr_col, na.rm = TRUE),
      logr_range = diff(range(logr_col, na.rm = TRUE)),
      n_segments = nrow(prof),
      mean_segment_length = mean_seg_length
    )
  })
  pb_metrics <- bind_rows(pb_metrics[!sapply(pb_metrics, is.null)])

  # Combine
  all_metrics <- bind_rows(sc_metrics, pb_metrics)

  # Create comparison plots
  p1 <- ggplot(all_metrics, aes(x = type, y = logr_variance, fill = type)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(title = "LogR Variance (noise)", x = "", y = "Variance") +
    theme(legend.position = "none")

  p2 <- ggplot(all_metrics, aes(x = type, y = logr_mad, fill = type)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(title = "LogR MAD (robustness)", x = "", y = "MAD") +
    theme(legend.position = "none")

  p3 <- ggplot(all_metrics, aes(x = type, y = n_segments, fill = type)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(title = "Segment Count (resolution)", x = "", y = "Count") +
    theme(legend.position = "none")

  p4 <- ggplot(all_metrics, aes(x = type, y = mean_segment_length, fill = type)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(title = "Mean Segment Length", x = "", y = "Length (bp)") +
    theme(legend.position = "none")

  # Combine plots
  library(gridExtra)
  combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2,
                                top = textGrob("Pseudobulk vs Single-Cell Quality Metrics\n(Lower variance/MAD = better; More segments = higher resolution)",
                                               gp = gpar(fontsize = 12, font = 2)))

  ggsave(output_file, plot = combined_plot, width = 12, height = 10)

  message(sprintf("Saved quality comparison: %s", output_file))

  return(all_metrics)
}



#' Plot cluster correlation heatmap
#'
#' @param logr_matrix LogR matrix
#' @param clusters Cluster assignments
#' @param output_file Output PDF path
#' @importFrom ComplexHeatmap Heatmap draw
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar grid.text
#' @importFrom grDevices pdf dev.off
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline annotate theme_minimal labs ggsave
#' @export
#' @examples
#' \dontrun{
#' plot_cluster_correlation(logr_matrix, clusters, "correlation.pdf")
#' }
plot_cluster_correlation <- function(logr_matrix, clusters, output_file) {

  message("Generating cluster correlation heatmap...")

  # Calculate cluster centroids
  cluster_ids <- sort(unique(clusters))

  centroids <- do.call(rbind, lapply(cluster_ids, function(cid) {
    cluster_samples <- logr_matrix[clusters == cid, , drop = FALSE]
    colMeans(cluster_samples, na.rm = TRUE)
  }))

  rownames(centroids) <- paste0("Cluster ", cluster_ids)

  # Calculate correlations
  cor_matrix <- cor(t(centroids), use = "pairwise.complete.obs")

  # Create heatmap with correlation values displayed
  ht <- Heatmap(
    cor_matrix,
    name = "Correlation",
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    use_raster = FALSE,  # Don't rasterize - this is a small matrix with text
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    cell_fun = function(j, i, x, y, width, height, fill) {
      # Add correlation value as text
      grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y,
                gp = gpar(fontsize = 7, col = "black"))
    },
    heatmap_legend_param = list(
      title = "Correlation"
    )
  )

  pdf(output_file, width = 8, height = 6)
  draw(ht)
  dev.off()

  message(sprintf("Saved correlation heatmap: %s", output_file))

  # Also create a histogram of correlations (excluding diagonal)
  cor_values <- cor_matrix[upper.tri(cor_matrix)]  # Get upper triangle only

  hist_file <- sub("\\.pdf$", "_histogram.pdf", output_file)

  p_hist <- ggplot(data.frame(correlation = cor_values), aes(x = correlation)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = median(cor_values, na.rm = TRUE),
               linetype = "dashed", color = "red", linewidth = 1) +
    geom_vline(xintercept = 0.9, linetype = "dotted", color = "orange", linewidth = 1) +
    annotate("text", x = median(cor_values, na.rm = TRUE), y = Inf,
             label = sprintf("Median: %.2f", median(cor_values, na.rm = TRUE)),
             hjust = -0.1, vjust = 2, color = "red") +
    annotate("text", x = 0.9, y = Inf,
             label = "Threshold: 0.90",
             hjust = -0.1, vjust = 4, color = "orange") +
    theme_minimal() +
    labs(
      title = "Distribution of Pairwise Cluster Correlations",
      subtitle = sprintf("n = %d cluster pairs", length(cor_values)),
      x = "Pearson Correlation",
      y = "Count"
    )

  ggsave(hist_file, plot = p_hist, width = 8, height = 6)
  message(sprintf("Saved correlation histogram: %s", hist_file))
}





#' Plot global spatial coherence using Moran's I statistics
#'
#' @param morans_results Results from calculate_morans_i
#' @param output_file Output PDF path
#' @importFrom ggplot2 ggplot aes geom_bar geom_point geom_errorbar geom_text scale_fill_manual theme_minimal theme labs ggsave element_blank
#' @export
#' @examples
#' \dontrun{
#' plot_spatial_coherence(morans_results, "coherence.pdf")
#' }
plot_spatial_coherence <- function(morans_results, output_file) {

  message("Generating spatial coherence plot...")

  # Handle new list format (global + local) or old data frame format
  if (is.list(morans_results) && !is.data.frame(morans_results)) {
    morans_df <- morans_results$global
  } else {
    morans_df <- morans_results
  }

  if (is.null(morans_df) || nrow(morans_df) == 0) {
    message("No Moran's I results to plot")
    return(NULL)
  }

  # Create bar plot of Moran's I values with random baseline
  p <- ggplot(morans_df, aes(x = factor(cluster), y = morans_i)) +
    geom_bar(stat = "identity", aes(fill = morans_i > random_mean), alpha = 0.8) +
    geom_point(aes(y = random_mean), color = "red", size = 3, shape = 18) +
    geom_errorbar(aes(ymin = random_mean - random_sd, ymax = random_mean + random_sd),
                  width = 0.3, color = "red", alpha = 0.6) +
    geom_text(aes(label = sprintf("%.3f\n%s", morans_i, interpretation)),
              vjust = ifelse(morans_df$morans_i > 0, -0.3, 1.3),
              size = 3) +
    scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "orange"),
                      labels = c("TRUE" = "Spatially clustered", "FALSE" = "Spatially dispersed"),
                      name = "") +
    theme_minimal() +
    theme(
      legend.position = "top",
      panel.grid.major.x = element_blank()
    ) +
    labs(
      title = "Global Spatial Coherence (Moran's I)",
      subtitle = "Bars = observed; Red diamonds = random baseline (mean ± SD from permutations)",
      x = "Cluster",
      y = "Moran's I"
    )

  ggsave(output_file, plot = p, width = 10, height = 6)

  message(sprintf("Saved spatial coherence plot: %s", output_file))
}
