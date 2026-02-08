#!/usr/bin/env Rscript
# =============================================================================
# analysis/abundance_plotting.R — Cluster Abundance Plots
# =============================================================================
# Refactored from: R_scripts/05_analysis/Plot_Abundance.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))

#' Create abundance barplots from cluster tracking matrix
#'
#' @param tracking_matrix_path Path to cluster_tracking_matrix.csv
#' @param conditions_csv Path to file_name -> condition mapping
#' @param output_dir Output directory for plots
#' @param plot_type "stacked", "grouped", or "proportional"
#' @return ggplot object (invisibly)
plot_abundance <- function(tracking_matrix_path,
                            conditions_csv = NULL,
                            output_dir,
                            plot_type = "stacked") {
  require(data.table)
  require(ggplot2)

  make_outdir(output_dir)
  tracking <- fread(tracking_matrix_path)
  cat("Loaded tracking matrix:", nrow(tracking), "samples x",
      ncol(tracking) - 1, "clusters\n")

  # Pivot to long format
  id_col <- names(tracking)[1]  # First col is file identifier
  cluster_cols <- setdiff(names(tracking), id_col)

  long_df <- tidyr::pivot_longer(tracking, cols = all_of(cluster_cols),
                                  names_to = "cluster", values_to = "count")
  long_df <- as.data.frame(long_df)

  # Add conditions if available
  if (!is.null(conditions_csv) && file.exists(conditions_csv)) {
    conditions <- fread(conditions_csv)
    long_df <- merge(long_df, conditions, by.x = id_col, by.y = "file_name",
                     all.x = TRUE)
  }

  # Create plot
  position <- switch(plot_type,
    "stacked"      = "stack",
    "grouped"      = "dodge",
    "proportional" = "fill",
    "stack"
  )

  p <- ggplot(long_df, aes(x = cluster, y = count, fill = .data[[id_col]])) +
    geom_bar(stat = "identity", position = position) +
    theme_minimal() +
    labs(title = "Cluster Abundance", x = "Cluster", y = "Cell Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(output_dir, paste0("abundance_", plot_type, ".png")),
         p, width = 12, height = 8)

  # Per-condition plot if conditions available
  if ("condition" %in% names(long_df)) {
    p2 <- ggplot(long_df, aes(x = cluster, y = count, fill = condition)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      labs(title = "Abundance by Condition", x = "Cluster", y = "Cell Count") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(output_dir, "abundance_by_condition.png"),
           p2, width = 12, height = 8)
  }

  cat("Plots saved to:", output_dir, "\n")
  invisible(p)
}
