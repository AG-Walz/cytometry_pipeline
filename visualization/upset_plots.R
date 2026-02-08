#!/usr/bin/env Rscript
# =============================================================================
# visualization/upset_plots.R — UpSet Plot Functions
# =============================================================================
# Merged from: R_scripts/05_analysis/Simple_UpSet_Plot.R
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "R_shared")
source(file.path(shared_dir, "io_utils.R"))

#' Create an UpSet plot from binary membership data
#'
#' @param data Data frame with binary (0/1) columns for set membership
#' @param set_columns Character vector of column names representing sets
#' @param output_path Output file path
#' @param title Plot title
#' @param min_size Minimum intersection size to show
plot_upset <- function(data, set_columns, output_path,
                        title = "UpSet Plot", min_size = 1) {
  require(UpSetR)

  png(output_path, width = 1200, height = 800, res = 100)
  tryCatch({
    UpSetR::upset(
      as.data.frame(data),
      sets = set_columns,
      order.by = "freq",
      decreasing = TRUE,
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size",
      text.scale = 1.2,
      main.bar.color = "steelblue",
      sets.bar.color = "grey50",
      mb.ratio = c(0.6, 0.4)
    )
    title(main = title, line = 2.5)
  }, error = function(e) {
    message("UpSet plot error: ", e$message)
  }, finally = dev.off())

  cat("UpSet plot saved to:", output_path, "\n")
}

#' Create UpSet plot from a list of sets
#'
#' @param sets Named list of vectors (set elements)
#' @param output_path Output file path
#' @param title Plot title
plot_upset_from_list <- function(sets, output_path,
                                  title = "UpSet Plot") {
  require(UpSetR)

  png(output_path, width = 1200, height = 800, res = 100)
  tryCatch({
    UpSetR::upset(
      UpSetR::fromList(sets),
      order.by = "freq",
      decreasing = TRUE,
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size",
      text.scale = 1.2
    )
    title(main = title, line = 2.5)
  }, error = function(e) {
    message("UpSet plot error: ", e$message)
  }, finally = dev.off())

  cat("UpSet plot saved to:", output_path, "\n")
}
