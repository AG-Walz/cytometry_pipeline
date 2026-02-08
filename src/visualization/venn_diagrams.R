#!/usr/bin/env Rscript
# =============================================================================
# visualization/venn_diagrams.R — Venn/Euler Diagram Functions
# =============================================================================
# Extracted from: R_scripts/Functions.R (plot_venn, plot_euler, plot_venneuler)
# =============================================================================

shared_dir <- file.path(dirname(getwd()), "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))

#' Plot Euler diagram from set data
#'
#' @param sets Named list of index vectors
#' @param output_dir Output directory
#' @param filename Base filename (without extension)
#' @param custom_name Optional custom set names
plot_euler <- function(sets, output_dir, filename, custom_name = NULL) {
  require(eulerr)
  n_sets <- length(sets)
  if (n_sets < 2) stop("At least two sets required")
  if (!is.null(custom_name) && length(custom_name) == n_sets) {
    names(sets) <- custom_name
  }
  fit <- eulerr::euler(sets)
  set_colors <- scales::hue_pal()(n_sets)

  png(file.path(output_dir, paste0("euler_", filename, ".png")),
      width = 600, height = 400)
  tryCatch({
    p <- plot(fit, main = paste(filename, "population"),
              fills = set_colors, quantities = TRUE,
              legend = names(sets))
    print(p)
  }, error = function(e) {
    message("Error plotting '", filename, "': ", e$message)
  }, finally = dev.off())
}

#' Plot Venn diagram
#'
#' @param sets Named list of index vectors
#' @param filename Base filename
#' @param output_dir Output directory
plot_venn <- function(sets, filename, output_dir) {
  require(venn)
  p <- venn::venn(
    x = sets, borders = FALSE, ggplot = TRUE,
    zcolor = "style", box = FALSE, plotsize = 10,
    ilcs = 1, sncs = 1, ilabels = "counts"
  ) + ggplot2::ggtitle(paste(filename, "population"))
  ggplot2::ggsave(
    file.path(output_dir, paste0("venn_", filename, ".png")),
    plot = p, width = 8, height = 6
  )
}

#' Plot Venn/Euler for each population in a binary matrix
#'
#' @param lists Binary matrix (rows = cells, cols = samples/workspaces)
#' @param population Population label to filter
#' @param filename Base filename
#' @param output_dir Output directory
#' @param scale_size Use Euler (TRUE) or Venn (FALSE)
#' @param custom_name Optional custom names for sets
plot_and_save_venn <- function(lists, population, filename, output_dir,
                                scale_size = FALSE, custom_name = NULL) {
  sets <- lapply(seq_len(ncol(lists)),
                 function(i) which(lists[, i] == population))
  non_empty <- sapply(sets, length) > 0
  sets <- sets[non_empty]
  if (length(sets) < 2) return(invisible(NULL))

  if (!is.null(custom_name) && length(custom_name) == ncol(lists)) {
    set_names <- custom_name[non_empty]
  } else {
    set_names <- colnames(lists)[non_empty]
  }
  names(sets) <- set_names

  if (scale_size) {
    plot_euler(sets, output_dir, filename, custom_name = set_names)
  } else {
    plot_venn(sets, filename, output_dir)
  }
}
