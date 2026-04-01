# =============================================================================
# R_shared/plot_utils.R — Plotting Utility Functions
# =============================================================================
# Deduplicated from: Functions.R, visualization_functions.R
# =============================================================================

#' Save a ggplot in multiple formats
#' @param plot ggplot object
#' @param filename_base Base filename (without extension)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param formats Named logical list, e.g., list(pdf=TRUE, png=TRUE, svg=FALSE)
#' @param png_dpi DPI for PNG output
#' @return Character vector of formats saved
save_plot_formats <- function(plot, filename_base,
                              width = 14, height = 10,
                              formats = list(pdf = TRUE, png = TRUE, svg = TRUE),
                              png_dpi = 200) {
  format_configs <- list(
    pdf = list(ext = ".pdf", device = "pdf", extra_args = list()),
    png = list(ext = ".png", device = "png", extra_args = list(dpi = png_dpi)),
    svg = list(ext = ".svg", device = "svg", extra_args = list())
  )
  formats_saved <- c()
  for (format_name in names(formats)) {
    if (isTRUE(formats[[format_name]])) {
      config <- format_configs[[format_name]]
      args <- c(
        list(filename = paste0(filename_base, config$ext),
             plot = plot, width = width, height = height,
             device = config$device),
        config$extra_args
      )
      do.call(ggplot2::ggsave, args)
      formats_saved <- c(formats_saved, toupper(format_name))
    }
  }
  formats_saved
}

#' Save a base R plot in PNG + SVG formats
#' @param filename_base Base filename (without extension)
#' @param plot_fn Function that draws the plot (called twice: once for PNG, once for SVG)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param png_dpi DPI for PNG output
save_base_plot_formats <- function(filename_base, plot_fn,
                                   width = 14, height = 10,
                                   png_dpi = 200) {
  # PNG
  grDevices::png(paste0(filename_base, ".png"),
                 width = width, height = height, units = "in", res = png_dpi)
  plot_fn()
  grDevices::dev.off()
  # SVG
  grDevices::svg(paste0(filename_base, ".svg"), width = width, height = height)
  plot_fn()
  grDevices::dev.off()
  cat("Saved:", paste0(filename_base, ".{png,svg}"), "\n")
}

#' Create a grouped barplot of population counts
#' @param population_counts Data frame with Population, Count, Sample columns
#' @param out_file Output file path
#' @param title Plot title
#' @param y_label Y-axis label
population_barplot <- function(population_counts, out_file, title, y_label) {
  population_df <- as.data.frame(population_counts)
  p <- ggplot2::ggplot(population_df,
                       ggplot2::aes(x = Population, y = Count, fill = Sample)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Population", y = y_label) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  # Strip extension from out_file to get base name for save_plot_formats
  filename_base <- sub("\\.[^.]+$", "", out_file)
  save_plot_formats(p, filename_base, width = 16, height = 10)
  cat("Barplot saved to:", out_file, "\n")
}

#' Plot barplots per sample (faceted by sample)
#' @param df Data frame with Population, Sample, Workspace columns
#' @param value_col Name of the value column
#' @param out_dir Output directory
#' @param out_prefix Filename prefix
#' @param y_label Y-axis label
#' @param title Plot title
plot_per_sample <- function(df, value_col, out_dir, out_prefix, y_label, title) {
  unique_samples <- unique(df$Sample)
  for (sample in unique_samples) {
    df_sample <- df[df$Sample == sample, ]
    p <- ggplot2::ggplot(df_sample,
                         ggplot2::aes(x = Population,
                                      y = !!rlang::sym(value_col),
                                      fill = Workspace)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = paste(title, sample), x = "Population", y = y_label) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    save_plot_formats(p,
      file.path(out_dir, paste0(out_prefix, "_", sample)),
      width = 16, height = 10
    )
  }
}

#' Create a filesystem-safe filename
#' @param name Character string to sanitize
#' @param max_length Maximum filename length
#' @return Sanitized filename string
create_safe_filename <- function(name, max_length = 100) {
  safe <- gsub("[^A-Za-z0-9_.-]", "_", name)
  safe <- gsub("_+", "_", safe)
  safe <- gsub("^_|_$", "", safe)
  if (nchar(safe) > max_length) safe <- substr(safe, 1, max_length)
  safe
}

#' Create compact display names from long strings
#' @param names Character vector of names
#' @param max_length Maximum display length
#' @return Character vector of compact names
create_compact_names <- function(names, max_length = 25) {
  sapply(names, function(n) {
    if (nchar(n) <= max_length) return(n)
    paste0(substr(n, 1, max_length - 3), "...")
  })
}
