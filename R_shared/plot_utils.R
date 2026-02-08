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
                              width = 10, height = 8,
                              formats = list(pdf = TRUE, png = TRUE, svg = FALSE),
                              png_dpi = 300) {
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
  ggplot2::ggsave(filename = out_file, plot = p, width = 12, height = 8)
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
    ggplot2::ggsave(
      filename = file.path(out_dir, paste0(out_prefix, "_", sample, ".png")),
      plot = p, width = 12, height = 8
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
