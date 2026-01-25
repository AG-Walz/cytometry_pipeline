### ========================= VISUALIZATION FUNCTIONS ============================== ###
# Functions for plotting and data visualization
#
# TODOs:
#   - [ ] Add support for custom color palettes in UMAP plots
#   - [ ] Implement faceted UMAP plots (by sample or marker)
#   - [ ] Add density contour options for UMAP visualization
#   - [ ] Support for 3D UMAP visualization
#   - [ ] Add publication-ready theme options
# =============================================================================

# ========================= UMAP HELPER FUNCTIONS =========================

#' Concatenate expression matrices from flow frames
#' @param flowframe_list List of flowFrame objects
#' @param file_names Optional vector of file names for origin tracking
#' @param col_indices Optional column indices to select
#' @return Combined expression matrix with file_origin column
concatenate_exprs_mtx <- function(flowframe_list, file_names = NULL, col_indices = NULL){
  if (is.null(file_names)) file_names <- seq_along(flowframe_list)
  expression_matrices <- lapply(seq_along(flowframe_list), function(i) {
    expr_data <- flowCore::exprs(flowframe_list[[i]])
    if (!is.null(col_indices)) expr_data <- expr_data[, col_indices, drop = FALSE]
    cbind(expr_data, file_origin = rep(file_names[i], nrow(expr_data)))
  })
  concat_ffr <- do.call(rbind, expression_matrices)
  return(concat_ffr)
}

#' Create a file-backed matrix (FBM) for large datasets
#' @param ff_list List of flowFrame objects
#' @param backingfile Path for the backing file
#' @param col_indices Optional column indices to select
#' @param add_origin_col Whether to add file origin column
#' @param overwrite Whether to overwrite existing files
#' @param set_names Whether to set column names
#' @return FBM object
concatenate_to_file <- function(ff_list, backingfile,
                                col_indices = NULL, add_origin_col = TRUE,
                                overwrite = FALSE, set_names = FALSE){
  dir.create(dirname(backingfile), recursive = TRUE, showWarnings = FALSE)
  bk <- paste0(backingfile, ".bk")
  rds <- paste0(backingfile, ".rds")
  if ((file.exists(bk) | file.exists(rds)) && !overwrite) {
    stop("Backing file exists; set overwrite=TRUE.")
  }
  if (overwrite) {
    suppressWarnings({unlink(bk); unlink(rds)})
  }
  first <- flowCore::exprs(ff_list[[1]])
  sel <- if (is.null(col_indices)) seq_len(ncol(first)) else col_indices
  fbm <- bigstatsr::FBM(
    sum(vapply(ff_list, nrow, 0L)),
    length(sel) + as.integer(add_origin_col),
    type = "double", backingfile = backingfile, init = NA_real_
  )
  if (set_names) {
    op <- getOption("bigmemory.allow.dimnames")
    options(bigmemory.allow.dimnames = TRUE)
    on.exit(options(bigmemory.allow.dimnames = op), add = TRUE)
    try(colnames(fbm$bm()) <- c(colnames(first)[sel], 
                                if(add_origin_col) "file_origin"), silent = TRUE)
  }
  cur <- 0L
  for (i in seq_along(ff_list)) {
    X <- flowCore::exprs(ff_list[[i]])[, sel, drop = FALSE]
    k <- nrow(X)
    if (k == 0L) next
    rng <- (cur + 1L):(cur + k)
    fbm[rng, seq_len(ncol(X))] <- X
    if (add_origin_col) fbm[rng, ncol(fbm)] <- i
    cur <- cur + k
  }
  ok <- tryCatch({fbm$save(); TRUE}, error = function(e) FALSE)
  if (!ok) saveRDS(fbm, rds)
  invisible(fbm)
}

#' Create compact sample names by removing common parts
#' @param file_paths Vector of file paths or names
#' @return Vector of compact names
create_compact_sample_names <- function(file_paths) {
  base_names <- tools::file_path_sans_ext(basename(file_paths))
  if (length(base_names) <= 1) return(base_names)
  
  name_components <- strsplit(base_names, "[._-]")
  min_length <- min(lengths(name_components))
  
  # Find common prefix
  common_prefix_length <- 0
  for (i in 1:min_length) {
    parts <- sapply(name_components, function(x) x[i])
    if (length(unique(parts)) == 1) {
      common_prefix_length <- i
    } else {
      break
    }
  }
  
  # Find common suffix
  common_suffix_length <- 0
  for (i in 1:min_length) {
    parts <- sapply(name_components, function(x) x[length(x) - i + 1])
    if (length(unique(parts)) == 1) {
      common_suffix_length <- i
    } else {
      break
    }
  }
  
  # Create compact names
  compact_names <- sapply(name_components, function(parts) {
    if (length(parts) <= (common_prefix_length + common_suffix_length)) {
      return(paste(parts, collapse = "_"))
    }
    start_idx <- common_prefix_length + 1
    end_idx <- length(parts) - common_suffix_length
    if (start_idx <= end_idx) {
      return(paste(parts[start_idx:end_idx], collapse = "_"))
    } else {
      return(paste(parts, collapse = "_"))
    }
  })
  
  # Ensure uniqueness
  if (length(unique(compact_names)) < length(compact_names)) {
    compact_names <- make.unique(compact_names, sep = "_")
  }
  return(compact_names)
}

#' Enhanced UMAP plotting for individual marker-sample combinations
#' @param umap_coords UMAP coordinates matrix
#' @param marker_data Vector of marker expression values
#' @param sample_name Name of the sample
#' @param marker_name Name of the marker
#' @param output_path Base output path (without extension)
#' @param plot_formats Formats to save ("png", "pdf", "svg")
#' @param plot_width Width in inches
#' @param plot_height Height in inches
#' @param dpi Resolution for PNG
#' @return Vector of saved format names
umap_marker_individual <- function(umap_coords, marker_data, 
                                   sample_name, marker_name, output_path, 
                                   plot_formats = c("png"), plot_width = 8, 
                                   plot_height = 6, dpi = 300) {
  plot_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    expression = marker_data
  )
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = UMAP1, y = UMAP2, color = expression)) +
    ggplot2::geom_point(size = 0.1, alpha = 0.7) +
    ggplot2::scale_color_viridis_c(name = marker_name, option = "viridis") +
    ggplot2::labs(title = paste(marker_name, "Expression"),
                  subtitle = paste("Sample:", sample_name),
                  x = "UMAP1", y = "UMAP2") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "right")
  
  saved_formats <- c()
  for (format in plot_formats) {
    file_path <- paste0(output_path, ".", format)
    if (format == "pdf") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, device = "pdf")
    } else if (format == "png") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, dpi = dpi, device = "png")
    } else if (format == "svg") {
      ggplot2::ggsave(file_path, p, width = plot_width, height = plot_height, device = "svg")
    }
    saved_formats <- c(saved_formats, format)
  }
  return(saved_formats)
}

#' Batch process UMAP for multiple samples
#' @param flowframe_list List of flowFrame objects
#' @param sample_names Vector of sample names
#' @param marker_indices Column indices for markers
#' @param output_dir Output directory
#' @param umap_params List of UMAP parameters
#' @return Data frame with UMAP results
umap_batch_process <- function(flowframe_list, sample_names, marker_indices, 
                               output_dir, umap_params = list()) {
  default_params <- list(n_neighbors = 15, min_dist = 0.1, verbose = TRUE)
  umap_params <- modifyList(default_params, umap_params)
  
  all_expression <- purrr::map(seq_along(flowframe_list), function(i) {
    expr_data <- flowCore::exprs(flowframe_list[[i]])[, marker_indices, drop = FALSE]
    data.frame(expr_data, sample = rep(sample_names[i], nrow(expr_data)))
  })
  
  combined_expression <- do.call(rbind, all_expression)
  expr_matrix <- as.matrix(combined_expression[, -ncol(combined_expression)])
  
  umap_result <- uwot::umap(expr_matrix, 
                            n_neighbors = umap_params$n_neighbors,
                            min_dist = umap_params$min_dist,
                            verbose = umap_params$verbose)
  
  results <- data.frame(
    UMAP1 = umap_result[, 1],
    UMAP2 = umap_result[, 2],
    sample = combined_expression$sample,
    combined_expression[, -ncol(combined_expression)]
  )
  return(results)
}

#' Create safe filename from marker and sample names
#' @param base_name Base name for the file
#' @param marker_name Marker name
#' @param sample_name Sample name
#' @return Safe filename string
create_safe_filename <- function(base_name, marker_name, sample_name) {
  clean_marker <- gsub("[^A-Za-z0-9._-]", "_", marker_name)
  clean_marker <- gsub("_{2,}", "_", clean_marker)
  clean_marker <- gsub("^_|_$", "", clean_marker)
  
  clean_sample <- gsub("[^A-Za-z0-9._-]", "_", sample_name)
  clean_sample <- gsub("_{2,}", "_", clean_sample)
  clean_sample <- gsub("^_|_$", "", clean_sample)
  
  filename <- paste(base_name, clean_marker, clean_sample, sep = "_")
  return(filename)
}

# Legacy UMAP plotting functions for backward compatibility
umap_marker <- function(umap, matrix, file_name, df_col, marker_name, plot_path) {
  as.data.frame(umap) %>%
    ggplot2::ggplot(ggplot2::aes(x = UMAP1, y = UMAP2)) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(ggplot2::aes(color = matrix[, df_col]), shape = 20, size = .1) +
    ggplot2::scale_color_gradientn(colors = c("#0000FFFF", "#008080", "#00FF00FF",
                                              "#80ff00", "#FFFF00FF", "#ff8000",
                                              "#FF0000FF")) +
    ggplot2::labs(x = "UMAP1", y = "UMAP2", subtitle = file_name,
                  color = paste0(marker_name, "          ", sep = ""))
  ggplot2::ggsave(paste(plot_path, file_name, ".pdf", sep = ""), dpi = 400)
}

umap_vec <- function(umap, vector, file_name, palette, plot_path) {
  as.data.frame(umap) %>%
    ggplot2::ggplot(ggplot2::aes(x = UMAP1, y = UMAP2)) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(ggplot2::aes(color = vector), shape = 1, size = 0.8, stroke = 0.1) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::labs(x = "UMAP1", y = "UMAP2", subtitle = file_name, color = file_name)
  ggplot2::ggsave(paste(plot_path, file_name, ".pdf", sep = ""), width = 7, height = 5, dpi = 800)
}

tsne_vector <- function(tsne, vector, file_name, palette, plot_path) {
  as.data.frame(tsne) %>%
    ggplot2::ggplot(ggplot2::aes(x = dimred_1, y = dimred_2)) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(ggplot2::aes(color = vector), shape = 1, size = 0.8, stroke = 0.2) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::labs(x = "tSNE1", y = "tSNE2", subtitle = file_name, color = file_name)
  ggplot2::ggsave(paste(plot_path, file_name, ".pdf", sep = ""), width = 7, height = 5, dpi = 800)
}

uwot_plot <- function(df, uwot_embedding, file_name, file_location, df_col) {
  cat("Computing uwot_embedding for file", file_name)
  uwot_plot <- as.data.frame(uwot_embedding) %>%
    ggplot2::ggplot(ggplot2::aes(x = UMAP1, y = UMAP2)) +
    ggplot2::geom_point(ggplot2::aes(color = df_col), shape = 20, size = .1) +
    ggplot2::scale_color_gradientn(colors = c("#0000FFFF", "#008080", "#00FF00FF",
                                              "#80ff00", "#FFFF00FF", "#ff8000",
                                              "#FF0000FF")) +
    ggplot2::labs(x = "UMAP1", y = "UMAP2", subtitle = file_name)
  ggplot2::ggsave(paste(file_location, file_name, ".png", sep = ""), width = 5, height = 3, dpi = 300)
}

# ========================= VENN/EULER DIAGRAM FUNCTIONS =========================

# Custom function to format seconds as s, min, h
time_formatter <- function(x) {
  sapply(x, function(val) {
    if (val < 60) {
      sprintf("%.0fs", val)
    } else if (val < 3600) {
      sprintf("%.0fmin", val/60)
    } else {
      sprintf("%.1fh", val/3600)
    }
  })
}

plot_venneuler <- function(sets, output_dir, filename) {
  n_sets <- length(sets)
  if (n_sets < 2) {
    stop("At least two sets are required for an Euler diagram.")
  }
  nms <- names(sets)
  all_indices <- sort(unique(unlist(sets)))
  data <- as.data.frame(sapply(sets, function(idx) all_indices %in% idx))
  colnames(data) <- nms
  fit <- venneuler::venneuler(data)
  
  set_colors <- scales::hue_pal()(n_sets)
  names(set_colors) <- nms
  
  png(filename = file.path(output_dir, paste0("venneuler_", filename, ".png")), 
      width = 800, height = 600)
  tryCatch({
    plot(fit, main = paste(filename, "population"), col = set_colors)
    # Annotate overlaps if possible
    if (is.numeric(fit$weights) && !is.null(fit$centers)) {
      overlap_counts <- round(fit$weights)
      for (i in seq_along(overlap_counts)) {
        region <- fit$centers[i, ]
        text(region[1], region[2], labels = overlap_counts[i], cex = 1.2, font = 2)
      }
    }
    legend("topright", legend = nms, fill = set_colors, border = NA, cex = 1.2, box.lwd = 0)
  }, error = function(e) {
    message(sprintf("Error plotting '%s': %s", filename, e$message))
  }, finally = {
    dev.off()
  })
}

plot_euler <- function(sets, output_dir, filename, custom_name = NULL) {
  n_sets <- length(sets)
  if (n_sets < 2) {
    stop("At least two sets are required for an Euler diagram.")
  }
  # Use custom_name if provided and valid, else fallback to current names
  if (!is.null(custom_name) && length(custom_name) == n_sets && all(!is.na(custom_name)) && all(custom_name != "")) {
    names(sets) <- custom_name
  }
  fit <- eulerr::euler(sets)
  set_colors <- scales::hue_pal()(n_sets)
  names(set_colors) <- custom_name

  png(filename = file.path(output_dir, paste0("euler_", filename, ".png")),
      width = 600, height = 400)
  tryCatch({
    p <- plot(
      fit,
      main = paste(filename, "population"),
      fills = set_colors,
      quantities = TRUE,
      legend = custom_name
    )
    print(p)
  }, error = function(e) {
    message(sprintf("Error plotting '%s': %s", filename, e$message))
  }, finally = {
    dev.off()
  })
}

plot_venn <- function(sets, filename, output_dir){
  p <- venn::venn(
    x = sets,
    borders = FALSE,
    ggplot = TRUE,
    zcolor = "style",
    box = FALSE,
    plotsize = 10,
    ilcs = 1,
    sncs = 1,
    ilabels = "counts"
  ) +
    ggplot2::ggtitle(paste(filename, "population"))
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0("venn_", filename, ".png")),
    plot = p, width = 8, height = 6
  )
}

# Function to plot and save venn for each string
plot_and_save_venn <- function(lists, population, filename, output_dir, 
                               scale_size = FALSE, custom_name = NULL) {
  sets <- lapply(1:ncol(lists), function(i) which(lists[, i] == population))
  # Remove empty sets
  non_empty <- sapply(sets, length) > 0
  sets <- sets[non_empty]
  if (length(sets) < 2) return(invisible(NULL)) # Need at least 2 sets for a Venn/Euler diagram
  
  # Use custom_name if provided and valid, else fallback to colnames
  if (!is.null(custom_name) && length(custom_name) == ncol(lists) &&
      all(!is.na(custom_name)) && all(custom_name != "")) {
    set_names <- custom_name[non_empty]
  } else {
    set_names <- colnames(lists)[non_empty]
  }
  names(sets) <- set_names
  if (any(is.na(names(sets))) || any(names(sets) == "")) {
    stop("All elements of `sets` must be named and non-empty for eulerr.")
  }
  # TODO: Replace scale_size for plot type and add upset plot and venneuler
  if (scale_size) {
    plot_euler(sets = sets, output_dir = output_dir, filename = filename, 
               custom_name = set_names)
  } else {
    plot_venn(sets, output_dir, filename)
  }
}