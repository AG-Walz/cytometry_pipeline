concatenate_exprs_mtx <- function(flowframe_list, file_names = NULL, col_indices = NULL){
  if (is.null(file_names)) file_names <- seq_along(flowframe_list)
  expression_matrices <- lapply(seq_along(flowframe_list), function(i) {
    expr_data <- exprs(flowframe_list[[i]])
    if (!is.null(col_indices)) expr_data <- expr_data[, col_indices,
                                                      drop = FALSE]
    cbind(expr_data, file_origin = rep(file_names[i], nrow(expr_data)))
  })

  concat_ffr <- do.call(rbind, expression_matrices)
  return(concat_ffr)
}

# Create an FBM on disk and fill it row-wise from a list of flowFrames
concatenate_to_file <- function(ff_list, backingfile,
                                col_indices=NULL, add_origin_col=TRUE,
                                overwrite=FALSE, set_names=FALSE){
  dir.create(dirname(backingfile), recursive=TRUE, showWarnings=FALSE)
  bk <- paste0(backingfile,".bk"); rds <- paste0(backingfile,".rds")
  if ((file.exists(bk)|file.exists(rds)) && !overwrite) stop("Exists; set overwrite=TRUE.")
  if (overwrite){suppressWarnings({unlink(bk); unlink(rds)})}
  first <- flowCore::exprs(ff_list[[1]])
  sel <- if (is.null(col_indices)) seq_len(ncol(first)) else col_indices
  fbm <- bigstatsr::FBM(sum(vapply(ff_list, nrow, 0L)),
                        length(sel)+as.integer(add_origin_col),
                        type="double", backingfile=backingfile, init=NA_real_)
  if (set_names){
    op <- getOption("bigmemory.allow.dimnames"); options(bigmemory.allow.dimnames=TRUE)
    on.exit(options(bigmemory.allow.dimnames=op), add=TRUE)
    try(colnames(fbm$bm()) <- c(colnames(first)[sel], if(add_origin_col)"file_origin"), silent=TRUE)
  }
  cur <- 0L
  for(i in seq_along(ff_list)){
    X <- flowCore::exprs(ff_list[[i]])[, sel, drop=FALSE]; k <- nrow(X); if(k==0L) next
    rng <- (cur+1L):(cur+k); fbm[rng, seq_len(ncol(X))] <- X
    if(add_origin_col) fbm[rng, ncol(fbm)] <- i; cur <- cur + k
  }
  ok <- tryCatch({fbm$save(); TRUE}, error=function(e) FALSE); if(!ok) saveRDS(fbm, rds)
  invisible(fbm)
}

# Enhanced UMAP plotting function for individual marker-sample combinations
umap_marker_individual <- function(umap_coords, expression_data, marker_data, 
                                  sample_name, marker_name, output_path, 
                                  plot_formats = c("png"), plot_width = 8, 
                                  plot_height = 6, dpi = 300) {
  # Create plot data
  plot_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    expression = marker_data
  )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = expression)) +
    geom_point(size = 0.1, alpha = 0.7) +
    scale_color_viridis_c(name = marker_name, option = "viridis") +
    labs(title = paste(marker_name, "Expression"),
         subtitle = paste("Sample:", sample_name),
         x = "UMAP1", y = "UMAP2") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Save in specified formats
  saved_formats <- c()
  for (format in plot_formats) {
    file_path <- paste0(output_path, ".", format)
    if (format == "pdf") {
      ggsave(file_path, p, width = plot_width, height = plot_height, device = "pdf")
    } else if (format == "png") {
      ggsave(file_path, p, width = plot_width, height = plot_height, dpi = dpi, device = "png")
    } else if (format == "svg") {
      ggsave(file_path, p, width = plot_width, height = plot_height, device = "svg")
    }
    saved_formats <- c(saved_formats, format)
  }
  
  return(saved_formats)
}

# Enhanced function for batch processing multiple samples
umap_batch_process <- function(flowframe_list, sample_names, marker_indices, 
                              output_dir, umap_params = list()) {
  # Set default UMAP parameters
  default_params <- list(n_neighbors = 15, min_dist = 0.1, verbose = TRUE)
  umap_params <- modifyList(default_params, umap_params)
  
  # Extract and combine expression data
  all_expression <- map(seq_along(flowframe_list), function(i) {
    expr_data <- flowCore::exprs(flowframe_list[[i]])[, marker_indices, drop = FALSE]
    data.frame(expr_data, sample = rep(sample_names[i], nrow(expr_data)))
  })
  
  combined_expression <- do.call(rbind, all_expression)
  
  # Compute UMAP on combined data (excluding sample column)
  expr_matrix <- as.matrix(combined_expression[, -ncol(combined_expression)])
  umap_result <- uwot::umap(expr_matrix, 
                           n_neighbors = umap_params$n_neighbors,
                           min_dist = umap_params$min_dist,
                           verbose = umap_params$verbose)
  
  # Create results dataframe
  results <- data.frame(
    UMAP1 = umap_result[, 1],
    UMAP2 = umap_result[, 2],
    sample = combined_expression$sample,
    combined_expression[, -ncol(combined_expression)]
  )
  
  return(results)
}

# Function to create compact sample names by removing common parts
create_compact_sample_names <- function(file_paths) {
  # Extract base names without extensions
  base_names <- tools::file_path_sans_ext(basename(file_paths))
  
  if (length(base_names) <= 1) return(base_names)
  
  # Split names into components
  name_components <- strsplit(base_names, "[._-]")
  
  # Find common prefix
  min_length <- min(lengths(name_components))
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

# Legacy functions for backward compatibility (updated to use proper namespacing)
umap_marker <- function(umap, matrix, file_name, df_col, marker_name, plot_path) {
  as.data.frame(umap) %>%
    ggplot2::ggplot(ggplot2::aes(x = UMAP1, y = UMAP2)) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(ggplot2::aes(color = matrix[, df_col]),
                        shape = 20, size = .1) +
    ggplot2::scale_color_gradientn(colors = c("#0000FFFF","#008080","#00FF00FF",
                                             "#80ff00","#FFFF00FF", "#ff8000",
                                             "#FF0000FF")) +
    ggplot2::labs(x = "UMAP1",
                  y = "UMAP2",
                  subtitle = file_name,
                  color = paste0(marker_name, "          ", sep = ""))
  ggplot2::ggsave(paste(plot_path, file_name, ".pdf", sep=""),
                  dpi = 400)
}

umap_vector <- function(umap1, umap2, vector, file_name, palette, plot_path){
  ggplot2::ggplot(data = data.frame(x = umap1, y = umap2),
                  ggplot2::aes(x = umap1, y = umap2)) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(ggplot2::aes(color = c(vector)), shape = 1, size = 0.8,
                        stroke = 0.1) +
    ggplot2::labs(x = "UMAP1",
                  y = "UMAP2",
                  subtitle = file_name,
                  color = file_name)
  ggplot2::ggsave(paste(plot_path, file_name, ".pdf", sep = ""), dpi = 400)
}

tsne_vector <- function(tsne, vector, file_name, palette, plot_path) {
  as.data.frame(tsne) %>%
    ggplot2::ggplot(ggplot2::aes(x = dimred_1, y = dimred_2)) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(ggplot2::aes(color = vector), shape = 1, size = 0.8, stroke = 0.2) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::labs(x = "tSNE1",
                  y = "tSNE2",
                  subtitle = file_name,
                  color = file_name)
  ggplot2::ggsave(paste(plot_path, file_name, ".pdf", sep = ""),
                  width = 7, height = 5, dpi = 800)
}

umap_vec <- function(umap, vector, file_name, palette, plot_path) {
  as.data.frame(umap) %>%
    ggplot2::ggplot(ggplot2::aes(x = UMAP1, y = UMAP2)) +
    ggplot2::theme_bw() +
    ggplot2::geom_point(ggplot2::aes(color = vector), shape = 1, size = 0.8, stroke = 0.1) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::labs(x = "UMAP1",
                  y = "UMAP2",
                  subtitle = file_name,
                  color = file_name)
  ggplot2::ggsave(paste(plot_path, file_name, ".pdf", sep = ""),
                  width = 7, height = 5, dpi = 800)
}

uwot_plot <- function(df, uwot_embeding, file_name, file_location, df_col){
  cat("Computing uwot_embedding for file", file_name)
  uwot_plot <- as.data.frame(uwot_embeding) %>%
    ggplot2::ggplot(ggplot2::aes(x = UMAP1, y = UMAP2)) +
    ggplot2::geom_point(ggplot2::aes(color = df_col), shape = 20, size = .1) +
    ggplot2::scale_color_gradientn(colors = c("#0000FFFF", "#008080", "#00FF00FF",
                                             "#80ff00","#FFFF00FF", "#ff8000",
                                             "#FF0000FF")) +
    ggplot2::labs(x = "UMAP1",
                  y = "UMAP2",
                  subtitle = file_name)
  ggplot2::ggsave(paste(file_location, file_name, ".png", sep = ""),
                  width = 5, height = 3, dpi = 300)
}
