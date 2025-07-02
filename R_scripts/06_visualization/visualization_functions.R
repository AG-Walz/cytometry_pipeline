### ========================= VISUALIZATION FUNCTIONS ============================== ###
# Functions for plotting and data visualization

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