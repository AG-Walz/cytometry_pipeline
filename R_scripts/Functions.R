#=============================== Functions.R ==============================

# This file contains a collection of helper functions for analyzing
# flow cytometry data in R.
# The functions included in this file are designed to assist with various tasks such as:
# - Checking if a path is a directory and creating directories if they do not exist.
# - Extracting information from flowFrame objects, such as dye and marker names.
# - Comparing panels of flow cytometry data to group flow frames with identical marker sets.
# - Reading and parsing the text segment of FCS (Flow Cytometry Standard) files.
# - Retrieving specified parameters from FCS files.
# - Checking the total number of events in a flow frame against the header information.
#
# VERSION NOTES:
# This version (R_scripts/Functions.R) contains the full function implementations.
# The archive version (archive/R_scripts_20251121/Functions.R) is a modular loader
# that sources functions from thematic subfolders. The actual function files have
# been copied to the appropriate subfolders in R_scripts/.
#
# TODO (from archive modular structure):
# - [ ] Consider implementing load_pipeline_functions() for automatic sourcing
#       of all function files from subfolders (see archive version for reference)
# - [ ] Add display_pipeline_overview() function for user guidance
# - [ ] The archive version organizes functions into:
#       - 01_data_loading/data_loading_functions.R
#       - 02_preprocessing/preprocessing_functions.R (exists as Preprocessing_functions.R)
#       - 03_normalization/homogenization_functions.R (copied)
#       - 05_analysis/analysis_functions.R (copied)
#       - 06_visualization/visualization_functions.R (copied)
#       - 04_gating/gating_functions.R (copied)
#       - 00_utilities/file_management.R, flowframe_functions.R, general_utilities.R (copied)
# ============================================================================

setwd_cwd <- function(fallback_dir) {
  #' Set Working Directory to Script Location
  #' This function sets the working directory to the location of the script being executed.
  #'
  #' @param fallback_dir The directory to use if the script location cannot be determined.
  #' 
  if (Sys.getenv("RSTUDIO")=="1" &&
      requireNamespace("rstudioapi", quietly = TRUE)) {
    # From RStudio
    try({
      script_path <- rstudioapi::getSourceEditorContext()$path
    }, silent = TRUE)
  } else if(is.na(script_path) || script_path == "") {
    # From Command line
    ca <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", ca[grep("^--file=", ca)])
    if (length(file_arg) > 0) script_path <- file_arg
  }else if ((is.na(script_path) || script_path == "") &&
              !is.null(sys.frames()[[1]]$ofile)) {
    # Via source()
    script_path <- sys.frames()[[1]]$ofile
  }else if ((is.na(script_path) || script_path == "") &&
              requireNamespace("knitr", quietly = TRUE)) {
    # Via RMarkdown
    ci <- try(knitr::current_input(), silent = TRUE)
    if (!inherits(ci, "try-error") && !is.null(ci)) script_path <- ci
  }else{
    # Otherwise default to fallback
    script_dir <-  fallback_dir
  }
  setwd(script_dir)
}

# Helper function to save plots in selected formats
save_plot_formats <- function(plot, filename_base,
                              width = PLOT_WIDTH, height = PLOT_HEIGHT) {
  formats_saved <- c()
  
  format_configs <- list(
    pdf = list(ext = ".pdf", device = "pdf", extra_args = list()),
    png = list(ext = ".png", device = "png", extra_args = list(dpi = PNG_DPI)),
    svg = list(ext = ".svg", device = "svg", extra_args = list())
  )
  
  for (format_name in names(OUTPUT_FORMATS)) {
    if (OUTPUT_FORMATS[[format_name]]) {
      config <- format_configs[[format_name]]
      args <- c(
        list(
          filename = paste0(filename_base, config$ext),
          plot = plot,
          width = width,
          height = height,
          device = config$device
        ),
        config$extra_args
      )
      do.call(ggplot2::ggsave, args)
      formats_saved <- c(formats_saved, toupper(format_name))
    }
  }
  return(formats_saved)
}

# Function to read CSV files and remove matching columns
remove_columns_from_csvs <- function(folder_path, colname_pattern, output_folder = NULL) {
  require(data.table)
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No CSV files found in the specified folder")
  
  output_folder <- output_folder %||% file.path(folder_path, "filtered")
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  
  results <- map(csv_files, ~ {
    cat("Processing:", basename(.x), "...")
    data <- fread(.x, sep = ";", header = TRUE)
    matching_cols <- grep(colname_pattern, names(data), value = TRUE)
    
    if (length(matching_cols) > 0) {
      data_filtered <- data[, !names(data) %like% colname_pattern, with = FALSE]
      cat(" removed", length(matching_cols), "columns")
    } else {
      data_filtered <- data
      cat(" no matches found")
    }
    
    output_file <- file.path(output_folder, basename(.x))
    fwrite(data_filtered, output_file, sep = ";")
    cat(" ->", ncol(data_filtered), "cols saved\n")
    
    list(original_cols = ncol(data), filtered_cols = ncol(data_filtered), 
         removed_cols = matching_cols)
  }) %>% setNames(basename(csv_files))
  
  cat("Summary: Processed", length(csv_files), "files -", output_folder, "\n")
  invisible(results)
}


# ========================== GUI section ========================================

# Function for advanced multi-selection (compact)
tk_multiselect <- function(choices, title = "Select items") {
  if (!require("tcltk", quietly = TRUE)) stop("tcltk package is required")
  tt <- tktoplevel(); tkwm.title(tt, title); tkwm.geometry(tt, "400x300")
  main_frame <- tkframe(tt); tkpack(main_frame, expand = TRUE, fill = "both")
  listbox <- tklistbox(main_frame, selectmode = "extended", height = 10, width = 30)
  scrollbar <- tkscrollbar(main_frame, orient = "vertical")
  tkconfigure(listbox, yscrollcommand = function(...) tkset(scrollbar, ...))
  tkconfigure(scrollbar, command = function(...) tkyview(listbox, ...))
  tkpack(scrollbar, side = "right", fill = "y"); tkpack(listbox, side = "left", expand = TRUE, fill = "both")
  for (item in choices) tkinsert(listbox, "end", item)
  selected_items <- NULL
  ok_function <- function() {
    idx <- as.integer(tkcurselection(listbox)) + 1
    selected_items <<- if (length(idx) > 0) choices[idx] else character(0)
    tkdestroy(tt)
  }
  cancel_function <- function() { selected_items <<- NULL; tkdestroy(tt) }
  ok_btn <- tkbutton(main_frame, text = "OK", command = ok_function)
  cancel_btn <- tkbutton(main_frame, text = "Cancel", command = cancel_function)
  tkpack(ok_btn, side = "left", padx = 5, pady = 5)
  tkpack(cancel_btn, side = "left", padx = 5, pady = 5)
  tkbind(listbox, "<Double-Button-1>", function() { ok_function() })
  tkbind(tt, "<Return>", ok_function)
  tkbind(tt, "<Escape>", cancel_function)
  tkfocus(listbox); tkgrab.set(tt); tkwait.window(tt)
  selected_items
}

# Helper for marker checkbox creation
.create_marker_checkboxes <- function(parent, marker_info) {
  vars <- vector("list", length(marker_info$indices))
  for (i in seq_along(marker_info$indices)) {
    vars[[i]] <- tclVar("1")
    txt <- paste0("Column ", marker_info$indices[i], ": ", marker_info$names[i])
    cb <- tkcheckbutton(parent, text = txt, variable = vars[[i]], anchor = "w")
    tkpack(cb, side = "top", fill = "x")
  }
  vars
}

# Helper for marker checkbox creation with scrollbar
.create_marker_checkboxes_scroll <- function(parent, marker_info) {
  canvas <- tkcanvas(parent, bg = "white")
  scrollbar <- tkscrollbar(parent, orient = "vertical", command = function(...) tkyview(canvas, ...))
  tkconfigure(canvas, yscrollcommand = function(...) tkset(scrollbar, ...))
  tkpack(scrollbar, side = "right", fill = "y")
  tkpack(canvas, side = "left", expand = TRUE, fill = "both")
  frame <- tkframe(canvas, bg = "white")
  window_id <- tkcreate(canvas, "window", 0, 0, anchor = "nw", window = frame)
  vars <- vector("list", length(marker_info$indices))
  for (i in seq_along(marker_info$indices)) {
    vars[[i]] <- tclVar("1")
    txt <- paste0("Column ", marker_info$indices[i], ": ", marker_info$names[i])
    cb <- tkcheckbutton(frame, text = txt, variable = vars[[i]], anchor = "w", bg = "white")
    tkpack(cb, side = "top", fill = "x")
  }
  tkbind(frame, "<Configure>", function() {
    bbox <- tkbbox(canvas, "all"); tkconfigure(canvas, scrollregion = bbox)
  })
  vars
}

# Marker selection popup
select_markers <- function(marker_info, title = "Select Markers") {
  tryCatch({
    if (!require("tcltk", quietly = TRUE)) stop("tcltk package is required")
    tt <- tktoplevel(); tkwm.title(tt, title); tkwm.geometry(tt, "400x350")
    main_frame <- tkframe(tt); tkpack(main_frame, expand = TRUE, fill = "both", padx = 10, pady = 10)
    label <- tklabel(main_frame, text = title, font = tkfont.create(size = 12, weight = "bold"))
    tkpack(label, side = "top", pady = 5)
    cb_frame <- tkframe(main_frame); tkpack(cb_frame, expand = TRUE, fill = "both")
    vars <- .create_marker_checkboxes_scroll(cb_frame, marker_info)
    selected_markers <- NULL
    ok_function <- function() {
      selected_markers <<- marker_info$indices[which(sapply(vars, tclvalue) == "1")]
      tkdestroy(tt)
    }
    cancel_function <- function() { selected_markers <<- NULL; tkdestroy(tt) }
    ok_btn <- tkbutton(main_frame, text = "OK", command = ok_function)
    cancel_btn <- tkbutton(main_frame, text = "Cancel", command = cancel_function)
    tkpack(ok_btn, side = "left", padx = 5, pady = 5)
    tkpack(cancel_btn, side = "left", padx = 5, pady = 5)
    tkbind(tt, "<Return>", ok_function)
    tkbind(tt, "<Escape>", cancel_function)
    tkfocus(tt); tkgrab.set(tt); tkwait.window(tt)
    selected_markers
  }, error = function(e) {
    warning(paste("Error creating marker selection popup:", e$message)); NULL
  })
}

# Simplified marker selection function
select_markers_by_mode <- function(mode = c("all", "popup", "terminal"),
                                   filtered_markers, filtered_indices) {
  mode <- match.arg(mode)
  switch(mode,
    "all" = filtered_indices,
    "popup" = {
      cat("\nFiltered markers available:\n")
      iwalk(filtered_markers, ~ cat(sprintf("%d. Column %d: %s\n", .y, filtered_indices[.y], .x)))
      tryCatch({
        selected <- select_markers(list(indices = filtered_indices, names = filtered_markers))
        if (is.null(selected) || length(selected) == 0) filtered_indices else selected
      }, error = function(e) { cat("Popup failed, using all markers\n"); filtered_indices })
    },
    "terminal" = {
      cat("\nAvailable markers:\n")
      iwalk(filtered_markers, ~ cat(sprintf("%d. Column %d: %s\n", .y, filtered_indices[.y], .x)))
      cat("Enter marker numbers (comma-separated, or 'all'): ")
      input <- readline()
      if (input == "all") filtered_indices else {
        nums <- as.numeric(strsplit(input, ",")[[1]])
        filtered_indices[nums[nums %in% seq_along(filtered_indices)]]
      }
    },
    # Default: use as-is if numeric vector, or fallback to all
    if (is.numeric(mode)) mode else filtered_indices
  )
}

# Load list of FCS files into flowframes
load_fcs_list <- function(fcs_file_list) {
  fcs_file_list %>%
    map2(seq_along(.), ~ {
      cat("Loading", basename(.x), "...")
      ff <- flowWorkspace::load_cytoframe_from_fcs(.x, transformation = FALSE,
                                                   truncate_max_range = FALSE)
      cat(" Sample", .y, "has", nrow(ff), "cells and", ncol(ff), "markers\n")
      return(ff)
  })
}

load_fcs_folder <- function(input_file_path) {
  # Load FCS files and extract sample names
  cat("Reading FCS-files from", input_file_path, "\n")
  input_fcs_files <- list.files(path = input_file_path,
                                pattern = "\\.fcs$", full.names = TRUE)
  sample_names <- input_fcs_files %>%
    basename() %>%
    gsub("\\.fcs$", "", .) %>%
    sapply(function(x) strsplit(x, "_")[[1]][1:2] %>% paste(collapse = "_"))

  cat("Found", length(input_fcs_files), "FCS files in", input_file_path, "\n")
  cat("Sample names:", paste(sample_names, collapse = ", "), "\n")

  # Load flow frames into a list
  flow_frames <- load_fcs_list(input_fcs_files)

  return(flow_frames)
}

# =================== Load FlowJo wsp functions ================================

population_barplot <- function(population_counts, out_file, title, y_label) {
  # Convert the population statistics to a data frame
  population_df <- as.data.frame(population_counts)
  
  # Plot the barplot
  population_plot <- ggplot(population_df, aes(x = Population,
                                               y = Count, fill = Sample)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(title = title, x = "Population", y = y_label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the plot to a file
  ggsave(filename = out_file, plot = population_plot, width = 12, height = 8)
  
  # Print completion message
  cat("Barplot saved to:", out_file, "\n")
}

get_population_stats <- function(gating_set, type = NULL) {
  sample_names <- sampleNames(gating_set)
  stats_list <- list()
  for (sample in sample_names) {
    if (is.null(type)) {
      stats <- tibble(gh_pop_get_stats(gating_set[[sample]]))
      stats_list[[sample]] <- data.frame(
        Population = stats$pop,
        Sample = sample,
        Value = stats$count
      )
    } else {
      stats <- tibble(gh_pop_get_stats(gating_set[[sample]], type = type))
      stats_list[[sample]] <- data.frame(
        Population = stats$pop,
        Sample = sample,
        Value = stats$percent
      )
    }
  }
  do.call(rbind, stats_list)
}

get_population_strings_list <- function(gating_set, sample_names) {
  population_strings_list <- list()
  for (i in seq_along(sample_names)) {
    gating_hierarchy <- gating_set[[i]]
    population_names <- gs_get_pop_paths(gating_hierarchy, path = "full")
    cytoframe <- gh_pop_get_data(gating_hierarchy)
    num_events <- nrow(exprs(cytoframe))
    population_strings <- rep(NA, num_events)
    for (pop_name in population_names) {
      indices <- gh_pop_get_indices(gating_hierarchy, pop_name)
      population_strings[indices] <- pop_name
    }
    population_strings_list[[sample_names[i]]] <- population_strings
  }
  population_strings_list
}

summarize_remainder_counts <- function(population_strings_list) {
  population_data <- do.call(rbind, lapply(names(population_strings_list),
                                           function(sample_name) {
                                             data.frame(
                                               Population = population_strings_list[[sample_name]],
                                               Sample = sample_name
                                             )
                                           }))
  population_data %>%
    group_by(Sample, Population) %>%
    summarise(Count = n(), .groups = "drop")
}

save_population_strings <- function(population_strings_list, out_dir) {
  saveRDS(population_strings_list,
          file = paste0(out_dir, "population_strings_list.rds"))
  cat(paste0("Population strings list saved to '",
             out_dir, "population_strings_list.rds'\n"))
}

plot_per_sample <- function(df, value_col, out_prefix, y_label, title) {
  unique_samples <- unique(df$Sample)
  for (sample in unique_samples) {
    df_sample <- df[df$Sample == sample, ]
    p <- ggplot(df_sample, aes(x = Population,
                               y = !!sym(value_col), fill = Workspace)) +
      geom_bar(stat = "identity", position = "dodge") +
      theme_minimal() +
      labs(title = paste(title, sample), x = "Population", y = y_label) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(filename = paste0(out_dir, out_prefix, "_", sample, ".png"),
           plot = p, width = 12, height = 8)
  }
}

# TODO: Function to run loading at once centralized

#=================== Cytolution gating ========================================

flatten_label_list <- function(label_lists, gating_tree) {
  sapply(label_lists, function(labels) {
    if (length(labels) == 0) {
      return("root")
    } else if (length(labels) == 1) {
      return(labels[[1]])
    } else {
      depths <- gating_tree[labels]
      depths[is.na(depths)] <- -Inf
      return(labels[which.max(depths)])
    }
  }, USE.NAMES = FALSE)  # Ensure we get a simple vector
}

detect_csv_sep <- function(file) {
  first_line <- readLines(file, n = 1)
  if (grepl(";", first_line)) {
    return(";")
  } else {
    return(",")
  }
}

labels_from_binary_df <- function(df) {
  apply(df, 1, function(row) {
    names(row)[which(row == 1)]
  })
}
# =========================== Plotting Functions ===============================
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

#============================= Helper Functions ================================

ensure_and_load_packages <- function(req_packages) {
#' Ensure and Load Required Packages
#'
#' This function checks if the required packages are installed and loads them.
#' If any of the required packages are not installed, it installs them first.
#'
#' @param req_packages A character vector of package names that are required.
#'
#' @return None. The function is used for its side effects of ensuring that
#' the required packages are installed and loaded.
#'
  for (pkg in req_packages) {
    # Check if the package is not already installed
    if (!requireNamespace(pkg, quietly = TRUE)) {
      # Install the package
      install.packages(pkg)
    }
    # Load the package
    library(pkg, character.only = TRUE)
  }
}

is_directory <- function(path) {
  #' Check if a Path is a Directory
#'
#' This function checks whether a given path exists and if it is a directory.
#'
#' @param path A character string representing the path to check.
#' @return A logical value: \code{TRUE} if the path exists and is a directory, \code{FALSE} otherwise.
#' @examples
#' is_directory("/path/to/directory")
#' is_directory("/path/to/file.txt")
#' is_directory("/non/existent/path")
#' 
  if (!file.exists(path)) {return(FALSE)}
    info <- file.info(path)
    return(info$isdir)
}


dir_create <- function(dir_path_list){
  #' Create Directories from a List of Paths
#'
#' This function takes a list of directory paths and creates the directories if they do not already exist.
#' It ensures that each path ends with a forward slash before attempting to create the directory.
#'
#' @param dir_path_list A list of directory paths to be created.
#'
#' @details
#' The function iterates over each path in the provided list. If a path does not end with a forward slash,
#' it appends one. Then, it checks if the directory exists using the `is_directory` function. If the directory
#' does not exist, it creates the directory using `dir.create`.
#' 
  for(path in dir_path_list){
    if (!endsWith(path, "/")){
      path <- paste0(path, "/")
      }
    if (!is_directory(path)){
      suppressWarnings({
        # Create directory if it does not exist
        dir.create(path, showWarnings = FALSE) 
      })
    }
  }
} # lapply(paths, dir.create, recursive = TRUE)

#======================= Compare_gating Functions ==============================

translate_labels <- function(string_list, dictionary) {
  # dictionary: named vector, names are original, values are standardized
  idx <- match(string_list, names(dictionary))
  # If found, use dictionary value, else keep original
  out <- ifelse(!is.na(idx), dictionary[idx], string_list)
  return(out)
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
#==================== Extract information from FF ==============================

explore_data <- function(flowframe){
#' Explore Flow Cytometry Data
#'
#' This function takes a flowFrame object and extracts information about the 
#' dyes and markers used in the flow cytometry experiment. It returns a list 
#' containing the names of the dyes, the names of the markers, and the positions 
#' of the features with marker names.
#'
#' @param flowframe A flowFrame object containing flow cytometry data.
#' @return A list with three elements:
#' \describe{
#'   \item{dyes}{A character vector of dye names.}
#'   \item{markers}{A character vector of marker names.}
#'   \item{mcols}{An integer vector indicating the positions of features with marker names.}
#' }
#' @importFrom flowCore exprs featureNames markernames
#' 
  marker_cols <- c()
  dye_names <- c()
  marker_names <- c()
  for (i in 1:ncol(flowCore::exprs(flowframe))){
    dye_names[i] <- toString(flowCore::featureNames(flowframe[,i]))
    marker_names[i] <- toString(flowCore::markernames(flowframe[,i]))
  }
  # List position of features with marker names
  for (j in 1:length(marker_names)){
    if(marker_names[j]!=""){marker_cols <- c(marker_cols,j)}
  }
  return_list <- list("dyes" = dye_names,"markers" = marker_names,
                      "mcols" = marker_cols)
  return(return_list)
}


compare_panels <- function(flowframes_list) {
  # Check whether 2 markernames are equal or not
#  all.equal(markernames(flow_frames[[1]]),markernames(flow_frames[[length(flow_frames)]]))

#' Compare Panels Function
#'
#' This function compares the markers in a list of flow frames and groups them based on identical marker sets.
#'
#' @param flowframes_list A list of flow frames to be compared.
#' @return A list where each element contains the indices of flow frames with identical marker sets.
#' 
  marker_dict <- list()
  marker_dict[[1]] <- explore_data(flowframes_list[[1]])$markers
  cat("\nComparing panels. \n")
  for (i in 2:length(flowframes_list)) {
    current_markers <- explore_data(flowframes_list[[i]])$markers
    marker_found <- FALSE
    
    for (j in seq_along(marker_dict)) {
      if (identical(marker_dict[[j]], current_markers)) {
        marker_dict[[j]] <- c(marker_dict[[j]], i)
        marker_found <- TRUE
        break
      }
    }
    
    if (!marker_found) {
      marker_dict[[length(marker_dict) + 1]] <- current_markers
    }
  }
  
  return(marker_dict)
}



readFCStext <- function(con, offsets, emptyValue = TRUE, cpp = TRUE, ...) {
  #' Read FCS Text Segment
#'
#' This function reads the text segment of an FCS (Flow Cytometry Standard) file.
#'
#' @param con A connection object to the FCS file.
#' @param offsets A named vector containing the start and end positions of the text segment and the FCS version.
#' @param emptyValue A logical value indicating whether to include empty values in the parsed output. Default is TRUE.
#' @param cpp A logical value indicating whether to use the C++ implementation for parsing. Default is TRUE.
#' @param ... Additional arguments passed to the parsing functions.
#'
#' @return A named vector containing the parsed text segment of the FCS file.
#'
#' @details
#' The function reads the text segment of an FCS file by seeking to the specified start position and reading the raw bytes until the end position.
#' It then converts the raw bytes to a character string using the specified encoding and trims any trailing whitespace.
#' Depending on the FCS version, it either splits the text segment using a delimiter or uses a parsing function to extract key-value pairs.
#' The parsed output is returned as a named vector with the FCS version as the first element.
#'
  seek(con, offsets["textstart"])
  txt <- iconv(rawToChar(readBin(con, "raw", offsets["textend"] - offsets["textstart"] + 1)), "", "latin1", sub = "byte")
  txt <- trimws(txt, "right")
  
  if (offsets["FCSversion"] <= 2) {
    delimiter <- substr(txt, 1, 1)
    sp <- strsplit(substr(txt, 2, nchar(txt)), split = delimiter, fixed = TRUE)[[1]]
    rv <- c(offsets["FCSversion"], sp[seq(2, length(sp), by = 2)])
    names(rv) <- gsub("^ *| *$", "", c("FCSversion", sp[seq(1, length(sp) - 1, by = 2)]))
  } else {
    rv <- if (cpp) fcsTextParse(txt, emptyValue = emptyValue) else fcs_text_parse(txt, emptyValue = emptyValue)
    rv <- c(offsets["FCSversion"], rv)
    names(rv) <- gsub("^ *| *$", "", names(rv))
  }
  
  return(rv)
}


readFCSgetPar <- function(x, pnam, strict = TRUE) {
  #' Read FCS Parameters
#'
#' This function retrieves specified parameters from an FCS file.
#'
#' @param x A character vector representing the names of the parameters in the FCS file.
#' @param pnam A character vector of parameter names to be retrieved.
#' @param strict A logical value indicating whether to enforce strict matching of parameter names.
#'   If TRUE (default), the function will stop with an error if any parameter names are not found.
#'   If FALSE, the function will return NA for any parameter names not found.
#'
#' @return A character vector of the retrieved parameters. If `strict` is FALSE, the returned vector
#'   will contain NA for any parameter names not found.
#'
  stopifnot(is.character(x), is.character(pnam))
  i <- match(pnam, names(x))
  if (any(is.na(i)) && strict) stop(paste("Parameter(s)", pnam[is.na(i)], "not contained in 'x'\n"))
  if (!strict) {
    i[!is.na(i)] <- x[i[!is.na(i)]]
    names(i) <- pnam
    return(i)
  }
  return(x[i])
}


check_tot <- function(flowframe){
  #' Check Total Events in Flow Frame
#'
#' This function compares the total number of events recorded in the header 
#' of a flow frame with the actual number of events in the flow frame.
#'
#' @param flowframe A flow frame object containing cytometry data.
#'
#' @return Prints a message indicating whether the total number of events 
#' in the header matches the actual number of events in the flow frame.
#'
#' @note The function assumes that the header of the flow frame contains 
#' a keyword named "$TOT" which indicates the total number of events.
#'
  header <- str(flowframe)
  tot_keyword <- header[["  .. ..$ $TOT"]]
  n_events <- nrow(flowframe)
  if (tot_keyword == n_events) {
    print(paste("The value of the '$TOT' keyword and the actual number of",
                "events match for flow frame",i))
  }else {
    print(paste("The value of the '$TOT' keyword and the actual number of",
                "events do not match for flow frame",i))
  }
}
# ==================  Load csv gating functions=============================
# Memory-efficient CSV to FCS mapping
create_csv_fcs_mapping <- function(csv_files, fcs_files) {
  csv_to_fcs <- character(length(csv_files))
  names(csv_to_fcs) <- csv_files
  
  fcs_bases_norm <- gsub("[-_]", "", gsub("\\.fcs$", "", basename(fcs_files)))
  
  for (i in seq_along(csv_files)) {
    csv_file <- csv_files[i]
    csv_base_norm <- gsub("[-_]", "", gsub("\\.csv$", "", basename(csv_file)))
    
    match_idx <- which(sapply(fcs_bases_norm, function(fcs_norm) {
      grepl(fcs_norm, csv_base_norm, fixed = TRUE)
    }))
    
    if (length(match_idx) != 1) {
      stop(sprintf("Could not uniquely match FCS file for CSV '%s' (found %d matches)",
                   basename(csv_file), length(match_idx)))
    }
    csv_to_fcs[csv_file] <- fcs_files[match_idx]
  }
  
  return(csv_to_fcs)
}

# Get FCS event count without loading full file
get_fcs_event_count <- function(fcs_file) {
  tryCatch({
    ff <- flowCore::read.FCS(fcs_file, truncate_max_range = FALSE, 
                            transformation = FALSE, which.lines = 1)
    as.numeric(flowCore::keyword(ff, "$TOT")[[1]])
  }, error = function(e) {
    flowcore_ff <- flowCore::read.FCS(fcs_file, truncate_max_range = FALSE)
    nrow(flowcore_ff)
  })
}

# Check if CSV has manual gating columns
has_manual_gating <- function(csv_file, pattern = "Manual gating") {
  tryCatch({
    separator <- detect_csv_sep(csv_file)
    header <- fread(csv_file, sep = separator, nrows = 0)
    any(grepl(pattern, names(header)))
  }, error = function(e) FALSE)
}

# Filter CSV files to only those with manual gating
filter_valid_csv_files <- function(csv_files, pattern = "Manual gating") {
  valid_files <- csv_files[sapply(csv_files, has_manual_gating, pattern = pattern)]
  cat("Found", length(valid_files), "valid CSV files with pattern:", pattern, "\n")
  valid_files
}

# Read and filter CSV to essential columns only
read_csv_filtered <- function(csv_file, pattern = "Manual gating") {
  require(data.table)
  separator <- detect_csv_sep(csv_file)
  csv_data <- fread(csv_file, sep = separator, header = TRUE)
  
  # Filter to manual gating columns
  manual_cols <- grep(pattern, names(csv_data), value = TRUE)
  if (length(manual_cols) == 0) return(NULL)
  
  essential_cols <- c("index", manual_cols)
  available_cols <- intersect(essential_cols, names(csv_data))
  csv_data[, ..available_cols]
}

# Process binary labels efficiently with data.table
labels_from_binary_dt <- function(dt) {
  mat <- as.matrix(dt[, -1])  # Exclude index column
  col_names <- names(dt)[-1]
  # Return a list where each element contains the selected column names for that row
  apply(mat, 1, function(row) {
    selected <- col_names[which(row == 1)]
    if (length(selected) == 0) selected <- character(0)
    return(selected)
  }, simplify = FALSE)  # Keep as list to avoid recycling issues
}

# Process small CSV data chunk
process_csv_chunk <- function(csv_data, csv_label_dict, gating_tree, num_events) {
  require(data.table)
  
  idx_list <- csv_data[[1]]  # First column should be index
  setnames(csv_data, names(csv_data), translate_labels(names(csv_data), csv_label_dict))
  
  label_lists <- labels_from_binary_dt(csv_data)
  flat_list <- flatten_label_list(label_lists, gating_tree)
  
  # Create result vector with proper length
  idx_labels_vec <- rep("", num_events)
  
  # Ensure indices are within bounds and translate labels
  valid_indices <- idx_list + 1
  valid_indices <- valid_indices[valid_indices > 0 & valid_indices <= num_events]
  
  if (length(valid_indices) > 0) {
    translated_labels <- translate_labels(flat_list[seq_along(valid_indices)], csv_label_dict)
    idx_labels_vec[valid_indices] <- translated_labels
  }
  
  # Return as data.table with explicit column creation
  result_dt <- data.table(
    cell_id = seq_len(num_events),
    label = idx_labels_vec
  )
  
  # Validate result has correct number of rows
  if (nrow(result_dt) != num_events) {
    stop(sprintf("Row count mismatch in process_csv_chunk: expected %d, got %d", 
                 num_events, nrow(result_dt)))
  }
  
  return(result_dt)
}

# Process large CSV in chunks
process_csv_chunked <- function(csv_data, csv_label_dict, gating_tree, 
                                num_events, chunk_size = 50000) {
  require(data.table)
  idx_labels_vec <- rep("", num_events)
  n_chunks <- ceiling(nrow(csv_data) / chunk_size)
  
  for (i in seq_len(n_chunks)) {
    start_row <- (i - 1) * chunk_size + 1
    end_row <- min(i * chunk_size, nrow(csv_data))
    
    chunk_data <- csv_data[start_row:end_row]
    chunk_result <- process_csv_chunk(chunk_data, csv_label_dict, gating_tree, num_events)
    
    # Only update non-empty labels
    non_empty_idx <- which(chunk_result$label != "")
    if (length(non_empty_idx) > 0) {
      idx_labels_vec[non_empty_idx] <- chunk_result$label[non_empty_idx]
    }
    
    # Cleanup chunk data
    rm(chunk_data, chunk_result)
    if (i %% 10 == 0) {
      gc()
      cat(".")  # Progress indicator
    }
  }
  
  # Return as proper data.table
  result_dt <- data.table(
    cell_id = seq_len(num_events),
    label = idx_labels_vec
  )
  
  # Validate result has correct number of rows
  if (nrow(result_dt) != num_events) {
    stop(sprintf("Row count mismatch in process_csv_chunked: expected %d, got %d", 
                 num_events, nrow(result_dt)))
  }
  
  return(result_dt)
}

# Get file size in MB
get_file_size_mb <- function(file_path) {
  file.info(file_path)$size / 1024^2
}

# Save results and clean memory
save_and_cleanup <- function(result, output_file, cleanup_vars = NULL) {
  require(data.table)
  fwrite(result, output_file, sep = ";")
  if (!is.null(cleanup_vars)) {
    for (var in cleanup_vars) {
      if (exists(var, envir = parent.frame())) {
        rm(list = var, envir = parent.frame())
      }
    }
  }
  gc()
  output_file
}

# Process single CSV file (main function)
process_single_csv <- function(csv_file, pattern, csv_label_dict, gating_tree,
                               csv_to_fcs, output_path, chunk_size = 50000) {
  require(data.table)
  
  file_size_mb <- get_file_size_mb(csv_file)
  cat("Processing:", basename(csv_file), sprintf("(%.1f MB)...", file_size_mb))
  
  # Read and filter CSV
  csv_data <- read_csv_filtered(csv_file, pattern)
  if (is.null(csv_data)) {
    cat(" - No manual gating columns, skipping\n")
    return(NULL)
  }
  
  cat(sprintf(" -> %d cols", ncol(csv_data)))
  
  # Get FCS info
  fcs_file <- csv_to_fcs[[csv_file]]
  num_events <- get_fcs_event_count(fcs_file)
  
  # Process based on size
  if (nrow(csv_data) > chunk_size) {
    cat(sprintf(" (chunked: %d rows)", nrow(csv_data)))
    result <- process_csv_chunked(csv_data, csv_label_dict, gating_tree, 
                                 num_events, chunk_size)
  } else {
    result <- process_csv_chunk(csv_data, csv_label_dict, gating_tree, num_events)
  }
  
  # Validate row count consistency
  if (nrow(result) != num_events) {
    cat(sprintf(" WARNING: Row count mismatch! Expected %d, got %d", num_events, nrow(result)))
  }
  
  # Save and cleanup
  output_file <- file.path(output_path, paste0(basename(csv_file), "_labels.csv"))
  save_and_cleanup(result, output_file, c("csv_data", "result"))
  cat(" -> saved")
  
  output_file
}

validate_csv_results <- function(results) {
  
  # Row count validation summary for saved files
  if (length(results) > 0) {
    cat("\nRow count validation summary:\n")
    for (file_name in names(results)) {
      result_file <- results[[file_name]]
      if (file.exists(result_file)) {
        result_data <- fread(result_file, sep = ";")
        csv_file <- valid_csvs[basename(valid_csvs) == file_name]
        if (length(csv_file) > 0) {
          fcs_file <- csv_to_fcs[[csv_file]]
          expected_rows <- get_fcs_event_count(fcs_file)
          actual_rows <- nrow(result_data)
          status <- if (actual_rows == expected_rows) "✓" else "✗"
          cat(sprintf("  %s %s: %d/%d expected rows \n", 
                     status, file_name, actual_rows, expected_rows))
        }
      }
    }
  }
}
  # Load and combine results (if needed later)
load_csv_results <- function(result_files) {
  require(data.table)
  results <- list()
  for (file in result_files) {
    if (file.exists(file)) {
      sample_name <- gsub("_labels\\.csv$", "", basename(file))
      results[[sample_name]] <- fread(file, sep = ";")
    }
  }
  results
}

# Memory status checker
check_memory_usage <- function(warn_threshold = 0.8) {
  if (Sys.info()["sysname"] == "Linux") {
    mem_info <- system("free -m", intern = TRUE)
    mem_line <- mem_info[2]
    mem_values <- as.numeric(strsplit(mem_line, "\\s+")[[1]][-1])
    used_pct <- mem_values[2] / mem_values[1]

    if (used_pct > warn_threshold) {
      cat("Warning: Memory usage at", round(used_pct * 100, 1), "%\n")
      gc()
    }
  }
}