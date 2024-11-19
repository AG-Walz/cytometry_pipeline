#=============================== COMPARE.R ==============================

# This file contains a collection of helper functions for analyzing flow cytometry data in R.
# The functions included in this file are designed to assist with various tasks such as:
# - Checking if a path is a directory and creating directories if they do not exist.
# - Extracting information from flowFrame objects, such as dye and marker names.
# - Comparing panels of flow cytometry data to group flow frames with identical marker sets.
# - Reading and parsing the text segment of FCS (Flow Cytometry Standard) files.
# - Retrieving specified parameters from FCS files.
# - Checking the total number of events in a flow frame against the header information.

#==================== Helper Functions ==============================

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

# Remove unused channels
# for (i in 1:length(flow_frames)){flow_frames[[i]] <-cyto_channels_restrict(flow_frames[[i]])}
# for(i in 1:length(flow_frames)){print(ncol(flow_frames[[i]]@exprs))}
# all.equal(explore_data(flow_frames[[1]])$mcols,explore_data(flow_frames[[17]])$mcols)
