is_directory <- function(path) {
  if (!file.exists(path)) {return(FALSE)}
    info <- file.info(path)
    return(info$isdir)
}

dir_create <- function(dir_path_list){
  for(path in dir_path_list){
    if (!endsWith(path, "/")){
      path <- paste0(path, "/")
      }
    if (!is_directory(path)){
      dir.create(path, showWarnings = FALSE) # Create directory 
    }
  }
}


#==================== Extract information from FF ==============================
explore_data <- function(flowframe){
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

# Check whether 2 markernames are equal or not
#  all.equal(markernames(flow_frames[[1]]),markernames(flow_frames[[length(flow_frames)]]))

compare_panels <- function(flowframes_list) {
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
