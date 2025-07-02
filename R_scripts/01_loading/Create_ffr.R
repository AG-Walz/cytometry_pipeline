### ============================ GENERATE FLOWFRAMES =========================== ###
# This script allows to generate new flowframes including their metadata
#
# Credit to Yann Abraham for the original implementation
# Source: https://gist.github.com/yannabraham/c1f9de9b23fb94105ca5
#
extract_exprmtx_list <- function(list_of_ffr){
  #' Extract Expression Matrices from a List of FlowFrames
#'
#' This function takes a list of flowFrame objects and extracts the expression matrices from each flowFrame.
#' The expression matrices are returned as a list of data frames.
#'
#' @param list_of_ffr A list of flowFrame objects from which to extract expression matrices.
#' 
#' @return A list of data frames, each containing the expression matrix of the corresponding flowFrame.
#' @importFrom flowCore exprs
#' 
  expression_mtx_list <- lapply(list_of_ffr, function(flowframe) {
    as.data.frame(flowCore::exprs(flowframe))
    })
  return(expression_mtx_list)
}


generate_cellID <- function(list_of_ffr){
  #' Generate Cell IDs from a List of FlowFrame Objects
#'
#' This function takes a list of FlowFrame objects and generates unique cell IDs for each cell in the FlowFrames.
#'
#' @param list_of_ffr A list of FlowFrame objects from which cell IDs will be generated.
#' @return A list of FlowFrame objects with an additional column for cell IDs.
#'
#' @details
#' The function iterates over each FlowFrame in the provided list and assigns a unique cell ID to each cell.
#' This is useful for tracking individual cells across different FlowFrames.
  expression_mtx_list <- extract_exprmtx_list(list_of_ffr)
    
  # Initialize a variable to keep track of the cell ID
  total_cells <- 0
  cell_ID <- c()
  # Loop through the list of dataframes
  for (i in seq_along(expression_mtx_list)) {
    # Get the current dataframe
    current_df <- expression_mtx_list[[i]]
    
    # Add a "cell_ID" column with unique numbering
    cell_ID[[i]] <- c((total_cells + 1):(total_cells + nrow(current_df)))
    
    # Update the total rows count
    total_cells <- total_cells + nrow(current_df)
    
    # Update the dataframe in the list
    expression_mtx_list[[i]] <- current_df
  }
  
  # Print the modified list of dataframes
  return(cell_ID)
}


from_mtx_to_ffr <- function(data_mtx, ordered_markernames, sample_name, cellID){
  #' Convert Matrix to FlowFrame
#'
#' This function converts a given data matrix into a FlowFrame object.
#'
#' @param data_mtx A numeric matrix where rows represent cells and columns represent markers.
#' @param ordered_markernames A character vector of marker names in the order they appear in the matrix.
#' @param sample_name A character string representing the name of the sample.
#' @param cellID A character vector or numeric vector representing the cell IDs.
#'
#' @return A FlowFrame object containing the data from the input matrix.
#'

  print(paste0("Generating a reordered and renamed flowframe for ", sample_name,sep = ""))
  # Appending additional columns to the expression matrix
  rownames(data_mtx) <- cellID
  # Replace markernames of matrix for those provided
  #dimnames(data_mtx)[[2]] <- ordered_markernames
  dimnames(data_mtx)[[2]] <-c(ordered_markernames)
  meta <- data.frame(name=dimnames(data_mtx)[[2]], desc=dimnames(data_mtx)[[2]])
  meta$range <- apply(apply(data_mtx,2,range),2,diff)
  meta$minRange <- apply(data_mtx,2,min)
  meta$maxRange <- apply(data_mtx,2,max)
  
  #meta <- data.frame(name=dimnames(data_mtx)[[2]], desc=dimnames(data_mtx)[[2]])
  #meta$minRange <- apply(data_mtx[,1:length(ordered_markernames)],2,min)
  #meta$maxRange <- apply(data_mtx[,1:length(ordered_markernames)],2,max)
  #meta$range <- as.numeric(meta$maxRange) - as.numeric(meta$minRange)
  write.csv(x = meta,file = paste0(rw_files_path,"/",sample_name,"_metadata.csv"))
  ff <- new("flowFrame",
            exprs=data_mtx,
            parameters=Biobase::AnnotatedDataFrame(meta)
  )
  cat("New flowframe created. \n")
  markernames(ff)<-setNames(ordered_markernames, ordered_markernames)
  return(ff)
}


from_ffl_to_fSet <- function(ffl){
  #' Convert FlowFrame List to FlowSet
#'
#' This function takes a list of FlowFrames (ffl) and converts it into a FlowSet.
#'
#' @param ffl A list of FlowFrame objects.
#' @return A FlowSet object containing the FlowFrames from the input list.
#' 
  flowSet(ffl)
}


extract_exprmtx_list <- function(list_of_ffr){
  #' Extract Expression Matrices from a List of FlowFrames
#'
#' This function takes a list of flowFrame objects and extracts the expression matrices from each flowFrame.
#' The expression matrices are returned as a list of data frames.
#'
#' @param list_of_ffr A list of flowFrame objects from which to extract expression matrices.
#' 
#' @return A list of data frames, each containing the expression matrix of the corresponding flowFrame.
#' @importFrom flowCore exprs
#' 
  expression_mtx_list <- lapply(list_of_ffr, function(flowframe) {
    as.data.frame(flowCore::exprs(flowframe))
    })
  return(expression_mtx_list)
}






