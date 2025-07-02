### ========================= HOMOGENIZE ============================== ###
# This script includes functions to remove panel differences from flow frames


translate_vector <- function(vector, dictionary){
#' Translates elements of a vector based on a provided dictionary.
#'
#' @param vector A vector of elements to be translated.
#' @param dictionary A named list where names are original elements and values are their translations.
#' @return A vector with translated elements. If an element is not found in the dictionary, "Unknown" is returned.
#'
  translated_vector <- sapply(vector, function(letter) {
    if (letter %in% names(dictionary)) {
      return(dictionary[[letter]])
    } else {
      return("Unknown")  # Handle cases where the letter is not in the dictionary
    }
  })
}


replace_markernames <- function(flowframe, replace_map){
#' Replace Marker Names in a FlowFrame
#'
#' This function replaces marker names in a given flowFrame object based on a provided mapping.
#'
#' @param flowframe A flowFrame object containing cytometry data.
#' @param replace_map A named list or vector where the names are the original marker names and the values are the new marker names.
#'
#' @return A flowFrame object with updated marker names.
#' 
  marker_list <-markernames(flowframe)
  marker_pos <- explore_data(flowframe)$mcol
  for (i in 1:length(marker_list)) {
    if (marker_list[i] %in% names(replace_map)) {
      new_marker <- replace_map[[marker_list[i]]]
      print(paste0("Marker ",marker_list[i], " was replaced for ", new_marker))
      names(new_marker) <- featureNames(flowframe)[marker_pos[i]]
      markernames(flowframe)[i] <- new_marker
    }
  }
  return(flowframe)
}



reorder_exprmtx <-function(ffr, marker_name_list){
#' Reorder Expression Matrix
#'
#' This function reorders the expression matrix based on a given list of marker names.
#'
#' @param ffr A flowFrame object containing the expression matrix to be reordered.
#' @param marker_name_list A character vector containing the names of the markers in the desired order.
#'
#' @return A reordered expression matrix.
#' 
  exprs_matr <- data.frame(matrix(NA, nrow = nrow(flowCore::exprs(ffr)),
                                              ncol = length(marker_name_list)))
  colnames(exprs_matr)<-marker_name_list
  for(markernam in 1:length(marker_name_list)){
    column_index <- which(explore_data(ffr)$markers 
                          == marker_name_list[markernam])
    exprs_matr[,markernam] <- flowCore::exprs(ffr)[, column_index]
  }
  # reordered_ffr[[ffr]] <-flowFrame(exprs = as.matrix(exprs_matr))
  reordered_mtx <- as.matrix(exprs_matr)
  return(reordered_mtx)
}



replace_reorder <- function(ffr, markers, marker_replacing, sample_name, cellID){
#' Replace and Reorder Markers in Flow Frame
#'
#' This function replaces specified markers in a flow frame and reorders them.
#'
#' @param ffr A flow frame object containing cytometry data.
#' @param markers A character vector of marker names to be reordered.
#' @param marker_replacing A named character vector where names are the markers to be replaced and values are the new marker names.
#' @param sample_name A character string specifying the sample name.
#' @param cellID A character string or numeric value specifying the cell ID.
#'
#' @return A modified flow frame object with markers replaced and reordered.
#'
  # Replace markernames
  rp_ff <- replace_markernames(ffr, marker_replacing)
  
  # Create an expression matrix with reordered markernames
  exprs_mtx <- reorder_exprmtx(rp_ff, markers)
  
  # Integrate the new expression matrix into a new flowframe
  new_ff <- from_mtx_to_ffr(exprs_mtx, markers, sample_name, cellID)
  return(new_ff)
}



reorder_flowframes <-function(ffr_list, marker_name_list, replace_mapping){
#' Reorder FlowFrames
#'
#' This function reorders a list of FlowFrames based on a given list of marker names and a replacement mapping.
#'
#' @param ffr_list A list of FlowFrame objects to be reordered.
#' @param marker_name_list A list of marker names used for reordering the FlowFrames.
#' @param replace_mapping A named list or vector where names are the original marker names and values are the new marker names.
#'
#' @return A list of reordered FlowFrame objects.
#'
  reordered_ffr <- vector(mode = "list", length = length(ffr_list))
  for(ffr in 1:length(ffr_list)){
    print(paste0("Reordering flowframe ",ffr))
    reordered_ffr[[ffr]] <- replace_reorder(ffr_list[[ffr]], marker_name_list, 
                                            replace_mapping)
  }
  return(reordered_ffr)
}



homogenize_ffl <- function(ff_list, sample_nams, expression_fcs, expression_csv,
                           marker_name_list, replace_mapping, output_path){
#' Homogenize FlowFrame List
#'
#' This function homogenizes a list of FlowFrames by applying a series of transformations and mappings.
#'
#' @param ff_list A list of FlowFrame objects to be homogenized.
#' @param sample_nams A vector of sample names corresponding to the FlowFrames.
#' @param expression_fcs A vector of file paths to the FCS files containing expression data.
#' @param expression_csv A vector of file paths to the CSV files containing expression data.
#' @param marker_name_list A list of marker names to be used for homogenization.
#' @param replace_mapping A named list where names are original marker names and values are the new marker names to replace them.
#' @param output_path A string specifying the directory path where the output files will be saved.
#'
#' @return A list of homogenized FlowFrame objects.
#' 
  #Generate cellID for each flowframe
  cell_IDs <- generate_cellID(ff_list)
  
  # Homogenize all flowframes if necessary
  homogene_ffL <- list()
  if(length(ff_list)<= 1){
    cat("Flowframe homogenization skipped automatically")
  } else{
    for(i in 1:length(ff_list)){
      
      if(i==1){cat("\nHomogenizing all flowframes to have the same panel according to the markers of interest:\n")}
      cat("Replacing markernames of", sample_nams[i], ":\n")
      homogene_ffL[[i]] <-replace_reorder(ff_list[[i]],marker_name_list,
                                          replace_mapping, sample_nams[i],
                                          cell_IDs[[i]])
      
      if(expression_fcs){
        cat("Saving reordered and renamed", sample_nams[i],"flowframe as a FCS file.\n")
        write.FCS(homogene_ffL[[i]], paste0(output_path, "/", sample_nams[i],"_af_ro.fcs"))
      }
      if(expression_csv){
        cat("Saving reordered and renamed ", sample_nams[i]," flowframe as a CSV file.\n \n")
        write.csv(homogene_ffL[[i]]@exprs, paste0(output_path, "/", sample_nams[i], "_af_ro.csv"))
      }
    } 
  }
  return(homogene_ffL)
}


