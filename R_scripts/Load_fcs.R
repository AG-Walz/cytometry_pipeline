# Path to the directory containing associated FCS files
fcs_dir <- paste0("/home/ngr24/Desktop/Maerklin/cytometry_pipeline/",
                  "Gating_comparisons/raw_data/")  # Replace with actual path

min_limit <- NULL
truncation <- FALSE
extract_rw_csv <- TRUE
rw_files_path <- out_dir

# List all files in input directory ending with .fcs
input_fcs_files <- list.files(path = fcs_dir, pattern = ".fcs$",  
                              full.names = TRUE)

# Extract the file names from the paths
input_fcs_names <-gsub(fcs_dir, "", input_fcs_files)

# Check that there are .fcs files in the input_file_directory
cat("Reading FCS-files from", fcs_dir,"\n")
if (length(input_fcs_files) == 0) {
  cat("There are no FCS files in the specified input file directory,",
      "please select another folder to proceed.\n")
} else {
  # Read all fcs files into a list of flowFrame objects
  flow_frames <- c()
  for (i in 1:length(input_fcs_files)) {
    flow_frames[[i]] <- flowWorkspace::load_cytoframe_from_fcs(input_fcs_files[i],
                                                               transformation = FALSE, 
                                                               truncate_max_range = truncation,
                                                               min.limit = min_limit)
    # Print number of cells before preprocessing
    cat("Sample",  input_fcs_names[i], "had a total of", 
        nrow(flow_frames[[i]]), "cells.\n")
    # Extract expression matrix and plot density profile per input file
    if(extract_rw_csv){
      cat("Extracting expression matrix of input sample", input_fcs_names[i],"into a CSV file.\n \n")
      write.csv(flow_frames[[i]]@exprs, paste0(rw_files_path, "/", input_fcs_names[i],"_bf_ro.csv"))
    }
  }
}

