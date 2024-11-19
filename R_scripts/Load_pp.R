### ============================ LOAD INPUT DATA =========================== ###

# Import loading functions
source_files <- c(parameter_file,"Compare.R", "Create_ffr.R", "Homogenize.R", 
                  "Density_plots.R", "Downsampling.R") # UMAP.R

lapply(source_files, source, chdir = TRUE)

## ====================== Create output directory paths ===================== ##
results_path <- ifelse(is.na(file.info(project_name)$isdir) || !file.info(project_name)$isdir,
                       file.path(getwd(), "Results", project_name), project_name)

paths <- c(
  results_path,
  paste0(results_path, "/raw_files/"),
  paste0(results_path, "/raw_plots/"),
  paste0(results_path, "/downsampled_files/")
)

dir_create(paths)

## ========================= Read & Check FCS files ==========================##

# List all files in input directory ending with .fcs
input_fcs_files <- list.files(path = input_file_path,pattern = ".fcs$",  full.names = TRUE)

# Extract the file names from the paths
input_fcs_names <-gsub(input_file_path, "", input_fcs_files)

# Check that there are .fcs files in the input_file_directory
cat("Reading FCS-files from", input_file_path,"\n")
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
    cat("Sample", sample_names[i],"=", input_fcs_names[i], "had a total of", 
        nrow(flow_frames[[i]]), "cells.\n")
    # Extract expression matrix and plot density profile per input file
    if(extract_rw_csv){
      cat("Extracting expression matrix of input sample", sample_names[i],"into a CSV file.\n \n")
      write.csv(flow_frames[[i]]@exprs, paste0(rw_files_path, "/", sample_names[i],"_bf_ro.csv"))
    }
  }
}

