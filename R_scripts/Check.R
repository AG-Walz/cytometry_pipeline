### ============================ CHECK RAW DATA =========================== ###
# This script is responsible for comparing panels in the raw cytometry data
# and in case homogenize = TRUE, homogenizing the data before downsampling.

# ============================ Source scripts  ============================ #

source_files <- c(parameter_file,"Compare.R", "Create_ffr.R", "Homogenize.R", 
                 "Downsampling_functions.R") # UMAP.R

lapply(source_files, source, chdir = TRUE)

rw_files_path <- paste0(results_path, "/raw_files/")

# ============================ Compare panels ============================ #

# Compare the different panels within the list of flowframes
## Iterate over all pairwise combinations of flow frames and check if their 
## markernames are equal, adding them to an existing subgroup if they are or
## create a new one if they aren't. Finally, print the subgroups and 
## markernames within each subgroup.
if (length(flow_frames) > 1) {
  if(check_panels){
    panel_dict <- compare_panels(flow_frames)
    cat(sprintf("Among the %d FCS-files, there are a total of %d different", 
                length(flow_frames),length(panel_dict)), 
        "panels. \n You can access them",
        sprintf("by typing <panel_dict[1:%d]>.\n",length(panel_dict)))
    message("Consensus of all markers present in all panels:")
    print(Reduce(intersect, panel_dict))
    print("Unique markers present in any list:")
    print(unique(unlist(panel_dict)))
  }
} else {print("You only inputted one file, there is nothing to compare. \n \n")}

# ============================ Homogenize data ============================ #
# Create metadata.csv for all files
metdat <- data.frame(filename = paste0( sample_names,"_af_ro.fcs")
                     #                      ,batch_names, condition_names
)
metdat_filepath <- paste0(rw_files_path, "reordered_metadata.csv")
write.csv(write.csv(metdat, file = metdat_filepath, row.names = FALSE))

# Homogenize all flowframes if necessary
if(homogenize == "before"){
  flow_frames <- homogenize_ffl(flow_frames,
                                sample_names, 
                                extract_ro_fcs,
                                extract_ro_csv,
                                markers_of_interest,
                                markername_replacing,
                                rw_files_path)
}

# ============================ Downsample data ============================ #
# Downsampling
ss_ffr <- subsample_flowframes(flow_frames, 
                               target_ncell_subsampling, 
                               subsampling_method, 
                               paste0(results_path, "/downsampled_files/"),
                               sample_names)


