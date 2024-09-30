source("Homogenize.R", chdir = TRUE)

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