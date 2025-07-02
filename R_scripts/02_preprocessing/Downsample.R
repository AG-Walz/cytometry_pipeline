source("Downsampling_functions.R", chdir = TRUE)

# Downsampling
ss_ffr <- subsample_flowframes(flow_frames, 
                               target_ncell_subsampling, 
                               subsampling_method, 
                               paste0(results_path, "/downsampled_files/"),
                               sample_names)