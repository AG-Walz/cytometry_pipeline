# Load the Compare.R script to access concatenate_names function
source("Compare.R")

# Concatenate expression matrices
cc_ff <- concatenate_exprs_mtx(flow_frames)

# Append sample names to the concatenated expression matrix 
cc_sample_names <- concatenate_names(flow_frames)

#cc_ff$sample_names <- cc_sample_names
  