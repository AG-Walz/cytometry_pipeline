# PREPROCESSING

# Load required functions
source(parameter_file, chdir=TRUE)
source("Preprocessing_functions.R", chdir = TRUE)

# Create folder to store output preprocessed files & another for quality control results
pp_output_path <-paste0(results_path, "/preprocessed_files/")
qc_output_path <- switch(qc_type,
                         flowAI = file.path(results_path, "quality_control"),
                         PeacoQC = file.path(results_path, "PeacoQC_results", "fcs_files"),
                         results_path)
dir_create(c(pp_output_path,qc_output_path))

# ========================= Preprocessing ======================================

# if(optimize_per_marker){
#   list_of_cofactors <- get_opt_b_cofactor(qn_new_ffL,
#                                                   target_bimodality_coefficient,
#                                                   markers_of_interest)
# }else{
#   print(paste0("The cofactor for all markers will be by default: ",
#                cofactor_default))
#   list_of_cofactors <-rep(cofactor_default, times = length(markers_of_interest))
#   }

pp_flowframes <- list()

#function_order <- if (pp_order == "QT") {
#  list(quality_control,transform_flowframe)
#} else if (pp_order == "TQ") {
#  list(transform_flowframe,quality_control)
#} else {
#  stop("Please give preprocessing order as either 'QT' (Quality -> Transform) or 'TQ' (Transform -> Quality).")
#}

for (frame_index in 1:length(flow_frames)) {
  cat("\nPreprocessing frame:", frame_index, "\n")
  
  # List of marker positions
  marker_indices_of_interest <- explore_data(flow_frames[[frame_index]])$mcol #Check if true
  
  # Execute functions in the specified order
  # quality_control(flow_frames[[frame_index]], sample_names[frame_index])
  # Apply the first function in function_order
  pp_flowframes[[frame_index]] <- quality_control(flow_frames[[frame_index]], sample_names[frame_index])
  
  # Apply the second function in function_order to the result of the first
  pp_flowframes[[frame_index]] <- transform_flowframe(pp_flowframes[[frame_index]], sample_names[frame_index],
                                                      markers_of_interest, offset_list, width_list)
  
  if (quant_norm){
    cat("Quantile normalizing", sample_names[frame_index], "\n")
    pp_flowframes[[frame_index]] <- quantile_norm(pp_flowframes[[frame_index]],
                                          marker_indices_of_interest)
  }
  if (scaling){
    cat("Scaling", sample_names[frame_index], "\n")
    pp_flowframes[[frame_index]] <- z_score_scaling(pp_flowframes[[frame_index]],
                                            marker_indices_of_interest)
  }
  if (save_pp_fcs){
    # cyto_channels_restrict(pp_flowframes[[frame_index]])
    cat("Saving preprocessed FCS file of sample", sample_names[frame_index],". \n")
    flowCore::write.FCS(pp_flowframes[[frame_index]], 
                        paste0(pp_output_path, "pp_", 
                               sample_names[frame_index],".fcs"))
  }
  # marker_data <- pp_flowframes[[frame_index]]@exprs[,marker_indices_of_interest]
  # colnames(marker_data)<- sort(markernames(pp_flowframes[[frame_index]][,marker_indices_of_interest]))
  if (save_pp_csv){  
    cat("Saving preprocessed expression matrix of file", sample_names[frame_index],"as a CSV file. \n\n")
    write.csv(pp_flowframes[[frame_index]]@exprs, 
              paste0(pp_output_path, "pp_",
                     sample_names[frame_index], ".csv"))
  }
}
if (homogenize == "after"){
  hg_pp_flowframes <- homogenize_ffl(pp_flowframes, 
                                     sample_names, 
                                     extract_ro_fcs,
                                     extract_ro_csv,
                                     markers_of_interest,
                                     markername_replacing,
                                     pp_output_path)
}
# Plot per sample
if(plot_density_per_samplenames){
  overlap_per_marker(pp_flowframes, paste0(plot_path,"pp_density_per_sample/"),
                     sample_names,target_ncell_subsampling)
}
# Plot per batch
if(plot_density_per_batch){
  overlap_per_marker(pp_flowframes, paste0(plot_path,"pp_batches/"),
                     batch_names,target_ncell_subsampling)
}
# Plot per condition
if(plot_density_per_condition){
  overlap_per_marker(pp_flowframes, paste0(plot_path,"pp_conditions/"),
                     condition_names,target_ncell_subsampling)
}
