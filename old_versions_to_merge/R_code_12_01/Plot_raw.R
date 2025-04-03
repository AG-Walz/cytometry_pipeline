### ========================= PLOT RAW ============================== ###
# Plot different density plots for raw data to see the profile of the data 
# before any preprocessing.

# Plot separate density plots per marker for each file
if(plot_density_raw_per_file){
  cat("\nPlotting separate density plots per marker for file", sample_names[i])
  plotMarkerDensity("raw_density_", flow_frames[[i]],
                    paste(plot_path, "raw_", sample_names[i], "/", sep = ""))
}

# Plot per sample
if(plot_density_per_samplenames){
  overlap_per_marker(ss_ffr$ffr, paste0(plot_path,"density_per_sample/"),
                     sample_names,target_ncell_subsampling)
}
# Plot per batch
if(plot_density_per_batch){
  overlap_per_marker(ss_ffr$ffr, paste0(plot_path,"batches/"),
                     batch_names,target_ncell_subsampling)
}
# Plot per condition
if(plot_density_per_condition){
  overlap_per_marker(ss_ffr$ffr, paste0(plot_path,"conditions/"),
                     condition_names,target_ncell_subsampling)
}
# ========================= Test Other plots ============================== #
# ggcyto(flow_frames[[1]], aes(x = CD4, y = CD8 )) + 
#   geom_hex(bins = 100) + 
#   ylim(c(-100,4e3)) + 
#   xlim(c(-100,3e3))+              # Set x-axis limits
#   #scale_fill_continuous(name = "Density") +  # Add a legend with a label
#   theme(legend.position = "right")  # Adjust legend position (optional)
# 
# ggcyto(flow_frames[[1]], aes(x = CD4)) +  geom_density(fill = "green")
# 
# autoplot(flow_frames[[1]],"CD3")
# ggcyto(flow_frames[[1]], aes(CD4, CD8)) + geom_hex(bins = 200)
