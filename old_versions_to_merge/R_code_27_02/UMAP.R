source_files <- c(parameter_file,"Compare.R", "UMAP_functions.R", "Downsampling_functions.R") # UMAP.R

lapply(source_files, source, chdir = TRUE)

results_path <- ifelse(is.na(file.info(project_name)$isdir) || !file.info(project_name)$isdir,
                       file.path(getwd(), "Results", project_name), project_name)
umap_dir <- paste0(results_path, "/umap_plots/")
dir_create(umap_dir)

# df <- data.table::fread(paste0(umap_dir, "/df.csv"))
# umap_df<- data.table::fread(paste0(umap_dir,"/umap_df.csv"))
# umapdf_with_labels<- data.table::fread(paste0(umap_dir,"/umapdf_with_labels.csv"))

marker_indices_of_interest <-explore_data(flow_frames[[1]])$mcols[c(1:8, 10:33)]

list_of_ffr <- flow_frames #subsample_flowframes(flow_frames,target_nrow = 1000)
samples <- sample_names
condition <- batch_names
patient <- patient_names


#===============================================================================

# Create a dictionary mapping values from vector1 to vector2
cond_dictionary <- setNames(condition, samples)
pati_dictionary <- setNames(patient, samples)

# Extracting the number of rows from each flowframe's expression matrix
num_rows <- sapply(list_of_ffr, function(flowframe) nrow(exprs(flowframe)))

# Replicating each letter based on the number of rows in the corresponding flowframe
sample_procedence <- unlist(mapply(rep, samples, num_rows))

# Condition procedence using condition dictionary
condition_procedence <- cond_dictionary[sample_procedence]
patient_procedence <- pati_dictionary[sample_procedence]

#Concatenate all expression matrices in the list of flowframes
df <-concatenate_exprs_mtx(list_of_ffr)

# Remove all instances (rows) with NA values
na_rows <- which(apply(df, 1, anyNA))  # Rows with at least one NA
na_cols <- which(apply(df, 2, anyNA))  # Columns with at least one NA
#condition_procedence  <- condition_procedence[-na_rows]
#patient_procedence  <- patient_procedence[-na_rows]
df <- na.omit(df)

# Subset the dataframe
#df_subset <- as.data.frame(df[c(sample(nrow(df), 10000)),]) %>% 
#  select(all_of(marker_indices_of_interest))
df_subset <- as.data.frame(df) %>% 
    select(all_of(marker_indices_of_interest))
data.table::fwrite(df_subset, paste0(umap_dir,"/df.csv"))

# Compute the UMAP embedding for the dataframe 
umap_df <- uwot::umap(X=df_subset,n_neighbors = 15, min_dist = 0.1)
colnames(umap_df) <- c("UMAP1", "UMAP2")
data.table::fwrite(umap_df, paste0(umap_dir,"/umap_df.csv"))


umapdf_with_labels <- cbind(umap_df, sample_procedence, condition_procedence, patient_procedence)
data.table::fwrite(umapdf_with_labels, paste0(umap_dir,"/umapdf_with_labels.csv"))


#=============================================================================
subset <- umapdf_with_labels[sample(nrow(umapdf_with_labels), 100000),]

custom_palette <- c("darkred","red","darkorange","gold", "greenyellow","limegreen","darkgreen",  "#1F9698",
                    "navy","cyan2", "#009FFF", "deeppink", "purple",  "chocolate4", "black")

# Assign colors explicitly to the levels of Patient_ID
color_sample <- setNames(custom_palette, levels(sample_names))
color_patient <- setNames(custom_palette, levels(patient_names))


ggplot(data.frame(subset), aes(x=as.numeric(UMAP1), y= as.numeric(UMAP2))) +
  geom_point(shape = 19, size = 0.1, stroke = 0.1) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme_bw()

samp_plot <- ggplot(data.frame(subset), 
               aes(as.numeric(UMAP1), y= as.numeric(UMAP2),color = sample_procedence)) +
  geom_point(shape = 19, size = 0.1, stroke = 0.1) +
  labs(x = "UMAP1", y = "UMAP2",color = "Sample") +
  #scale_color_manual(values = c("deepskyblue", "brown1"))+
  guides(colour = guide_legend(override.aes = list(size=5),  ncol = 2))+
  ggplot2::scale_color_manual(values = color_sample) +
  theme_bw()

# Save the plot as a PNG file
ggsave(plot=samp_plot,filename=paste0(umap_dir,"/sample_procedence.pdf"),width = 6.05, height = 3.9, units = "in")

cond_plot <- ggplot(data.frame(subset), 
               aes(as.numeric(UMAP1), y= as.numeric(UMAP2), 
                   color = condition_procedence)) +
  geom_point(shape = 19, size = 0.1, stroke = 0.1) +
  labs(x = "UMAP1", y = "UMAP2",color = "Timepoint") +
  #scale_color_manual(values = c("deepskyblue", "brown1"))+
  guides(fill=guide_legend(ncol=1, byrow=TRUE))+
  guides(colour = guide_legend(override.aes = list(size=5),  ncol = 1))+
  theme_bw()

# Save the plot as a PNG file
ggsave(plot=cond_plot,filename=paste0(umap_dir,"/condition_procedence.pdf"), width = 6.05, height = 3.9, units = "in")


pat_plot <- ggplot(data.frame(subset), 
               aes(as.numeric(UMAP1), y= as.numeric(UMAP2), 
                   color = patient_procedence)) +
  geom_point(shape = 19, size = 0.1, stroke = 0.1) +
  labs(x = "UMAP1", y = "UMAP2",color = "Patient") +
  #scale_color_manual(values = c("deepskyblue", "brown1"))+
  guides(fill=guide_legend(ncol=1, byrow=TRUE))+
  guides(colour = guide_legend(override.aes = list(size=5),  ncol = 2))+
  ggplot2::scale_color_manual(values = color_patient) +
  theme_bw()

# Save the plot as a PNG file
ggsave(plot=pat_plot,filename=paste0(umap_dir,"/patient_procedence.pdf"), width = 6.05, height = 3.9, units = "in")


#for (i in 1:length(markers_of_interest)){
#  print(markers_of_interest[i])
#  umap_marker(umap_df, concat_subsample, paste("UMAP coloured by ", markers_of_interest[i], " expression", sep=""), i, markers_of_interest[i])
#}

#umap_marker(umap_df, concat_subsample, "test_1", 1, "CD49a")

file_precedence <- rep(sample_names,each=downsample_size)
umap_vector(umap_df, file_precedence,"file_procedence", rainbow(16)) 

patient <- rep(c("E","F","G", "H"), each=4*downsample_size)
umap_vector(umap_df, patient,"Patient       ", c("red", "yellow", "green", "blue"))


umap_vector(umap_df, salt_concentration,"Salt condition", c("aquamarine", "lightsalmon"))
tsne_vector(dimred_res$layout[,c(1,2)], salt_concentration,"Salt condition", c("aquamarine", "lightsalmon"))
umap_vec(tsne_dim$Y, salt_concentration,"Salt condition", c("aquamarine", "lightsalmon"))

#umap_vector(umap_df, cluster_origin,"cluster origin")


umap_vector(umap_df, udn_vector,"Cluster group ", c("Blue", "Grey", "Red"))

metacluster_vector <- GetMetaclusters(fsom)
metacluster_origin <- metacluster_vector[subsampling_idx]
umap_vector(umap_df, metacluster_origin,"metacluster origin", rainbow(22))
print(table(udn_vector))

