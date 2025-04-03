# Input directory path (either absolute or relative to the current working directory path)
input_file_path <- "T:/User/Nacho/Data/Cytometry/Endometriosis_data/"


# List all files in input directory ending with .fcs
input_fcs_files <- list.files(path = input_file_path,pattern = ".fcs$",  full.names = TRUE)

# Output directory path(either absolute or relative to the current working directory path)
project_name <- "Endometriosis_24_08_06"


sample_names <- c( "E61","H03","E67", "E68", "E69",  "H04", "H05", "H06", "H07",
                   "H08", "E70", "E71", "E72", "H09", "H10", "H12","E60", "E64", 
                   "E65", "E66", "H01", "H02", "E55", "E58", "E62", "E63"
)

batch_names <- c("Endometriosis","Healthy","Endometriosis", "Endometriosis", "Endometriosis",  "Healthy", "Healthy", "Healthy", "Healthy",
             "Healthy", "Endometriosis", "Endometriosis", "Endometriosis", "Healthy", "Healthy", "Healthy","Endometriosis", "Endometriosis", 
             "Endometriosis", "Endometriosis", "Healthy", "Healthy", "Endometriosis", "Endometriosis", "Endometriosis", "Endometriosis"
)

# Create a dictionary mapping sample names to batches
sample_batch_dict <- setNames(batch_names, sample_names)

markers_of_interest <- c("CD45RO","CD16","CD62L","CD314","CD56","CD8","CD197",
                         "CD123","CD127","CD3", "CD10","CD28","CD196","CD279","CD57",  
                         "CD14","CD45","CD141","TCR.GD","CD195","CD159C","CD4","CD95",
                         "CD25","CD183","CD19",  
                         "CD39","CD38","CD185","CD326","CD11B","CD1A","HLA.DR","CD159A"
                         )

# Dictionary of possible markernames stored in the input FCS files (left column) 
# and the actual markernames they should be replaced with in the output (right)
markername_replacing <- list(
  #  `PD1 (CD279)` = 'PD1',
  #`Fas (CD95)` = 'CD95',
  #`CCR7 (CD197)` = 'CCR7',
  #`TCRgd` = 'gdTCR',
  #`HLA-DR` = 'HLADR'
  # Add more mappings as needed
)

#============================== Loading parameters ===================================
truncation = FALSE # Truncate the highest outliers according to the upper $THRESHOLD value
min_limit = NULL # Lower truncation limit

extract_rw_csv = FALSE # Takes up to 5min per file
check_panels = TRUE # Compare pannels with one another

# Plot the density plots of the raw files
plot_density_per_samplenames = FALSE # per marker for all sample names 
plot_density_per_batch = FALSE # per marker for all subgroup1 members
plot_density_per_condition = FALSE # per marker for all subgroup2 members

## Subsampling
# Target number of cells to subsample for dimensionality reduction
target_ncell_subsampling = 10000
# Subsampling type to choose from "Consistent" or "Maximal"
subsampling_method = "Maximal"
save_downsampling_fcs = FALSE
plot_umap_raw_per_file = FALSE

# Homogenize all flow_frames =  remove all columns that are not "of interest" and
# order those columns to match the order of the variable markers_of_interest
# CAVE: This step currently cannot be run before quality control 
homogenize = "before" # accepts either ["before" or "after" or "None"]
extract_ro_csv = FALSE # Save the reordered expression matrix in a CSV (May take a long time)
extract_ro_fcs = TRUE # Saves reordered FCS file (May take a long time)

#======================== Preprocessing parameters =============================

### Choose transformation type: either "arcsinh" "logicle" or "biexponential"
transform_type <- "logicle"
cofactor_default <- 6000
# Default list of cofactors = all markers will be transformed with the same cofactor
list_of_cofactors <- rep(cofactor_default, times = length(markers_of_interest))

# Choose quality control algorithm: "None","flowAI", "PeacoQC", (not yet implemented: "flowCut" or "flowClean")
qc_type <- "None"

# Choose preprocessing order, either first quality control, then transform "QT" 
# or the other way around "TQ"
pp_order <- "QT"

quant_norm = TRUE
scaling = TRUE

save_pp_fcs = TRUE
save_pp_csv = TRUE

normalize = FALSE
norm_covar = "batch"

