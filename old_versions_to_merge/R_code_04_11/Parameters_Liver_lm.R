# Input directory path (either absolute or relative to the current working directory path)
input_file_path <- "/home/ubuntu/CellCnn/Data/Liver/ungated_raw/"

# Output directory path(either absolute or relative to the current working directory path)
project_name <- "Liver_ungated"
# project_name <- "Salt_new_0613/"

sample_names <- c("01_DB_0w","01_DL_0w", "01_DS_0w","01_RB_0w","01_RB_3w","01_RB_4w","01_RL_0w", "01_RS_0w", "01_RS_4w",
                  "02_DB_0w","02_RB_0w", "02_RL_0w", "02_RS_0w",
                  "03_DB_0w","03_DL_0w", "03_DS_0w","03_RB_0w","03_RL_0w", "03_RS_0w",
                  "05_DB_0w","05_DL_0w", "05_DS_0w","05_RB_0w","05_RL_0w",
                  "06_DB_1w","06_DL_0w", "06_DS_0w","06_RB_0w","06_RB_1w","06_RB_1y","06_RB_2w","06_RB_3m","06_RB_3w","06_RB_4w","06_RB_5w","06_RB_6m","06_RB_6w","06_RB_7w","06RL_0w","06_RS_0w", "06_RS_1y","06_RS_3m","06_RS_6m",
                  "07_DB_0w","07_DL_0w", "07_DS_0w","07_RB_1w","07_RB_2w","07_RB_3w","07_RB_4w","07_RB_5w","07RB_9m","07_RL_0w", "07_RS_0w","07_RS_9m",
                  "08_DB_0w","08_DB_1w","08_DB_2w","08_DL_0w", "08_DS_0w","08_RB_0w","08_RB_1w","08_RB_1y","08_RB_2w","08_RB_3w","08RB_3m","08_RB_4w","08RB_6m","08RB_9m","08_RL_0w", "08_RS_0w","08_RS_1y","08_RS_3m","08_RS_6m","08_RS_9m",
                  "09_DB_0w","09_DB_1w","09_DL_0w","09_DS_0w","09_RB_0w","09_RB_1w","09_RB_2m","09_RB_2w","09_RB_3w","09_RB_4w","09RB_6m","09_RL_0w", "09_RS_0w","09_RS_2m","09_RS_3w","09_RS_6m",
                  "10_DB_0w","10_DB_1w","10_DL_0w","10_DS_0w","10_RB_0w","10_RB_1w","10RB_2w","10RB_3m","10_RB_3w","10_RB_6m","10_RG_6m","10_RL_0w", "10_RS_0w","10_RS_2w","10_RS_3m","10_RS_6m",
                  "11_DB_0w","11_DL_0w","11_DS_0w","11_RB_0w","11_RL_0w","11_RS_0w",
                  "12_DB_0w","12_DB_3w","12_DL_0w", "12_DS_0w","12_RB_0w","12_RB_1w","12RB_2w","12_RB_3w","12_RB_4m","12_RG_4m","12_RL_0w", "12_RS_0w","12_RS_2w","12_RS_4m"
)

markers_of_interest <- c("CCR7","CD3","CD11c","CD45RA","CD8","CCR6","CD5","CLA",
                         "CD27","CD56","CD69","CD28","CXCR5","CD103","CD45R0",
                         "CD14","PD1","CD11b","gdTCR","CD49a","CD4","CCR4","CD20",
                         "CD95","CD25","CXCR3","CXCR6","CD19","CD127","HLADR",
                         "CD38", "CD45")

# Dictionary of possible markernames stored in the input FCS files (left column) 
# and the actual markernames they should be replaced with in the output (right)
markername_replacing <- list(
  `PD1 (CD279)` = 'PD1',
  `Fas (CD95)` = 'CD95',
  `CCR7 (CD197)` = 'CCR7',
  `TCRgd` = 'gdTCR',
  `HLA-DR` = 'HLADR'
  # Add more mappings as needed
)

#============================== Loading parameters ===================================
truncation = FALSE # Truncate the highest outliers according to the upper $THRESHOLD value
min_limit = NULL # Lower truncation limit

extract_rw_csv = FALSE # Takes up to 5min per file
check_panels = FALSE # Compare pannels with one another

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
homogenize = "None" # accepts either ["before" or "after" or "None"]
extract_ro_csv = FALSE # Save the reordered expression matrix in a CSV (May take a long time)
extract_ro_fcs = TRUE # Saves reordered FCS file (May take a long time)

#======================== Preprocessing parameters =============================

### Choose transformation type: either "arcsinh" "logicle" or "biexponential"
transform_type <- "arcsinh"
cofactor_default <- 150
# Default list of cofactors = all markers will be transformed with the same cofactor
list_of_cofactors <- rep(cofactor_default, times = length(markers_of_interest))

# Choose quality control algorithm: "None","flowAI", "PeacoQC", (not yet implemented: "flowCut" or "flowClean")
qc_type <- "flowAI"

# Choose preprocessing order, either first quality control, then transform "QT" 
# or the other way around "TQ"
pp_order <- "QT"

quant_norm = TRUE
scaling = FALSE

save_pp_fcs = TRUE
save_pp_csv = TRUE

normalize = TRUE
norm_covar = "batch"
