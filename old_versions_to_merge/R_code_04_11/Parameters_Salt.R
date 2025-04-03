# Input directory path (either absolute or relative to the current working directory path)
input_file_path <- "/home/ubuntu/CellCnn/Data/Salt/salt_gated/raw/"


# Output directory path(either absolute or relative to the current working directory path)
project_name <- "Salt_24_02_26"


sample_names <- c( "B2+","B4+","B6+", "C2+","C4+","C6+","D2-","D4-","D6-","E2-a", 
                   "E2-b","E3-","E4+","E4-","E5+","E6-","F2-","F3-","F4+","F5+",
                   "G2-","G3-","G4+","G5+","H2-","H3-","H4+","H5+"
)

batches <- c("Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch1", "Batch2", "Batch1", "Batch2", "Batch2", "Batch1", 
             "Batch2", "Batch1", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2", "Batch2")


# Create a dictionary mapping sample names to batches
sample_batch_dict <- setNames(batches, sample_names)

batch_names <- c( "High NaCl","High NaCl","High NaCl", "High NaCl","High NaCl","High NaCl","Low NaCl","Low NaCl","Low NaCl","Low NaCl", 
                   "Low NaCl","Low NaCl","High NaCl","Low NaCl","High NaCl","Low NaCl","Low NaCl","Low NaCl","High NaCl","High NaCl",
                   "Low NaCl","Low NaCl","High NaCl","High NaCl","Low NaCl","Low NaCl","High NaCl","High NaCl"
)

markers_of_interest <- c("CXCR6","CD38", "CD127", "HLA-DR", "LAG3", "CD57" ,"CD45R0",
                         "CD49a", "CCR6", "CCR7" ,"CD3", "CD11c", "CD45RA" ,
                         "CXCR3", "CD27", "CD56", "CTLA4", "CD69" ,"CD28", "CXCR5" ,
                         "CD103", "CD62L", "CD25", "CD95", "CD31", "CCR4", "ICOS", 
                         "PD1", "CD11b", "gdTCR", "CD5", "CD8", "CD4") # CD4, CD45, CD14, CD20, CD3, CD8

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

