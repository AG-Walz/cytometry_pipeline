#============================= General parameters===============================

# Input directory path (either absolute or relative to the current working directory path)
input_file_path <- "T:/User/Nacho/Data/Cytometry/ICK_data/"


# List all files in input directory ending with .fcs
input_fcs_files <- list.files(path = input_file_path,pattern = ".fcs$",  full.names = TRUE)

# Output directory path(either absolute or relative to the current working directory path)
project_name <- "ICK_2"


sample_names <- c( "18_3","18_4","19_3","19_4","20_2","20_3","20_4","21_2","21_3",
                   "21_4a","21_4b","22_2","22_3","22_4","23_1","23_2","23_3","24_2",
                   "24_3","25_1","25_2","25_3a","25_3b","26_1","27_1","27_2","27_3a",
                   "27_3b","28_1","28_2","29_2","29_3","30_1","30_2","30_3","31_1",
                   "31_2","32_1","32_2a","32_2b","33_1","34_1","36_1a","36_1b"
)

batch_names <- c("3","4","3","4","2","3","4","2","3",
                 "4","4","2","3","4","1","2","3","2",
                 "3","1","2","3","3","1","1","2","3",
                 "3","1","2","2","3","1","2","3","1",
                 "2","1","2","2","1","1","1","1"
)

patient_names <- c( "18","18","19","19","20","20","20","21","21",
                    "21","21","22","22","22","23","23","23","24",
                    "24","25","25","25","25","26","27","27","27",
                    "27","28","28","29","29","30","30","30","31",
                    "31","32","32","32","33","34","36","36"
)

# Create a dictionary mapping sample names to batches
sample_batch_dict <- setNames(batch_names, sample_names)

markers_of_interest <- c("OX40","PDL-1" ,"DNAM-1","CD45RO",            
                         "CD16","CD62L (L-Selectin)","CD137 (4-1BB)","CD154 (CD40L)","CD56","CD8",              
                         "LOX-1","CD127 (IL-7Ra)","CD11b","TIGIT","CD28",              
                         "TCRgd","CD57","CD366","CD252 (OX40L)","CD152 (CTLA-4)",    
                         "CD233 _ LAG-3","4-1BBL","CD25","CD66b","KLRG-1","PD-1",              
                         "CD69","HLA-DR","CD19","CD278","CD4"
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

#============================== Loading parameters =============================
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
save_downsampling_fcs = TRUE
plot_umap_raw_per_file = FALSE

# Homogenize all flow_frames =  remove all columns that are not "of interest" and
# order those columns to match the order of the variable markers_of_interest
# CAVE: This step currently cannot be run before quality control 
homogenize = "None" # accepts either ["before" or "after" or "None"]
extract_ro_csv = FALSE # Save the reordered expression matrix in a CSV (May take a long time)
extract_ro_fcs = FALSE # Saves reordered FCS file (May take a long time)

#======================== Preprocessing parameters =============================

### Choose transformation type: either "arcsinh" "logicle" or "biexponential"
transform_type <- "arcsinh"

# Default list of cofactors = all markers will be transformed with the same cofactor
offset_list <- c(33.68, -746.6 ,-1889,-574.1,            
                 80.57, -120.3,753.9,831.1,-532,-432.6,              
                 -1322,923.5,1457,-593.5,315.9,              
                 -1805,-2400,-642.7,-8630,-3034,    
                 -923.5,-289.4,-944.7,-1383,-1414,-1238,              
                 -1058,-401.1,-1507,-2.11,-5253
                 
)#rep(0, times = length(markers_of_interest))

width_list <- c(1800, 2027,3210,1401,            
                4593,4645,2980,2468,1542,1605,              
                1670, 6334,4530,4129,3306,              
                2319,3352,1843,8525,2879,    
                2149,2820,2019,3314,2615,2243,              
                2096,2415,2961,1952,5307
                
)#rep(1000, times = length(markers_of_interest))

# Choose quality control algorithm: "None","flowAI", "PeacoQC", (not yet implemented: "flowCut" or "flowClean")
qc_type <- "None"

# Choose preprocessing order, either first quality control, then transform "QT" 
# or the other way around "TQ"
#pp_order <- "QT"

quant_norm = TRUE
scaling = TRUE

save_pp_fcs = TRUE
save_pp_csv = FALSE

normalize = FALSE
norm_covar = "batch"

