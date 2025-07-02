# Preprocessing functions
source("Homogenize.R", chdir = TRUE)
source("Create_ffr.R", chdir = TRUE)

#=========================== Quality control ===================================
# Quality control function
quality_control <- function(flowframe, sample_name) {
  qc_df <- switch(
    qc_type,
    "flowAI" = {
      cat("Sample",sample_name,"will be quality controlled with FlowAI. ")
      qc_df <- flowAI::flow_auto_qc(flowframe,
                                    folder_results = paste(qc_output_path,"/FlowAI_", sample_name, sep = ""))
      qc_output_dir <- dirname(qc_output_path)
      if (!dir.exists(qc_output_dir)) {
        dir_create(qc_output_dir)
      }
      flowCore::write.FCS(qc_df, paste0(qc_output_path,"/FlowAI_", sample_name, ".fcs"))
      qc_df
      cat(paste0("After QC ", nrow(qc_df), "/", nrow(flowframe)," cells remained."))
    },
    "PeacoQC" = {
      cat("Sample",sample_name,"will be quality controlled with PeacoQC. \n")
      qc_df <- PeacoQC::PeacoQC(flowframe, channels = marker_indices_of_interest, 
                                plot = TRUE, save_fcs = TRUE, 
                                determine_good_cells = "all", 
                                output_directory = qc_output_path)
      cat(paste0("After QC ", nrow(qc_df), "/", nrow(flowframe)," cells remained."))
      qc_df$FinalFF
    },
    "None" = {
      cat("Sample",sample_name,"will not be quality controlled. \n")
      flowframe
    },
    stop("The given quality control type is not valid. Try with PeacoQC or flowAI.")
  )
  return(qc_df)
}

#============================ Transformation ===================================
# Define a function for arcsinh transformation
transform_flowframe <- function(flowframe,  samplename, markers_of_interest, offste_list, width_list) {
  # Check if the length of markers and list_of_cofactors is the same
  if (length(markers_of_interest) != length(offset_list)) {
    stop("Marker and offset lists must have the same length.")
  }else if (length(markers_of_interest) != length(width_list)) {
    stop("Marker and width lists must have the same length.")
  }
  cat(paste(transform_type,"transform sample", samplename, "\n"))
  tf_flowframe <- flowframe
  # Iterate through markers and apply arcsinh transformation
  for (marker_index in 1:length(markers_of_interest)) {
    marker <- markers_of_interest[marker_index]
    offset <-offste_list[marker_index]
    width <- width_list[marker_index]
    cat(paste("Marker", marker,"transformed with offset", offset,"and width", width, "\n"))
    # Check if the marker exists in the flow frame
    if (!(marker %in% colnames(flowframe))) {
      warning(paste("Skipping marker", marker, "for it was not found in the flow frame."))
      next
    }else{cat(paste0("Transforming marker", marker, "with offset", offset, "and width", width))}
    
    # Transformation type
    trans <- switch(
      transform_type,
      "logicle" = flowCore::logicleTransform(),
      # f_asinh(x) = arcsinh((x-offset)/width)
      "arcsinh" = flowCore::arcsinhTransform(a = -offset, b = 1 / width),
      "biexponential" = flowCore::biexponentialTransform(),
      "None" = NULL,
      stop("The given transformation type is not valid. Try with arcsinh, biexponential, logicle, or None.")
    )
    
    if (!is.null(trans)) {
      tflist <-transformList(from = marker, tfun = trans)
      
      tf_flowframe <- flowCore::transform(tf_flowframe, tflist,truncate_max_range = FALSE)
    }
  }
  return(tf_flowframe)
}

# =========================== OPTIMIZE BIMODALITY ==============================

# Function to calculate the bimodality coefficient
calculate_bimodality_coefficient <- function(b_cofactor, data, target_coefficient) {
  transformed_data <- asinh(data * b_cofactor)
  bicoeff <- bimodality_coefficient(transformed_data)
  return(abs(bicoeff - target_coefficient))
}

# Function to find the optimal b cofactor
find_optimal_b_cofactor <- function(data, target_bimodality_coefficient = 0.9) {
  # Set the optimization boundaries for b (modify if needed)
  lower_bound <- 0.001
  upper_bound <- 10000
  
  # Use optimization to find the optimal b_cofactor
  optimal_b_result <- optimize(
    f = calculate_bimodality_coefficient,
    interval = c(lower_bound, upper_bound),
    data = data,
    target_coefficient = target_bimodality_coefficient,
    maximum = FALSE
  )
  
  return(optimal_b_result$minimum)
}

get_opt_b_cofactor <- function(flowframe_list, target_bimodality_coefficient = 0.9,
                              markers_of_interest){
  conc_eM <- concatenate_exprs_mtx(flowframe_list)
  cofactors <- list()
  for(j in 1:ncol(conc_eM)){
    column <- conc_eM[,j]
    cofactors[j] <-find_optimal_b_cofactor(column, target_bimodality_coefficient)
    print(paste("The optimal cofactor for marker", markers_of_interest[j], "is",
                cofactors[j]))
  }
  return(cofactors)
}

optimize_cofactor <- function(flowframe, transform_type, samplename, cofactors_list,
                              channels=markernames(flowframe), markernams){
  if(transform_type == "arcsinh"){
    cat("Arcsinh transform sample",samplename, "\n")
    cofactor_list <- list()
    opt_flowframe <- flowframe
    
    for(marker_pos in 1:length(channels)){
      for(cofactor in cofactors_list){
        test_flowframe <- transform_ff(flowframe,transform_type,samplename,channels, cofactor)
        if(LaplacesDemon::is.bimodal(test_flowframe@exprs[,channels[marker_pos]])){break}
      }
      cat("Transforming marker", markernams[marker_pos], "with cofactor", cofactor, "\n")
      opt_flowframe <- transform_ff(opt_flowframe,transform_type,samplename,channels[marker_pos], cofactor)
      cofactor_list[marker_pos] = cofactor
    }
  }else {
    opt_flowframe <- transform_ff(flowframe,transform_type,samplename,channels, cofactor)
  }
  return(opt_flowframe)
}
# opt_ff <- optimize_cofactor(flow_frames[[1]], transform_type,samplename, marker_indices_of_interest)
# opt_flowframe <- transform_ff(flowframe,transform_type,samplename,marker_indices_of_interest, cofactor)
# for(i in marker_indices_of_interest){print(is.bimodal(opt_ff@exprs[,i]))}

# ================================ Scaling =====================================
z_score_scaling <- function(flowframe, marker_indices){
  transform_list = transformList(featureNames(flowframe[,marker_indices]), scale)
  sc_flowframe <- flowCore::transform(flowframe,transform_list,
                                      truncate_max_range = FALSE) 
  return(sc_flowframe)
}

#=========================== Quantile normalization ============================
norm_quant <- function(matrix){
  qn_matrix <- preprocessCore::normalize.quantiles(as.matrix(matrix, 
                                                            nrow = length(markers_of_interest)))
  return(qn_matrix)
}
quantile_norm <- function(flowframe, marker_indices){
  transform_list = transformList(featureNames(flowframe[, marker_indices]), tfun=norm_quant)
  qn_flowframe <- flowCore::transform(flowframe,transform_list,
                                      truncate_max_range = FALSE) 
  return(qn_flowframe)
}

# cyto_plot(pp_flowframe[[12]], c("CD4"))
