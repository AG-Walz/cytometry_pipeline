### ========================= ANALYSIS FUNCTIONS ============================= ###
# Functions for dimensionality reduction, clustering, and statistical analysis

# PCA Analysis Functions
run_pca <- function(data_matrix, scale = TRUE, center = TRUE) {
  #' Perform Principal Component Analysis
  #'
  #' This function performs PCA on a data matrix
  #'
  #' @param data_matrix Numeric matrix or data frame with samples as rows
  #' @param scale Logical, whether to scale variables
  #' @param center Logical, whether to center variables
  #' @return PCA result object
  #'
  
  # Remove any non-numeric columns
  if (is.data.frame(data_matrix)) {
    numeric_cols <- sapply(data_matrix, is.numeric)
    data_matrix <- data_matrix[, numeric_cols, drop = FALSE]
  }
  
  # Perform PCA
  pca_result <- prcomp(data_matrix, scale. = scale, center = center)
  
  # Add variance explained information
  variance_explained <- summary(pca_result)$importance
  
  return(list(
    pca = pca_result,
    variance_explained = variance_explained,
    loadings = pca_result$rotation,
    scores = pca_result$x
  ))
}

# UMAP Analysis Functions
run_umap <- function(data_matrix, n_neighbors = 15, min_dist = 0.1, n_components = 2) {
  #' Perform UMAP Dimensionality Reduction
  #'
  #' This function performs UMAP on a data matrix
  #'
  #' @param data_matrix Numeric matrix with samples as rows
  #' @param n_neighbors Number of neighbors for UMAP
  #' @param min_dist Minimum distance parameter
  #' @param n_components Number of dimensions to reduce to
  #' @return UMAP result with coordinates
  #'
  
  if (!requireNamespace("umap", quietly = TRUE)) {
    stop("umap package is required for this function")
  }
  
  # Configure UMAP settings
  config <- umap::umap.defaults
  config$n_neighbors <- n_neighbors
  config$min_dist <- min_dist
  config$n_components <- n_components
  
  # Run UMAP
  umap_result <- umap::umap(data_matrix, config = config)
  
  return(list(
    layout = umap_result$layout,
    config = config,
    original_data = data_matrix
  ))
}

# Clustering Functions
leiden_clustering <- function(data_matrix, resolution = 0.5, n_neighbors = 10) {
  #' Perform Leiden Clustering
  #'
  #' This function performs Leiden clustering on a data matrix
  #'
  #' @param data_matrix Numeric matrix with samples as rows
  #' @param resolution Resolution parameter for clustering
  #' @param n_neighbors Number of neighbors for graph construction
  #' @return Cluster assignments
  #'
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required for this function")
  }
  
  # Create k-nearest neighbor graph
  dist_matrix <- dist(data_matrix)
  
  # For now, use hierarchical clustering as a placeholder
  # In a full implementation, this would use proper Leiden algorithm
  hclust_result <- hclust(dist_matrix, method = "ward.D2")
  
  # Cut tree to get clusters
  n_clusters <- max(2, round(nrow(data_matrix) / 50))
  clusters <- cutree(hclust_result, k = n_clusters)
  
  return(list(
    clusters = clusters,
    method = "hierarchical_ward",
    n_clusters = n_clusters
  ))
}

# Statistical Analysis Functions
calculate_statistics <- function(data_matrix, groups, test_type = "wilcox") {
  #' Calculate Statistical Comparisons
  #'
  #' This function calculates statistical tests between groups
  #'
  #' @param data_matrix Numeric matrix with features as columns
  #' @param groups Group assignments for each row
  #' @param test_type Type of statistical test ("wilcox", "t.test")
  #' @return Statistical test results
  #'
  
  if (length(groups) != nrow(data_matrix)) {
    stop("Number of group assignments must match number of rows in data matrix")
  }
  
  results <- list()
  unique_groups <- unique(groups)
  
  if (length(unique_groups) < 2) {
    warning("Need at least 2 groups for statistical comparison")
    return(NULL)
  }
  
  # Perform pairwise comparisons
  for (i in 1:(length(unique_groups) - 1)) {
    for (j in (i + 1):length(unique_groups)) {
      group1 <- unique_groups[i]
      group2 <- unique_groups[j]
      
      group1_data <- data_matrix[groups == group1, , drop = FALSE]
      group2_data <- data_matrix[groups == group2, , drop = FALSE]
      
      comparison_name <- paste(group1, "vs", group2)
      
      if (test_type == "wilcox") {
        # Wilcoxon rank sum test for each feature
        p_values <- apply(data_matrix, 2, function(x) {
          tryCatch({
            wilcox.test(x[groups == group1], x[groups == group2])$p.value
          }, error = function(e) NA)
        })
      } else if (test_type == "t.test") {
        # T-test for each feature
        p_values <- apply(data_matrix, 2, function(x) {
          tryCatch({
            t.test(x[groups == group1], x[groups == group2])$p.value
          }, error = function(e) NA)
        })
      }
      
      results[[comparison_name]] <- list(
        p_values = p_values,
        group1 = group1,
        group2 = group2,
        test_type = test_type
      )
    }
  }
  
  return(results)
}

# Utility analysis functions
get_top_features <- function(statistical_results, n_top = 10, p_threshold = 0.05) {
  #' Get Top Differentially Expressed Features
  #'
  #' Extract top features from statistical analysis results
  #'
  #' @param statistical_results Results from calculate_statistics
  #' @param n_top Number of top features to return
  #' @param p_threshold P-value threshold for significance
  #' @return List of top features per comparison
  #'
  
  top_features <- list()
  
  for (comparison in names(statistical_results)) {
    p_values <- statistical_results[[comparison]]$p_values
    
    # Remove NA values
    valid_pvals <- p_values[!is.na(p_values)]
    
    # Filter by significance threshold
    significant <- valid_pvals[valid_pvals < p_threshold]
    
    # Sort by p-value and take top features
    sorted_features <- sort(significant)
    top_n <- head(sorted_features, n_top)
    
    top_features[[comparison]] <- list(
      features = names(top_n),
      p_values = top_n,
      n_significant = length(significant)
    )
  }
  
  return(top_features)
}