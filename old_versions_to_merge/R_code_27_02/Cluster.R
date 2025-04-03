library(dbscan)

# Define columns of interest
cols_of_interest <- colnames(cc_ff)[9:41]
data <- cc_ff[, cols_of_interest]

# Perform HDBSCAN clustering
clustering_result <- hdbscan(data, minPts = 10)  # Adjust minPts as needed

# Extract cluster origin numbers
cluster_origins <- clustering_result$cluster

# Return the cluster origins
cluster_origins

# =============== Kmeans =================
# Define columns of interest
# TODOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
OOOOOOOOOOOOOOOOOOOOOOOOOOOO

install.packages("ClusterR")
library(ClusterR)

mrks_of_interest <- c("CD3", "CD4", "CD8")
dyes_of_interest <- translate_markers_to_dyes(mrks_of_interest, flow_frames[[1]])
data <- cc_ff[, dyes_of_interest]

# Perform k-means clustering
set.seed(123)  # Set seed for reproducibility
k <- 25  # Adjust k to the desired number of clusters
clustering_result <- kmeans(data, centers = k, nstart = 10, iter.max = 1000000)  # nstart for better initial centroids

# Extract cluster origin numbers
cluster_origins <- clustering_result$cluster

# Return the cluster origins
length(cluster_origins)
nrow(cc_ff)
# Normalize the data by column
normalized_data <- scale(data)

# Calculate mean expression for each marker by cluster
mean_expression <- aggregate(normalized_data, by = list(cluster_origins), FUN = mean)
rownames(mean_expression) <- mean_expression$Group.1
mean_expression <- mean_expression[, -1]

# Normalize by column
normalized_expression <- t(apply(mean_expression, 1, scale))

# Plot heatmap
# Perform hierarchical clustering on the normalized expression data
dendrogram <- hclust(dist(normalized_expression))

# Plot heatmap with dendrogram
pheatmap(normalized_expression, cluster_rows = dendrogram, cluster_cols = FALSE,
         main = "Heatmap of Scaled Mean Expression by Cluster",
         display_numbers = FALSE,
         labels_col = mrks_of_interest)
