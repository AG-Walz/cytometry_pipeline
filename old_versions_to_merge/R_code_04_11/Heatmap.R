
df <-read.csv("df.csv")
leiden_clusters <- read.csv("./clustered_harm.csv")$Leiden.cluster

expression_matrix <- do.call(rbind, df)

cluster_labels <- data.frame(Cluster = factor(leiden_clusters))

pheatmap(expression_matrix, 
         cluster_rows = FALSE,  # Cluster rows
         annotation_row = data.frame(Cluster = factor(leiden_clusters)),  # Annotation for rows
         cluster_cols = FALSE,
         main = "Expression Heatmap")  # Main title


