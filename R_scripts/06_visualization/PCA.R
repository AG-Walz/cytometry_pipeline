library(stats)

df <- read.csv("df.csv")
umapdf_with_labels<- read.csv("umapdf_with_labels.csv")
condition_procedence <- umapdf_with_labels$condition_procedence
sample_procedence <- umapdf_with_labels$sample_procedence

pca_result <- prcomp(df, scale. = FALSE, center = FALSE)

# Calculate proportion of variance explained
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Calculate Kaiser's rule threshold
kaiser_threshold <- sum(pca_result$sdev^2 > 1)

# Create dataframe for plotting
scree_data <- data.frame(
  PC = 1:length(variance_explained),
  Proportion_of_Variance_Explained = variance_explained
)

# Plot scree plot with Kaiser's rule threshold
ggplot(scree_data, aes(x = PC, y = Proportion_of_Variance_Explained)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = kaiser_threshold, linetype = "dashed", color = "blue") +
  labs(x = "Principal Component", y = "Proportion of Variance Explained") 
theme_minimal()
