#=============================== BATCH CORRECT ==============================
# This script is used to perform batch correction on the PCA embeddings using the Harmony algorithm.

# Load the required libraries
install.packages("harmony")
library(harmony)
library(stats)

# Load the expression matrix to batch correct
df <- read.csv("df.csv")

# Import the UMAP embedding with batches, and other labels
umapdf_with_labels<- read.csv("umapdf_with_labels.csv")
condition_procedence <- umapdf_with_labels$condition_procedence
sample_procedence <- umapdf_with_labels$sample_procedence
umap_harm <-read.csv("umap_harm.csv")

# Translate sample names to batches using the dictionary
translated_batches <- sample_batch_dict[sample_procedence]
umapdf_with_labels <- cbind(umapdf_with_labels, translated_batches)

my_harmony_embeddings <- RunHarmony(data_mat=pca_result$x[,1:kaiser_threshold], 
                                    meta_dat = umapdf_with_labels$translated_batches,
                                    theta=2, sigma = 0.05,
                                    nclust=14, max.iter= 2000, early_stop= FALSE)
write.csv(my_harmony_embeddings, "harm.csv", row.names = FALSE)

umap_harm <- uwot::umap(X=my_harmony_embeddings,n_neighbors = 50, min_dist = 0.05)
write.csv(umap_harm, "umap_harm.csv", row.names = FALSE)

colnames(umap_harm) <- c("UMAP1", "UMAP2")
umapharm_with_labels <- cbind(umap_harm, sample_procedence, condition_procedence)
write.csv(umapdf_with_labels, "umapharm_with_labels.csv", row.names = FALSE)

plot <- ggplot(data.frame(umapharm_with_labels), 
               aes(as.numeric(UMAP1), y= as.numeric(UMAP2), 
                   color = sample_procedence)) +
  geom_point(shape = 19, size = 0.1, stroke = 0.1) +
  labs(x = "UMAP1", y = "UMAP2",color = "Treatment") +
  #scale_color_manual(values = c("deepskyblue", "brown1"))+
  guides(colour = guide_legend(override.aes = list(size=5),  ncol = 4))+
  theme_bw()

# Save the plot as a PNG file
ggsave("harm_UMAP.png", plot,width = 6, height = 2.9, units = "in")

plot <- ggplot(data.frame(umapharm_with_labels), 
               aes(as.numeric(UMAP1), y= as.numeric(UMAP2), 
                   color = condition_procedence)) +
  geom_point(shape = 19, size = 0.1, stroke = 0.1) +
  labs(x = "UMAP1", y = "UMAP2",color = "Treatment") +
  #scale_color_manual(values = c("deepskyblue", "brown1"))+
  guides(colour = guide_legend(override.aes = list(size=5),  ncol = 1))+
  theme_bw()

# Save the plot as a PNG file
ggsave("harm_condition.png", plot,width = 4.4, height = 2.9, units = "in")

