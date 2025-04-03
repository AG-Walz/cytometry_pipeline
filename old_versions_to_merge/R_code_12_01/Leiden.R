#install.packages("leiden")
library(leiden)

# Replace 'dist_method' and 'dist_metric' with the desired distance method and metric
leiden_clusters <- read.csv("./clusters_harmony.csv")$Leiden.cluster
final_plots <- "./Salt_24_02_26/final_plots/"

umapharm_with_labels <- cbind(umap_harm, sample_procedence, condition_procedence)
umapharm_with_labels <- cbind(umapharm_with_labels, as.integer(leiden_clusters))

# Define a vector of 14 colors
cluster_colors <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00","blue1", # orange
  "gold1", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "maroon", "orchid1", "deeppink1",  
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "black"
)

plot <- ggplot(data.frame(umapharm_with_labels), 
               aes(as.numeric(UMAP1), y= as.numeric(UMAP2), 
                   color = factor(leiden_clusters))) +
  geom_point(shape = 19, size = 0.1, stroke = 0.1) +
  labs(x = "UMAP1", y = "UMAP2",color = "Leiden cluster") +
  scale_color_manual(values = cluster_colors) +
  guides(colour = guide_legend(override.aes = list(size=5),  nrow=5))+
  theme_bw()

# Save the plot as a PNG file
ggsave(paste0(final_plots,"leiden_harm.png"), plot,width = 4.55, height = 2.9, units = "in")


# Create a stacked barplot
plot <- ggplot(data.frame(umapharm_with_labels), aes(x = factor(leiden_clusters), fill = condition_procedence)) +
  geom_bar(position = "stack", color = "black") +
  labs(x = "Leiden clusters", y = "Number of cells", fill = "Treatment") +
  theme_bw()
# Save the plot as a PNG file
ggsave(paste0(final_plots,"barplot_treatment_harm.png"), plot,width = 6, height = 2.9, units = "in")

plot <- ggplot(data.frame(umapharm_with_labels), aes(x = factor(leiden_clusters), fill = translated_batches)) +
  geom_bar(position = "stack", color = "black") +
  labs(x = "Leiden clusters", y = "Number of cells", fill = "Batches") +
  theme_bw()
# Save the plot as a PNG file
ggsave(paste0(final_plots,"barplot_batches_harm.png"), plot,width = 6, height = 2.9, units = "in")

plot <- ggplot(data.frame(umapharm_with_labels), aes(x = factor(leiden_clusters), fill = sample_procedence)) +
  geom_bar(position = "stack") +
  labs(x = "Leiden clusters", y = "Number of cells", fill = "Sample") +
  theme_bw()
# Save the plot as a PNG file
ggsave(paste0(final_plots,"barplot_sample_harm.png"), plot+ guides(fill = guide_legend(ncol = 4)),width = 6, height = 2.9, units = "in")


# Create a stacked barplot
plot <- ggplot(data.frame(umapharm_with_labels), aes(x = factor(leiden_clusters), fill = condition_procedence)) +
  geom_bar(position = "fill", color = "black") +
  labs(x = "Leiden clusters", y = "Number of cells", fill = "Treatment") +
  theme_bw()
# Save the plot as a PNG file
ggsave(paste0(final_plots,"barplot_treatment_harm_fill.png"), plot,width = 6, height = 2.9, units = "in")

plot <- ggplot(data.frame(umapharm_with_labels), aes(x = factor(leiden_clusters), fill = translated_batches)) +
  geom_bar(position = "fill", color = "black") +
  labs(x = "Leiden clusters", y = "Number of cells", fill = "Batches") +
  theme_bw()
# Save the plot as a PNG file
ggsave(paste0(final_plots,"barplot_batches_harm_fill.png"), plot,width = 6, height = 2.9, units = "in")

plot <- ggplot(data.frame(umapharm_with_labels), aes(x = factor(leiden_clusters), fill = sample_procedence)) +
  geom_bar(position = "fill") +
  labs(x = "Leiden clusters", y = "Number of cells", fill = "Sample") +
  theme_bw()
# Save the plot as a PNG file
ggsave(paste0(final_plots,"barplot_sample_harm_fill.png"), plot+ guides(fill = guide_legend(ncol = 4)),width = 6, height = 2.9, units = "in")
