#!/usr/bin/env Rscript
# =============================================================================
# ICK Cluster Abundance Analysis
# =============================================================================
# Loads preprocessed CSVs with binary clustering columns, builds a sample-by-
# cluster tracking matrix (parent-level populations only), and generates
# abundance barplots.
#
# Usage:
#   cd pipelines
#   Rscript run_ICK_abundance.R
# =============================================================================

cat("=== ICK Cluster Abundance Analysis ===\n")
start_time <- Sys.time()

# --- Load config ---
base_dir <- dirname(getwd())
source(file.path(base_dir, "config", "study_params", "ICK_abundance.R"))

# --- Source shared utilities ---
shared_dir <- file.path(base_dir, "src", "R_shared")
source(file.path(shared_dir, "io_utils.R"))
source(file.path(shared_dir, "cluster_utils.R"))
source(file.path(shared_dir, "plot_utils.R"))

# --- Step 1: Build tracking matrix ---
cat("\n--- Step 1: Building tracking matrix ---\n")

require(data.table)

make_outdir(output_dir)

csv_info <- list_csvs(input_dir)
cat("Found", nrow(csv_info), "CSV files in", input_dir, "\n")

# Read first file to discover columns
first_header <- names(fread(csv_info$file_path[1], nrows = 0))
all_cluster_cols <- grep(cluster_regex, first_header, value = TRUE, perl = TRUE)

# Apply exclusion filters to keep only parent-level clusters
parent_cols <- all_cluster_cols
for (pat in exclude_patterns) {
  parent_cols <- parent_cols[!grepl(pat, parent_cols)]
}
cat("Total clustering columns:", length(all_cluster_cols), "\n")
cat("Parent-level clusters:  ", length(parent_cols), "\n")
cat("Clusters to track:\n")
for (col in parent_cols) cat("  -", clean_cluster_name(col), "\n")

# Process each CSV and count cells per parent cluster
tracking_list <- vector("list", nrow(csv_info))
n_files <- nrow(csv_info)
batch_time <- Sys.time()

for (i in seq_len(n_files)) {
  dt <- fread(csv_info$file_path[i], select = parent_cols, showProgress = FALSE)
  counts <- colSums(dt, na.rm = TRUE)
  names(counts) <- sapply(names(counts), clean_cluster_name)
  tracking_list[[i]] <- c(file = csv_info$file_name[i], counts)

  # Log every 10 files or on the last file
  if (i %% 10 == 0 || i == n_files) {
    elapsed_s <- as.numeric(difftime(Sys.time(), batch_time, units = "secs"))
    rate <- i / elapsed_s
    eta <- (n_files - i) / rate
    cat(sprintf("  [%d/%d] %s | %.1f files/sec | ETA: %.0fs\n",
                i, n_files, csv_info$file_name[i], rate, eta))
    flush(stdout())
  }
}

tracking_matrix <- as.data.frame(do.call(rbind, tracking_list), stringsAsFactors = FALSE)
# First column is file name; remaining are counts
for (col in names(tracking_matrix)[-1]) {
  tracking_matrix[[col]] <- as.integer(tracking_matrix[[col]])
}

tracking_path <- file.path(output_dir, "cluster_tracking_matrix.csv")
fwrite(tracking_matrix, tracking_path, sep = ";")
cat("\nTracking matrix saved:", tracking_path, "\n")
cat("  ", nrow(tracking_matrix), "samples x", ncol(tracking_matrix) - 1, "clusters\n")

# --- Step 2: Abundance plots ---
cat("\n--- Step 2: Generating abundance plots ---\n")

require(ggplot2)
require(tidyr)

# Pivot to long format: sample | cluster | count
long_df <- tidyr::pivot_longer(tracking_matrix, cols = -file,
                                names_to = "cluster", values_to = "count")
long_df <- as.data.frame(long_df)

# 2a. Stacked barplot (samples on x, clusters as fill)
cat("  Plotting: stacked barplot\n")
p_stacked <- ggplot(long_df, aes(x = file, y = count, fill = cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Cluster Abundance per Sample (Stacked)",
       x = NULL, y = "Cell Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm")) +
  scale_fill_viridis_d(option = "turbo")
ggsave(file.path(output_dir, "abundance_stacked.png"),
       p_stacked, width = 24, height = 10, dpi = 150)

# 2b. Proportional barplot (samples on x, clusters as fill, position = fill)
cat("  Plotting: proportional barplot\n")
p_prop <- ggplot(long_df, aes(x = file, y = count, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(title = "Cluster Proportions per Sample",
       x = NULL, y = "Proportion") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm")) +
  scale_fill_viridis_d(option = "turbo")
ggsave(file.path(output_dir, "abundance_proportional.png"),
       p_prop, width = 24, height = 10, dpi = 150)

# 2c. Aggregate: total cells per cluster across all samples
cat("  Plotting: total cells per cluster\n")
cluster_totals <- aggregate(count ~ cluster, data = long_df, sum)
cluster_totals <- cluster_totals[order(cluster_totals$count, decreasing = TRUE), ]

p_cluster <- ggplot(cluster_totals,
                     aes(x = reorder(cluster, -count), y = count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Total Cells per Cluster (All Samples)",
       x = NULL, y = "Cell Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  scale_y_continuous(labels = scales::comma)
ggsave(file.path(output_dir, "total_cells_per_cluster.png"),
       p_cluster, width = 14, height = 8, dpi = 150)

# --- Step 3: Per-sample cell count summary plot ---
cat("\n--- Step 3: Per-sample total cell count ---\n")

require(ggplot2)

total_counts <- data.frame(
  sample = tracking_matrix$file,
  total  = rowSums(tracking_matrix[, -1])
)
total_counts <- total_counts[order(total_counts$total, decreasing = TRUE), ]

p_total <- ggplot(total_counts, aes(x = reorder(sample, -total), y = total)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Total Cells per Sample", x = NULL, y = "Cell Count") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))
ggsave(file.path(output_dir, "total_cells_per_sample.png"),
       p_total, width = 20, height = 8, dpi = 150)

# --- Step 4: Heatmap of cluster proportions per sample ---
cat("\n--- Step 4: Cluster proportion heatmap ---\n")

mat <- as.matrix(tracking_matrix[, -1])
rownames(mat) <- tracking_matrix$file

# Normalize each row to proportions
prop_mat <- mat / rowSums(mat)
prop_mat[is.nan(prop_mat)] <- 0

# Order by hierarchical clustering
if (nrow(prop_mat) > 2 && ncol(prop_mat) > 1) {
  row_order <- hclust(dist(prop_mat))$order
  col_order <- hclust(dist(t(prop_mat)))$order
  prop_mat <- prop_mat[row_order, col_order]
}

long_hm <- reshape2::melt(prop_mat)
names(long_hm) <- c("sample", "cluster", "proportion")

p_hm <- ggplot(long_hm, aes(x = cluster, y = sample, fill = proportion)) +
  geom_tile() +
  scale_fill_viridis_c(option = "inferno", name = "Proportion") +
  theme_minimal() +
  labs(title = "Cluster Proportion per Sample", x = "Cluster", y = "Sample") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 3)
  )
ggsave(file.path(output_dir, "cluster_proportion_heatmap.png"),
       p_hm, width = 16, height = 24, dpi = 150)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat(sprintf("\n=== ICK analysis complete (%.1f minutes) ===\n", elapsed))
cat("Output directory:", output_dir, "\n")
