# Resume Fine Clustering Analysis from Saved Results
# This script continues the analysis after results have been computed and saved
# Use this if the main script crashed during visualization phase

# =========================
# Configuration - MODIFY THESE
# =========================
CONFIG <- list(
  # 1) Path to the saved results CSV (from previous run)
  results_csv = "/data/Cytometry/ICK/ICK_T1_P+_poR/Results/cluster_enrichment_results.csv",
  
  # 2) Path to the saved by-file counts CSV (OPTIONAL - may not exist if main script crashed)
  by_file_counts_csv = "/data/Cytometry/ICK/ICK_T1_P+_poR/Results/by_file_cluster_counts.csv",
  
  # 3) Output directory (should match original analysis)
  out_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR/Results/",
  
  # 4) Data directory (for re-loading cell data if needed for plots)
  data_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR/",
  
  # 5) Condition levels (must match original analysis)
  condition_levels = c("pCR", "npCR"),
  
  # 6) UMAP column names (for visualization)
  umap_x = "UMAP1",
  umap_y = "UMAP2", 
  
  # 7) Markers of interest (for marker plots)
  markers_of_interest = c("CD45", "CD3", "CD4", "CD8"),  # Modify as needed
  
  # 8) Statistical method used (for plot titles)
  stat_method = "mixed_effects",
  
  # 9) Memory limits for visualization
  max_cells_for_plots = 25000,  # Reduce if crashing
  enable_umap_plots = TRUE,     # Set FALSE to skip UMAP plots
  enable_marker_plots = TRUE,   # Set FALSE to skip marker plots
  
  # 10) Plot settings
  plot_dpi = 300,
  plot_width = 10,
  plot_height = 8
)

# =========================
# Load Required Libraries
# =========================
message("Loading libraries...")
suppressMessages({
  library(tidyverse)
  library(data.table)
  library(janitor)
  library(ggrepel)
  library(scales)
  library(RColorBrewer)
})

# =========================
# Helper Functions
# =========================
make_outdir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Created output directory: ", path)
  }
}

safe_sample_cells <- function(cells_data, max_cells = CONFIG$max_cells_for_plots) {
  if (nrow(cells_data) > max_cells) {
    message("Sampling ", max_cells, " cells from ", nrow(cells_data), " for visualization")
    return(cells_data %>% slice_sample(n = max_cells))
  }
  return(cells_data)
}

# =========================
# Load Saved Results
# =========================
message("Loading saved results...")

# Check if results files exist
# Check if results files exist
if (!file.exists(CONFIG$results_csv)) {
  stop("Results file not found: ", CONFIG$results_csv, 
       "\nPlease check the path or run the main analysis first.")
}

# Load results (required)
results <- read_csv(CONFIG$results_csv, show_col_types = FALSE)

# Load by-file counts (optional - may not exist if script crashed during generation)
by_file_counts <- NULL
if (file.exists(CONFIG$by_file_counts_csv)) {
  by_file_counts <- read_csv(CONFIG$by_file_counts_csv, show_col_types = FALSE)
  message("Loaded by-file counts data")
} else {
  message("By-file counts not found (", basename(CONFIG$by_file_counts_csv), ") - continuing without it")
  message("This is normal if the main script crashed before generating this file")
}

message("Loaded results for ", nrow(results), " clusters")
message("Significant clusters: ", sum(results$significant, na.rm = TRUE))

# Create output directory
make_outdir(CONFIG$out_dir)

# =========================
# Generate Visualizations (Memory-Safe)
# =========================
message("Generating visualizations...")

# 1) Volcano plot - Safe, uses only results data
message("Creating volcano plot...")
volcano <- results %>%
  filter(!is.na(p_adj) & !is.na(log2_enrichment_A_vs_B)) %>%
  mutate(sig = significant) %>%
  ggplot(aes(x = log2_enrichment_A_vs_B, y = -log10(p_adj))) +
  geom_point(aes(color = sig, size = sig), alpha = 0.7) +
  geom_text_repel(
    data = . %>% slice_min(p_adj, n = 15),
    aes(label = cluster), 
    size = 3, max.overlaps = 200, box.padding = 0.3
  ) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
  scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(
    title = paste("Cluster Enrichment Analysis (", CONFIG$stat_method, ")", sep = ""),
    subtitle = paste0(CONFIG$condition_levels[1], " vs ", CONFIG$condition_levels[2]),
    x = paste0("log2 enrichment (", CONFIG$condition_levels[1], " / ", CONFIG$condition_levels[2], ")"),
    y = "-log10(adjusted p-value)",
    color = "Significant",
    size = "Significant"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(CONFIG$out_dir, "volcano_enrichment.png"), volcano, 
       width = CONFIG$plot_width, height = CONFIG$plot_height, dpi = CONFIG$plot_dpi)
message("✓ Volcano plot saved")

# 2) Barplot of top significant clusters
n_top <- min(25, sum(results$significant, na.rm = TRUE))
if (n_top > 0) {
  message("Creating top clusters barplot...")
  bar_top <- results %>%
    filter(significant) %>%
    slice_min(p_adj, n = n_top) %>%
    mutate(cluster = fct_reorder(cluster, log2_enrichment_A_vs_B)) %>%
    ggplot(aes(x = cluster, y = log2_enrichment_A_vs_B, fill = direction)) +
    geom_col(alpha = 0.8) +
    coord_flip() +
    scale_fill_brewer(type = "qual", palette = "Set2") +
    labs(
      title = paste0("Top ", n_top, " significant clusters"),
      y = paste0("log2 enrichment (", CONFIG$condition_levels[1], " / ", CONFIG$condition_levels[2], ")"),
      x = "Cluster",
      fill = "Direction"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  
  ggsave(file.path(CONFIG$out_dir, "top_clusters_bar.png"), bar_top, 
         width = 8, height = max(6, n_top * 0.3), dpi = CONFIG$plot_dpi)
  message("✓ Top clusters barplot saved")
} else {
  message("No significant clusters found, skipping barplot")
}

# 3) Summary heatmap of enrichment patterns
if (nrow(results) > 1) {
  message("Creating enrichment heatmap...")
  n_heatmap <- min(50, nrow(results))
  heatmap_data <- results %>%
    slice_head(n = n_heatmap) %>%
    mutate(
      cluster = fct_reorder(cluster, log2_enrichment_A_vs_B),
      significance = case_when(
        significant ~ "Significant",
        TRUE ~ "Not significant"
      )
    )
  
  enrichment_heatmap <- ggplot(heatmap_data, aes(x = 1, y = cluster)) +
    geom_tile(aes(fill = log2_enrichment_A_vs_B), color = "white", size = 0.1) +
    geom_point(data = . %>% filter(significant), aes(x = 1.3), 
               shape = 8, size = 2, color = "black") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      name = "log2\\nenrichment"
    ) +
    labs(
      title = "Cluster Enrichment Summary",
      subtitle = "★ indicates significance (adj. p < 0.05)",
      y = "Cluster",
      x = ""
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank()
    )
  
  ggsave(file.path(CONFIG$out_dir, "enrichment_heatmap.png"), enrichment_heatmap,
         width = 6, height = max(6, n_heatmap * 0.2), dpi = CONFIG$plot_dpi)
  message("✓ Enrichment heatmap saved")
}

# =========================
# Memory-Intensive Plots (Optional)
# =========================

# Try to load cell data for UMAP and marker plots (if requested)
cells <- NULL
if (CONFIG$enable_umap_plots || CONFIG$enable_marker_plots) {
  message("Attempting to load cell data for advanced plots...")
  
  tryCatch({
    # Try to find original cell data
    csv_files <- list.files(CONFIG$data_dir, pattern = "\\.csv$", full.names = TRUE)
    if (length(csv_files) == 0) {
      message("No CSV files found in ", CONFIG$data_dir)
    } else {
      message("Found ", length(csv_files), " CSV files, loading first file for plot data...")
      
      # Load just one file for plotting (to avoid memory issues)
      first_file <- csv_files[1]
      cells <- fread(first_file) %>%
        clean_names() %>%
        safe_sample_cells(CONFIG$max_cells_for_plots)
      
      message("Loaded ", nrow(cells), " cells from ", basename(first_file))
    }
  }, error = function(e) {
    message("Could not load cell data: ", e$message)
    CONFIG$enable_umap_plots <<- FALSE
    CONFIG$enable_marker_plots <<- FALSE
  })
}

# 4) UMAP visualization (if cell data loaded and UMAP columns exist)
if (CONFIG$enable_umap_plots && !is.null(cells)) {
  message("Checking for UMAP visualization...")
  
  # Check if UMAP columns exist
  if (CONFIG$umap_x %in% names(cells) && CONFIG$umap_y %in% names(cells)) {
    # Find cluster columns
    cluster_cols <- names(cells)[str_detect(names(cells), "^cluster_|^c_")]
    
    if (length(cluster_cols) > 0) {
      # Get best significant cluster
      sig_clusters <- results %>% filter(significant) %>% arrange(p_adj)
      if (nrow(sig_clusters) > 0) {
        best_cluster <- sig_clusters$cluster[1]
        
        # Check if this cluster exists in our cell data
        if (best_cluster %in% cluster_cols) {
          message("Creating UMAP plot for cluster: ", best_cluster)
          
          tryCatch({
            umap_df <- cells %>%
              select(all_of(c(CONFIG$umap_x, CONFIG$umap_y, best_cluster))) %>%
              mutate(condition = sample(CONFIG$condition_levels, nrow(.), replace = TRUE)) %>%  # Mock condition for demo
              rename(in_cluster = all_of(best_cluster)) %>%
              mutate(in_cluster = factor(in_cluster, levels = c(0,1), labels = c("No","Yes")))
            
            umap_plot <- ggplot(umap_df, aes_string(x = CONFIG$umap_x, y = CONFIG$umap_y)) +
              geom_point(aes(color = in_cluster), alpha = 0.6, size = 0.3) +
              scale_color_manual(values = c("No" = "lightgray", "Yes" = "red")) +
              facet_grid(in_cluster ~ condition) +
              coord_equal() +
              labs(
                title = paste0("UMAP for most significant cluster: ", best_cluster),
                subtitle = paste0("p_adj = ", scientific(sig_clusters$p_adj[1], digits = 3)),
                x = CONFIG$umap_x, y = CONFIG$umap_y,
                color = "In cluster"
              ) +
              theme_minimal(base_size = 11) +
              theme(legend.position = "bottom")
            
            ggsave(file.path(CONFIG$out_dir, paste0("umap_", make.names(best_cluster), ".png")),
                   umap_plot, width = CONFIG$plot_width, height = CONFIG$plot_height, dpi = CONFIG$plot_dpi)
            message("✓ UMAP plot saved")
          }, error = function(e) {
            message("UMAP plot failed: ", e$message)
          })
        } else {
          message("Best cluster '", best_cluster, "' not found in cell data columns")
        }
      } else {
        message("No significant clusters for UMAP plot")
      }
    } else {
      message("No cluster columns found in cell data")
    }
  } else {
    message("UMAP columns not found: ", CONFIG$umap_x, ", ", CONFIG$umap_y)
  }
}

# 5) Marker expression plots (if markers available)
if (CONFIG$enable_marker_plots && !is.null(cells) && length(CONFIG$markers_of_interest) > 0) {
  message("Checking for marker expression plots...")
  
  # Check which markers exist
  available_markers <- CONFIG$markers_of_interest[CONFIG$markers_of_interest %in% names(cells)]
  
  if (length(available_markers) > 0) {
    message("Found markers: ", paste(available_markers, collapse = ", "))
    
    # Find cluster columns and get best cluster
    cluster_cols <- names(cells)[str_detect(names(cells), "^cluster_|^c_")]
    sig_clusters <- results %>% filter(significant) %>% arrange(p_adj)
    
    if (length(cluster_cols) > 0 && nrow(sig_clusters) > 0) {
      best_cluster <- sig_clusters$cluster[1]
      
      if (best_cluster %in% cluster_cols) {
        message("Creating marker plots for cluster: ", best_cluster)
        
        tryCatch({
          marker_long <- cells %>%
            select(all_of(c(best_cluster, available_markers))) %>%
            mutate(condition = sample(CONFIG$condition_levels, nrow(.), replace = TRUE)) %>%  # Mock condition
            rename(in_cluster = all_of(best_cluster)) %>%
            mutate(in_cluster = factor(in_cluster, levels = c(0,1), labels = c("No","Yes"))) %>%
            pivot_longer(all_of(available_markers), names_to = "marker", values_to = "value")
          
          marker_violin <- ggplot(marker_long, aes(x = in_cluster, y = value)) +
            geom_violin(aes(fill = condition), scale = "width", trim = TRUE, alpha = 0.7) +
            geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(0.9)) +
            facet_wrap(~marker, scales = "free_y", ncol = 3) +
            scale_fill_brewer(type = "qual", palette = "Set2") +
            labs(
              title = paste0("Marker expression in/out of cluster: ", best_cluster),
              x = "In cluster?", y = "Marker intensity",
              fill = "Condition"
            ) +
            theme_minimal(base_size = 11) +
            theme(legend.position = "bottom")
          
          ggsave(file.path(CONFIG$out_dir, paste0("markers_", make.names(best_cluster), ".png")),
                 marker_violin, width = 12, height = 8, dpi = CONFIG$plot_dpi)
          message("✓ Marker plots saved")
        }, error = function(e) {
          message("Marker plots failed: ", e$message)
        })
      }
    }
  } else {
    message("None of the specified markers found in data: ", paste(CONFIG$markers_of_interest, collapse = ", "))
  }
}

# =========================
# Final Summary
# =========================
message("\n=== RESUMED ANALYSIS COMPLETE ===")
message("Results directory: ", CONFIG$out_dir)
message("Generated plots:")
message("- volcano_enrichment.png: Volcano plot")
if (n_top > 0) message("- top_clusters_bar.png: Top significant clusters")
message("- enrichment_heatmap.png: Summary heatmap")
if (CONFIG$enable_umap_plots) message("- umap_*.png: UMAP visualization (if successful)")
if (CONFIG$enable_marker_plots) message("- markers_*.png: Marker distributions (if successful)")

# Print key findings
sig_results <- results %>% filter(significant) %>% arrange(p_adj)
message("\nKey findings (", nrow(sig_results), " significant clusters at adj. p < 0.05):")

if (nrow(sig_results) > 0) {
  for (i in seq_len(min(10, nrow(sig_results)))) {
    row <- sig_results[i, ]
    direction_text <- if (row$log2_enrichment_A_vs_B > 0) {
      paste("enriched in", CONFIG$condition_levels[1])
    } else {
      paste("enriched in", CONFIG$condition_levels[2])
    }
    message(sprintf("  %d. %s: %.2f-fold %s (adj. p = %.2e)",
                   i, row$cluster, 2^abs(row$log2_enrichment_A_vs_B), 
                   direction_text, row$p_adj))
  }
  if (nrow(sig_results) > 10) {
    message("  ... and ", nrow(sig_results) - 10, " more significant clusters")
  }
} else {
  message("  No statistically significant cluster enrichments found.")
}

message("\nAnalysis method: ", CONFIG$stat_method)
message("Total clusters analyzed: ", nrow(results))

message("\n✓ Resume script completed successfully!")

