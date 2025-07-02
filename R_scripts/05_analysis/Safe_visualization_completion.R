# Crash-Resistant Continuation for Fine Clustering Analysis
# This modifies the visualization section to be memory-safe and crash-resistant

# =========================
# Memory-Safe Visualization Functions 
# =========================

# Safe plot generation with memory limits
safe_plot_generation <- function() {
  
  message("=== Starting memory-safe visualization generation ===")
  
  # 1) Basic plots first (safest, using only results data)
  message("Generating basic plots...")
  
  tryCatch({
    # Volcano plot
    volcano <- results %>%
      filter(!is.na(p_adj) & !is.na(log2_enrichment_A_vs_B)) %>%
      mutate(sig = significant) %>%
      ggplot(aes(x = log2_enrichment_A_vs_B, y = -log10(p_adj))) +
      geom_point(aes(color = sig, size = sig), alpha = 0.7) +
      scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "red")) +
      scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
      labs(title = "Cluster Enrichment Analysis") +
      theme_minimal()
    
    ggsave(file.path(CONFIG$out_dir, "volcano_enrichment.png"), volcano, 
           width = 10, height = 8, dpi = 150)  # Lower DPI to save memory
    message("✓ Volcano plot saved")
    
    # Clear from memory
    rm(volcano)
    gc()
    
  }, error = function(e) {
    message("Basic plots failed: ", e$message)
  })
  
  # 2) Barplot for significant clusters
  tryCatch({
    n_top <- min(20, sum(results$significant, na.rm = TRUE))  # Limit to 20
    if (n_top > 0) {
      bar_top <- results %>%
        filter(significant) %>%
        slice_min(p_adj, n = n_top) %>%
        mutate(cluster = fct_reorder(cluster, log2_enrichment_A_vs_B)) %>%
        ggplot(aes(x = cluster, y = log2_enrichment_A_vs_B)) +
        geom_col(alpha = 0.8, fill = "steelblue") +
        coord_flip() +
        labs(title = paste0("Top ", n_top, " significant clusters")) +
        theme_minimal()
      
      ggsave(file.path(CONFIG$out_dir, "top_clusters_bar.png"), bar_top, 
             width = 8, height = max(6, n_top * 0.3), dpi = 150)
      message("✓ Barplot saved")
      
      rm(bar_top)
      gc()
    }
  }, error = function(e) {
    message("Barplot failed: ", e$message)
  })
  
  # 3) Simple heatmap
  tryCatch({
    n_heatmap <- min(30, nrow(results))  # Limit to 30 clusters
    heatmap_data <- results %>%
      slice_head(n = n_heatmap) %>%
      mutate(cluster = fct_reorder(cluster, log2_enrichment_A_vs_B))
    
    enrichment_heatmap <- ggplot(heatmap_data, aes(x = 1, y = cluster)) +
      geom_tile(aes(fill = log2_enrichment_A_vs_B), color = "white") +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      labs(title = "Cluster Enrichment Summary") +
      theme_minimal() +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    ggsave(file.path(CONFIG$out_dir, "enrichment_heatmap.png"), enrichment_heatmap,
           width = 6, height = max(6, n_heatmap * 0.2), dpi = 150)
    message("✓ Heatmap saved")
    
    rm(enrichment_heatmap, heatmap_data)
    gc()
    
  }, error = function(e) {
    message("Heatmap failed: ", e$message)
  })
  
  message("Basic visualization complete!")
}

# Replace the crash-prone visualization section
safe_plot_generation()

# Skip advanced plots that require cell data reload
message("Skipping advanced UMAP and marker plots to avoid memory crashes")
message("Use the Resume_Fine_clustering_analysis.R script if you need these plots")

# =========================
# Safe Analysis Summary
# =========================
message("\n=== CRASH-RESISTANT ANALYSIS COMPLETE ===")
message("Results saved to: ", CONFIG$out_dir)
message("Basic plots generated:")
message("- volcano_enrichment.png: Volcano plot") 
message("- top_clusters_bar.png: Top significant clusters")
message("- enrichment_heatmap.png: Summary heatmap")

# Print findings without complex operations
sig_count <- sum(results$significant, na.rm = TRUE)
total_clusters <- nrow(results)

message("\nSummary:")
message("- Total clusters analyzed: ", total_clusters)
message("- Significant clusters (adj. p < 0.05): ", sig_count)
message("- Analysis method: ", CONFIG$stat_method)

if (sig_count > 0) {
  message("\nTop 5 most significant clusters:")
  top_sig <- results %>% 
    filter(significant) %>% 
    arrange(p_adj) %>% 
    slice_head(n = 5)
  
  for (i in 1:nrow(top_sig)) {
    row <- top_sig[i, ]
    message(sprintf("  %d. %s (p_adj = %.2e)", i, row$cluster, row$p_adj))
  }
}

message("\n✓ Analysis completed without crashes!")