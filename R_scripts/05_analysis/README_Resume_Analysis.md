# Instructions for Resuming Fine Clustering Analysis

## Problem
Your main Fine_clustering_cytolution.R script ran successfully until line 947 (saved results) but crashed during the visualization phase. This happens because:

1. **Memory issues**: UMAP and marker plots require reloading large cell datasets
2. **Complex ggplot operations**: Multiple large plots generated simultaneously
3. **Data manipulation**: Pivoting and sampling large datasets for visualization

## Solutions

### Option 1: Use the Resume Script (Recommended)
Use `Resume_Fine_clustering_analysis.R` to generate visualizations from saved results:

1. **Modify the CONFIG section** in `Resume_Fine_clustering_analysis.R`:
   ```r
   CONFIG <- list(
     # Update these paths to match your setup
     results_csv = "cluster_enrichment_out/cluster_enrichment_results.csv",
     by_file_counts_csv = "cluster_enrichment_out/by_file_cluster_counts.csv", 
     data_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR/",
     out_dir = "cluster_enrichment_out",
     
     # Your UMAP and marker columns
     umap_x = "UMAP1",
     umap_y = "UMAP2",
     markers_of_interest = c("CD45", "CD3", "CD4", "CD8"),
     
     # Memory limits (reduce if still crashing)
     max_cells_for_plots = 25000,
     enable_umap_plots = TRUE,
     enable_marker_plots = TRUE
   )
   ```

2. **Run the resume script**:
   ```r
   source("Resume_Fine_clustering_analysis.R")
   ```

### Option 2: Insert Safe Completion Code
Add the code from `Safe_visualization_completion.R` to your main script just after line 947 (where results are saved), replacing the existing visualization section.

### Option 3: Modify Original Script for Safety
Add these memory-safe modifications to your original script:

1. **Reduce plot complexity**:
   - Lower DPI from 300 to 150
   - Limit number of clusters in heatmaps (max 30)
   - Sample fewer cells for UMAP (max 25000)

2. **Add garbage collection**:
   ```r
   # After each plot
   rm(plot_object)
   gc()
   ```

3. **Skip memory-intensive plots**:
   ```r
   # Set these in CONFIG to disable problematic plots
   CONFIG$enable_umap_plots <- FALSE
   CONFIG$enable_marker_plots <- FALSE
   ```

## What Each Script Does

### Resume_Fine_clustering_analysis.R
- **Input**: Saved results CSV files
- **Output**: All visualizations
- **Safety**: Samples data, handles errors gracefully
- **Memory**: Optimized for large datasets

### Safe_visualization_completion.R  
- **Input**: Existing results object in memory
- **Output**: Basic visualizations only
- **Safety**: Skips advanced plots that cause crashes
- **Memory**: Minimal memory usage

## Recommended Workflow

1. **First try**: Run `Resume_Fine_clustering_analysis.R` with default settings
2. **If it crashes**: Reduce `max_cells_for_plots` to 10000 or 5000
3. **If still crashing**: Set `enable_umap_plots = FALSE` and `enable_marker_plots = FALSE`
4. **Minimal safe**: Use `Safe_visualization_completion.R` for basic plots only

## Files You Should Have
After the main script ran to line 947, you should have:
- `cluster_enrichment_out/cluster_enrichment_results.csv` (main results) ✅ REQUIRED
- `cluster_enrichment_out/analysis_summary.csv` (summary stats) ✅ REQUIRED

You may NOT have:
- `cluster_enrichment_out/by_file_cluster_counts.csv` ❌ MISSING if crashed at lines 949-957

**This is normal!** The resume script will work without the by_file_cluster_counts.csv file.
The crash likely happened during the memory-intensive pivot/group operations that create this file.

The resume scripts will use these to generate all visualizations without recomputing the statistical analysis.