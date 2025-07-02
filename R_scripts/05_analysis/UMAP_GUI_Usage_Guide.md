# Enhanced UMAP Enrichment Visualization - Usage Guide

### 1. Automatic Interactive Selection (Default)
```r
CONFIG <- list(
  results_csv = "/path/to/cluster_enrichment_results.csv",
  data_dir = "/path/to/data/",
  out_dir = "/path/to/results/",
  
  # Set to NULL to trigger GUI selection
  umap_x = NULL,  # Will open GUI to select X coordinate
  umap_y = NULL,  # Will open GUI to select Y coordinate
  
  # ... other settings
)
```

When you run the script:
1. 📁 Reads column names from your first CSV file
2. 🖱️ Opens GUI: "Select UMAP X coordinate" 
3. 🖱️ Opens GUI: "Select UMAP Y coordinate" (excludes already selected X)
4. ✅ Proceeds with analysis using your selections

### 2. Manual Specification (Skip GUI)
```r
CONFIG <- list(
  # ... paths ...
  
  # Specify directly to skip GUI
  umap_x = "UMAP1",
  umap_y = "UMAP2",
  
  # ... other settings
)
```

## GUI Features
- **Scrollable list**: See all available columns
- **Search-friendly**: Column names displayed clearly
- **Single selection**: One column at a time for X and Y
- **Smart filtering**: Y coordinate options exclude already selected X
- **Cancel option**: Exit gracefully if needed

## Example Workflow

1. **Edit CONFIG section** with your file paths:
   ```r
   results_csv = "/data/Cytometry/ICK/ICK_T1_P+_poR/Results/cluster_enrichment_results.csv"
   data_dir = "/data/Cytometry/ICK/ICK_T1_P+_poR/"
   ```

2. **Set UMAP coordinates to NULL** (for interactive selection):
   ```r
   umap_x = NULL,
   umap_y = NULL,
   ```

3. **Run the script**:
   ```r
   source("UMAP_enrichment_visualization.R")
   ```

4. **Select coordinates in GUI**:
   - First dialog: Choose X coordinate (e.g., "UMAP1", "tsne_1", "pc1")
   - Second dialog: Choose Y coordinate (e.g., "UMAP2", "tsne_2", "pc2")

5. **Get beautiful plots**:
   - 🔴 Red: pCR-enriched clusters
   - 🔵 Blue: npCR-enriched clusters  
   - ⚪ Gray: Non-significant clusters

## Advantages of Interactive Selection

- **Flexible**: Works with any coordinate system (UMAP, t-SNE, PCA, etc.)
- **No guessing**: See exactly what columns are available
- **Error-proof**: Can't select non-existent columns
- **User-friendly**: Same interface as main analysis script
- **Reusable**: Save your selections in CONFIG for future runs

## Tips

1. **For batch processing**: Set specific column names in CONFIG to skip GUI
2. **For exploration**: Use NULL to interactively explore different coordinate systems
3. **For new datasets**: Always use interactive selection first to see available options
4. **For automation**: Once you know your column names, hardcode them in CONFIG

## Output
Same beautiful UMAP visualizations as before, but now with complete flexibility in coordinate selection!

- `UMAP_cluster_enrichment.png`: Main enrichment visualization
- `UMAP_cluster_enrichment_by_condition.png`: Faceted by condition (if available)

Enjoy exploring your cytometry data with complete coordinate flexibility! 🎯