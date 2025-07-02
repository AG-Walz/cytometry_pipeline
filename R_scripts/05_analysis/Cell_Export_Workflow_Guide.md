# Updated Workflow: Cell-Level Export and UMAP Visualization

## What's New

### 1. **Fine_clustering_cytolution.R** - Enhanced Export
The main clustering script now exports cell-level data for downstream analysis:

**New Output File:**
- `cell_cluster_assignments.csv` - Contains every cell with:
  - `cell_index`: Unique identifier for each cell
  - `file_name`: Origin file
  - `condition`: pCR/npCR assignment
  - `UMAP1`, `UMAP2` (if available): Coordinate data
  - `cluster_*`: All cluster assignments (0/1 for each cluster)

### 2. **UMAP_enrichment_visualization.R** - Streamlined Data Loading
The visualization script now uses the exported cell data instead of reloading raw files:

**Benefits:**
- ✅ **Faster**: No need to reload large CSV files
- ✅ **Consistent**: Uses exact same data as statistical analysis
- ✅ **Memory efficient**: Pre-processed and cleaned data
- ✅ **Reliable**: No risk of data processing differences

## Complete Workflow

### Step 1: Run Main Analysis
```r
# Edit CONFIG in Fine_clustering_cytolution.R
source("Fine_clustering_cytolution.R")
```

**Outputs:**
- `cluster_enrichment_results.csv` (statistical results)
- `cell_cluster_assignments.csv` ← **NEW!** (cell-level data)
- Various plots and summaries

### Step 2: Create UMAP Enrichment Visualization
```r
# Edit CONFIG in UMAP_enrichment_visualization.R
CONFIG <- list(
  results_csv = "path/to/cluster_enrichment_results.csv",
  cell_data_csv = "path/to/cell_cluster_assignments.csv",  # NEW!
  out_dir = "path/to/results/",
  umap_x = NULL,  # Interactive selection
  umap_y = NULL   # Interactive selection
)

source("UMAP_enrichment_visualization.R")
```

**What happens:**
1. Loads enrichment results and cell data
2. GUI prompts to select UMAP coordinates
3. Creates beautiful colored UMAP plots

## Key Advantages

### **Data Consistency**
- UMAP script uses **identical** cell data from main analysis
- No risk of different preprocessing or filtering
- Cell indices maintain traceability

### **Performance**
- **Before**: Load many large CSV files → process → visualize
- **After**: Load single processed file → visualize
- Significantly faster and more memory efficient

### **Flexibility**
- Cell-level data can be used for other downstream analyses
- Easy to create custom visualizations
- Full traceability with cell indices

### **Reliability**
- No file loading errors during visualization
- Consistent data formats
- Robust error handling

## File Structure
```
Results/
├── cluster_enrichment_results.csv     # Statistical results
├── cell_cluster_assignments.csv       # NEW: Cell-level data
├── UMAP_cluster_enrichment.png        # Main UMAP plot
├── UMAP_cluster_enrichment_by_condition.png
└── other plots...
```

## Cell Data Format
The new `cell_cluster_assignments.csv` contains:

| cell_index | file_name | condition | UMAP1 | UMAP2 | cluster_1 | cluster_2 | ... |
|------------|-----------|-----------|-------|-------|-----------|-----------|-----|
| 1 | file1.csv | pCR | -2.1 | 1.4 | 1 | 0 | ... |
| 2 | file1.csv | pCR | 0.5 | -1.2 | 0 | 1 | ... |
| 3 | file2.csv | npCR | 1.8 | 0.3 | 0 | 0 | ... |

## Migration Guide

### If you have the old workflow:
1. **Update** `Fine_clustering_cytolution.R` (already done)
2. **Run it once** to generate `cell_cluster_assignments.csv`
3. **Update** `UMAP_enrichment_visualization.R` CONFIG:
   ```r
   # REMOVE these lines:
   # data_dir = "path/to/raw/data/"
   # file_pattern = "\\.csv$"
   # max_cells_per_file = 50000
   
   # ADD this line:
   cell_data_csv = "path/to/cell_cluster_assignments.csv"
   ```
4. **Run visualization** - much faster now!

### For new analyses:
Just follow the complete workflow above - everything is optimized!

## Troubleshooting

**Error: "Cell data file not found"**
→ Run the updated `Fine_clustering_cytolution.R` first

**Error: "UMAP columns not found"**
→ Your original data didn't have UMAP coordinates, or they were filtered out

**Memory issues:**
→ Reduce `max_total_cells` in CONFIG

**GUI not working:**
→ Set specific column names: `umap_x = "UMAP1", umap_y = "UMAP2"`

This new workflow is much more robust and efficient! 🚀