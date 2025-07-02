# Create Labelled Expression Matrix

This script extracts and concatenates cytolution CSV files into a simplified format with cluster assignments.

## Purpose

Transform multiple cytometry CSV files with binary cluster columns into a single dataset where each cell has one cluster assignment, along with UMAP coordinates and selected expression markers.

## Usage

```bash
Rscript Create_labelled_exprmat.R
```

## Configuration

Edit the CONFIG block at the top of the script:

- `data_dir`: Path to directory containing CSV files
- `cluster_regex`: Pattern to match cluster columns (default: "Clustering annotated")
- `out_dir`: Output directory for results
- `cluster_assignment`: Method for handling multiple cluster memberships
  - `"highest_membership"`: Use first cluster found per cell
  - `"concatenate"`: Join multiple clusters with "|"
  - `"first_found"`: Use first matching cluster column

## Interactive Selection

The script will open GUI windows to select:

1. **UMAP X coordinate column** (single selection)
2. **UMAP Y coordinate column** (single selection)  
3. **Markers of interest** (multiple selection allowed)

## Output Files

### `labelled_exprmat.csv`
Main output with columns:
- `cell_index`: Unique cell identifier
- `cluster_name`: Single cluster assignment per cell
- `file_origin`: Original CSV filename 
- `umap_x`, `umap_y`: Selected UMAP coordinates
- Expression markers: Selected marker columns

### `labelled_exprmat_summary.txt`
Summary report containing:
- Processing statistics
- Cluster distribution
- File distribution
- Configuration used

## Memory Management

The script includes several memory optimization features:

- **File size assessment**: Large files are sampled or chunked
- **Column filtering**: Only essential columns are loaded
- **Memory monitoring**: Stops processing if memory usage exceeds 85%
- **Progress tracking**: Visual feedback during processing

## Example Output Structure

```csv
cell_index,cluster_name,file_origin,umap_x,umap_y,CD3,CD4,CD8
1,Clustering annotated_cluster_01,sample_001,12.3,8.7,1.2,0.8,0.1
2,Clustering annotated_cluster_02,sample_001,15.6,9.2,0.9,0.2,1.8
3,unassigned,sample_001,10.1,7.4,0.5,0.3,0.4
...
```

## Requirements

- R packages: tidyverse, data.table, janitor, stringr, progress, tcltk
- CSV files with:
  - Binary cluster columns (0/1 values)
  - UMAP/dimensionality reduction coordinates
  - Expression marker columns
- GUI environment (for tkinter column selection)

## Notes

- Cells not assigned to any cluster get `cluster_name = "unassigned"`
- Multiple cluster assignments are handled based on `cluster_assignment` method
- Large files are automatically sampled to prevent memory issues
- Missing columns are handled gracefully with informative warnings