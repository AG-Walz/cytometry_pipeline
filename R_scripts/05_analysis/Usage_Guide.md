# File-by-File Cluster Column Condenser - Usage Guide

## Overview
This script processes cytometry CSV files **one at a time** to condense multiple cluster columns into a single `cluster_name` column, creating space-efficient versions that save significant disk space.

## What it does
1. **Processes files individually** (no memory crashes!)
2. **Finds cluster columns** matching your regex pattern
3. **Condenses them** into a single `cluster_name` column
4. **Saves space-efficient files** to a separate output folder
5. **Tracks space savings** and provides detailed reports

## Configuration (Edit the CONFIG section)

### Required Settings:
- `input_dir`: Folder containing your CSV files
- `output_dir`: Where to save condensed files
- `cluster_regex`: Pattern to find cluster columns (e.g., "Clustering annotated")

### Cluster Assignment Methods:
- `"highest_membership"`: Single cluster or concatenate multiples with "|"
- `"concatenate"`: Always concatenate multiple clusters with "|" 
- `"first_found"`: Use first cluster found (column order)

### Memory Settings:
- `max_memory_gb`: Memory limit (default: 2GB)
- `max_cells_per_file`: Sample large files (default: 500k cells)
- `sample_large_files`: TRUE to sample, FALSE to process all

### Processing Options:
- `overwrite_existing`: FALSE to skip existing files
- `create_summary`: TRUE to generate detailed report

## Usage

### Basic Usage:
```bash
Rscript Create_labelled_exprmat.R
```

### Example Output:
```
Processing: sample_file.csv
  Reading 150000 rows (45.2 MB)
  Found 12 cluster columns: Clustering annotated_cluster_1, ...
  ✓ Processed successfully!
    Original: 150000 rows × 45 cols (45.2 MB)
    Condensed: 150000 rows × 34 cols (38.1 MB)
    Space saved: 7.1 MB (15.7%)
```

## Output Files

### Condensed Files:
- Location: `CONFIG$output_dir`
- Naming: `original_filename_condensed.csv`
- Content: All original columns EXCEPT cluster columns, plus new `cluster_name` column

### Summary Report:
- File: `condensation_summary.txt`
- Contains: Processing statistics, space savings, individual file results

## Key Features

### Memory Efficient:
- Processes **one file at a time** 
- **No memory crashes** from loading multiple files
- Automatic chunked reading for very large files
- Garbage collection every 5 files

### Space Savings:
- Typical savings: **10-30%** depending on number of cluster columns
- Removes redundant binary cluster columns
- Preserves all other data intact

### Robust Processing:
- Skips files that already exist (unless `overwrite_existing = TRUE`)
- Continues processing if individual files fail
- Detailed error reporting and progress tracking

## Troubleshooting

### If No Cluster Columns Found:
- Check your `cluster_regex` pattern
- Script will show available column names
- Files will be copied without condensation

### For Very Large Files:
- Set `sample_large_files = TRUE` to sample instead of processing all
- Reduce `max_cells_per_file` for smaller samples
- Increase `max_memory_gb` if you have more RAM

### Memory Issues:
- Reduce `max_memory_gb` 
- Enable `sample_large_files`
- Process files in smaller batches by splitting input directory

## Example Results

### Before (Original):
- 50 columns including 15 cluster columns (0s and 1s)
- File size: 120 MB

### After (Condensed):
- 36 columns with single `cluster_name` column
- File size: 95 MB  
- **Space saved: 25 MB (21%)**

The condensed files maintain all your original data but are much more storage-efficient and easier to work with!