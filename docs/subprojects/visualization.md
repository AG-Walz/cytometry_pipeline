[< Back to Documentation Index](../README.md)

# visualization/ -- Plotting

7 scripts for generating publication-quality visualizations of cytometry data, including composition barplots, density overlays, heatmaps, PCA, UMAP, UpSet plots, and Venn/Euler diagrams.

---

## barplots.R

Creates composition barplots showing cell distribution across clusters and samples.

### Main Function

```r
plot_composition_barplot(data, output_dir, x_col = "cluster",
                         fill_col = "sample", position = "stack",
                         title = "Cell Composition")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data` | data.frame | *required* | Data with at least `x_col` and `fill_col` columns |
| `output_dir` | character | *required* | Output directory |
| `x_col` | character | `"cluster"` | Column for x-axis categories |
| `fill_col` | character | `"sample"` | Column for fill (stacking) variable |
| `position` | character | `"stack"` | ggplot2 position: `"stack"`, `"dodge"`, or `"fill"` |
| `title` | character | `"Cell Composition"` | Plot title |

### Output

`barplot_{x_col}_{fill_col}.png`

---

## density_plots.R

Plots per-marker density overlays or ridge plots, grouped by sample or condition.

### Main Function

```r
plot_density(flow_frames, sample_names, markers, output_dir,
             plot_type = "density", group_by = "sample",
             group_labels = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `sample_names` | character vector | *required* | Sample identifiers |
| `markers` | character vector | *required* | Markers to plot |
| `output_dir` | character | *required* | Output directory |
| `plot_type` | character | `"density"` | `"density"`, `"ridge"`, or `"both"` |
| `group_by` | character | `"sample"` | Grouping variable for overlays |
| `group_labels` | character vector | `NULL` | Custom labels for groups |

### Output

| Plot Type | Filename Pattern |
|-----------|-----------------|
| `"density"` | `density_{marker}.png` |
| `"ridge"` | `ridge_{marker}.png` |
| `"both"` | Both of the above |

Ridge plots use the `ggridges` package.

### R_shared Dependencies

`io_utils`, `flowframe_utils` (`explore_data`)

---

## heatmaps.R

Generates clustered heatmaps of expression data using `pheatmap`.

### Functions

#### plot_heatmap

```r
plot_heatmap(expr_matrix, output_path,
             title = "Expression Heatmap",
             cluster_rows = TRUE, cluster_cols = TRUE,
             scale = "row", color_palette = NULL,
             width = 10, height = 8)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Numeric matrix (rows = observations, columns = markers) |
| `output_path` | character | *required* | Full path for the output image |
| `title` | character | `"Expression Heatmap"` | Plot title |
| `cluster_rows` | logical | `TRUE` | Hierarchically cluster rows |
| `cluster_cols` | logical | `TRUE` | Hierarchically cluster columns |
| `scale` | character | `"row"` | Scale data by `"row"`, `"column"`, or `"none"` |
| `color_palette` | character vector | `NULL` | Custom color palette; `NULL` uses default |
| `width` | numeric | `10` | Output width in inches |
| `height` | numeric | `8` | Output height in inches |

#### plot_cluster_heatmap

```r
plot_cluster_heatmap(expr_data, cluster_col, marker_cols,
                     output_path, title = "Cluster Expression")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_data` | data.frame | *required* | Expression data with a cluster column |
| `cluster_col` | character | *required* | Name of the cluster column |
| `marker_cols` | character vector | *required* | Names of marker columns |
| `output_path` | character | *required* | Full path for the output image |
| `title` | character | `"Cluster Expression"` | Plot title |

Calculates the **median expression** per cluster for each marker, then calls `plot_heatmap()`.

### Packages

`pheatmap`

---

## pca.R

Performs PCA and generates scatter plots of components plus loadings plots.

### Main Function

```r
plot_pca(expr_matrix, sample_labels = NULL, output_dir,
         n_components = 2, color_by = "sample", title = "PCA")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells-by-markers expression matrix |
| `sample_labels` | character vector | `NULL` | Sample label per cell for coloring |
| `output_dir` | character | *required* | Output directory |
| `n_components` | integer | `2` | Number of principal components |
| `color_by` | character | `"sample"` | Variable to color points by |
| `title` | character | `"PCA"` | Plot title |

### Output

| File | Description |
|------|-------------|
| `pca_plot.{pdf,png,svg}` | PCA scatter plot (PC1 vs PC2) |
| `pca_loadings.{pdf,png,svg}` | Variable loadings plot |

### Return Value

```r
list(
  pca_result    = prcomp(),   # PCA result object
  scatter_plot  = ggplot(),   # PCA scatter plot
  loadings_plot = ggplot()    # Loadings plot
)
```

### R_shared Dependencies

`io_utils`, `stats_utils` (`run_pca`), `plot_utils` (`save_plot_formats`)

---

## umap.R

Computes UMAP embeddings and generates per-sample and per-marker overlay plots.

### Main Function

```r
plot_umap(flow_frames, sample_names, markers_of_interest = NULL,
          output_dir, n_neighbors = 15, min_dist = 0.1,
          marker_selection_mode = "all",
          output_formats = list(pdf = FALSE, png = TRUE, svg = FALSE),
          plot_width = 8, plot_height = 6, seed = 42)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `sample_names` | character vector | *required* | Sample identifiers |
| `markers_of_interest` | character vector | `NULL` | Markers for UMAP computation; `NULL` uses all |
| `output_dir` | character | *required* | Output directory |
| `n_neighbors` | integer | `15` | UMAP nearest neighbors parameter |
| `min_dist` | numeric | `0.1` | UMAP minimum distance parameter |
| `marker_selection_mode` | character | `"all"` | `"all"` or `"interactive"` |
| `output_formats` | list | `list(pdf=FALSE, png=TRUE, svg=FALSE)` | Which formats to save |
| `plot_width` | numeric | `8` | Plot width in inches |
| `plot_height` | numeric | `6` | Plot height in inches |
| `seed` | integer | `42` | Random seed for reproducibility |

### Output

| File | Description |
|------|-------------|
| `umap_sample.png` | UMAP colored by sample identity |
| `umap_{marker}.png` | UMAP colored by marker expression (viridis scale), one per marker |

### Return Value

```r
list(
  umap_coords   = matrix(),  # N x 2 UMAP coordinates
  sample_labels = character(), # Sample label per cell
  plots         = list()      # Named list of ggplot objects
)
```

### R_shared Dependencies

`io_utils`, `gui_utils`, `plot_utils`, `stats_utils` (`run_umap`), `flowframe_utils` (`explore_data`)

---

## upset_plots.R

Creates UpSet plots for visualizing intersections of multiple sets.

### Functions

#### plot_upset

```r
plot_upset(data, set_columns, output_path,
           title = "UpSet Plot", min_size = 1)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data` | data.frame | *required* | Binary matrix (rows = elements, columns = sets) |
| `set_columns` | character vector | *required* | Column names representing sets |
| `output_path` | character | *required* | Full path for the output image |
| `title` | character | `"UpSet Plot"` | Plot title |
| `min_size` | integer | `1` | Minimum intersection size to display |

#### plot_upset_from_list

```r
plot_upset_from_list(sets, output_path, title = "UpSet Plot")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sets` | named list | *required* | Named list of character vectors (set elements) |
| `output_path` | character | *required* | Full path for the output image |
| `title` | character | `"UpSet Plot"` | Plot title |

### Packages

`UpSetR`

---

## venn_diagrams.R

Creates Venn diagrams and area-proportional Euler diagrams for comparing sets.

### Functions

#### plot_euler

```r
plot_euler(sets, output_dir, filename, custom_name = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sets` | named list | *required* | Named list of character vectors |
| `output_dir` | character | *required* | Output directory |
| `filename` | character | *required* | Base filename (without extension) |
| `custom_name` | character | `NULL` | Custom title; `NULL` uses `filename` |

Generates an area-proportional Euler diagram using the `eulerr` package.

#### plot_venn

```r
plot_venn(sets, filename, output_dir)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sets` | named list | *required* | Named list of character vectors |
| `filename` | character | *required* | Base filename |
| `output_dir` | character | *required* | Output directory |

Generates a non-proportional Venn diagram using the `venn` package.

#### plot_and_save_venn

```r
plot_and_save_venn(lists, population, filename, output_dir,
                   scale_size = FALSE, custom_name = NULL)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `lists` | data.frame | *required* | Binary matrix (columns = sets, rows = elements) |
| `population` | character | *required* | Population name to filter by |
| `filename` | character | *required* | Base filename |
| `output_dir` | character | *required* | Output directory |
| `scale_size` | logical | `FALSE` | If `TRUE`, use Euler (proportional); if `FALSE`, use Venn |
| `custom_name` | character | `NULL` | Custom plot title |

Filters the binary matrix by the specified population, then dispatches to `plot_euler()` or `plot_venn()` based on `scale_size`.

### Packages

`eulerr`, `venn`
