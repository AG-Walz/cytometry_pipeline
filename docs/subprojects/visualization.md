[< Back to Documentation Index](../README.md)

# src/visualization/ — Plotting

Seven scripts for generating publication-quality visualizations: composition barplots, density overlays, clustered heatmaps, PCA scatter plots, UMAP embeddings, UpSet intersection plots, and Venn/Euler diagrams. All plots are saved in configurable formats (PDF, PNG, SVG) via `save_plot_formats()`.

---

## umap.R

Computes UMAP embeddings from flowFrame objects and generates per-sample and per-marker overlay scatter plots. Used by Pipelines A and E.

### Function

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
| `markers_of_interest` | character vector | `NULL` | Markers for UMAP; `NULL` uses all |
| `output_dir` | character | *required* | Output directory |
| `n_neighbors` | integer | `15` | UMAP neighbors parameter |
| `min_dist` | numeric | `0.1` | UMAP minimum distance |
| `output_formats` | list | see above | Which formats to save |

### Returns

```r
list(
  umap_coords   = matrix(),    # N x 2 coordinates
  sample_labels = character(),  # Sample per cell
  plots         = list()        # Named list of ggplot objects
)
```

### Generated Files

```
output_dir/
├── umap_sample.png        # Points colored by sample identity (categorical)
├── umap_CD3.png           # Points colored by CD3 expression (viridis)
├── umap_CD4.png           # Points colored by CD4 expression
├── umap_CD8.png
└── ...                    # One plot per marker
```

### Example

```r
source("umap.R")
result <- plot_umap(
  flow_frames         = my_frames,
  sample_names        = c("S1", "S2", "S3"),
  markers_of_interest = c("CD3", "CD4", "CD8"),
  output_dir          = "/output/umap",
  n_neighbors         = 15,
  seed                = 42
)
# Creates: umap_sample.png, umap_CD3.png, umap_CD4.png, umap_CD8.png
```

### Dependencies

`io_utils`, `gui_utils`, `plot_utils`, `stats_utils` (`run_umap`), `flowframe_utils`

---

## barplots.R

Creates composition barplots showing cell distribution across clusters and samples.

### Function

```r
plot_composition_barplot(data, output_dir, x_col = "cluster",
                         fill_col = "sample", position = "stack",
                         title = "Cell Composition")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data` | data.frame | *required* | Long-format data with category and count columns |
| `output_dir` | character | *required* | Output directory |
| `x_col` | character | `"cluster"` | Column for x-axis categories |
| `fill_col` | character | `"sample"` | Column for fill (stacking) variable |
| `position` | character | `"stack"` | `"stack"`, `"dodge"`, or `"fill"` |
| `title` | character | `"Cell Composition"` | Plot title |

### Example Input Data

```r
data <- data.frame(
  cluster = c("CD4+", "CD4+", "CD8+", "CD8+"),
  file    = c("sample_001", "sample_002", "sample_001", "sample_002"),
  count   = c(1234, 987, 567, 1234)
)
```

### Generated File

`barplot_{x_col}_{fill_col}.png` — a stacked/dodged/proportional barplot.

### Example

```r
source("barplots.R")
plot_composition_barplot(
  data       = long_data,
  output_dir = "/output/plots",
  x_col      = "cluster",
  fill_col   = "file",
  position   = "stack",
  title      = "Cell Composition by Cluster"
)
```

---

## density_plots.R

Plots per-marker density overlays or ridge plots, grouped by sample or condition.

### Function

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
| `group_by` | character | `"sample"` | Grouping variable |

### Generated Files

| Plot Type | Filename Pattern |
|-----------|-----------------|
| `"density"` | `density_{marker}.png` |
| `"ridge"` | `ridge_{marker}.png` |
| `"both"` | Both of the above |

### Example

```r
source("density_plots.R")
plot_density(
  flow_frames  = my_frames,
  sample_names = c("S1", "S2", "S3"),
  markers      = c("CD3", "CD4", "CD8"),
  output_dir   = "/output/density",
  plot_type    = "both"
)
# Creates: density_CD3.png, ridge_CD3.png, density_CD4.png, ridge_CD4.png, ...
```

### Dependencies

`io_utils`, `flowframe_utils`. Ridge plots use `ggridges`.

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
| `expr_matrix` | matrix | *required* | Numeric matrix (rows = features, columns = observations) |
| `output_path` | character | *required* | Full path for output image |
| `scale` | character | `"row"` | Scale by `"row"`, `"column"`, or `"none"` |

#### plot_cluster_heatmap

```r
plot_cluster_heatmap(expr_data, cluster_col, marker_cols,
                     output_path, title = "Cluster Expression")
```

Computes **median expression** per cluster for each marker, then calls `plot_heatmap()`.

### Example

```r
source("heatmaps.R")

# Raw expression heatmap (sampled cells)
plot_heatmap(
  expr_matrix = t(expr_matrix[1:500, ]),
  output_path = "/output/expression_heatmap.png",
  title       = "Expression Heatmap",
  scale       = "row"
)

# Cluster-level summary heatmap
plot_cluster_heatmap(
  expr_data   = clustered_df,
  cluster_col = "cluster",
  marker_cols = c("CD3", "CD4", "CD8", "CD19"),
  output_path = "/output/cluster_heatmap.png"
)
```

### Packages

`pheatmap`

---

## pca.R

Performs PCA and generates scatter plots of components plus loadings plots.

### Function

```r
plot_pca(expr_matrix, sample_labels = NULL, output_dir,
         n_components = 2, color_by = "sample", title = "PCA")
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells × markers matrix |
| `sample_labels` | character vector | `NULL` | Sample label per cell |
| `output_dir` | character | *required* | Output directory |
| `n_components` | integer | `2` | Number of PCs |

### Returns

```r
list(
  pca_result    = prcomp(),
  scatter_plot  = ggplot(),
  loadings_plot = ggplot()
)
```

### Generated Files

| File | Description |
|------|-------------|
| `pca_plot.png` | PCA scatter plot (PC1 vs PC2) |
| `pca_loadings.png` | Variable loadings plot |

### Dependencies

`io_utils`, `stats_utils` (`run_pca`), `plot_utils` (`save_plot_formats`)

---

## upset_plots.R

Creates UpSet plots for visualizing intersections of multiple sets. Used by Pipeline D.

### Functions

```r
plot_upset(data, set_columns, output_path,
           title = "UpSet Plot", min_size = 1)
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `data` | data.frame | Binary matrix (rows = elements, columns = sets) |
| `set_columns` | character vector | Column names representing sets |
| `output_path` | character | Output image path |

```r
plot_upset_from_list(sets, output_path, title = "UpSet Plot")
```

| Parameter | Type | Description |
|-----------|------|-------------|
| `sets` | named list | Named list of character vectors (set elements) |

### Packages

`UpSetR`

---

## venn_diagrams.R

Creates Venn diagrams and area-proportional Euler diagrams for comparing sets. Used by Pipeline D.

### Functions

```r
plot_euler(sets, output_dir, filename, custom_name = NULL)
```

Area-proportional Euler diagram using `eulerr`.

```r
plot_venn(sets, filename, output_dir)
```

Non-proportional Venn diagram using `venn`.

```r
plot_and_save_venn(lists, population, filename, output_dir,
                   scale_size = FALSE, custom_name = NULL)
```

Convenience wrapper: filters binary matrix by population, dispatches to `plot_euler()` or `plot_venn()` based on `scale_size`.

### Example

```r
source("venn_diagrams.R")

sets <- list(
  "Analyst A" = c("cell_1", "cell_2", "cell_3", "cell_4"),
  "Analyst B" = c("cell_2", "cell_3", "cell_5", "cell_6"),
  "Analyst C" = c("cell_1", "cell_3", "cell_5", "cell_7")
)

plot_euler(sets, output_dir = "/output/venn", filename = "CD4_overlap")
# Creates: CD4_overlap.png — area-proportional Euler diagram
```

### Packages

`eulerr`, `venn`
