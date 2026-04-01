[< Back to Documentation Index](../README.md)

# src/visualization/ -- Plotting

Thirteen scripts for generating publication-quality visualizations: UMAP embeddings, t-SNE plots, diffusion maps, heatmaps, density overlays, biaxial scatter plots, violin/box/sina plots, radar/spider charts, contour density plots, composition barplots, PCA scatter plots, UpSet intersection plots, and Venn/Euler diagrams. All plots are saved in configurable formats (PDF, PNG, SVG) via `save_plot_formats()`.

## Summary Table

| Script | Functions | Description |
|--------|-----------|-------------|
| `umap.R` | `plot_umap()`, `plot_umap_from_matrix()` | UMAP with sample/cluster/marker/condition overlays |
| `tsne.R` | `plot_tsne()`, `plot_tsne_from_matrix()` | t-SNE with sample/cluster/marker overlays |
| `diffusion_map.R` | `plot_diffusion_map()` | Diffusion map via destiny package |
| `heatmaps.R` | `plot_heatmap()`, `plot_cluster_heatmap()` | Clustered heatmaps via pheatmap |
| `density_plots.R` | `plot_density()`, `plot_density_from_matrix()` | Density and ridge plots per marker |
| `biaxial_plots.R` | `plot_biaxial()` | Classic biaxial scatter / pseudocolor plots |
| `violin_plots.R` | `plot_violin()` | Violin, box, sina, and combined plots |
| `radar_plots.R` | `plot_radar()` | Radar / spider charts for cluster profiles |
| `contour_plots.R` | `plot_contour()` | 2D density contour overlays |
| `barplots.R` | `plot_composition_barplot()` | Composition barplots (stacked/dodged/proportional) |
| `pca.R` | `plot_pca()` | PCA scatter + loadings plots |
| `upset_plots.R` | `plot_upset()`, `plot_upset_from_list()` | UpSet intersection plots |
| `venn_diagrams.R` | `plot_euler()`, `plot_venn()`, `plot_and_save_venn()` | Venn and Euler diagrams |

---

## umap.R

Computes UMAP embeddings and generates per-sample, per-cluster, per-condition, and per-marker overlay scatter plots. Used by Pipelines A and E.

### `plot_umap()`

```r
plot_umap(flow_frames, sample_names, markers_of_interest = NULL,
          output_dir, n_neighbors = 15, min_dist = 0.1,
          marker_selection_mode = "all",
          output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
          plot_width = 10, plot_height = 8, seed = 42)
```

flowFrame-based UMAP. Concatenates expression data from flowFrame objects, computes UMAP via `run_umap()`, and generates sample-colored and marker-colored plots.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `sample_names` | character vector | *required* | Sample identifiers |
| `markers_of_interest` | character vector | `NULL` | Markers for UMAP; `NULL` uses all (or interactive) |
| `output_dir` | character | *required* | Output directory |
| `n_neighbors` | integer | `15` | UMAP neighbors parameter |
| `min_dist` | numeric | `0.1` | UMAP minimum distance |
| `marker_selection_mode` | character | `"all"` | `"all"`, `"popup"`, or `"terminal"` |
| `output_formats` | list | see above | Which formats to save |
| `plot_width` | numeric | `10` | Plot width in inches |
| `plot_height` | numeric | `8` | Plot height in inches |
| `seed` | integer | `42` | Random seed |

**Returns:**

```r
list(
  umap_coords   = matrix(N x 2),
  sample_labels = character(),
  plots         = list(sample = ggplot, marker1 = ggplot, ...)
)
```

**Generated Files:** `umap_sample.png`, `umap_{marker}.png` (one per marker).

**Input from:** `load_fcs_raw()`, `downsample_flowframes()`, `homogenize_flowframes()`, `create_flowframes_from_matrices()`.
**Output to:** Standalone visualization output.

### `plot_umap_from_matrix()`

```r
plot_umap_from_matrix(expr_matrix, sample_labels = NULL,
                       markers_of_interest = NULL,
                       cluster_labels = NULL,
                       condition_labels = NULL,
                       output_dir, n_neighbors = 15, min_dist = 0.1,
                       umap_coords = NULL,
                       output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                       plot_width = 10, plot_height = 8, seed = 42)
```

Matrix-based UMAP. Accepts a numeric matrix directly, without requiring flowFrame objects. Supports pre-computed UMAP coordinates, cluster overlay, and condition faceting.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells x markers matrix |
| `sample_labels` | character vector | `NULL` | Sample label per cell |
| `markers_of_interest` | character vector | `NULL` | Markers to plot; `NULL` = all |
| `cluster_labels` | vector | `NULL` | Cluster assignments for overlay |
| `condition_labels` | character vector | `NULL` | Condition labels for faceted plot |
| `output_dir` | character | *required* | Output directory |
| `umap_coords` | matrix | `NULL` | Pre-computed 2-column UMAP coordinates |

**Returns:** `list(umap_coords, sample_labels, plots)`.

**Generated Files:** `umap_sample.png`, `umap_cluster.png` (if cluster_labels), `umap_condition.png` (if condition_labels), `umap_{marker}.png`.

**Input from:** `build_combined_expr_matrix()`, `combine_expression_matrices()`, `run_leiden_with_umap()`.
**Output to:** Standalone visualization output.

**Example:**

```r
source("umap.R")
result <- plot_umap_from_matrix(
  expr_matrix     = combined_matrix,
  sample_labels   = labels,
  cluster_labels  = clusters,
  output_dir      = "/output/umap",
  n_neighbors     = 15
)
```

**Dependencies:** `io_utils`, `gui_utils`, `plot_utils`, `stats_utils`, `flowframe_utils`

---

## tsne.R

t-SNE dimensionality reduction with dual input modes.

### `plot_tsne()`

```r
plot_tsne(flow_frames, sample_names, markers_of_interest = NULL,
          output_dir, perplexity = 30, max_iter = 1000,
          marker_selection_mode = "all",
          output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
          plot_width = 10, plot_height = 8, seed = 42)
```

flowFrame-based t-SNE. Concatenates expression data, runs `Rtsne::Rtsne()`, and generates sample-colored and marker-colored plots.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `sample_names` | character vector | *required* | Sample identifiers |
| `markers_of_interest` | character vector | `NULL` | Markers for t-SNE |
| `output_dir` | character | *required* | Output directory |
| `perplexity` | numeric | `30` | t-SNE perplexity parameter |
| `max_iter` | integer | `1000` | Maximum t-SNE iterations |
| `seed` | integer | `42` | Random seed |

**Returns:** `list(tsne_coords, sample_labels, plots)`.

### `plot_tsne_from_matrix()`

```r
plot_tsne_from_matrix(expr_matrix, sample_labels = NULL,
                       cluster_labels = NULL, output_dir,
                       perplexity = 30, max_iter = 1000,
                       tsne_coords = NULL,
                       output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                       plot_width = 10, plot_height = 8, seed = 42)
```

Matrix-based t-SNE. Accepts pre-computed t-SNE coordinates via `tsne_coords`. Generates sample, cluster, and per-marker plots.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells x markers matrix |
| `sample_labels` | character vector | `NULL` | Sample label per cell |
| `cluster_labels` | vector | `NULL` | Cluster assignments for overlay |
| `output_dir` | character | *required* | Output directory |
| `perplexity` | numeric | `30` | t-SNE perplexity |
| `max_iter` | integer | `1000` | Maximum iterations |
| `tsne_coords` | matrix | `NULL` | Pre-computed 2-column t-SNE coordinates |
| `seed` | integer | `42` | Random seed |

**Returns:** `list(tsne_coords, sample_labels, plots)`.

**Packages:** `Rtsne`, `ggplot2`, `viridis`

---

## diffusion_map.R

Diffusion map dimensionality reduction using the destiny package.

### `plot_diffusion_map()`

```r
plot_diffusion_map(expr_matrix, sample_labels = NULL,
                    cluster_labels = NULL, output_dir,
                    n_components = 3, sigma = "local", seed = 42,
                    output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
                    plot_width = 10, plot_height = 8)
```

Computes diffusion maps via `destiny::DiffusionMap()` and generates plots of DC1 vs DC2 and DC1 vs DC3, colored by sample, cluster, and each marker.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells x markers matrix |
| `sample_labels` | character vector | `NULL` | Sample labels |
| `cluster_labels` | vector | `NULL` | Cluster assignments for overlay |
| `output_dir` | character | *required* | Output directory |
| `n_components` | integer | `3` | Number of diffusion components |
| `sigma` | character/numeric | `"local"` | Kernel width: `"local"` for adaptive |
| `seed` | integer | `42` | Random seed |

**Returns:**

```r
list(
  dm_coords    = matrix,      # N x n_components
  eigenvalues  = numeric,     # Top eigenvalues
  plots        = list()       # Named list of ggplot objects
)
```

**Generated Files:** `diffmap_sample_dc1_vs_dc2.png`, `diffmap_cluster_dc1_vs_dc2.png`, `diffmap_{marker}_dc1_vs_dc2.png`, plus dc1_vs_dc3 variants.

**Packages:** `destiny`, `ggplot2`, `viridis`

---

## heatmaps.R

Generates clustered heatmaps of expression data using `pheatmap`.

### `plot_heatmap()`

```r
plot_heatmap(expr_matrix, output_path, title = "Expression Heatmap",
             cluster_rows = TRUE, cluster_cols = TRUE,
             scale = "row", color_palette = NULL,
             width = 10, height = 8)
```

General-purpose heatmap from a numeric matrix.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Numeric matrix (rows = features, columns = observations) |
| `output_path` | character | *required* | Full path for output image |
| `title` | character | `"Expression Heatmap"` | Plot title |
| `cluster_rows` | logical | `TRUE` | Cluster rows |
| `cluster_cols` | logical | `TRUE` | Cluster columns |
| `scale` | character | `"row"` | Scale by `"row"`, `"column"`, or `"none"` |
| `color_palette` | vector | `NULL` | Custom color palette |
| `width` | numeric | `10` | Plot width in inches |
| `height` | numeric | `8` | Plot height in inches |

**Returns:** Invisible pheatmap object.

### `plot_cluster_heatmap()`

```r
plot_cluster_heatmap(expr_data, cluster_col, marker_cols,
                     output_path, title = "Cluster Expression")
```

Computes median expression per cluster for each marker, then calls `plot_heatmap()`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_data` | data.frame | *required* | Expression data with cluster column |
| `cluster_col` | character | *required* | Name of cluster column |
| `marker_cols` | character vector | *required* | Marker column names |
| `output_path` | character | *required* | Output path |
| `title` | character | `"Cluster Expression"` | Plot title |

**Returns:** Invisible pheatmap object.

**Packages:** `pheatmap`

---

## density_plots.R

Plots per-marker density overlays or ridge plots, grouped by sample or condition.

### `plot_density()`

```r
plot_density(flow_frames, sample_names, markers, output_dir,
             plot_type = "density", group_by = "sample",
             group_labels = NULL)
```

flowFrame-based density plots. Generates per-marker density overlays or ridge plots.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flow_frames` | list | *required* | Named list of flowFrame objects |
| `sample_names` | character vector | *required* | Sample identifiers |
| `markers` | character vector | *required* | Markers to plot |
| `output_dir` | character | *required* | Output directory |
| `plot_type` | character | `"density"` | `"density"`, `"ridge"`, or `"both"` |
| `group_by` | character | `"sample"` | Grouping variable |
| `group_labels` | named vector | `NULL` | Maps samples to group labels |

**Returns:** Invisible `NULL`.

**Generated Files:** `density_{marker}.png`, `ridge_{marker}.png`.

### `plot_density_from_matrix()`

```r
plot_density_from_matrix(expr_matrix, sample_labels = NULL,
                          markers = NULL, output_dir,
                          plot_type = "density",
                          group_labels = NULL)
```

Matrix-based density plots. Accepts a numeric matrix directly.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells x markers matrix |
| `sample_labels` | character vector | `NULL` | Sample labels per cell |
| `markers` | character vector | `NULL` | Markers to plot; `NULL` = all columns |
| `output_dir` | character | *required* | Output directory |
| `plot_type` | character | `"density"` | `"density"`, `"ridge"`, or `"both"` |
| `group_labels` | named vector | `NULL` | Maps samples to group labels |

**Returns:** Invisible `NULL`.

**Dependencies:** `io_utils`, `flowframe_utils`, `plot_utils`. Ridge plots use `ggridges`.

---

## biaxial_plots.R

Classic flow cytometry biaxial scatter/pseudocolor plots with density, marker, or cluster coloring.

### `plot_biaxial()`

```r
plot_biaxial(flowframe = NULL, expr_matrix = NULL,
              x_marker, y_marker, color_by = "density",
              sample_label = NULL, cluster_labels = NULL,
              gate_coords = NULL, output_dir,
              output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
              plot_width = 10, plot_height = 8)
```

Dual-input biaxial plot (flowFrame or matrix). Supports density pseudocolor (via `MASS::kde2d`), marker expression overlay, and cluster assignment coloring. Optionally draws polygon gate boundaries.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `flowframe` | flowFrame | `NULL` | Optional flowFrame input |
| `expr_matrix` | matrix | `NULL` | Optional matrix input (one of the two required) |
| `x_marker` | character | *required* | Marker for x axis |
| `y_marker` | character | *required* | Marker for y axis |
| `color_by` | character | `"density"` | `"density"`, a marker name, or `"cluster"` |
| `sample_label` | character | `NULL` | Sample name for title |
| `cluster_labels` | vector | `NULL` | Cluster assignments (required when `color_by = "cluster"`) |
| `gate_coords` | list | `NULL` | `list(x = numeric, y = numeric)` for polygon gate overlay |
| `output_dir` | character | *required* | Output directory |

**Returns:** Invisible ggplot object.

**Generated File:** `biaxial_{x}_vs_{y}_{color}.png`.

**Example:**

```r
source("biaxial_plots.R")
plot_biaxial(
  flowframe    = my_ff,
  x_marker     = "CD3",
  y_marker     = "CD4",
  color_by     = "density",
  output_dir   = "/output/biaxial"
)
```

---

## violin_plots.R

Expression distribution plots grouped by cluster or condition.

### `plot_violin()`

```r
plot_violin(expr_data, markers, group_col = "cluster",
             plot_type = "violin_box", output_dir,
             output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
             plot_width = 12, plot_height = 8)
```

Generates one plot per marker showing expression distributions across groups.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_data` | data.frame | *required* | Data with marker columns and a grouping column |
| `markers` | character vector | *required* | Marker column names |
| `group_col` | character | `"cluster"` | Grouping column name |
| `plot_type` | character | `"violin_box"` | `"violin"`, `"box"`, `"sina"`, `"violin_box"` |
| `output_dir` | character | *required* | Output directory |

**Returns:** Invisible `NULL`.

**Generated Files:** `{plot_type}_{marker}.png` (one per marker).

**Packages:** `ggplot2`, `ggforce` (for sina plots)

---

## radar_plots.R

Radar (spider) charts for visualizing cluster marker expression profiles.

### `plot_radar()`

```r
plot_radar(cluster_medians, markers, output_dir,
            normalize = TRUE,
            output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
            plot_width = 10, plot_height = 10)
```

Generates one radar plot per cluster plus one combined overlay plot. When `normalize = TRUE`, each marker is scaled to 0-1 across all clusters.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cluster_medians` | data.frame | *required* | Rows = clusters, columns include markers (+ optional `cluster` column) |
| `markers` | character vector | *required* | Marker column names |
| `output_dir` | character | *required* | Output directory |
| `normalize` | logical | `TRUE` | Scale each marker to 0-1 |

**Returns:** Invisible `NULL`.

**Generated Files:** `radar_{cluster}.png` (one per cluster), `radar_combined.png`.

---

## contour_plots.R

2D contour density overlay plots for dimensionality reduction coordinates.

### `plot_contour()`

```r
plot_contour(coords, group_labels = NULL, x_col = "x", y_col = "y",
              method = "kde2d", output_dir,
              output_formats = list(pdf = FALSE, png = TRUE, svg = TRUE),
              plot_width = 10, plot_height = 8)
```

Generates contour plots from 2D coordinates with optional per-group coloring.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `coords` | data.frame | *required* | Must have `x` and `y` columns |
| `group_labels` | vector | `NULL` | Group labels for colored contours |
| `x_col` | character | `"x"` | Name of x coordinate column |
| `y_col` | character | `"y"` | Name of y coordinate column |
| `method` | character | `"kde2d"` | `"kde2d"` (filled via MASS) or `"geom_density2d"` (line contours) |
| `output_dir` | character | *required* | Output directory |

**Returns:** Invisible ggplot object.

**Generated File:** `contour_{method}[_grouped].png`.

**Packages:** `MASS` (for kde2d), `ggplot2`

---

## barplots.R

Creates composition barplots showing cell distribution across clusters and samples.

### `plot_composition_barplot()`

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

**Returns:** Invisible ggplot object.

**Generated File:** `barplot_{x_col}_{fill_col}.png`.

---

## pca.R

Performs PCA and generates scatter plots of components plus loadings plots.

### `plot_pca()`

```r
plot_pca(expr_matrix, sample_labels = NULL, output_dir,
         n_components = 2, color_by = "sample", title = "PCA")
```

Runs `prcomp()` on the expression matrix and generates a scatter plot of PC1 vs PC2 plus a variable loadings plot.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expr_matrix` | matrix | *required* | Cells x markers matrix |
| `sample_labels` | character vector | `NULL` | Sample label per cell |
| `output_dir` | character | *required* | Output directory |
| `n_components` | integer | `2` | Number of PCs |
| `color_by` | character | `"sample"` | Coloring variable |
| `title` | character | `"PCA"` | Plot title |

**Returns:**

```r
list(
  pca_result    = prcomp(),
  scatter_plot  = ggplot(),
  loadings_plot = ggplot()
)
```

**Generated Files:** `pca_plot.png`, `pca_loadings.png`.

**Dependencies:** `io_utils`, `stats_utils`, `plot_utils`

---

## upset_plots.R

UpSet plots for visualizing intersections of multiple sets. Used by Pipeline D.

### `plot_upset()`

```r
plot_upset(data, set_columns, output_path,
           title = "UpSet Plot", min_size = 1)
```

Creates an UpSet plot from a binary matrix.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data` | data.frame | *required* | Binary matrix (rows = elements, columns = sets) |
| `set_columns` | character vector | *required* | Column names representing sets |
| `output_path` | character | *required* | Output image path |
| `title` | character | `"UpSet Plot"` | Plot title |
| `min_size` | integer | `1` | Minimum intersection size to display |

**Returns:** Invisible UpSetR object.

### `plot_upset_from_list()`

```r
plot_upset_from_list(sets, output_path, title = "UpSet Plot")
```

Creates an UpSet plot from a named list of character vectors.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sets` | named list | *required* | Named list of character vectors (set elements) |
| `output_path` | character | *required* | Output image path |

**Returns:** Invisible UpSetR object.

**Packages:** `UpSetR`

---

## venn_diagrams.R

Venn diagrams and area-proportional Euler diagrams. Used by Pipeline D.

### `plot_euler()`

```r
plot_euler(sets, output_dir, filename, custom_name = NULL)
```

Area-proportional Euler diagram using `eulerr`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `sets` | named list | *required* | Named list of character vectors |
| `output_dir` | character | *required* | Output directory |
| `filename` | character | *required* | Output filename (without extension) |
| `custom_name` | character | `NULL` | Custom set names |

**Returns:** Invisible euler plot object.

### `plot_venn()`

```r
plot_venn(sets, filename, output_dir)
```

Non-proportional Venn diagram using `venn`.

### `plot_and_save_venn()`

```r
plot_and_save_venn(lists, population, filename, output_dir,
                   scale_size = FALSE, custom_name = NULL)
```

Convenience wrapper: filters binary matrix by population, dispatches to `plot_euler()` or `plot_venn()` based on `scale_size`.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `lists` | data.frame | *required* | Label lists (columns = raters/analysts) |
| `population` | character | *required* | Population label to compare |
| `filename` | character | *required* | Output filename |
| `output_dir` | character | *required* | Output directory |
| `scale_size` | logical | `FALSE` | `TRUE` = Euler (proportional), `FALSE` = Venn |
| `custom_name` | character | `NULL` | Custom names for each set |

**Example:**

```r
source("venn_diagrams.R")

sets <- list(
  "Analyst A" = c("cell_1", "cell_2", "cell_3", "cell_4"),
  "Analyst B" = c("cell_2", "cell_3", "cell_5", "cell_6"),
  "Analyst C" = c("cell_1", "cell_3", "cell_5", "cell_7")
)

plot_euler(sets, output_dir = "/output/venn", filename = "CD4_overlap")
```

**Packages:** `eulerr`, `venn`
