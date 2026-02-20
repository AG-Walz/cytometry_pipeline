# Set Visualization: Venn Diagrams, Euler Diagrams, and UpSet Plots

## Summary

Set visualization methods display the overlap and intersection structure among two or more categorical groupings of cells. In this pipeline, they are used primarily to compare gating strategies across FlowJo workspaces, automated clustering algorithms, or analysts, answering the question "which cells are classified identically and which are classified differently?" Three complementary techniques are available: fixed-circle Venn diagrams for simple 2--3 set comparisons, area-proportional Euler diagrams for visually accurate overlap representation, and matrix-based UpSet plots for precise quantification of complex multi-set intersections.

---

## Venn Diagrams

### Historical Background

The Venn diagram was [introduced by John Venn in 1880 as a scheme for representing logical propositions and syllogisms using overlapping closed curves](https://doi.org/10.1080/14786448008626877). Each region formed by the intersection of curves corresponds to a possible combination of set memberships, and the presence or absence of elements in that region indicates whether that combination is realized in the data. The original formulation used circles, which [can represent all possible intersections for up to three sets but not for four or more](https://doi.org/10.1080/14786448008626877), a limitation that motivated later extensions using ellipses and irregular shapes.

### How Venn Diagrams Work

A Venn diagram draws one closed curve per set, positioned so that every possible intersection region exists as a distinct, non-empty area in the diagram. For two sets, this produces three regions (A only, B only, A and B). For three sets, it produces seven regions (A only, B only, C only, AB, AC, BC, ABC). The number of regions grows as 2^n - 1 for n sets.

In a standard Venn diagram, the circles (or curves) are drawn at **fixed size and position** regardless of the actual cardinalities of the sets or their intersections. The count of elements in each region is typically annotated as a label overlaid on the corresponding area. This means the visual area of a region [does not correspond to the number of elements it contains](https://doi.org/10.1371/journal.pone.0101717), which can mislead viewers into overestimating or underestimating overlaps based on the apparent size of regions.

### Implementation in This Pipeline

The pipeline's `plot_venn()` function (`src/visualization/venn_diagrams.R`) uses the [`venn` R package](https://cran.r-project.org/package=venn) with `ggplot = TRUE` to render Venn diagrams as ggplot2 objects:

```r
plot_venn <- function(sets, filename, output_dir) {
  require(venn)
  p <- venn::venn(
    x = sets, borders = FALSE, ggplot = TRUE,
    zcolor = "style", box = FALSE, plotsize = 10,
    ilcs = 1, sncs = 1, ilabels = "counts"
  ) + ggplot2::ggtitle(paste(filename, "population"))
  save_plot_formats(p, file.path(output_dir, paste0("venn_", filename)),
                    width = 10, height = 8)
}
```

The function accepts a named list of index vectors (one per set), draws fixed-size circles with region counts, and saves output in PNG and SVG formats via the shared `save_plot_formats()` utility.

### Strengths

- **Universally recognized.** The Venn diagram is [one of the most widely understood data visualizations in biology and the sciences generally](https://doi.org/10.1371/journal.pone.0101717), requiring no special training to interpret.
- **Complete intersection display.** Every possible intersection region is shown, even if empty, making it immediately clear which combinations exist and which do not.
- **Simplicity.** For two or three sets, the layout is clean and unambiguous.

### Weaknesses

- **Not area-proportional.** The visual size of each region does not reflect its cardinality, which [can create misleading impressions about the relative sizes of overlaps](https://doi.org/10.1371/journal.pone.0101717).
- **Unreadable beyond three sets.** While [Venn diagrams can theoretically be drawn for any number of sets using non-circular curves](https://doi.org/10.1080/14786448008626877), diagrams with four or more sets become visually complex, with many small, oddly shaped regions that are difficult to label and interpret. [Empirical studies show that readers make significantly more errors interpreting Venn diagrams as the number of sets increases beyond three](https://doi.org/10.1371/journal.pone.0101717).
- **Fixed geometry.** The arrangement of circles is the same regardless of the data, providing no visual signal about which sets are more or less similar.

---

## Euler Diagrams

### Concept and Theory

An Euler diagram is a generalization of the Venn diagram in which the areas of regions are **proportional to the number of elements** they contain. Unlike Venn diagrams, Euler diagrams [omit regions that correspond to empty intersections](https://doi.org/10.1371/journal.pone.0101717), producing a cleaner layout that directly reflects the data structure. If two sets have no overlap, their circles do not intersect; if one set is a subset of another, its circle is drawn entirely within the larger one.

[Micallef and Rodgers (2014)](https://doi.org/10.1371/journal.pone.0101717) conducted a systematic empirical evaluation comparing area-proportional Euler diagrams against non-proportional Venn diagrams and found that [area-proportional diagrams led to significantly faster and more accurate judgments of relative set sizes and overlap magnitudes](https://doi.org/10.1371/journal.pone.0101717). This advantage was particularly pronounced for tasks that required comparing the sizes of two regions, where the visual area encoding provided an [intuitive analog representation that participants could estimate without reading numeric labels](https://doi.org/10.1371/journal.pone.0101717).

### The eulerr Algorithm

This pipeline uses the [`eulerr` R package](https://cran.r-project.org/package=eulerr) by [Larsson (2022)](https://doi.org/10.21105/joss.04986), which computes area-proportional Euler diagrams through a two-stage optimization process:

1. **Initial layout.** The package computes an initial configuration of circles (or ellipses) by finding positions and radii that approximate the target areas using a loss function based on the [sum of squared differences between the actual and ideal region areas](https://doi.org/10.21105/joss.04986).
2. **Nonlinear optimization.** A second stage refines the layout using nonlinear least squares optimization ([Nelder-Mead or NLOPT algorithms](https://doi.org/10.21105/joss.04986)) to minimize residual error. The package reports a `diagError` statistic that quantifies the maximum deviation between any region's drawn area and its target area, along with a `stress` statistic analogous to [stress in multidimensional scaling](https://doi.org/10.21105/joss.04986).

For some input data, a perfect area-proportional solution using circles alone does not exist. When [three or more sets have intersection structures that cannot be represented with circles, `eulerr` can use ellipses](https://doi.org/10.21105/joss.04986), which provide additional degrees of freedom (orientation and eccentricity) to improve the fit.

### Implementation in This Pipeline

The `plot_euler()` function (`src/visualization/venn_diagrams.R`) wraps `eulerr::euler()`:

```r
plot_euler <- function(sets, output_dir, filename, custom_name = NULL) {
  require(eulerr)
  n_sets <- length(sets)
  if (n_sets < 2) stop("At least two sets required")
  if (!is.null(custom_name) && length(custom_name) == n_sets) {
    names(sets) <- custom_name
  }
  fit <- eulerr::euler(sets)
  set_colors <- scales::hue_pal()(n_sets)

  filename_base <- file.path(output_dir, paste0("euler_", filename))

  draw_euler <- function() {
    tryCatch({
      p <- plot(fit, main = paste(filename, "population"),
                fills = set_colors, quantities = TRUE,
                legend = names(sets))
      print(p)
    }, error = function(e) {
      message("Error plotting '", filename, "': ", e$message)
    })
  }

  save_base_plot_formats(filename_base, draw_euler, width = 10, height = 8)
}
```

The function accepts a named list of index vectors, computes the area-proportional layout via `eulerr::euler()`, and renders the diagram with fill colors from `scales::hue_pal()`, region quantity labels, and a legend. Output is saved in both PNG and SVG formats through `save_base_plot_formats()`.

### Strengths

- **Area-proportional.** The visual size of each region [directly encodes the number of elements](https://doi.org/10.1371/journal.pone.0101717), enabling intuitive magnitude comparisons without relying on numeric labels.
- **Omits empty regions.** Intersections that contain no elements are not drawn, reducing visual clutter and [making the non-empty structure of the data immediately apparent](https://doi.org/10.1371/journal.pone.0101717).
- **Effective for 2--4 sets.** Euler diagrams [maintain readability for up to approximately four sets](https://doi.org/10.21105/joss.04986) when the intersection structure is not too complex.
- **Diagnostic statistics.** The `eulerr` package provides `diagError` and `stress` metrics that [quantify how well the drawn diagram matches the target areas](https://doi.org/10.21105/joss.04986), allowing the user to assess whether the visual representation is trustworthy.

### Weaknesses

- **No exact solution guaranteed.** For some combinations of set sizes and overlaps, [no arrangement of circles or even ellipses can perfectly reproduce the target areas](https://doi.org/10.21105/joss.04986). In such cases, the diagram is an approximation, and the `diagError` should be checked.
- **Less familiar.** Readers accustomed to standard Venn diagrams may [misinterpret an Euler diagram by assuming that all possible intersection regions exist, when in fact empty regions are omitted by design](https://doi.org/10.1371/journal.pone.0101717).
- **Scalability limit.** Beyond four or five sets, the optimization problem becomes increasingly constrained and the visual output increasingly difficult to parse, even with ellipses.

---

## UpSet Plots

### Motivation

As the number of sets increases beyond three, both Venn and Euler diagrams become difficult to read. [Lex et al. (2014)](https://doi.org/10.1109/TVCG.2014.2346248) introduced the UpSet plot as a fundamentally different approach to set visualization that [scales effectively to dozens of sets and provides precise quantitative information about every intersection](https://doi.org/10.1109/TVCG.2014.2346248). Rather than encoding set membership through overlapping geometric shapes, UpSet uses a matrix layout combined with bar charts.

### How UpSet Plots Work

An UpSet plot consists of three coordinated components:

1. **Set size bars (left).** A horizontal bar chart showing the total number of elements in each set, providing context for interpreting intersection sizes.

2. **Intersection matrix (bottom).** A binary matrix where rows correspond to sets and columns correspond to specific intersections. A filled dot in the matrix indicates that the corresponding set participates in that intersection; connected dots (linked by a vertical line) indicate a multi-set intersection. An empty dot means the set does not participate.

3. **Intersection size bars (top).** A vertical bar chart directly above the matrix, where each bar shows the number of elements in the corresponding intersection. Bars are [typically sorted by frequency (decreasing)](https://doi.org/10.1109/TVCG.2014.2346248), making it immediately obvious which intersections are largest.

This encoding separates the identity of an intersection (which sets participate, shown in the matrix) from its magnitude (how many elements, shown in the bar), enabling [precise quantitative comparisons that are difficult or impossible with area-based encodings](https://doi.org/10.1109/TVCG.2014.2346248). In their evaluation, [Lex et al. demonstrated that UpSet plots enabled significantly more accurate identification of intersections and their sizes than either Venn or Euler diagrams for datasets with more than three sets](https://doi.org/10.1109/TVCG.2014.2346248).

### Implementation in This Pipeline

The pipeline provides two UpSet plot functions in `src/visualization/upset_plots.R`, both using the [`UpSetR` R package](https://cran.r-project.org/package=UpSetR):

**From binary membership data:**

```r
plot_upset <- function(data, set_columns, output_path,
                       title = "UpSet Plot", min_size = 1) {
  require(UpSetR)
  filename_base <- sub("\\.[^.]+$", "", output_path)
  draw_upset <- function() {
    tryCatch({
      UpSetR::upset(
        as.data.frame(data),
        sets = set_columns,
        order.by = "freq",
        decreasing = TRUE,
        mainbar.y.label = "Intersection Size",
        sets.x.label = "Set Size",
        text.scale = 1.2,
        main.bar.color = "steelblue",
        sets.bar.color = "grey50",
        mb.ratio = c(0.6, 0.4)
      )
      title(main = title, line = 2.5)
    }, error = function(e) {
      message("UpSet plot error: ", e$message)
    })
  }
  save_base_plot_formats(filename_base, draw_upset, width = 16, height = 10)
}
```

**From a named list of sets:**

```r
plot_upset_from_list <- function(sets, output_path, title = "UpSet Plot") {
  require(UpSetR)
  filename_base <- sub("\\.[^.]+$", "", output_path)
  draw_upset <- function() {
    tryCatch({
      UpSetR::upset(
        UpSetR::fromList(sets),
        order.by = "freq",
        decreasing = TRUE,
        mainbar.y.label = "Intersection Size",
        sets.x.label = "Set Size",
        text.scale = 1.2
      )
      title(main = title, line = 2.5)
    }, error = function(e) {
      message("UpSet plot error: ", e$message)
    })
  }
  save_base_plot_formats(filename_base, draw_upset, width = 16, height = 10)
}
```

Both functions order intersections by frequency (`order.by = "freq"`) so that the most populated intersections appear first (leftmost bars). The `plot_upset()` variant accepts a data frame with binary (0/1) columns for each set, while `plot_upset_from_list()` converts a named list of element vectors using `UpSetR::fromList()`.

### Strengths

- **Scales to many sets.** UpSet plots [remain readable and interpretable with 10, 20, or even 40+ sets](https://doi.org/10.1109/TVCG.2014.2346248), whereas Venn and Euler diagrams become unusable beyond 3--4 sets.
- **Precise quantification.** Intersection sizes are [encoded as bar heights on a common linear scale](https://doi.org/10.1109/TVCG.2014.2346248), enabling accurate magnitude comparisons. Readers do not need to estimate areas of irregular shapes.
- **Sortable and filterable.** Intersections can be sorted by frequency, degree (number of participating sets), or custom criteria. Small intersections below a threshold (`min_size`) can be hidden to reduce clutter.
- **Explicit intersection identity.** The matrix encoding makes it [unambiguous which combination of sets defines each intersection](https://doi.org/10.1109/TVCG.2014.2346248), unlike Venn/Euler diagrams where complex region shapes can be difficult to associate with specific set combinations.

### Weaknesses

- **Less intuitive for simple cases.** For two or three sets, a Venn or Euler diagram provides a more immediate, spatially intuitive overview. The UpSet matrix encoding carries an overhead of visual decoding that is unnecessary when the number of sets is small.
- **No area-proportional encoding.** UpSet does not visually show containment or subset relationships in the spatial sense. Two sets with 90% overlap look the same in the matrix as two sets with 10% overlap; only the bar heights differ.
- **Requires more plot space.** The combined matrix-plus-bars layout needs more horizontal and vertical space than a compact Venn or Euler diagram, especially when many intersections are shown.
- **Base R graphics.** The `UpSetR` package uses base R plotting rather than ggplot2, which makes integration with ggplot2-based figure composition (e.g., `patchwork`) more involved. The pipeline wraps UpSet rendering in `save_base_plot_formats()` to handle this.

---

## When to Use Each Method

The three set visualization methods serve complementary roles. The choice depends on the number of sets, the importance of area proportionality, and the analytical question.

### Venn Diagrams: 2--3 Sets, Quick Overview

Use `plot_venn()` when:

- You are comparing **exactly two or three** gating strategies, workspaces, or analysts.
- You want a **familiar, instantly recognizable** figure that any collaborator can interpret.
- The **absolute counts** in each region (displayed as labels) matter more than visual proportionality.
- You are producing a figure for a presentation or report where audience familiarity with the diagram type is important.

Do not use Venn diagrams when you have four or more sets, or when the [visual proportionality of regions is important for correct interpretation](https://doi.org/10.1371/journal.pone.0101717).

### Euler Diagrams: 2--4 Sets, Proportional Representation

Use `plot_euler()` when:

- You have **two to four** sets and want the **visual size of each region to reflect its cardinality**.
- The relative magnitude of overlaps is the key message (e.g., showing that 80% of cells are classified identically by two analysts).
- You want to convey subset relationships visually (e.g., one gating strategy classifies a strict subset of the cells classified by another).
- You are willing to check `diagError` to confirm the diagram is a faithful representation.

Do not use Euler diagrams when the intersection structure [cannot be well approximated by circles or ellipses](https://doi.org/10.21105/joss.04986) (check `diagError`), or when you have more than four or five sets.

### UpSet Plots: 4+ Sets, Complex Intersections

Use `plot_upset()` or `plot_upset_from_list()` when:

- You have **four or more** sets to compare (e.g., comparing gating from five or more FlowJo workspaces).
- You need **precise quantification** of intersection sizes and want to sort or filter intersections.
- You want to identify which specific combination of analysts agree or disagree on cell classification.
- The question is "what are the largest areas of agreement and disagreement?" rather than "what is the overall overlap structure?"

[Lex et al. (2014)](https://doi.org/10.1109/TVCG.2014.2346248) specifically recommended UpSet as the default visualization for any set analysis involving more than three sets.

### Summary Decision Table

| Criterion | Venn | Euler | UpSet |
|-----------|------|-------|-------|
| Number of sets | 2--3 | 2--4 | 4+ (any) |
| Area proportionality | No | Yes | No (bar charts) |
| Shows empty intersections | Yes | No | Optional |
| Precise quantification | Labels only | Labels + area | Bar chart (linear scale) |
| Audience familiarity | High | Medium | Lower |
| Scalability | Poor (>3 sets) | Moderate (>4 sets) | Excellent |
| Best for | Quick overview | Proportional comparison | Complex multi-set analysis |

---

## Cytometry Applications

### Comparing Gating Strategies Across Analysts (Pipeline D)

The primary use case for set visualization in this pipeline is **inter-analyst gating comparison**. In a typical flow cytometry workflow, [different analysts may apply gating strategies that produce substantially different cell population assignments for the same data](https://doi.org/10.1002/cyto.a.22625), even when following the same protocol. This inter-operator variability is a [well-documented challenge in clinical and research cytometry](https://doi.org/10.1002/cyto.a.22764).

Pipeline D (`pipelines/pipeline_D_flowjo_comparison.R`) loads multiple FlowJo workspaces containing gates applied to the same FCS file, extracts per-cell population labels from each workspace, and then generates Venn/Euler/UpSet diagrams for each population. For each cell population (e.g., "CD4+ T cells"), the diagram shows how many cells are classified as that population by each analyst, and where the analysts agree or disagree.

The `compare_gatings()` function (`src/metadata_qc/compare_gatings.R`) orchestrates this workflow:

1. Extracts label vectors from each FlowJo workspace for a target FCS file.
2. Translates raw gate paths to standardized labels using a user-provided dictionary.
3. Optionally merges CSV-derived labels (e.g., from automated clustering) for cross-method comparison.
4. Iterates over all unique population labels and generates a Venn or Euler diagram for each.

### Comparing Automated vs. Manual Gating

When both manual gates (from FlowJo) and automated cluster assignments (from Leiden clustering) are available for the same sample, set visualization reveals the concordance between the two approaches. This is valuable for [validating unsupervised clustering results against expert-defined populations](https://doi.org/10.1002/cyto.a.22764). The `idx_labels` parameter of `compare_gatings()` allows merging CSV-derived automated labels alongside FlowJo-derived manual labels for direct visual comparison.

### Comparing Cluster Membership Across Methods

Set visualization can also be applied to compare the output of different clustering algorithms or parameter settings on the same data. For example, comparing Leiden clustering at two resolution values, or comparing Leiden against FlowSOM metaclustering, using the cluster assignment of each cell as set membership.

### Quantifying Agreement with Fleiss' Kappa

Set visualization shows **where** analysts agree and disagree, but does not directly quantify the **overall level of agreement**. For a statistical measure of inter-rater reliability, the pipeline computes [Fleiss' kappa](fleiss-kappa.md) via `compute_fleiss_kappa()` in `src/metadata_qc/overlap_statistics.R`. Fleiss' kappa [extends Cohen's kappa to multiple raters and provides a single statistic summarizing agreement above chance](https://doi.org/10.1037/h0031619), complementing the visual per-population breakdown provided by set diagrams.

The recommended workflow is to:

1. Generate set visualizations (Venn, Euler, or UpSet) to **identify which populations have the most disagreement**.
2. Compute Fleiss' kappa to **quantify the overall level of agreement** and per-population agreement.
3. Use both results together to report inter-analyst reproducibility.

---

## How to Use in the Pipeline

### Euler Diagram (Area-Proportional)

For two to four sets where visual proportionality matters:

```r
source("src/visualization/venn_diagrams.R")

# Prepare sets as a named list of index vectors
sets <- list(
  "Analyst_1" = which(labels_analyst1 == "CD4+ T cells"),
  "Analyst_2" = which(labels_analyst2 == "CD4+ T cells"),
  "Analyst_3" = which(labels_analyst3 == "CD4+ T cells")
)

plot_euler(
  sets       = sets,
  output_dir = "output/gating_comparison/",
  filename   = "CD4_T_cells",
  custom_name = c("Analyst 1", "Analyst 2", "Analyst 3")
)
```

### Venn Diagram (Fixed-Size)

For two to three sets when a familiar, non-proportional layout is preferred:

```r
source("src/visualization/venn_diagrams.R")

plot_venn(
  sets       = sets,   # Same named list as above
  filename   = "CD4_T_cells",
  output_dir = "output/gating_comparison/"
)
```

### Choosing Between Euler and Venn via `plot_and_save_venn()`

The convenience function `plot_and_save_venn()` accepts a binary label matrix and a `scale_size` flag:

```r
source("src/visualization/venn_diagrams.R")

# label_lists: data frame where each column is a rater's labels for each cell
plot_and_save_venn(
  lists       = label_lists,
  population  = "CD4+ T cells",
  filename    = "CD4_comparison",
  output_dir  = "output/gating_comparison/",
  scale_size  = TRUE,       # TRUE = Euler, FALSE = Venn
  custom_name = c("WS1", "WS2", "WS3")
)
```

### UpSet Plot (For Many Workspaces)

For four or more sets, or when precise intersection quantification is needed:

```r
source("src/visualization/upset_plots.R")

# From a named list of element vectors
sets <- list(
  "Analyst_1" = cell_ids_analyst1,
  "Analyst_2" = cell_ids_analyst2,
  "Analyst_3" = cell_ids_analyst3,
  "Analyst_4" = cell_ids_analyst4,
  "Leiden"    = cell_ids_leiden
)

plot_upset_from_list(
  sets        = sets,
  output_path = "output/gating_comparison/upset_CD4.png",
  title       = "CD4+ T cell Classification Agreement"
)
```

```r
# From a binary membership data frame
plot_upset(
  data        = binary_membership_df,
  set_columns = c("Analyst_1", "Analyst_2", "Analyst_3", "Analyst_4", "Leiden"),
  output_path = "output/gating_comparison/upset_all.png",
  title       = "Multi-Analyst Gating Comparison"
)
```

### Full Gating Comparison Workflow (Pipeline D)

The `compare_gatings()` function automates the entire comparison workflow across FlowJo workspaces:

```r
source("src/metadata_qc/compare_gatings.R")

compare_gatings(
  pop_strings_list = pop_strings,        # Named list from load_flowjo_wsp()
  target_fcs       = "Sample_001.fcs",   # FCS file to compare across workspaces
  label_dict       = label_dictionary,   # Gate path -> simplified label mapping
  custom_name      = c("Analyst A", "Analyst B", "Analyst C"),
  output_dir       = "output/gating_comparison/",
  idx_labels       = csv_labels,         # Optional: CSV-derived labels for comparison
  scale_size       = TRUE                # TRUE = Euler, FALSE = Venn
)
```

This function:

1. Validates that the target FCS file exists in all workspaces.
2. Extracts per-cell labels from each workspace for that FCS file.
3. Translates gate paths to standardized population names using `label_dict`.
4. Optionally merges CSV-derived automated labels (`idx_labels`).
5. Saves label lists as CSV files for downstream analysis.
6. Generates one Venn or Euler diagram per unique population.

To add UpSet plots to this workflow for populations with many raters, source `upset_plots.R` separately:

```r
source("src/visualization/upset_plots.R")

# After compare_gatings() has produced label_lists
label_lists <- data.table::fread("output/gating_comparison/label_lists_Sample_001.fcs.csv")

# Build binary membership for a specific population
pop <- "CD4+ T cells"
binary_df <- as.data.frame(lapply(label_lists, function(col) as.integer(col == pop)))

plot_upset(
  data        = binary_df,
  set_columns = names(binary_df),
  output_path = "output/gating_comparison/upset_CD4.png",
  title       = paste("UpSet:", pop)
)
```

### Combining with Fleiss' Kappa

After visual comparison, quantify agreement statistically:

```r
source("src/metadata_qc/overlap_statistics.R")

label_lists <- data.table::fread("output/gating_comparison/label_lists_Sample_001.fcs.csv")

kappa_results <- compute_fleiss_kappa(
  label_lists = label_lists,
  output_dir  = "output/gating_comparison/"
)

# Inspect overall kappa
print(kappa_results[["Full data"]])
```

---

## Limitations

### Venn Diagram Limitations

- **Maximum effective sets: 3.** While the `venn` package can draw diagrams for higher set counts using non-circular shapes, the resulting figures are [visually complex and error-prone to interpret beyond three sets](https://doi.org/10.1371/journal.pone.0101717).
- **Misleading area perception.** Because [circle sizes are fixed and not proportional to set cardinalities](https://doi.org/10.1371/journal.pone.0101717), viewers may unconsciously infer proportionality from area, leading to systematic misinterpretation.
- **No sorting or filtering.** All regions are always shown at fixed positions, with no ability to highlight the most important intersections.

### Euler Diagram Limitations

- **Approximate solutions.** For certain intersection structures, [no configuration of circles or ellipses can exactly represent the target areas](https://doi.org/10.21105/joss.04986). The `eulerr` package minimizes the residual error but cannot eliminate it when the mathematical constraints are incompatible. Always check `diagError` and `stress` in the fitted object.
- **Ambiguity with missing regions.** When a possible intersection is empty, the corresponding region is absent from the diagram. Readers unfamiliar with Euler diagrams may [incorrectly assume the diagram is a Venn diagram and misinterpret the absent region as occupied](https://doi.org/10.1371/journal.pone.0101717).
- **Computation time.** For four or more sets, the nonlinear optimization can take several seconds, though this is negligible compared to the upstream data processing in cytometry workflows.

### UpSet Plot Limitations

- **Less intuitive for novice audiences.** The matrix-plus-bars encoding requires [a brief learning period that is not needed for Venn diagrams](https://doi.org/10.1109/TVCG.2014.2346248). Include a legend or caption explaining the plot structure when presenting to non-specialist audiences.
- **Large number of intersections.** With n sets, there are 2^n - 1 possible intersections. For n = 10, this is 1,023 intersections, most of which may be empty or very small. The `min_size` parameter in `plot_upset()` filters out intersections below a threshold, but careful threshold selection is needed to balance completeness against visual clutter.
- **No spatial overlap encoding.** UpSet plots [do not visually convey containment or subset relationships between sets](https://doi.org/10.1109/TVCG.2014.2346248) in the way that nested circles in an Euler diagram do. The plot shows that two sets share elements, but does not show this through spatial nesting.

### General Considerations for Cytometry

- **Cell-level vs. population-level comparison.** All three methods operate on sets of cells. They show whether the same cells are classified into the same population by different methods or analysts, but they do not capture finer-grained differences in how populations are defined (e.g., differences in gate boundary placement that do not change the majority of assignments).
- **Multiple populations.** Each diagram compares set membership for a single population. When comparing gating strategies that define many populations, the pipeline generates one diagram per population, which can produce a large number of output files. Review key populations of interest first.
- **Complementary to statistical agreement.** Set visualizations are descriptive. For inferential statements about inter-rater agreement, pair them with [Fleiss' kappa](fleiss-kappa.md) or other agreement statistics from `src/metadata_qc/overlap_statistics.R`.

---

## Related Documentation

- [Fleiss' Kappa](fleiss-kappa.md) -- Statistical measure of inter-rater agreement, complementary to visual set comparison
- [Leiden Clustering](leiden-clustering.md) -- Unsupervised clustering method whose output can be compared via set visualization
- [UMAP](umap.md) -- Dimensionality reduction for visualizing cell populations
- [Pipeline D: FlowJo Comparison](../pipelines/pipeline-D-flowjo-comparison.md) -- Full pipeline using set visualization for gating comparison
- [Metadata QC Subproject](../subprojects/metadata-qc.md) -- Parent documentation for the gating comparison and overlap statistics modules
