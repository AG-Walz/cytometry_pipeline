# Gate-Anchored UMAP Layout

## Summary

The gate-anchored layout builds the UMAP features from the FlowJo gates instead of from the raw marker values. For every gate below the chosen population it reads which marker the gate decides and where the cut lies in that sample, and records the gater's call for every cell. Each decided marker is then re-expressed as a bounded continuous term, centred on the sample's own cells of the same combination of calls, plus a fixed jump whose sign is the call. Cells called positive and cells called negative on any marker are then guaranteed to be apart in feature space, so every combination of calls (for example `CD4+ | IFNy+ TNFa+`) can form its own island, while the continuous term keeps a natural cloud shape inside each island. Rare combinations are kept by stratified sampling, the embedding is fitted on a set balanced across combinations, and every frequency shown is computed from the full, unsampled counts. It is run with `scripts/umap_from_flowjo.R --layout gates`; see [the script guide](../scripts/umap-from-flowjo.md#gate-separated-layout) for usage.

Because the gates decide the layout, the separation of islands restates the gating. It is not independent evidence that the populations differ biologically (see [Reading the Plots](#reading-the-plots)).

---

## Motivation

A plain [UMAP](umap.md) of FlowJo values embeds every marker on the same footing. A cytokine-positive T cell differs from its cytokine-negative neighbours in one channel, often by only a few units on the 0 to 256 FlowJo scale, while it shares every other channel with them. In the k-nearest-neighbour graph it is surrounded by negatives, so it is easily drawn among them, and the polyfunctional cells a gater has carefully defined need not form islands of their own.

Two further problems come from real cohorts:

- **Gates move between acquisitions.** Analysts follow the detector voltages, so the same gate sits at a different position in each visit. In the MTB cohort every tube of a visit shares one set of gate coordinates, but the cytokine cuts differ between visits and between the CD4 and CD8 parents of the same visit. A fixed global threshold would disagree with the gater in many samples.
- **The interesting cells are rare.** Cytokine-positive cells are typically well under 10 % of a lineage. Uniform downsampling keeps a handful per tube, too few to form a visible island.

The gate-anchored layout addresses all three: it reads the gater's own decisions per sample, turns each decision into a guaranteed gap in feature space, and samples so that rare combinations keep enough cells.

---

## The Method

### Decisions Read from the Gates

The values come from FlowJo's working space: compensated and transformed with FlowJo's biexponential transform (a member of the [logicle family, Parks et al. 2006](https://doi.org/10.1002/cyto.a.20258)) onto a 0 to 256 channel scale. Gates are stored on the same scale, so a gate's bounds can be compared with the cell values directly.

For each sample, every non-Boolean gate strictly below the population is examined. A gate dimension counts as a **decision** on marker *m* when its channel maps to a named marker (scatter, `Time` and unnamed channels never do, and marker names that match `--viability_pattern`, by default the common viability dyes, are never decisions) and the gate excludes a real share of its parent along that axis:

- **Rectangle gates.** A finite bound binds when on its own it excludes at least 0.5 % of the parent's cells. A gate that binds only at its lower bound is a positive decision, with the cut at that bound; one that binds only at its upper bound is a negative decision. When both bounds bind (a band such as "CD4 dim") the dimension makes no threshold decision, unless the smaller side excludes under 5 % of the parent: that is a far edge clipping a tail (a box not dragged to the end of the axis), and the larger side is the decision. Bounds that exclude almost nothing, such as the top of a cytokine gate drawn to the end of its axis, drop out under the same rule. The FSC-A range of a cytokine-versus-FSC gate never counts, since FSC-A is a scatter channel.
- **Polygon, ellipse and other gates.** These are decided from membership. With at least 10 members and 10 non-members in the parent, the area under the ROC curve of the marker for members against non-members is computed; an AUC of at least 0.8 makes a positive decision and one of at most 0.2 a negative decision. The cut is the one-dimensional threshold that best reproduces membership (maximum Youden's J). Otherwise the gate makes no decision in that sample.
- **Boolean gates** are ignored, since they are derived from the gates they combine.

**Ancestors decide first.** A gate never decides a marker that one of its ancestor gates below the population already decides. A cytokine box drawn on CD4 x IFNy under `CD4+` may clip the CD4-dim tail; that edge trims cells inside the side the `CD4+` gate chose and is not a second CD4 threshold. Gates in sibling subtrees are not ancestors: an IFNy gate directly under the population and another under `CD4+` are both kept, and the call rules below let the deeper one speak for CD4+ cells.

Each decision is recorded in `gate_decisions.csv` with its gate, parent, depth below the population, channel, marker, polarity, cut, method, the share of the parent it excludes, the parent and member cell counts, and `n_axes`: the number of marker axes the gate binds on after the ancestor rule, bands included. Bands themselves have no row.

### Call Rules

For cell *i* and marker *m*, the **applicable** decisions are those whose gate's parent contains *i*. Of these, only the deepest count: when a marker is gated both under a lineage-independent parent and again under `CD4+`, the gate under `CD4+` is the decision for a CD4+ cell, and the shallower gate only speaks for cells that no deeper gate covers. Decisions at the same depth are used together.

- If *i* is inside one or more of the deciding gates, it is positive when any of them decides *m* positively, otherwise negative.
- If *i* is inside none of them and every one of them binds on *m* alone (`n_axes` of 1), being outside is the decision: negative when every one keeps cells above its cut (a cell outside the IFNg+ gate), positive when every one keeps cells below its cut.
- Otherwise the cell's own value decides: positive when `x >= cut`, with `cut` the median cut of the deciding gates. This covers gates that also bind on another marker (a cell can be outside a CD4+CD8- box because of its CD8, which says nothing about its CD4) and gates of opposite polarity.
- With no applicable decision, *m* is **not decided** for *i* (for example, cytokines for CD4-CD8- cells when the cytokine gates sit only under CD4 and CD8).

The call is the gater's decision, membership included. A cell outside a cytokine gate because of a dimension that decides nothing, such as the FSC-A edge of a cytokine-versus-FSC box, is negative even when its value on *m* lies above the cut.

### Cuts and Anchors

All of the above is done per sample, so each decision's cut is the position of that gate **in that sample**. The cut enters the calls only where they are taken from the value (above).

The continuous part of the feature is measured from an **anchor**, not from the cut: for each marker, the median value over the sample's own cells that share the cell's combination of calls, taken over all of the sample's cells before sampling. In each tube, the CD4+ IFNy+ cells are centred on their own median IFNy, CD4 and so on, and the CD8+ cells on theirs.

A voltage shift between visits moves the cells and their median together, so anchoring removes it as centring on the cut would. Unlike the cut, the anchor does not depend on where the gater drew the line in that tube, and gaters do not always draw the same line on identical cells. Where they did not, centring on the cut split one population into sample-pure halves inside its island, which looks like a per-sample effect but is a gating artefact ([plan section 11.1](../plans/umap_gate_layout_plan.md#111-after-the-adversarial-review-same-day)). Centring per combination rather than per side of a gate also keeps the anchor independent of which other populations share that side: CD8 T cells and CD4-CD8- cells are both CD4-negative, in proportions that differ between tubes, so a median over the whole negative side would put the CD8 T cells of different tubes at different CD4 positions.

What anchoring gives up is the location of a whole combination between tubes: if the IFNy+ cells of one tube are all brighter than those of another, both tubes' IFNy+ cells still sit in the same place in the island. In the MTB cohort such differences between tubes are mostly instrument setting; the marker plots keep the raw values.

A marker a sample has no decision for at all (a gate missing from its workspace) is centred in the same way, per combination, on the sample's kept cells, and its reported cut is the marker's median cut over the whole cohort.

### The Feature Transform

Each decided marker *m* becomes

```
f = tanh((x - anchor) / (softness * sigma)) + gap_m * s

s     = +1 for a positive call, -1 for a negative or undecided call
gap_m = gate_gap * gate_tier_ratio ^ (deepest tier - tier of m)
```

where `x` is the cell's FlowJo value, `anchor` the anchor of its sample and combination, `sigma` the marker's scale (next section), `softness` is `--gate_softness`, and `gap_m` is `--gate_gap` for the deepest tier, multiplied by `--gate_tier_ratio` for each tier above it. With the defaults (2 and 2) the MTB cytokines get a jump of 2 and CD4 and CD8 a jump of 4.

**Why this guarantees separation.** The within-side term `tanh(...)` lies in (-1, 1) whatever the value. Cells called positive on *m* therefore have `f` in (gap_m - 1, gap_m + 1) and the others in (-gap_m - 1, -gap_m + 1). The smallest possible distance between the two on that coordinate is greater than

```
(gap_m - 1) - (-gap_m + 1) = 2 * (gap_m - 1)
```

so with `gap_m > 1` the two sides of every gate can never overlap or interleave, whatever voltages the samples were acquired with and wherever their gates were drawn. `--gate_tier_ratio` must be at least 1, so every `gap_m` is at least `--gate_gap`: two cells whose calls on some marker differ as positive against negative or undecided are more than `2 * (gate_gap - 1)` apart in Euclidean distance, 2 at the defaults. This holds for the markers that are features: a decided marker that some sample does not measure, or that `--markers` leaves out, still labels the cells but adds no coordinate.

Negative and undecided calls share a side, so the bound does not separate them: CD4+ cells with no cytokine decision (label `CD4+`) and CD4+ cells called negative for every cytokine (`CD4+ | none`) can share a region. Cells with no decision on any marker are handled separately (below).

The call, not the value, picks the side. A cell called negative although its value lies above the cut (outside a gate because of a non-marker edge, say) can have a tanh term near +1 but stays on the negative side; the bound covers exactly these disagreements.

The bound is a statement about feature space. Whether UMAP draws two regions as visibly separate islands also depends on how large the gap is relative to the spread inside each island and on the UMAP parameters, which is what the [separation metrics](#separation-metrics) measure.

**Why a larger jump higher up the tree.** With the lineage jump larger than the cytokine jump, the first split of the gating tree is the largest distance in feature space. The islands of one lineage then gather in one region, with the cytokine combinations arranged inside it from cytokine-negative to polyfunctional, instead of being scattered over the whole embedding.

**Why tanh.** A plain shift keeps the gate's difference as small as it was. Binarising the calls throws away all within-island structure and leaves UMAP a few thousand identical points per combination. A bounded, smooth term keeps the order of cells by distance from the anchor (a bright IFNy+ cell sits on the far side of the IFNy+ island from a dim one) while capping how far any cell can stray. `softness` sets how much of that spread survives: a large value keeps the tanh term small and nearly linear, so cells crowd around the island's centre; a small one saturates it, so most cells sit near the island's edges at ±1. In the tuning runs ([plan section 11.1](../plans/umap_gate_layout_plan.md#111-after-the-adversarial-review-same-day)) softness 2 left the polyfunctional islands striped, with several axes saturated at both ends, and 4 gave smooth islands with the same separation.

**Cells with no decision.** When some cells of a run have no decision on any marker (with a population above the live gate, the cells that gate excluded), they would sit on the negative side of every marker, on top of the all-negative combination. They get one more coordinate, `no_gate_decision`, equal to the largest `gap_m` for them and to its negative for every other cell, and the label `(no decision)`. The log reports how many cells this concerns.

Markers that no gate decides but that the user asks for with `--markers` enter as `tanh((x - median) / (softness * MAD))` without a jump, so they shape the islands without splitting them.

### Scale: the Negative Population's MAD

`sigma` is computed once per marker after all samples are read:

```
sigma_m = 1.4826 * weighted MAD of (x - anchor) over the cells called negative for m
```

pooled across samples, with each sampled cell weighted back to the full counts, and floored at 1 % of the channel range (2.56 units on the 0 to 256 scale). When fewer than 10 cells are negative, all cells with a finite value are used instead.

The negative population is the bulk of most channels, and its spread is mostly instrument noise and compensation spread. It therefore defines what "far from the centre" means for that channel in a way that does not depend on how many positive cells a cohort happens to have, or on how heterogeneous they are. The factor 1.4826 makes the MAD equal the standard deviation for normally distributed values; the floor prevents a nearly constant channel from blowing tiny differences up to the full tanh range.

### Combination Labels and Tiers

A marker's **tier** is the depth below the population of the gates that decide it (the rounded median depth when it varies between samples). In the MTB tree under `Live cells`, CD4 and CD8 are decided by the quadrant gates directly below the population (tier 1), and the cytokines by gates under CD4 and CD8 (tier 2).

A combination label lists, per tier, the positive markers in panel order, with tiers joined by ` | `:

| Calls | Label |
|-------|-------|
| CD4+, CD8-, IFNy+, TNFa+, others negative | `CD4+ \| IFNy+ TNFa+` |
| CD8+, CD4-, no cytokine positive | `CD8+ \| none` |
| CD4-, CD8-, cytokines undecided | `CD4- CD8-` |
| nothing decided | `(no decision)` |

A tier with no positive marker reads `none`, except the first tier, which spells its markers out as negative. Tiers in which every marker is undecided are left out. Labels are built from marker names, not gate names, so gates spelled `IFNy` in one workspace and `INFy` in another give one label once the marker names agree (see `--marker_aliases`).

Combinations with fewer than `--min_combination_cells` cells in total (counted with the sampling weights, so this is the full count over all included samples) are pooled as `other (rare)` in the colours, fit balancing, supervision, separation metrics, `combination_counts.csv` and the frequency and change tables. Pooling changes only the label: the cells keep their own features and stay in their own region of feature space, and their calls stay in `cells.csv.gz`. `--min_combination_cells 0` turns pooling off.

**Missing decisions.** Within each first-tier group (all cells with the same tier-1 calls, for example all CD4+ cells), a deeper marker that is decided for more than half of the group's cells is set negative for the group's cells that have no decision on it; the log reports the fills and `run_info.json` lists them per sample. The share is of cells, counted with the sampling weights, not of kept cells or of tubes. This covers a sample whose workspace lacks one gate that the others have; without it those cells would form a spurious combination of their own.

---

## Sampling and Embedding

### Counts First

Combination counts per sample are taken from **all** cells of the population before any sampling and written to `combination_counts.csv` (sample, combination, cells, population cells), with a zero for every combination a sample lacks. Every frequency shown later, in plot legends, `combination_frequencies_<col>.csv` and `combination_change_<col>.csv`, is computed from these counts.

### Stratified Sampling and Weights

Per sample, cells are grouped by their exact call pattern and at most `--cells_per_combination` cells of each pattern are kept at random; patterns with fewer cells are kept whole. Each kept cell carries a weight

```
weight = cells of its pattern in that sample / cells of that pattern kept
```

so a weighted sum over the kept cells reproduces the full count of every pattern in every sample. A sample in which no gate decides anything has a single pattern and is sampled as such, not dropped. This is the [stratified downsampling](downsampling.md#stratified-downsampling) idea with the gater's calls as strata instead of clusters.

The sampling is deliberately not proportional: a sample with 5,000 cytokine-negative CD4 cells and 40 polyfunctional ones keeps 2,000 and 40 (at the default cap). Any frequency read off the plotted points would therefore overstate the rare combinations, which is why frequencies come from the full counts and why the weights are carried wherever frequencies are rebuilt from plotted cells.

### Balanced Fit and Projection

The UMAP is fitted on a set in which no combination label has more than `--umap_subsample` / (number of labels) cells (100,000 / number of labels when `--umap_subsample` is 0 in this layout). Smaller labels are fitted whole, and the share they leave unused is not handed to the larger ones, so the fit set can be much smaller than `--umap_subsample`, and a run with fewer cells than the budget is balanced too. `--umap_init auto` chooses spectral or PCA initialisation by the number of fitted cells.

The remaining sampled cells are then placed onto that embedding with `uwot::umap_transform`. Its neighbour search runs on every plotting core; it finds the same neighbours at any thread count and its optimiser runs serially, so the coordinates do not depend on the thread count.

UMAP tends to give a population room in proportion to how many points it has in the fit. Balancing gives every combination comparable room, which is what produces evenly sized, readable islands; projection then puts the abundant cells into the islands their features already point to. The consequence for reading the plot is that **island area says nothing about frequency**.

### Supervised UMAP

`--supervise w` with `w > 0` passes the combination label to `uwot::umap` as the target `y` with `target_weight = w`. In uwot's terms, [a target weight of 0 weights entirely on the data and 1 entirely on the target](https://jlmelville.github.io/uwot/reference/umap.html), with 0.5 balancing the two. The default `0` fits without a target; the gap transform alone already separates combinations in feature space.

Supervision shapes only the fitted cells. [`umap_transform` takes no target](https://jlmelville.github.io/uwot/reference/umap_transform.html), so projected cells are placed from their features alone.

**Interpretation risk.** With supervision the embedding is pulled towards the labels by construction: two combinations land apart because they were labelled differently, whether or not their marker values differ much. Supervision makes a tidier picture of the gating, not a test of it.

---

## Separation Metrics

Each run scores the embedding against the combination labels and writes `separation_metrics.csv` (this needs the `RANN` or `FNN` package; without either, the metrics are skipped with a message). By default the cells scored are at most 50,000, balanced across labels as for the fit, and neighbours are the 15 nearest cells in the 2-D embedding.

| Column | Meaning |
|--------|---------|
| `label` | The combination, or `(macro average)` for the first row |
| `cells_scored` | Cells of the label among the scored cells, not the run's total |
| `knn_purity` | Share of a cell's 15 nearest neighbours that carry its label, averaged over the label's cells |
| `isolated` | Share of the label's cells whose 15 nearest neighbours all carry its label |
| `silhouette` | Simplified silhouette: `(b - a) / max(a, b)`, with `a` the distance to the cell's own label centroid and `b` the distance to the nearest other centroid; 1 is far from every other island, values near 0 or below mean overlap |
| `group_mixing` | Within a label, the share of same-label neighbours that come from a different workspace, divided by what random mixing among that label's workspaces would give: 1 is as mixed as chance, 0 means every neighbour comes from the cell's own workspace. NA for a label found in one workspace only |
| `sample_mixing` | The same over samples: 0 means each sample sits apart inside the island. NA for a label found in one sample only |

The first row, `(macro average)`, averages each column over labels so that every combination counts equally, whatever its size. High purity and isolation say the layout keeps combinations apart; given the construction, that is expected and mainly confirms that the gap and UMAP parameters are adequate. The two mixing columns are the more informative: a low value means the same combination from different workspaces or samples lands in different places, which points at residual acquisition or gating differences, or at a real per-sample shift. The log lists every island of at least 200 scored cells with `sample_mixing` or `group_mixing` below 0.5; look at those before reading a split or fold-change plot of them as biology. It also lists every island with `isolated` below 0.9: islands that are apart in feature space can still be drawn touching, since UMAP places small disconnected pieces of the graph with little to keep them apart, and two touching colours are then not one population.

---

## Split Plots and Fold-Change Maps

For each split column (`--split_by`), levels are in natural order unless `--level_order` sets one, and the first level is the reference unless `--split_reference` picks another.

- **Frequencies.** For each sample, the percentage of the population in each combination, from `combination_counts.csv`. The mean over a level's samples is the level's frequency of that combination, so every sample counts equally: a tube with 200,000 cells does not outweigh one with 600. The means, their SD and the number of samples are in `combination_frequencies_<col>.csv`, and the legends of the split plots show the means.
- **Split plots** draw a level's plotted cells, coloured by combination, over a grey silhouette of the whole embedding, all on the shared coordinates.
- **Fold-change maps.** For each level except the reference, every combination gets

  ```
  log2_fold_change = log2((mean_pct_level + p) / (mean_pct_reference + p))
  ```

  where `p` is one cell in the median sample of the two levels, as a percentage (100 / the median population size of their samples). It gives a combination absent at one level a finite value and keeps changes of well under one cell per sample pale. A combination absent at both levels has no change and is drawn grey. Each island of the embedding is filled with its combination's colour on a diverging scale that stops at 8-fold either way (±3 in log2), inside a grey outline of the embedding. The numbers are in `combination_change_<col>.csv`: `column`, `level`, `reference`, `combination`, `mean_pct_level`, `mean_pct_reference`, `log2_fold_change`, `n_samples_level`, `n_samples_reference`, `pseudo_pct`.

Colouring whole islands by fold change replaced a first version that mapped the difference in cell density between the levels. That map was set by the abundant lineage islands, so every cytokine island stayed near-white however much it changed: an 11-fold rise in `CD4+ | IFNy+ IL-2+ TNFa+ CD107a+` read as no change.

Per-sample counts for testing are in `combination_counts.csv` (cells and population cells per sample and combination); `combination_frequencies_<col>.csv` and `combination_change_<col>.csv` only summarise them per level.

`--split_by` also works in the markers layout. There the plots use the subgate labels of `--subpop_mode` instead of gate combinations, and the percentages are shares of the plotted, possibly downsampled, cells.

---

## Parameters

| Option | Default | Effect |
|--------|---------|--------|
| `--gate_gap` | `2` | Jump at each cut of the deepest tier, in feature units. Must exceed 1 for the separation guarantee; sides are more than `2 * (gap_m - 1)` apart |
| `--gate_tier_ratio` | `2` | Factor by which the jump grows per tier up the tree; must be 1 or more. Above 1, the islands of one lineage gather in one region |
| `--gate_softness` | `4` | Width of the tanh term in units of `sigma`; larger flattens islands, smaller pushes cells to the island edges |
| `--viability_pattern` | common viability dyes | Case-insensitive regular expression of marker names that never decide; `""` turns it off |
| `--supervise` | `0` | uwot `target_weight` with the combination label as target; 0 = unsupervised |
| `--cells_per_combination` | `2000` | Per sample, cap on cells kept per call pattern |
| `--min_combination_cells` | `50` | Combinations with fewer cells in total are pooled as `other (rare)`; 0 turns pooling off |
| `--umap_subsample` | `0` (= 100,000) | Fit budget: each combination contributes at most this / (number of combinations) cells |
| `--n_neighbors`, `--min_dist` | `15`, `0.25` (`0.1` in the markers layout) | Usual [UMAP parameters](umap.md#parameter-tuning) |

Defaults were chosen on the MTB cohort; see [docs/plans/umap_gate_layout_plan.md](../plans/umap_gate_layout_plan.md) sections 11 and 11.1.

---

## Reading the Plots

The figures are easy to over-read, because the layout is built from the gates. In plain terms:

**The islands.** Each island is one combination of gate calls, for example CD4+ T cells that the gates call IFNy+ and TNFa+. A cell is in an island because the gates put it there: the layout encodes the gater's calls with a guaranteed gap before UMAP sees the data. That islands are separate, pure and far apart is therefore expected by construction and restates the gating. It is not independent evidence that the combinations are biologically distinct, and a gate drawn in the wrong place gives a clean island of wrongly called cells. What an island does show is the spread of its cells: a cell's position inside it follows its marker values relative to the other cells of the same combination in its own sample, so a bright subset or a smear towards a gate is visible. A shift of a whole combination between tubes or timepoints (every IFNy+ cell brighter after vaccination, say) does not move the cells; read intensities from the marker plots. An island's size says nothing about how common the combination is, because the fit gives every combination comparable room. The distance between two islands means only that their calls differ, as in any UMAP ([see the UMAP pitfalls](umap.md#over-interpreting-distances-between-clusters)); the islands of one lineage sit near each other only because the lineage jump is made larger on purpose.

**The split plots.** A split plot shows the cells of one level (one timepoint, say) in colour on the shared map, over a grey silhouette of all cells, so an island that is empty at that level is still visible. The number next to each combination in the legend is its frequency at that level: the percentage of the population in that combination, averaged over the level's samples, from the full counts. The density of coloured points is not the frequency. Each tube contributes at most 2,000 cells per call pattern (`--cells_per_combination`), so rare combinations are over-represented among the points on purpose, and a level with more tubes shows more points.

**The fold-change maps.** A fold-change map colours each island by how much more or less common its combination is at one level than at the reference level: red for more, blue for less, near-white for no change, grey for a combination absent at both. It compares the two levels' average percentages per sample, so it says nothing about how the frequency changed within any one patient, and it has no test or confidence interval. The patients behind two levels can differ: in MTB, pre-vaccination (`V0`) holds 8 stimulated tubes from 4 patients, none of them previously peptide-vaccinated, while the level after the first vaccination holds 60 tubes from 12 patients, so a change against V0 mixes time with differences between patients. One patient with a strong response can drive a mean, and a rare combination can show a large fold change from a few cells. The colours stop at 8-fold. Use the map to see where to look, then test the per-sample counts in `combination_counts.csv` with a model that accounts for patients.

**Further cautions.**

- **Undecided is placed with negative.** A marker with no gate for a cell (cytokines in CD4-CD8- cells, for example) puts the cell on the negative side, and its label leaves that tier out. Bright cells in such a population are not called positive by anything.
- **Supervision** (`--supervise`, off by default) pulls combinations apart by their labels, which makes the picture tidier and the separation even less informative.
- **Low mixing.** An island with low `sample_mixing` or `group_mixing` holds samples or workspaces that sit apart inside it; check the gating of those samples before reading a change in that island as biology.
- **Touching islands.** Two colours drawn touching (the log names islands with `isolated` below 0.9) are still two combinations; the contact is where UMAP placed them in 2-D, not an overlap of their cells.

---

## How to Use in the Pipeline

```bash
Rscript scripts/umap_from_flowjo.R \
  --wsp /mnt/volume/data/MTB_FACS/wsp \
  --fcs_dir /mnt/volume/data/MTB_FACS/fcs \
  --population "Lymphocytes/Single Cells/Live cells" \
  --marker_aliases "INFy=IFNy;IFN=IFNy;TNF=TNFa" \
  --layout gates \
  --sample_metadata mtb_umap_metadata.csv \
  --filter "stimulation=stimulated" \
  --split_by timepoint \
  --output_dir results/umap_mtb_gates
```

The full set of options, the metadata format and the output files are described in [UMAP from a FlowJo Population](../scripts/umap-from-flowjo.md#gate-separated-layout).

### Source Code Location

The pure functions (decisions, calls, anchors, tiers, labels, missing-decision fills, scales, features, stratified sampling, balanced fit index, sample filters, level order, separation metrics and the density grid behind the outline) are in [`src/visualization/umap_gate_layout.R`](../../src/visualization/umap_gate_layout.R). Everything there except `gate_layout_extract_sample()` and `gate_geometry()` is base R (the separation metrics also need `RANN` or `FNN`), so it can be tested without CytoML: `Rscript test_data/test_umap_gate_layout.R`. The frequency and change tables and the figures are in [`src/visualization/umap_split_plots.R`](../../src/visualization/umap_split_plots.R), their palette and theme in [`src/visualization/umap_style.R`](../../src/visualization/umap_style.R). The design and its rationale are in [the plan](../plans/umap_gate_layout_plan.md).

---

## References

- McInnes L, Healy J, Melville J (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv:1802.03426. [doi:10.48550/arXiv.1802.03426](https://doi.org/10.48550/arXiv.1802.03426)
- Melville J. uwot: The Uniform Manifold Approximation and Projection (UMAP) Method for Dimensionality Reduction. R package. [CRAN](https://cran.r-project.org/package=uwot); documentation of [`umap()`](https://jlmelville.github.io/uwot/reference/umap.html) (arguments `y` and `target_weight`) and [`umap_transform()`](https://jlmelville.github.io/uwot/reference/umap_transform.html).
- Parks DR, Roederer M, Moore WA (2006). A new "Logicle" display method avoids deceptive effects of logarithmic scaling for low signals and compensated data. Cytometry A 69(6):541-551. [doi:10.1002/cyto.a.20258](https://doi.org/10.1002/cyto.a.20258)

---

## Related Documentation

- [UMAP from a FlowJo Population](../scripts/umap-from-flowjo.md): The script that runs this layout, its options and outputs
- [UMAP](umap.md): The embedding itself, parameters and common pitfalls
- [Downsampling](downsampling.md): Uniform and stratified sampling strategies
- [Logicle Transform](logicle-transform.md): The transform family behind FlowJo's biexponential scale
