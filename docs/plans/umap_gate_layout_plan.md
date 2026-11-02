[< Back to Documentation Index](../README.md)

# Gate-separated UMAP, split by sample metadata

Status: plan written 2026-10-13. Implementation, tests and the real-data tuning run follow in the same session; section 11 records the outcome.

## 1. What the user asked for

- One UMAP over every patient, built once, then drawn several times: pre-vaccination, then each later timepoint, each on the same coordinates with only that timepoint's cells shown.
- Every gater-defined population, including every cytokine combination (CD4/CD8 x IFNg, IL-2, TNFa, CD107a), should land in its own, ideally visually separate, island.
- Get there by dropping channels nobody is interested in (scatter, live/dead) and by transforming the FlowJo values so that the differences the gates decide become as large as possible.
- Reproducible, part of the Dash app, and generic: any panel, any gating tree, any metadata. The user should get the plots by supplying a metadata table for their files.
- Plots should look like `/mnt/volume/data/pasted-image.png`: muted pastel palette, soft translucent point clouds, compact separate islands, no axes, a small UMAP1/UMAP2 arrow in a corner, and no leader lines pointing at clusters.

## 2. What the data and code look like today (facts this plan relies on)

- MTB cohort: 328 FCS, 20 workspaces, about 9.8 million Live cells. Identical copies on the dev box (`/mnt/volume/data/MTB_FACS`) and on highmem1 (`/mnt/volume/patrick/UMAP`).
- Values from `gh_pop_get_data()` are compensated and FlowJo-biexp transformed onto a 0 to 256 scale. Gates sit on the same scale.
- Under `Live cells` every gate is a `rectangleGate`: CD4 is APC-Cy7 above a cut and PE-Cy7 below one, CD8 the opposite, and the two blank-named gates are the CD4+CD8+ and CD4-CD8- quadrants. Each of CD4 and CD8 has four single-cytokine gates drawn as cytokine x FSC-A rectangles. There are no Boolean combination gates.
- Gate positions change per visit (the analysts followed the detector voltages), and the IFNg gate is spelled `IFNy` in some workspaces and `INFy` in others.
- Cytokine-positive cells are rare (0.1 to 6.5 % of a lineage in stimulated tubes, often zero in neg tubes). Uniform downsampling to 1000 cells per tube, as in the 2026-09-17 run, keeps a handful per tube.
- The current module labels a polyfunctional cell by one single-cytokine gate, splits IFNy/INFy into two labels, keeps the live/dead channel unless `--markers` removes it, never reads gate geometry, and embeds z-scored FlowJo values. In the 2026-09-17 run the cytokine-positive cells are sprinkled through the CD4 and CD8 islands instead of forming their own.
- Metadata (`--sample_metadata`) already resolves rows to samples by file name and draws one single-colour plot per level, but cannot order levels, cannot filter samples, and does not colour a level's plot by population.

## 3. Approach in one paragraph

Read the gates themselves. For every gate under the chosen population, find which marker it decides and where its cut lies in that sample. For every cell, record the gater's call per marker (positive, negative, not decided). Re-express each decision marker relative to the gate cut of the cell's own sample, compress the spread within each side, and add a fixed jump at the cut so positives and negatives sit far apart. Cells that differ in any gate call then differ by at least the jump in at least one coordinate, so every combination of calls is its own region of feature space, while the continuous part keeps a natural cloud shape inside each island. Sample the cells so that rare combinations are kept, fit the UMAP on a set balanced across combinations, project the rest, and draw it in the reference style, once overall and once per metadata level.

## 4. Gate-anchored features (`--layout gates`)

Per sample, for the cells of the selected population:

1. **Decision gates.** Every non-Boolean gate below the population. A gate dimension is a decision on marker *m* when the channel maps to a marker (via the existing channel dictionary and aliases; scatter and Time never do) and the gate excludes a real share of its parent along that axis:
   - rectangle: a finite lower bound that alone excludes at least 0.5 % of the parent's cells is a positive decision with cut = lower bound; a finite upper bound doing the same is a negative decision with cut = upper bound; both binding (a band) is not used as a call and gets no row in `gate_decisions.csv` (superseded, see 11.1: it only counts towards `n_axes`, and a far edge excluding under 5 % no longer makes a band). FSC-A bounds on the cytokine gates exclude nothing and FSC is not a marker, so they drop out by both rules.
   - polygon, ellipse and anything else: decide from membership. With at least 10 members and 10 non-members in the parent, the marker is a decision when the AUC of member against non-member is at least 0.8 or at most 0.2; polarity is the side of the AUC and the cut is the Youden-optimal 1-D threshold. Otherwise that gate makes no decision in that sample.
   - Boolean gates are ignored; they are derived from the others.
2. **Calls.** For cell *i* and marker *m*, the applicable decisions are those whose gate's parent contains *i*. Only the deepest of them count: when a marker is gated both under `Viable` and again under `CD4+` (as in `test_data/flowjo_data/Killing_AML.wsp`), the CD4+-level gate decides for CD4+ cells and the shallower gate only for cells no deeper gate covers. (Superseded in part, see 11.1: a gate also never decides a marker that an ancestor gate already decides.)
   - If *i* is inside one or more applicable gates: positive if any of them decides *m* positively, otherwise negative.
   - If *i* is inside none: negative when every applicable decision is positive (a cell outside the IFNg+ gate), positive when every one is negative, and `x >= cut` when they are mixed. (Superseded, see 11.1: the opposite call is inferred only when every deciding gate binds on that marker alone; otherwise `x >= cut` decides.)
   - No applicable decision: not decided (NA). Example: CD4-CD8- cells have no cytokine gate.
   The call is the gater's decision, membership included, so a cell outside a rectangle because of its FSC-A is negative even if its cytokine value is above the cut.
3. **Cut per cell** = the median cut of its applicable decisions on *m* (the four quadrant gates share one cut; a cytokine has one gate per lineage). Undecided cells use the sample's median cut for *m*. (Superseded, see 11.1: the cut now only decides `x >= cut` calls and is reported; features are centred on anchors.)
4. **Scale.** After all samples are read, sigma_m = weighted MAD of (x - cut) (superseded, see 11.1: x - anchor) over negative cells, pooled across samples, floored at 1 % of the channel range.
5. **Feature.** `f = tanh((x - cut) / (softness * sigma_m)) + gap_m * (call ? +1 : -1)`, NA calls treated as negative, with `gap_m = gap * tier_ratio^(deepest tier - tier of m)`, so the first split of the tree (CD4 against CD8) is the largest distance and each lineage's islands gather in one region. The tanh term is bounded by 1, so with `gap > 1` positive and negative cells of a marker can never overlap, whatever the sample's voltages were. Centring on the sample's own cut also removes the per-visit gate shifts, which is what makes the same combination from different patients and days land in the same island. (Superseded, see 11.1: the tanh term is centred on an anchor, not on the cut; in the final code the anchor is the sample's own median over the cells of the same combination of calls, `combination_anchors()`.)
6. **Markers used.** Default: the decision markers found under the population, which excludes scatter, Time, live/dead and anything gated only above the population. Marker names that look like a viability dye (`live`, `dead`, `viab`, `zombie`, `ALD`, `7-AAD`, `DAPI`, `PI`, ...) are never decisions, so choosing a population above the live gate does not give dead cells their own island. `--markers` restricts or extends this; a requested marker that no gate decides enters as `tanh(robust z / softness)` without a jump.

Labels: a marker's tier is the depth below the population of the gates that decide it (CD4 and CD8: tier 1; cytokines: tier 2). A combination label lists, per tier, the positive markers in panel order, for example `CD4+ | IFNy+ TNFa+`; a tier with none positive reads `none`, except tier 1, which lists its markers as negative (`CD4- CD8-`). Undecided tiers are omitted. The names come from markers, not gate names, so `IFNy` and `INFy` gates give one label. Combinations with fewer than `--min_combination_cells` cells in total are pooled as `other (rare)` for colouring and supervision. (In the final code the pooled label is used everywhere: fit balancing, separation metrics and every count, frequency and change table; the calls stay in `cells.csv.gz`.)

Consistency checks written to `gate_decisions.csv` and the log: per sample, which markers were decided by which gate, cut, polarity and method, plus warnings when a sample lacks a decision that most samples of the same tier-1 group have (that sample's cells are then treated as negative for that marker). (Superseded, see 11.1: the majority is over the group's cells, counted with the sampling weights, not over samples.)

## 5. Sampling and embedding

- **Counts first.** Combination counts per sample are taken from all cells before any sampling and written to `combination_counts.csv`; every frequency shown later comes from these, not from the plotted cells.
- **Stratified sampling.** Per sample, keep up to `--cells_per_combination` cells of each combination (default 2000; rare combinations are kept whole). Each kept cell carries a weight = cells of its combination in that sample / cells kept.
- **Balanced fit.** Fit the UMAP on at most `--umap_subsample` cells (default 100000 in this layout) drawn evenly across combinations, then project the remaining sampled cells with `uwot::umap_transform`. (Superseded, see 11.1: every combination is capped at budget / number of combinations, and what small ones leave unused is not redistributed.) Balancing is what gives every combination comparable room, which is the even island layout of the reference figure; projection puts the abundant cells into the islands their features already point to.
- **Optional supervision.** `--supervise` sets uwot's `target_weight` with the combination label as target (0 = off). The gap alone should already separate combinations; supervision is tuned in section 8 and documented as making the layout gate-driven by construction.
- `run_umap()` gains `y`, `target_weight` and `fit_idx` arguments (defaults reproduce today's behaviour). Threads stay `embed_cores` for the fit; the run records the thread count because the batch optimiser differs from the serial one. (Added later: `transform_threads`, so the projection's neighbour search uses every plotting core; it gives the same coordinates at any thread count.)

## 6. Metadata-driven plots

- `--sample_metadata` stays the single metadata input (CSV/TSV/XLSX upload in the app). Any column can drive the plots.
- `--filter "col=a,b;col2!=c"` keeps only matching samples before anything is extracted (for MTB: `stimulation=stimulated`). Unassigned samples fail any `=` clause.
- `--split_by col[,col2]`: per column, one plot per level showing that level's cells coloured by combination over a light grey silhouette of all cells; one panel figure with every level side by side; one map per level of each combination's change in frequency against a reference level (superseded the density-difference map first planned here, see 11.1).
- Level order: natural sort by default (numbers compare as numbers, so `W5 < W12 < W110`); `--level_order "col=lvA,lvB,...;col2=..."` overrides; `--split_reference` picks the reference level (default: the first).
- Change maps (superseded the density difference first planned here, see 11.1): each island is coloured by the log2 fold change of its combination's mean % per sample at the level against the reference, with one cell in the median sample of the two levels added to both means, capped at 8-fold, grey where the combination is absent at both levels; symmetric muted red/blue scale like panel E of the reference, grey outline of the embedding. The numbers go to `combination_change_<col>.csv`.
- Per-level frequency tables `combination_frequencies_<col>.csv`: mean over samples of each combination's share of the population, with its SD and the number of samples per level.
- App: a "Download metadata template" button writes a CSV with one row per loaded sample (`filename` plus columns parsed by the optional regex), so a user only has to fill in columns.

## 7. Style

A shared `umap_style.R`:

- muted palette modelled on the reference: dusty rose, sage, dark green, teal, lavender, ochre, brick, slate blue, navy, sand; extended in HCL at the same chroma when more colours are needed (superseded: each extra colour is the muted HCL candidate furthest in CIE-Lab from the colours already chosen);
- combinations coloured by tier-1 group so lineage reads at a glance (one hue family per group, the all-negative combination of a group in its lightest tint, CD4-CD8- in neutral grey);
- `theme_void`, small sans title, legend with small round keys and per-label counts or percentages, a short L-shaped UMAP1/UMAP2 arrow in the lower-left corner, no leader lines or centroid dots;
- small, low-alpha, rasterised points (ragg PNG; `ggrastr` for SVG/PDF when installed), point size and alpha scaled to the number of cells unless set, drawn from the largest combination to the smallest so rare islands stay visible.

## 8. Tuning on the real cohort (highmem1)

Extract once with the new layout, save the per-cell features, then sweep gap, softness, supervision weight, `n_neighbors` and `min_dist` on the balanced fit set. Score each embedding with:

- macro-averaged 15-NN label purity in 2-D (every combination counts equally);
- share of cells whose 15 nearest 2-D neighbours all share its label;
- simplified silhouette (distance to own centroid against the nearest other centroid);
- for alignment: within each combination, how evenly patients mix among 15-NN, relative to their shares.

The defaults become the best-scoring setting after a visual check of the top candidates against the reference figure. Every normal run writes the same metrics to `separation_metrics.csv` (with mixing over workspaces and over samples, `group_mixing` and `sample_mixing`, not over patients).

## 9. Reproducibility and outputs

- `run_info.json`: argv, git commit and dirty state, R and package versions, seed, thread count, feature scales, decisions summary, and a copy of the metadata table (`sample_metadata.csv`).
- `cells.csv.gz`: coordinates, sample, workspace, combination, calls, weight, fit flag and annotation columns for every plotted cell.
- `gate_decisions.csv`, `combination_counts.csv`, `combination_frequencies_<col>.csv`, `separation_metrics.csv`.
- Plots: `umap_<pop>_by_combination`, `umap_<pop>_split_<col>_<level>`, `umap_<pop>_split_<col>_panel`, `umap_<pop>_diff_<col>_<level>`, plus the existing marker, file and workspace plots.

## 10. Integration and delivery

- R: new `src/visualization/umap_gate_layout.R` (pure functions, testable without CytoML), `src/visualization/umap_style.R`, `src/visualization/umap_split_plots.R`; small hooks in `umap_from_flowjo.R`, `scripts/umap_from_flowjo.R`, `stats_utils.R` (`dimred_dispatch.R` was left unchanged: the gates layout calls `run_umap()` directly). The default `--layout markers` keeps every existing output byte-compatible in behaviour.
- App: layout switch and its parameters in step 3; split columns, level order and filter in the metadata section of step 2; template download; new gallery section for split and difference plots; new progress phases.
- Tests: `test_data/test_umap_gate_layout.R` (base-R harness like `test_umap_multi_wsp.R`) for decisions, calls, features, labels, sampling, filters, level order, density grid and palette; the app driven in Chromium against the stand-in; a real run on highmem1.
- `environment.yml`: declare `r-rann` (kNN for the metrics) and `r-ggrastr`, both already present on highmem1.
- highmem1: development runs from a staging copy under `/mnt/volume/tmp_claude_umap/`; the launcher checkout is left alone, since its `git pull origin main` would otherwise merge into local edits. The code reaches the launcher when the user commits and pushes. `cytometry_env` there lacks pandas and openpyxl, which the existing metadata upload already needs.

## 11. Outcome of the tuning (2026-10-13, highmem1)

Data: all 272 stimulated MTB tubes (`--filter stimulation=stimulated`), population `Lymphocytes/Single Cells/Live cells`, marker aliases `INF=IFNy;TNF=TNFa;IFN=IFNy;INFy=IFNy`. Extraction of 328 tubes with the gate calls took 1.4 min after the workspaces were loaded. Stratified at 2000 cells per combination per tube: 1.53 million cells, 32 combinations plus 90 cells pooled as `other (rare)` (every lineage x cytokine combination occurs; the smallest kept ones have 81 to 164 cells). Every tube has the same six decision markers (CD4, CD8 at tier 1; IFNy, IL-2, TNFa, CD107a at tier 2), so no decision had to be filled in. Feature scales (negative-population MAD, biex units): CD4 11.6, CD8 24.2, IFNy 6.9, IL-2 17.6, TNFa 8.2, CD107a 14.2.

Each configuration was fitted on a 60000-cell set balanced across combinations, 50000 further cells were projected, and both were scored together (scripts and per-configuration outputs in `/mnt/volume/tmp_claude_umap/tune/` on highmem1).

| Setting | 15-NN purity | cells with all 15 neighbours in their own combination | silhouette | patient mixing in islands | worst combination purity |
|---|---|---|---|---|---|
| z-scored FlowJo values (as in the 2026-09-17 run), n_neighbors 30 or 50 | 0.54 to 0.57 | 18 to 19 % | -0.04 | 0.68 | 0.08 |
| gate layout, any of 18 settings of gap 1.5 to 3, softness 1 to 4, supervision 0 or 0.5 | 1.000 | 100 % | 0.64 to 0.78 | 0.87 to 0.91 | 1.000 |

- Separation is not sensitive to gap, softness, n_neighbors (15 to 100) or min_dist (0.1 to 0.4). Supervision (`target_weight` 0.5) changed nothing measurable, so the default is unsupervised and the layout comes from the features alone.
- Lineage coherence (mean distance between island centres of different first-tier groups over the same within a group): 1.32 to 1.38 with one gap for every tier, 1.46 to 1.49 with the lineage jump twice the cytokine jump (`tier_ratio` 2), 1.42 to 1.45 at 3. With ratio 2 the CD8 islands sit on one side, the CD4 islands on the other, CD4-CD8- and CD4+CD8+ between them, and within each lineage the islands run from cytokine-negative to polyfunctional.
- min_dist 0.25 gives rounder, softer islands than 0.1 (silhouette 0.62 against 0.74, both fully separated), which is closer to the reference figure.

Defaults first chosen (cut-centred features): gap 1.5, softness 2, tier ratio 2, supervision 0, n_neighbors 15, min_dist 0.25.

### 11.1 After the adversarial review (same day)

The review ran the layout on the repo's second dataset (`test_data/flowjo_data/Killing_AML.wsp`, Boolean gates, cytokine gates at two depths, 2-D lineage boxes) and on edited variants of it. Confirmed defects and what changed:

- A gate re-decided a marker its ancestor already decides (a CD4 x IFNg box under CD4+ clipping the CD4-dim tail made almost every CD4+ cell CD4-negative). Now a gate never decides a marker an ancestor gate decides.
- Outside a gate that binds on two markers (a lone CD4+CD8- box), the code inferred the opposite call on both axes, so double negatives became CD8+. Now polarity is only inferred for gates that bind on that marker alone (`n_axes` in `gate_decisions.csv`); otherwise the cell's own value against the cut decides. A far edge clipping under 5 % no longer turns a decision into a band.
- Centring the continuous term on each tube's cut split one population into sample-pure halves when the gater set different cuts on identical cells. The continuous term is now centred on the median of the tube's own cells on the same side of the same decision (the anchor), which still removes voltage shifts but no longer depends on where the line was drawn. Sample mixing within islands is now reported (`sample_mixing`).
- Cells no gate decides (a population above the live gate) get their own coordinate and island instead of landing on the all-negative combination; viability-dye names were extended and are configurable (`--viability_pattern`).
- Smaller fixes: a sample with no deciding gate lost all its cells in the sampling; the missing-decision fill counted kept rows instead of cells; `--gate_tier_ratio` below 1 is refused; split and filter columns are always resolved even with `--annotate`; the fit set now caps every combination at budget / combinations even when the run is smaller than the budget (a 57000-cell run was not balanced at all before).
- Difference maps showed absolute density differences dominated by the large lineage islands, so cytokine islands were always near-white (an 11-fold rise in CD4+ | IFNy+ IL-2+ TNFa+ CD107a+ read as no change). They now colour each island by the log2 fold change of its mean per-sample %, with the numbers in `combination_change_<col>.csv`.

Re-tuned with the anchors (same 272 tubes, 60000 fitted plus 50000 projected cells scored):

| Setting | 15-NN purity | isolated | silhouette | workspace mixing | sample mixing | worst combination | lineage coherence |
|---|---|---|---|---|---|---|---|
| cut-centred, gap 1.5, softness 2 | 1.000 | 1.000 | 0.605 | 0.910 | 0.985 | 1.000 | 1.46 |
| anchored, gap 1.5, softness 2 | 0.997 | 0.993 | 0.612 | 0.924 | 0.986 | 0.935 | 1.41 |
| anchored, gap 1.5, softness 1 | 0.915 | 0.798 | 0.225 | 0.912 | 0.968 | 0.421 | 1.54 |
| anchored, gap 2, softness 2 | 1.000 | 1.000 | 0.656 | 0.928 | 0.988 | 1.000 | 1.37 |
| anchored, gap 2, softness 4 | 1.000 | 1.000 | 0.671 | 0.919 | 0.986 | 1.000 | 1.41 |
| anchored, gap 2, softness 6 | 1.000 | 1.000 | 0.681 | 0.916 | 0.984 | 1.000 | 1.42 |

Anchoring centres each side, so the guaranteed distance between sides is the bound 2 (gap - 1) itself, where cut-centring had left most negatives below and most positives above the cut; gap 2 restores full separation. Softness 2 leaves striped polyfunctional islands (several within-side axes saturate at both ends); softness 4 gives smooth islands at no cost.

A second re-check on real data then showed that a side median still depends on which populations share a side: CD8 T cells and double negatives are both CD4-negative, in proportions that differ per tube, so on two MTB workspaces the 'CD8+ | CD107a+' island split by workspace (mixing 0.22). The anchor is therefore the median of the tube's cells **of the same combination** (`combination_anchors()`). Re-tuned once more on the 272 tubes:

| Setting (per-combination anchors) | 15-NN purity | isolated | silhouette | workspace mixing | sample mixing | worst combination | lineage coherence |
|---|---|---|---|---|---|---|---|
| gap 1.5, softness 4 | 1.000 | 1.000 | 0.674 | 0.966 | 0.997 | 1.000 | 1.39 |
| gap 2, softness 2 | 1.000 | 1.000 | 0.668 | 0.973 | 0.998 | 1.000 | 1.41 |
| gap 2, softness 4 | 1.000 | 1.000 | 0.661 | 0.969 | 0.997 | 1.000 | 1.43 |
| gap 2.5, softness 4 | 1.000 | 1.000 | 0.693 | 0.968 | 0.998 | 1.000 | 1.45 |

Workspace mixing within islands rose from about 0.92 to 0.97, and 'CD8+ | CD107a+' mixes at 0.97. On the Killing_AML fixture and on two MTB workspaces, gap 2.5 changed little against gap 2. One limitation remains on the small fixture: two small islands (101 and 202 cells) that differ on two markers are drawn crossing each other in 2-D (isolated 0.79 and 0.76), a placement artefact of small disconnected graph pieces, not overlap in feature space. The run log now names any island with isolated below 0.9.

Final defaults for `--layout gates`: gap 2, softness 4, tier ratio 2, supervision 0, n_neighbors 15, min_dist 0.25, fit set capped at 100000 / number of combinations per combination, 2000 cells per combination per tube, combinations under 50 cells pooled.

### 11.2 Final MTB run (2026-10-13, highmem1)

Command (from the staged working tree, `cytometry_env`): `--layout gates --population "/Lymphocytes/Single Cells/Live cells" --marker_aliases "INF=IFNy;TNF=TNFa;IFN=IFNy;INFy=IFNy" --sample_metadata /mnt/volume/patrick/UMAP/mtb_umap_metadata.csv --filter stimulation=stimulated --split_by timepoint --annotate timepoint --output_formats png`, all other settings at their defaults. Output: `/mnt/volume/patrick/UMAP/umap_gates_timepoint_20261013/` (1666 files, 287 MB, per-file marker plots included).

- 272 of 328 tubes kept by the filter; 1,526,691 sampled cells, 33 combinations; fit set 68,151 cells balanced over them, the rest projected.
- Separation (macro over combinations): 15-NN purity 1.000, isolated 1.000, silhouette 0.784, workspace mixing 0.966, sample mixing 0.997. No island flagged for low mixing or for touching another.
- 12 min 40 s wall clock, 3.3 GB peak memory. Two runs from the same code gave identical metrics.
- Example from `combination_change_timepoint.csv`, V2 (111 tubes) against V0 (8 tubes, 4 patients without prior peptide vaccination): CD4+ IFNy+ IL-2+ TNFa+ CD107a+ 0.04 % to 0.45 % of Live cells (log2 fold change 3.25), CD4+ IFNy+ TNFa+ 0.037 % to 0.14 % (1.7), CD8+ IFNy+ TNFa+ CD107a+ unchanged (0.083 % to 0.076 %). The patients differ between levels, so these are descriptive, not paired, comparisons.
