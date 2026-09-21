[< Back to Documentation Index](../README.md)

# Cytokine combinations per patient over time

`scripts/cytokine_combinations.R` turns FlowJo workspaces plus their FCS files
into one figure per patient: a stacked bar per FCS file showing the Boolean
combinations of the cytokine gates (CD107a, IFNy, IL-2, TNF) among CD4 and CD8
cells, grouped by timepoint and faceted by stimulation cocktail. It answers
"how did this patient's antigen-specific T-cell response drift across the
vaccinations?".

The script is a thin CLI wrapper around `cytokine_combinations()` in
[`src/analysis/cytokine_combinations.R`](../../src/analysis/cytokine_combinations.R).

## What it does

1. Discovers the workspaces under `--wsp`, pairs each with its FCS files by
   file name under `--fcs_dir` and builds one `GatingSet` per workspace, using
   the same loader as the [FlowJo UMAP](umap-from-flowjo.md).
2. For every sample and every parent gate (`--parents`, default `CD4,CD8`) it
   reads which cells fall in each cytokine gate beneath the parent
   (`--markers`, default `CD107a,IFNy,IL-2,TNF`) and counts the 2^k Boolean
   combinations: `none`, `IFNy+`, `IFNy+ IL-2+`, ... `CD107a+ IFNy+ IL-2+ TNF+`.
3. Parses patient, batch, vaccination number, week and cocktail from the sample
   name (`MTB<patient>_<batch>_V<vaccine>W<week>_<cocktail>`).
4. Writes the long table and, per patient, three figures.

## Usage

```bash
conda activate cytometry_env
# One patient
Rscript scripts/cytokine_combinations.R \
  --wsp /mnt/volume/data/MTB_FACS/wsp/MTB059_1.wsp \
  --fcs_dir /mnt/volume/data/MTB_FACS/fcs \
  --output_dir /mnt/volume/data/MTB_FACS

# Whole cohort, plots for two patients
Rscript scripts/cytokine_combinations.R --wsp /mnt/volume/data/MTB_FACS/wsp \
  --fcs_dir /mnt/volume/data/MTB_FACS/fcs --output_dir results/combos \
  --patients MTB059,MTB074

# Re-plot from an earlier table (no CytoML needed, only ggplot2)
Rscript scripts/cytokine_combinations.R \
  --from_table results/combos/cytokine_combinations.csv --output_dir results/combos
```

| Option | Default | Meaning |
|---|---|---|
| `--wsp` | | A `.wsp` file or a folder searched recursively |
| `--fcs_dir` | | Folder searched recursively for the FCS files |
| `--output_dir` | `results/cytokine_combinations` | Where PNGs and the CSV go |
| `--patients` | all | Comma-separated patients to plot (`MTB059,MTB074`) |
| `--parents` | `CD4,CD8` | Parent gates, one facet row each |
| `--markers` | `CD107a,IFNy,IL-2,TNF` | Cytokine gates under each parent |
| `--group` | `All Samples` | FlowJo sample group to load |
| `--output_formats` | `png` | `png`, `pdf`, `svg` |
| `--from_table` | | Re-plot from a `cytokine_combinations.csv` |

## Output

| File | Content |
|---|---|
| `cytokine_combinations.csv` | One row per workspace x sample x parent x combination: `workspace; sample; patient; batch; vaccine; week; timepoint; cocktail; parent; parent_count; combination; count; percent`. Semicolon separated. |
| `<patient>_cytokine_combinations.png` | Rows: parent gate. Columns: cocktail, `neg` first. x: timepoint (`V1 W5`) in chronological order on a discrete axis. Stack: the 16 combinations as % of the parent gate; `none` is grey, the other colours are grouped by the number of positive cytokines (blue 1, green 2, orange 3, red 4). |
| `<patient>_cytokine_combinations_minus_neg.png` | The same with the `neg` tube of the same timepoint subtracted from each combination (floored at 0, `none` recomputed so the bar closes at 100). Cocktails without a `neg` at their timepoint are left out with a warning. |
| `<patient>_cytokine_combinations_minus_neg_positive_only.png` | The neg-subtracted figure without the `none` segment and with a free y axis per row, so the few percent of responding cells are readable. |

## Name conventions and caveats

- Patients are zero-padded to three digits: `MTB59_1_...` (FCS) and
  `MTB059_1.wsp` (workspace) are the same patient.
- Anything ending in `neg` is the unstimulated control (`neg`, `1_neg`). A
  trailing `_NNN` on a cocktail is a re-acquisition index and is stripped
  (`194335_018` is cocktail `194335`). A leading `1` on a repeated pool
  (`1215008` for pool `215008` measured again later) is kept as is; treat it as
  the same cocktail when reading the plot.
- The gate spelled `INFy` in some workspaces is read as `IFNy`.
- A parent gate must occur exactly once in the tree; a cytokine gate missing
  under a parent counts every cell as negative for it, with a warning.
- The FCS data is taken as gated in FlowJo. No compensation, transform or gate
  is recomputed here.

## Sanity checks on a real run

For `MTB59_1_V0W0_224065` the workspace `MTB059_1.wsp` records CD4 = 14209
cells and, under CD4, CD107a = 124, IFNy = 3, IL-2 = 47, TNF = 47. In the CSV
the CD4 `parent_count` of that sample must be 14209, the counts of one sample
and parent must sum to `parent_count`, and the summed counts of all
combinations containing `IFNy+` must be 3.

## Test

```bash
Rscript test_data/test_cytokine_combinations.R
```

Base R only: name parsing, combination counting, background subtraction and
gate-path resolution. The `GatingSet` side and the rendering need an analysis
host with CytoML and ggplot2.
