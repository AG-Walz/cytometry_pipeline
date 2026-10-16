[< Back to Documentation Index](../README.md)

# Harmonise Marker Names Across a Cohort

`scripts/harmonize_marker_names.R` writes one canonical marker name into the
`$PnS` keyword of every FCS file in a folder, so a cohort stops spelling the
same antibody as `INFy`, `IFN` and `IFNg`, and files that never carried a
marker name get one from the channel they were measured on.

The script is a thin CLI wrapper. The work happens in
**`harmonize_marker_names()`**
([`src/metadata_qc/harmonize_marker_names.R`](../../src/metadata_qc/harmonize_marker_names.R)),
on top of the FCS keyword layer in
[`src/R_shared/fcs_keywords.R`](../../src/R_shared/fcs_keywords.R), the alias
rules in
[`src/harmonization/panel_homogenization.R`](../../src/harmonization/panel_homogenization.R)
and `assign_canonical()` from
[`src/harmonization/drift/drift_panel.R`](../../src/harmonization/drift/drift_panel.R).

## Prerequisites

**R 4.3+ and nothing else.** No conda, no Bioconductor, no flowCore. The CLI
uses `optparse` when it is installed and a minimal `--key value` parser when it
is not; `--canonical` needs the `yaml` package. The module and its tests need
neither.

---

## Read this before you use it

### 1. The problem

In the MTB cohort 162 of 328 files carry no `$PnS` at all, and the rest spell
the same marker as `INFy` / `IFN`, `TNFa` / `TNF`, `ALD` on two different
channels. Every downstream step patched this with its own per-cohort regexes.
This tool fixes the files once, reversibly, so later steps can address channels
by one marker name.

Only `$PnS` is written. `$PnN` (the detector name) and every structural
parameter keyword (`$PnB`, `$PnE`, `$PnR`) are never touched, so FlowJo gates
and compensation, which address channels by `$PnN`, are unaffected.

### 2. Where the new names come from

Three rules, applied in this order of precedence per parameter:

| Source | What it is | When to use it |
|---|---|---|
| `alias` | an explicit map, as written to to use: `--aliases "INFy=IFNg;TNF=TNFa"` or a two-column file | to confirm a low-confidence pair, or to override the other two rules |
| `canonical` | the `markers` block of a drift config: canonical name to regex patterns, matched case-insensitively against `$PnS` first and then `$PnN` | the normal case. Fills blank `$PnS` from the channel name, so it is the only rule that helps files with no marker names |
| `spelling` | `suggest_marker_aliases()` groups that are identical after case, punctuation and a trailing Latin Greek letter are dropped (`TNFa` / `TNF`), the most-used spelling winning | automatic, but only when neither rule above decided the parameter |

Groups that are identical only after the letters are **sorted** (`INFy` /
`IFN`) are *transpositions*: right nearly always and wrong spectacularly when
wrong, so they are never applied. The dry run lists them and you confirm with
`--aliases`. `CD4` and `CD40` never merge under either rule.

Scatter (`FSC*`, `SSC*`, including `SP SSC-A`) and `Time` are never renamed,
even when the canonical map has patterns for them.

### 3. Run the census before you write anything

`--probe` prints, for every channel, every `$PnS` spelling seen and how many
files use it. That table is the whole cohort's naming in one screen, and it is
how you check that the canonical map covers what is actually there.

---

## Usage

```bash
# 1. Census. Read-only.
Rscript scripts/harmonize_marker_names.R --fcs_dir data/fcs --probe

# 2. Dry run against the drift config. Prints the plan grouped by channel and
#    the low-confidence pairs. Still writes nothing.
Rscript scripts/harmonize_marker_names.R --fcs_dir data/fcs \
  --canonical config/drift_mtb.yaml

# 3. Apply in place, with an undo ledger next to the data.
Rscript scripts/harmonize_marker_names.R --fcs_dir data/fcs \
  --canonical config/drift_mtb.yaml --apply --in_place

# Or write patched copies and leave the originals alone.
Rscript scripts/harmonize_marker_names.R --fcs_dir data/fcs \
  --canonical config/drift_mtb.yaml --apply --out_dir data/fcs_harmonised

# Undo the last run. Repeat to walk back through earlier runs.
Rscript scripts/harmonize_marker_names.R --fcs_dir data/fcs --revert
```

Nothing is written without `--apply`, and `--apply` refuses to guess between
editing the originals and writing copies: say `--in_place` or `--out_dir`.
With `--out_dir`, files that need no change are copied verbatim so the folder
is a complete cohort.

A dry run **exits 0 even when it finds problems**; only `--apply` turns a file
that could not be written into a non-zero exit.

### Options

| Flag | Default | Meaning |
|---|---|---|
| `--fcs_dir` | `data` | Folder, searched recursively and case-insensitively |
| `--probe` | off | Census only: channel x spelling x files. Writes nothing |
| `--canonical` | none | YAML file whose top level is a `markers` block, or a drift config containing one |
| `--aliases` | none | `"INFy=IFNg;TNF=TNFa"` or a two-column file (as written, to use) |
| `--apply` | off | Actually write |
| `--in_place` | off | Patch the originals where they sit |
| `--out_dir` | none | Write patched copies here instead, mirroring subfolders. Must not be inside `--fcs_dir` |
| `--revert` | off | Undo one previous run from its ledger |
| `--ledger` | `marker_names_ledger.csv` | Undo ledger, written beside the written data |
| `--no_recursive` | off | Stay at the top level |
| `--src_dir` | auto | Path to the repo `src` directory |

### The canonical map

The same shape as the `markers` block of
[`config/drift_mtb.yaml`](../../config/drift_mtb.yaml), so the drift config
itself can be passed. The first pattern that hits wins; `$PnS` and `$PnN` are
both tried, area channels are preferred over height and width:

```yaml
markers:
  LiveDead: ["live", "dead", "^ALD$", "aqua", "^BV480-A$"]
  CD4:      ["^CD4$", "^APC-Cy7-A$"]
  IFNg:     ["IFN", "INF", "^PE-A$"]
  TNFa:     ["TNF", "^Pacific Blue-A$"]
```

After the files are harmonised, a downstream config only needs the canonical
names themselves: `IFNg: ["^IFNg$"]`.

---

## How it writes

Through `fcs_write_keywords()`, so the event data is **never decoded**: the new
TEXT segment is written in place when it fits and the file is re-emitted with
the data segment copied verbatim when it does not. Every write is verified
against the original (keyword sequence, segment length, every byte of the data
segment) and rolled back on any mismatch. See
[sync-fcs-names.md](sync-fcs-names.md) for the three strategies and the five
witnesses; they are the same code.

Symlinks are refused. A file that fails the keyword layer's safety checks
(`nextdata_nonzero`, `crc_present`, and so on) is reported and left
byte-identical.

### The ledger

Every replaced value is appended to a CSV next to the written data, one row per
keyword, with a run id (timestamp and process id) before the file is opened for
writing, and a status row after. `--revert` replays the most recent run that
has not been reverted yet: it restores every `$PnS` that run replaced and
removes the ones it created. Nothing is stamped into the files themselves; the
ledger is the provenance record.

---

## Calling it from R

```r
source("src/R_shared/fcs_keywords.R")
source("src/harmonization/panel_homogenization.R")
source("src/harmonization/drift/drift_panel.R")
source("src/metadata_qc/harmonize_marker_names.R")

files  <- list.files("data/fcs", "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
census <- marker_census(files, "data/fcs")
census_summary(census)                          # channel x spelling x n_files

pl <- plan_marker_names(census, canonical = yaml::read_yaml("config/drift_mtb.yaml")$markers,
                        aliases = c(INFy = "IFNg"))
pl$plan          # filename, index, channel, marker_old, marker_new, source
pl$suggestions   # transposition groups to confirm by hand

apply_marker_plan(pl$plan, "data/fcs")                              # dry run
res <- apply_marker_plan(pl$plan, "data/fcs", in_place = TRUE, dry_run = FALSE)
revert_marker_plan(attr(res, "ledger_path"))
```

## Tests

```bash
Rscript test_data/test_harmonize_marker_names.R
```

Builds its fixtures byte by byte with `test_data/fcs_fixtures.R`; needs neither
flowCore nor optparse nor yaml.

## Related

- [Sync FCS Embedded Names](sync-fcs-names.md), the sibling tool for `$FIL`
- [scripts/drift/README.md](../../scripts/drift/README.md), which expects this to have run before step 1
- [Metadata QC Reference](../subprojects/metadata-qc.md)
