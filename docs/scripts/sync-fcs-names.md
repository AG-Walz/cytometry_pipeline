[< Back to Documentation Index](../README.md)

# Sync FCS Embedded Names to the File Names on Disk

`scripts/sync_fcs_names.R` writes each FCS file's own file name into its `$FIL`
keyword, so a folder of files that were all acquired under one name stops
showing up in FlowJo as a column of identical samples.

The script is a thin CLI wrapper. The work happens in
**`sync_fcs_names()`** ([`src/metadata_qc/sync_fcs_names.R`](../../src/metadata_qc/sync_fcs_names.R)),
on top of the FCS keyword layer in
[`src/R_shared/fcs_keywords.R`](../../src/R_shared/fcs_keywords.R).

## Prerequisites

**R 4.3+ and nothing else.** No conda, no Bioconductor, no flowCore, no CytoML.
Only the CLI wrapper needs `optparse`; on a machine without it, source the two
modules and call `sync_fcs_names()` directly. This is deliberate — the tool's
job is to repair files *before* the cytometry stack will accept them, so it
cannot depend on that stack.

---

## Read this before you use it

### 1. The problem

A plate loader or a BD Diva tube template acquires a whole run under one file
name. You rename the files on disk afterwards, but the original name is *also*
recorded inside each file, in the `$FIL` keyword — and `$FIL` is what FlowJo
displays by default. So the workspace shows n rows all called `Specimen_001`.

FlowJo cannot fix this itself:

> "Keyword additions or changes made in FlowJo are implemented only in the
> workspace; the raw data is never altered."
> — [FlowJo: Keywords and Annotation](https://docs.flowjo.com/flowjo/workspaces-and-samples/keywords-and-annotation/)

Re-exporting from FlowJo would re-encode every event. Hence this tool.

### 2. Check the FlowJo preference first — you may not need this at all

Open **FlowJo → Preferences → Workspace → "How should data files be named"**.
There are four options:

| Setting | What FlowJo shows | What it means for you |
|---|---|---|
| **FCS Default** (the shipped default) | the `$FIL` keyword | this tool is the fix |
| **Use Short File Name** | the on-disk file name | **you do not need this tool** — flip the setting and you are done |
| **Use Full Path and File** | path + name | documented as being "for resolving many samples with the same name" |
| **Keyword Combination** | keywords you pick | find out *which* keywords and pass them to `--also`; the default keyword set will change nothing you can see |

"Keyword Combination" over `$SRC` + `Tube Name` is standard advice in BD Diva
shops — which is exactly the population that ends up with duplicate names. It
is worth the sixty seconds to look.

### 3. Run the census before you write anything

`--probe` is read-only and is the default posture of the tool. Besides the
per-file layout it prints two things worth reading carefully:

- **which identity keywords are actually duplicated** — `$FIL`, `GUID`, `$SRC`,
  `TUBE NAME` and the on-disk name are each checked separately, so you can see
  which one FlowJo is tripping over rather than assuming.
- **groups of files sharing both `$FIL` and `$TOT`.** This one outranks the
  naming problem. When the on-disk name does not match a workspace sample,
  CytoML falls back to matching on `$FIL` and then tie-breaks on event count;
  this repo passes `greedy_match = TRUE`
  ([`src/loading/flowjo_workspace.R:819`](../../src/loading/flowjo_workspace.R)),
  which takes the first `$TOT`-correct hit instead of erroring on ambiguity.
  Acquisitions stopped at a fixed event count have identical `$TOT`. If that
  group is non-empty, **gates may already have been applied to the wrong
  patient's events in past runs**, and that is worth investigating separately.

---

## Usage

```bash
# 1. Census. Read-only.
Rscript scripts/sync_fcs_names.R --fcs_dir data/renamed --probe

# 2. Dry run: exactly what would change, per file. Still writes nothing.
Rscript scripts/sync_fcs_names.R --fcs_dir data/renamed

# 3. Apply in place, recording an undo ledger next to the data.
Rscript scripts/sync_fcs_names.R --fcs_dir data/renamed --apply --in_place

# Undo the last run. Repeat it to walk back through earlier runs.
Rscript scripts/sync_fcs_names.R --fcs_dir data/renamed --revert

# See what every file has stored inside it, without changing anything.
Rscript scripts/sync_fcs_names.R --fcs_dir data/renamed --probe \
  --keyword_table keywords.csv
```

Nothing is written without `--apply`, and `--apply` still refuses to guess
between editing the originals and writing copies — you must say `--in_place` or
`--output_dir`.

### Seeing what is in there

The TEXT segment is a delimiter-packed blob, so `--keyword_table` is how you
read a whole cohort at once. It emits long format — one row per file/keyword —
because files need not carry the same keywords, and a consumer that wants a
grid can pivot far more easily than one that wants tidy rows can unpivot. Two
extra columns save the consumer from hard-coding anything:

- `keyword_rank` — 1..n for the identity keywords in a sensible order, 999 for
  everything else, so you can sort the interesting ones to the top.
- `name_like` — whether a value shared by two files is a *problem*. A cohort
  sharing `$TOT` or `$DATE` is ordinary; a cohort sharing `$FIL` or `GUID` is
  the thing this tool exists to fix. `ORIGINALFIL` is deliberately not
  name-like: it is supposed to be identical everywhere.

The [Dash app](../../app/README.md) renders exactly this, as a table with the
keyword column frozen.


### Options

| Flag | Default | Meaning |
|---|---|---|
| `--fcs_dir` | `data` | Folder, searched recursively and case-insensitively (`.fcs` and `.FCS`) |
| `--probe` | off | Census only. Writes nothing. |
| `--list_keywords` | off | Dump every keyword of the first few files, then exit |
| `--apply` | off | Actually write |
| `--in_place` | off | Patch the originals where they sit |
| `--output_dir` | — | Write patched copies here instead, mirroring subfolders |
| `--also` | — | Comma-separated extra keywords to sync, e.g. `'$SRC,TUBE NAME'` |
| `--set` | — | Semicolon-separated `KEY=TEMPLATE`; templates take `%name% %stem% %oldname% %relpath% %dir%` |
| `--unset` | — | Comma-separated keywords to remove |
| `--guid` | `never` | `never` / `auto` / `always` — see below |
| `--strip_extension` | off | Write `$FIL` without the `.fcs` suffix |
| `--keyword_table` | — | Write every keyword of every file to this CSV (`file, keyword, value, keyword_rank, name_like`). Read-only, and done before any write, so it works with `--probe` |
| `--no_original` | off | Do not preserve the replaced `$FIL` in `ORIGINALFIL` |
| `--revert` | off | Undo one previous run from its ledger |
| `--allow` | — | Refusal codes to downgrade, e.g. `crc_present` |
| `--ledger` | `fcs_name_sync_ledger.csv` | Undo ledger, written beside the data |
| `--no_recursive` | off | Stay at the top level |
| `--src_dir` | auto | Path to the repo `src` directory |

Run `--help` for the same list.

### Which keywords get written

Two keywords by default: **`$FIL`**, set to the file name *with* its `.fcs`
extension (which is what `$FIL` means in the standard and what acquisition
software writes), and **`ORIGINALFIL`**, which preserves the name being
replaced. Everything else is opt-in, on purpose:

| Keyword | Default | Why |
|---|---|---|
| `$FIL` | **written** | FlowJo's display name under the shipped preference. Created if absent. |
| `ORIGINALFIL` | **written** | The name being replaced, so the acquisition name survives the rename and can be read back. Write-once: a second run leaves it alone, which is both what "original" should mean and what keeps repeat runs idempotent. Nothing is written for a file that never had a `$FIL` — there is no original to preserve, and recording the new name as if it were the old one would be a lie. No `$` prefix: that namespace belongs to the standard (FCS 3.1 §3.2.12), which is why FlowJo's own provenance keyword is `ORIGINALGUID` and not `$ORIGINALGUID`. Turn it off with `--no_original`. |
| `FILENAME` | `--also` | A `flowCore::read.FCS` artefact recording the absolute path it read. FlowJo does not use it and it is not in FlowJo's keyword list. |
| `GUID` | `--guid` | CytoML's sample matcher never reads it and nothing in this repo derives a sample name from it. `auto` rewrites it **only** when its current value is already an echo of the old file name; a real vendor GUID is left alone, and a fresh UUID is never minted. Run `--guid auto` on the *first* pass — after a sync, the GUID no longer echoes anything recognisable. |
| `$SRC`, `TUBE NAME`, `$SMNO`, `WELL ID`, … | `--also` | Acquisition facts. Needed when FlowJo is set to "Keyword Combination". |
| `$TOT`, `$PAR`, `$PnB/N/S/E/R`, `$DATATYPE`, `$BYTEORD`, `$MODE`, `$CYT`, `$CYTSN`, `$DATE`, `$BTIM`, `ORIGINALGUID` | **never** | Structural, hardware and audit facts. Not writable even via `--set`. |
| `$BEGINDATA`, `$ENDDATA`, `$BEGINANALYSIS`, … | recomputed | Owned by the writer. |

Nothing is stamped into the files to record that the tool ran. A "renamed by /
renamed on" pair costs ~120 bytes, which forces the slow rewrite path on every
zero-slack file, re-injects the colliding value into a field a LIMS will index,
and makes a second run non-idempotent. **The ledger CSV is the provenance
record**, and it is what `--revert` replays.

---

## How it writes

The event data is **never decoded**. The DATA segment is an opaque byte range
addressed by ASCII offsets in the 58-byte header, so it only ever has to be
located and, at worst, moved. Which of three strategies applies is decided per
file by how much room sits between the start of TEXT and whatever follows it:

| Strategy | When | What happens |
|---|---|---|
| `neutral` | the new TEXT is no longer than the old | one write of ≤ a few KB. DATA never moves, the file never changes size, `$BEGINDATA`/`$ENDDATA` keep their exact original bytes. |
| `patch` | it is longer but fits the slack before DATA | identical write, identical guarantees |
| `rewrite` | it does not fit | the file is re-emitted: patched header, new TEXT, then `[body_start, EOF)` copied **verbatim** in 8 MiB chunks |

A file that already carries `$FIL` usually lands on `neutral` or `patch`, which
is O(1) no matter how big the file is. Adding a `$FIL` that was never there
means the segment grows, which usually means a rewrite.

Freed bytes are filled with spaces, which is what the standard requires for
unused space between segments — so no fragment of the old TEXT is ever left
where a lenient parser could reach it.

### What is checked after every write

Five independent witnesses, all fatal, all before any rename is committed:

1. the *original* TEXT round-tripped through the parser byte-for-byte (without
   this, every other check is self-referential)
2. the length change matches an independently computed per-keyword byte sum
3. the ordered `(key, value)` byte sequence is exactly what was planned
4. **every byte** of `[body_start, EOF)` is identical, compared in chunks — not
   a sample, because a 4 KB spot check cannot see an overrun at byte 2e9
5. the data segment's recorded extent is unchanged, header and TEXT offsets
   still agree, and the `$TOT × $PnB` invariant still holds if it held before

If any of them fails, the write is rolled back and the file is left exactly as
it was.

---

## Notes & troubleshooting

**"Refused: nextdata_nonzero"** — the file holds more than one data set. Shifting
anything in the first one would corrupt the rest, so it is left alone.

**"Refused: crc_present"** — the file carries a real CRC, which covers the header
and TEXT, so any edit invalidates it. `--allow crc_present` proceeds and zeroes
the CRC field (the standard's own opt-out), and says so.

**"Refused: text_roundtrip_failed"** — the TEXT segment does not survive a
parse/serialize round trip, usually a delimiter inside a value that was not
doubled as the standard requires. The file is malformed; nothing is written.

**"FATAL: these files would STILL collide after the sync"** — two files in
different subfolders share a basename, so syncing both would recreate exactly
the duplicate you came here to remove. There is no `--force`: rename one of them
on disk first.

**"a previous run left these behind"** — a `.fcsprefix` or `.fcstmp` file from a
run that died mid-write. Look at it before deleting: for an interrupted in-place
patch the `.fcsprefix` sidecar holds the only copy of the original header and
TEXT.

**"FlowJo still shows the old name"** — check the preference (above). Also note
that a sample already in a workspace may be showing FlowJo's *cached* copy of
the keyword; remove and re-add it.

**My gating disappeared** — a `.wsp` records a frozen copy of `$FIL` from the
moment the sample was added, and this repo's `resolve_fcs_names()` pairs a
workspace sample to a file **by on-disk basename**
([`src/loading/flowjo_workspace.R:192-234`](../../src/loading/flowjo_workspace.R)).
Renaming the files on disk already broke that, before this tool ran. Syncing
`$FIL` additionally removes CytoML's fallback of matching a sample to a file
whose `$FIL` still carries the old label. If you need those workspaces rather
than re-gating: re-open each one in FlowJo, re-point its samples, and re-save
*after* running this, so FlowJo refreshes its own copy. The tool prints this
warning whenever it finds a `.wsp` under the data folder — it warns, it does not
refuse.

**Revert does not restore the byte count.** It restores every keyword and every
event byte exactly. If the forward pass had to move DATA, undoing it shrinks
TEXT again, and the safe way to write that is the in-place patch — which leaves
DATA where it now is and fills the freed bytes with spaces. The file stays the
larger size with a little legal padding.

**Hard links are not detected.** Base R's `file.info()` does not expose
`st_nlink`. An in-place patch of a hard-linked file edits every name pointing at
that inode. Symlinks *are* detected and refused.

---

## Calling it from R

```r
source("src/R_shared/fcs_keywords.R")
source("src/metadata_qc/sync_fcs_names.R")

# census only
res <- sync_fcs_names("data/renamed", probe_only = TRUE)
res$probe        # one row per file
res$collisions   # exact / near / non_ascii

# apply
sync_fcs_names("data/renamed", apply = TRUE, in_place = TRUE)
```

The keyword layer is usable on its own, and is the cheapest way in this repo to
read a keyword without loading a `flowFrame`:

```r
source("src/R_shared/fcs_keywords.R")
fcs_keyword("sample.fcs", "$TOT")      # "50"
fcs_print_keywords("sample.fcs")       # everything
fcs_inspect("sample.fcs")$refusals     # is it safe to edit?
```

## Tests

```bash
Rscript test_data/test_fcs_keywords.R      # the byte-level layer
Rscript test_data/test_sync_fcs_names.R    # the folder-level tool
```

Both build their fixtures byte by byte and need neither flowCore nor CytoML nor
optparse.

## In the GUI

Everything above is also the **FCS names** tab of the
[Dash front-end](../../app/README.md) — same script, same flags, and it prints
the command it is about to run. The census is the mode it opens on, and the
write mode is deliberately not remembered between sessions.

## Related

- [Metadata QC Reference](../subprojects/metadata-qc.md)
- [File Formats](../file-formats.md)
- [R_shared API Reference](../r-shared-reference.md)
