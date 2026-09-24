[< Back to Documentation Index](../README.md)

# Drift analysis for a new cohort: from a folder of FCS and FlowJo files to the report

`scripts/drift/drift.R` is the single entry point of the acquisition-drift analysis
(`docs/plans/mtb_batch_drift_plan.md` for the method, `docs/plans/drift_usability_plan.md` for the
design of this command). It answers one question per cohort: is the change in peptide-specific T cell
responses over the weeks after vaccination measurable separately from the day and the cytometer on
which each tube was acquired, and does either technical clock need a correction?

This page walks through a new cohort. Every command is run from the repository root inside the R
environment (`conda activate <env>` or `conda run -n <env> ...`).

## 0. Environment

```bash
conda env create -f environment_drift.yml           # R packages, CmdStan, CytoML; once
bash scripts/drift/00_setup_env.sh drift_env bayesflow_env   # optional: also the BayesFlow env for steps 7 and 8
Rscript scripts/drift/drift.R doctor [--bf_env bayesflow_env]
```

`doctor` prints one line per component (R packages, CmdStan, the BayesFlow env) with OK, missing or
skipped, and an install hint for anything missing. It exits non-zero when something is missing.
Steps 1 and 2 need only base R and `yaml`; steps 3 to 5 need flowCore, CytoML and flowWorkspace; step
6a needs brms, 6b cmdstanr with a visible CmdStan, 6c and 6d lme4, step 10 rmarkdown and ggplot2.

## 1. Init: pre-fill the sample sheet

```bash
Rscript scripts/drift/drift.R init --fcs_dir data/cohortX/fcs --out results/cohortX --cohort "Cohort X"
```

This reads the header of every FCS file and writes two files into `--out`:

- `sample_sheet.csv`, one row per FCS file, and
- `drift_Cohort X.yaml`, a copy of `config/drift_template.yaml` with the cohort name filled in.

Pass `--config config/drift_mtb.yaml` when the file names follow a known pattern: the
`tubes.filename_regex` of that config is then used to guess subject, week and stimulation.

### Completing the sheet

Open `sample_sheet.csv` in a spreadsheet and fill the biological columns. Columns in order:

| Column | Who fills it | Meaning |
|---|---|---|
| `filename` | pipeline | FCS file name; do not edit |
| `acq_date`, `instrument`, `n_events`, `software` | pipeline, from the FCS header | the technical clock (`acq_date`) and the cytometer. They are regenerated on every run; if you edit them the pipeline warns and keeps the header value |
| `subject` | you (guessed when a regex applies) | the repeated-measures unit, typically the patient id |
| `visit` | you | the visit label as the clinic writes it, for example `V3W6` |
| `weeks_since_first_vaccination` | you | the biological clock, a number. This is the column the whole analysis is about: the response is modelled as a slope over these weeks, and the check in step 3 asks whether that slope can be told apart from `acq_date`. Fill it from the visit schedule, not from `acq_date`: if a visit's blood was drawn in week 6 but the tube was run months later, the row says 6. Two rows with the same `subject` and week but different `acq_date` are valuable: they pin the week effect at a fixed day |
| `stimulation` | you (guessed) | peptide cocktail id, or exactly `unstimulated` for the control tube of a visit |
| `control_for` | pipeline default, you may edit | which stimulated rows an `unstimulated` row is the background for. `subject+visit` (default): the control of the same subject and week. `subject+acq_date`: the control run on the same day for that subject, for cohorts where one control serves every visit acquired together. `subject`: one control per subject. Normally one value for the whole sheet; it becomes `controls.match_keys` in the config, and the `fallbacks` listed there cover rows whose control is missing |
| `group_*` | you, any number of columns | protected biological covariates, for example `group_prior_peptide_vaccination` with values 1 or 2. Any column starting with `group_` becomes a fixed effect that the technical terms are never allowed to absorb |
| `exclude` | you | `TRUE` to drop a tube by hand (a failed stain, a wrong file). Excluded rows are listed in `excluded_manual.csv` by step 1 |
| `note` | you | free text, kept next to the exclusion |

`init` prints how many rows have every biological field filled; aim for all of them. Step 1 rejects a
sheet with non-numeric weeks, an `exclude` value other than TRUE or FALSE, or an unknown `control_for` key.

### Editing the config

Open the generated `drift_<cohort>.yaml` and work through the fields marked `EDIT`:

- `instrument_labels`: short label per cytometer, matched against the FCS `$CYT` keyword.
- `markers`: canonical marker names with the regexes that recognise them in `$PnS` and `$PnN`.
- `neg_population`, `parent_populations`, `cytokine_gate_regex`: the FlowJo gate names. List the gate
  tree first with `Rscript scripts/drift/04_neg_summaries.R --wsp_dir ... --fcs_dir ... --out ... --config ... --list_populations`
  (after step 1 and 2 have run once).
- `reference_instrument`: the label a correction, if any, would map onto.

The `design` and `controls` blocks are written from the sheet and need no editing for a standard cohort.

## 2. Harmonise marker names

Several files of one cohort often spell the same marker differently in `$PnS` (`IFNy`, `INFy`, `IFN`)
or leave it blank. Step 2 tolerates this through the `markers` regexes, but the cleaner path is to
write one spelling into the files before the first run:

```bash
Rscript scripts/harmonize_marker_names.R --fcs_dir data/cohortX/fcs --probe                 # every spelling per channel
Rscript scripts/harmonize_marker_names.R --fcs_dir data/cohortX/fcs --canonical results/cohortX/drift_cohortX.yaml   # preview
Rscript scripts/harmonize_marker_names.R --fcs_dir data/cohortX/fcs --canonical results/cohortX/drift_cohortX.yaml --apply --in_place
```

The apply run writes an undo ledger; see `docs/scripts/harmonize-marker-names.md`. If FlowJo shows
duplicate samples after files were renamed, run `scripts/sync_fcs_names.R` as well (`docs/scripts/sync-fcs-names.md`).

## 3. Check: can week and acquisition day be separated?

Step 1 has to run once so that `tubes.csv` exists, then the check reads it:

```bash
Rscript scripts/drift/drift.R run --sample_sheet results/cohortX/sample_sheet.csv --fcs_dir data/cohortX/fcs \
  --out results/cohortX --config "results/cohortX/drift_Cohort X.yaml" --only 1
Rscript scripts/drift/drift.R check --out results/cohortX --config "results/cohortX/drift_Cohort X.yaml" \
  [--wsp_dir data/cohortX/wsp --fcs_dir data/cohortX/fcs]
```

It prints the bridge report (which subjects span days and instruments), the confounding index (the
share of the variance of the week within subject that `acq_date` explains; 1.0 means the two clocks are
the same variable) and a simulation with lme4: an endpoint with a known week slope and a known day shift
is simulated on the real sheet and refitted a few times, and the leakage from one into the other is
reported. With `--wsp_dir` it also checks that every tube is in a workspace and that the configured gates
exist. It writes `separability.csv` and `separability_summary.csv` and ends with a verdict line:

- `separable`: week and day are identifiable from the stimulated tubes alone. Run everything.
- `separable only through the unstimulated layer`: the endpoint model cannot tell them apart by itself; the
  day term is pinned by the unstimulated tubes, which should not depend on the week. Results for the day
  effect rest on that assumption; read the assumption check of step 6c.
- `not separable: add bridging acquisitions`: every visit of every subject sits on its own day. No model
  separates the clocks; acquire a bridging control (the same frozen aliquot or a repeat of one visit) on
  every acquisition day and add it to the sheet before running the models.

The same check runs as step 1b of `run`, so a full run records the verdict as well.

## 4. Run

```bash
screen -S cohortX      # the model steps take hours
Rscript scripts/drift/drift.R run --sample_sheet results/cohortX/sample_sheet.csv \
  --fcs_dir data/cohortX/fcs --wsp_dir data/cohortX/wsp \
  --out results/cohortX --config "results/cohortX/drift_Cohort X.yaml" --cores 4 --workers 8
```

Steps and what they produce under `--out`:

| Step | What | Outputs | Time |
|---|---|---|---|
| 1 | tube table from the sheet and the FCS headers | `tubes.csv`, `bridge_report.csv`, `excluded_manual.csv` | seconds |
| 1b | design check (section 3) | `separability.csv`, `separability_summary.csv` | minutes |
| 2 | panel reconciliation | `panel_map.csv`, `panel_signatures.csv`, `panel_gaps.csv` | seconds |
| 3 | tube QC | `qc.csv`, `excluded_tubes.csv` | 30 min |
| 4 | unstimulated-tube summaries from the FlowJo gates | `neg_summaries.csv`, `neg_emd.csv`, `gate_bounds.csv`, `neg_cells/` | 30 min |
| 5 | gated frequencies with matched controls | `gated_frequencies.csv`, `all_gate_counts.csv` | 30 min |
| 6a | brms yardstick | `effects_brms_{A,BG,B}.csv`, `covariates_brms.csv`, `fits_brms_summary.csv` | hours |
| 6b | hand-written Stan | `effects_stan_{A,BG,B}.csv`, `covariates_stan.csv` | 1 h |
| 6c | lme4 variance components and assumption check | `lme4_variance_components.csv`, `assumption_check.csv` | minutes |
| 6d | biological effect: week slope with and without technical terms, leave-one-day-out | `week_slope_sensitivity.csv`, `week_slope_lodo.csv` | minutes |
| 7, 8 | BayesFlow (only with `--with_bayesflow --bf_env NAME`) | `effects_bayesflow_v{1,2}.csv`, `bf_*` | hours |
| 9 | decision rules | `decisions.csv`, `controls_summary.csv`, `gate_drift.csv`, `fitter_agreement.csv` | a minute |
| 10 | report | `report.html` | minutes |

Everything printed also goes to `<out>/run_all.log`, and every finished step writes
`<out>/.done/<step>.json` with its start time, duration, git commit, config checksum and the exact
command line, so the report can be traced to what produced it. A failing step stops the run with its
exit status; fix the cause and run the same command again.

## 5. Status and reruns

```bash
Rscript scripts/drift/drift.R status --out results/cohortX --config "results/cohortX/drift_Cohort X.yaml"
```

prints one line per step: how many of its outputs exist, whether it is `ok`, `missing` or `stale`, and
when it last finished. A step is stale when one of its outputs is older than one of its inputs; the config
is an input of every step, so editing it makes everything stale. `run` uses the same rule and skips
every step that is `ok`, which makes a rerun after a failure or an edit cheap: only what changed runs again.

- Rerun one step: `drift.R run ... --only 9` (only if stale) or `drift.R run ... --force --only 9` (always). Step 10 becomes stale as a consequence and runs on the next plain `run`.
- Rerun from a step onwards: `drift.R run ... --from 6c` considers steps 6c and later.
- Several steps: `--only 6c,6d,9,10`. `6` is short for `6a`.
- Skip the FlowJo steps on a machine without CytoML: copy `tubes.csv`, `panel_*.csv`, `neg_summaries.csv`, `neg_emd.csv`, `gated_frequencies.csv` from a machine that has it and run `--from 6`.

## 6. The report

`<out>/report.html` opens in a browser and leads with the biological effect: the week slope per
endpoint with and without the technical terms and the leave-one-day-out sensitivity from step 6d,
followed by the instrument and day verdicts per layer (A technical, BG background, B endpoint) from
`decisions.csv`, the control matching, the gate-position drift and the agreement between fitters.
`decisions.csv` is the machine-readable version.

## Backwards compatibility

`bash scripts/drift/run_all.sh DATA_ROOT OUT CFG [r_env] [bf_env]` still works and now calls
`drift.R run`; the environment variables `FCS_DIR`, `WSP_DIR`, `SAMPLE_SHEET`, `SKIP_BF`, `CORES` and
`WORKERS` keep their meaning. The individual `NN_*.R` scripts under `scripts/drift/` are unchanged and
can still be called on their own with `--out` and `--config`.
