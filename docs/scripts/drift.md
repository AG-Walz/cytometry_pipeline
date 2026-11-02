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
| `group_*` | you, any number of columns | protected biological covariates (wanted variation), for example `group_prior_peptide_vaccination` with values 1 or 2. Any column starting with `group_` becomes a fixed effect that the technical terms are never allowed to absorb |
| `tech_*` | you, any number of columns | further technical factors (unwanted variation), for example `tech_operator` with the name of the person who acquired the tube. Each becomes a fixed contrast against its reference level, is tested by the decision rule like the instrument and appears in the variance partition. A blank cell becomes the level `unknown`. A `tech_*` column replaces an FCS header field of the same name (`tech_operator` replaces `$OP`); list the factor in `design.technical` after `acq_day` and `instrument` |
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
Two optional fields: `design.reference_levels` sets the contrast level of a `tech_*` factor (for example
`{operator: Patrick}`; default the first level in sort order), and `design.time_df` sets how flexible the
response over time is (default 5: a natural spline on log(1 + week) next to the linear slope; 1 keeps a
straight line). The spline matters: on a synthetic cohort with the MTB design, a response peak after a boost
and no day drift, the straight line reported 0.13 to 0.94 logit² of day variance; the spline 0.002 to 0.018.
A technical term listed in `design.technical` without a `tech_*` column is read from `tubes.csv`, so the
FCS operator keyword works without a sheet column (`technical: [acq_day, instrument, operator]`).

Three fields added on 2026-10-13 (`docs/plans/drift_methods_review.md`):

- `design.visit`: the columns of one sample, normally `[patient, week]`. Every model whose tubes repeat
  within a visit (endpoints, technical markers on all tubes, the variance partition) then gets a visit
  intercept. Without it, whatever is specific to one sample is booked as acquisition day whenever a day holds
  a single subject, and the tube residual used as rule 1 noise is too large. In the partition the visit
  share is its own category, neither biology nor batch, because it mixes both.
- `technical_tubes`: `all` fits the technical markers on every tube with a `stimulated` term and the visit
  intercept (step 4 always summarises all tubes); `neg` keeps the unstimulated tubes only.
- `relevance.count_cv`: the rule 1 margin for counts, a fold change of 1 + CV (default 0.35, the
  inter-laboratory CV of ICS responses, Jaimes et al. 2011). Rule 1 then has three states: detected,
  negligible and inconclusive.

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
| 4 | tube summaries from the FlowJo gates, and the parent and cytokine gate positions | `neg_summaries.csv`, `neg_emd.csv`, `gate_bounds.csv`, `neg_cells/` | 10 to 30 min |
| 5 | gated frequencies with matched controls | `gated_frequencies.csv`, `all_gate_counts.csv` | 30 min |
| 6a | brms yardstick | `effects_brms_{A,BG,B}.csv`, `covariates_brms.csv`, `fits_brms_summary.csv` | hours |
| 6b | hand-written Stan | `effects_stan_{A,BG,B}.csv`, `covariates_stan.csv` | 1 h |
| 6c | lme4 variance components and assumption check | `lme4_variance_components.csv`, `assumption_check.csv` | minutes |
| 6d | biological effect: week slope with and without technical terms, leave-one-day-out, and the drift-free slope from subject-days with several visits | `week_slope_sensitivity.csv`, `week_slope_lodo.csv`, `week_slope_within_day.csv` | under an hour |
| 7, 8 | BayesFlow (only with `--with_bayesflow --bf_env NAME`) | `effects_bayesflow_v{1,2}.csv`, `bf_*` | hours |
| 11 | MFI of every marker in the neg population and the parent populations, all tubes | `mfi.csv` | 30 min |
| 11b | acquisition audit: every FCS keyword, detector voltage per marker, compensation matrix, CS&T setup, FlowJo display transformation per sample and channel, gate positions per sample, sample quality (viability, live cells) and the overlap of all these with the design factors | `acquisition_settings.csv`, `sample_quality.csv`, `flowjo_transformations.csv`, `gate_positions.csv`, `factor_overlap.csv`, `fcs_keywords.csv.gz` | 10 min |
| 12a | variance partition per readout and representation, lme4 (plus refits without each extra technical factor) | `variance_partition_lme4.csv` | minutes |
| 12b | the same with brms, posterior intervals (`--workers` fits in parallel, 4 chains each) | `variance_partition_brms.csv`, `partition_fits_brms.csv`, `partition_brms/` | hours |
| 13 | responder calls by the lab rule, Fisher's exact test and MIMOSA | `calls_compare.csv`, `calls_agreement.csv` | minutes |
| 14 | figures and tables comparing the ten data settings (MFI linear / arcsinh, frequency % / arcsinh % / logit, each plain and stimulated minus unstimulated; the responder call is a decision, not a setting): biological against batch share. Figures, in reading order: `01_analysis_workflow`, `02_CD4_IFNg_per_tube_by_acquisition_day`, `03_combinations_ranked_biology_vs_batch`, `04_combinations_kept_and_dropped` (gating level x setting, with the keep/drop rule), `05_frequency_biology_vs_batch_by_source`, `06_mfi_biology_vs_batch_by_source`, `07_frequency_all_sources_of_variance`, `08_mfi_all_sources_of_variance`, `09_frequency_per_readout_biology_vs_batch`, `10_mfi_per_readout_biology_vs_batch_with_negative_controls`, `11_engines_brms_vs_lme4_variance_shares`, `12_engines_brms_vs_stan_correction_effects`, `13_bayesflow_vs_brms_run_<date>` (only when an earlier BayesFlow run is kept in `bayesflow_<date>/`). Gating levels: cytokine MFI on all events, live cells and the CD4/CD8 gate (step 11); frequencies in the manual FlowJo cytokine gate | `partition_grid.csv`, `partition_grid_summary.csv`, `partition_engine_agreement.csv`, `figures/*.png`, `figures/*.pdf` | a minute |
| 15 | sensitivity analyses: calendar-time smooth next to the day intercept, manual against fixed-threshold gates, call stability under counting noise and the minimum parent count, technical markers without the live gate (`--ungated_dir`: step 4 run with `neg_population: root` into that folder) | `sensitivity_calendar.csv`, `sensitivity_fixed_gates.csv`, `call_stability.csv`, `call_stability_summary.csv`, `min_parent_sensitivity.csv`, `sensitivity_ungated.csv` | minutes |
| 9 | decision rules | `decisions.csv`, `controls_summary.csv`, `gate_drift.csv` (parent and cytokine gates), `fitter_agreement.csv` | a few minutes |
| 10 | report | `report.html` | minutes |

Sources in the variance partition and why each is (or is not) a model term:

- **Batch:** acquisition day (random), the technical factors (instrument, operator; fixed) and **sample quality** (viability and live-cell count of the visit's unstimulated tube; the peptide does not change viability, so the unstimulated tube describes the thawed sample).
- **Biology:** patient, weeks since vaccination (log-time spline), protected covariates, cocktail (stimulated or not, and which pool).
- **Visit** (one patient at one time point, random): its own category. A visit determines the day, the week and the patient, so its variance holds both sample handling and visit-level biology and cannot be counted as either.
- **Counting noise:** for the percentage and arcsinh frequency settings, the variance a frequency has from counting k of n cells by chance, computed from the counts and split out of the unexplained rest (the logit count models exclude it by construction; MFI medians of thousands of cells carry almost none).
- **Not model terms:** detector voltages, compensation matrix and CS&T setup are almost fully determined by the acquisition day (`factor_overlap.csv`), so the model cannot separate them from it; they explain why days differ and are shown as audit figures. The gate position is determined by the day and set after looking at the data, so adjusting for it would count the day twice; it is a diagnostic only.
- **MFI on the FlowJo display scale** (`median_flowjo` in `mfi.csv`): the median through each sample's own FlowJo transformation as stored in the workspace.

The arcsinh frequency settings use `freq_asinh_cofactor` from the config, in percent (default 0.05, the
minimum subtracted frequency of the positivity rule): linear below it, logarithmic above.

Steps 4, 5 and 11 read the FlowJo workspaces through flowWorkspace, which writes HDF5 files to the
temporary folder; on a machine with a small root disk set `TMPDIR` to a folder on the data volume.
MIMOSA (step 13) is no longer in Bioconductor 3.22; install it with
`remotes::install_github("RGLab/MIMOSA")`. Without it step 13 compares the lab rule and Fisher's test only.

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

`<out>/report.html` opens in a browser and leads with three answers
(`docs/plans/drift_three_questions_plan.md`):

1. Do the data need correcting? Per technical term and layer, the contrasts detected, negligible or
   inconclusive against the rule 1 margin (count CV for counts, same-sample tube noise for technical
   markers), and the responder calls they flip over posterior draws against the flips that counting noise
   alone produces on the same tubes (rule 2).
2. Should the unstimulated tube be subtracted? The unwanted and wanted variance shares of the gated
   counts, the subtracted frequency, the frequency over the neg and the responder call, per readout, and
   the agreement of the calling rules.
3. Where do the differences come from? The variance partition per readout and representation, with
   the technical CV implied by the unwanted variance against the CLSI H62 / O'Hara targets and a
   Levey-Jennings view of the unstimulated tubes.

Then the biological effect: the week slope per
endpoint with and without the technical terms and the leave-one-day-out sensitivity from step 6d,
followed by the instrument and day verdicts per layer (A technical, BG background, B endpoint) from
`decisions.csv`, the control matching, the gate-position drift (parent and cytokine gates), the
sensitivity analyses of step 15 and the agreement between fitters.
`decisions.csv` is the machine-readable version.

## Backwards compatibility

`bash scripts/drift/run_all.sh DATA_ROOT OUT CFG [r_env] [bf_env]` still works and now calls
`drift.R run`; the environment variables `FCS_DIR`, `WSP_DIR`, `SAMPLE_SHEET`, `SKIP_BF`, `CORES` and
`WORKERS` keep their meaning. The individual `NN_*.R` scripts under `scripts/drift/` are unchanged and
can still be called on their own with `--out` and `--config`.
