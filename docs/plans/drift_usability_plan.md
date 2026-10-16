[< Back to Documentation Index](../README.md)

# Drift analysis: one sample sheet, one command, and a test that biology and acquisition time are separable

Status: plan written 2026-09-24 after the generalisation described in section 11 of `mtb_batch_drift_plan.md`; implementation of steps 1 to 5 started the same day (see section 8 for the outcome).

## 1. What the user wants to measure, and what gets in the way

The biological quantity of interest is the change in peptide-specific T cell responses over the weeks after the first vaccination. The technical nuisance is acquisition drift: the day a tube went through the cytometer, and which cytometer. Within one patient these two clocks move together by design, because every visit is a new acquisition day. The analysis can only separate them because

- patients are acquired on different days from each other, so the day effect is shared across patients while the week effect is not,
- several visits of one patient were sometimes acquired on the same day (15 patient-days carry two or three weeks), which pins the week effect at a fixed day,
- the unstimulated tubes give a signal that should not depend on the week at all.

Today the user has to know this. The pipeline should make it visible, test it, and refuse to run when the two clocks cannot be told apart.

## 2. One sample sheet as the single metadata input

Metadata currently comes from a file-name regex or an optional CSV, and the technical fields are read from the FCS header separately. Replace this with one sheet that the pipeline pre-fills and the user completes, so biological and technical time sit side by side in the same row.

`drift init --fcs_dir DIR --out DIR [--filename_regex ...]` writes `sample_sheet.csv`:

| Column | Filled by | Meaning |
|---|---|---|
| `filename` | pipeline | the FCS file |
| `acq_date`, `instrument`, `n_events`, `software` | pipeline, from the FCS header | technical clock and cytometer; read-only, shown so the user sees them next to the biology |
| `subject` | pipeline guess from the regex, user confirms | the repeated-measures unit |
| `visit` | user | visit label as the clinic writes it (`V3W6`) |
| `weeks_since_first_vaccination` | user, or derived from `visit` by a rule in the config | the biological clock, numeric |
| `stimulation` | pipeline guess, user confirms | peptide cocktail id, or `unstimulated` |
| `control_for` | pipeline default, user may edit | which rows this unstimulated tube is the background for: `subject+visit` by default |
| `group_*` | user, any number | protected biological covariates such as `group_prior_peptide_vaccination` (1 or 2) |
| `exclude`, `note` | user | manual exclusions with a reason |

Rules: a column prefixed `group_` is automatically a protected fixed covariate; `subject`, `visit`, `weeks_since_first_vaccination`, `stimulation` and `control_for` are fixed names; everything read from the header is regenerated on each run and a mismatch with the sheet is a warning, not an error. The config `design` block then needs no column names at all for the standard cohort, and `tubes.filename_regex` becomes optional, used only to pre-fill guesses.

The current MTB sheet is generated from the file names in one command, which is the migration path, and the file-name path keeps working for scripted use.

## 3. Checks that run before any model is fitted

`drift check --sample_sheet sample_sheet.csv --fcs_dir DIR --wsp_dir DIR --config CFG` runs in under a minute and stops the pipeline on hard failures:

1. Every FCS file is in the sheet and every sheet row has a file; every row has a subject, a visit and a stimulation; unstimulated rows exist for every `control_for` key or a fallback in `controls` covers them.
2. FlowJo workspaces reference every file, and the configured gate names exist in every workspace (the empty CD8 gate of patient 23 week 208 would be reported here, not discovered in step 9).
3. Marker names are consistent across files, or the harmonisation tool has been run.
4. Separability of the two clocks, the test the user asked for:
   - **Bridge report**, already produced: which subjects span days and instruments.
   - **Confounding index**: the share of the variance of `weeks_since_first_vaccination` explained by `acq_date` within subject. At 1.0 the two are the same variable for every subject and the day term is unidentifiable without the unstimulated layer; the MTB cohort sits below 1.0 only because of the 15 multi-week days. Report the index, the number of multi-week days and the number of subjects acquired on more than one day.
   - **Design power by simulation**: with the real sheet, simulate an endpoint with a known week slope and a known day shift (the synthetic generator already does this), fit the endpoint model with lme4, and report whether both are recovered with less than 25 percent leakage from one to the other. Also simulate a week effect with zero day shift and report the false detection rate of day effects, and a day shift with zero week effect and report the bias in the week slope. Ten repetitions each, a minute in total.
   - A verdict line: `separable`, `separable only through the unstimulated layer`, or `not separable: add bridging acquisitions`.

## 4. One command with resumable steps

`drift run --sample_sheet ... --out ...` replaces the bash runner. Each step declares its inputs and outputs; a step runs only when an output is missing or older than an input, and `--from 6` or `--only 9,10` overrides. `drift status --out DIR` prints what exists and what is stale. Heavy steps write a `.done` marker with the git commit and config hash so a report can state exactly what produced it. Environment comes from one lock file per platform (`environment_drift.yml` plus the BayesFlow env), and `drift doctor` checks the packages, CmdStan and the workspace reader before anything is run. The BayesFlow steps stay opt-in (`--with_bayesflow`).

Implementation: a thin R entry point `scripts/drift/drift.R` with sub-commands `init`, `check`, `run`, `status`, `doctor`, reusing the existing step scripts unchanged; the steps already take `--out` and `--config`, so resumability is a wrapper concern.

## 5. Making the biological effect the headline of the report

The week slope is currently reported among the protected covariates. It should lead: a first section "Biological effect over weeks since first vaccination" with the population slope and the per-subject slopes for every endpoint, with and without the technical terms in the model. The two estimates side by side answer the user's question directly: if removing `acq_date` and `instrument` does not move the week slope, technical drift is not what the response is measuring. Add a leave-one-day-out sensitivity: refit the endpoint model dropping each acquisition day in turn and plot how far the week slope moves. Both use the lme4 fitter so they run in minutes.

## 6. Order of work and effort

| Step | Deliverable | Effort |
|---|---|---|
| 1 | `drift init` and the sample-sheet schema, MTB sheet generated from file names, step 1 reads it | half a day |
| 2 | `drift check` with the separability tests and the simulation | one day |
| 3 | `drift run`, `status`, `doctor` | half a day |
| 4 | Report headline section and leave-one-day-out sensitivity | half a day |
| 5 | Documentation walk-through for a new cohort, from folder of FCS files to report | half a day |

Steps 1 and 2 can start now against the current results. Step 4 depends on the endpoint fits of the 2026-09-24 run only for its illustration.

## 7. What this does not solve

If a future cohort acquires every visit of every subject on its own day and never repeats a subject across days, no analysis separates week from day at the endpoint level. The check in section 3 will say so. The remedy is in the lab, not the code: acquire a bridging control (the same frozen PBMC aliquot, or a repeat of one patient's visit) on every acquisition day, and mark it `stimulation = bridge` in the sheet. The pipeline should treat such rows as pure technical replicates in layer A and as a direct estimate of the day effect.

## 8. Outcome of the implementation (2026-09-24)

All five steps of section 6 were implemented the same day. What exists now:

- `scripts/drift/00_init_sheet.R` and `src/harmonization/drift/drift_sheet.R`: the sample sheet of section 2, pre-filled from the FCS headers and the file-name regex, and a generated cohort config from `config/drift_template.yaml`. For the MTB cohort all 328 rows were guessed completely, and step 1 fed from the sheet reproduces the reference tube table exactly. The MTB sheet is at `/mnt/volume/data/MTB_FACS/mtb_drift/sample_sheet_mtb.csv`.
- `scripts/drift/check_design.R` and `drift_separability.R`: the checks of section 3. On the MTB design the confounding index is 0.98, yet the simulation recovers a week slope and a day SD injected together with 1 percent leakage each, a pure week effect produces no false day detection, and a pure day shift biases the week slope by 0.003 logit per year. The verdict lets the simulation outrank the index and reads `separable`, with the caveat that the identification rests on the 15 multi-week acquisition days. The workspace check finds the empty CD8 gate of patient 23 week 208 before any model is fitted.
- `scripts/drift/drift.R` with `init`, `doctor`, `check`, `run`, `status`; `run_all.sh` is now a wrapper around it. `environment_drift.yml` holds the R environment. Walkthrough in `docs/scripts/drift.md`.
- `scripts/drift/06_biology.R` and `drift_biology.R`: week slopes with and without the technical terms, leave-one-day-out refits and per-subject slopes, leading the report. On the MTB endpoints no week slope is distinguishable from zero in the full model, and removing the technical terms moves CD8 IFNg by 3.8 SE but none of the others by more than 1.3 SE.
- One model change fell out of this work: the count models had week only as a zero-mean per-patient random slope, so brms reported no population week effect while the Stan program did. `design_formula()` now adds the fixed week term for the background-count and endpoint classes; the Stan design matrix leaves it out because the beta-binomial program carries it as `mu_week`.

Not done: bridging-control rows (section 7) are not yet treated specially; the `exclude` column is honoured through the excluded-tube list but was not exercised on real data end to end.
