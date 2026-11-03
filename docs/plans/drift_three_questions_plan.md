[< Back to Documentation Index](../README.md)

# Drift analysis, round 3: three questions, backed by clinical and methods literature, for any cohort

## Context

The 2026-09-24 run asked whether any day or the instrument was shifted beyond tube noise, and whether that flipped responder calls. Two things have changed since:

- The cohort now has operator metadata, and the user wants the analysis answerable for any cohort.
- A logic audit found that the time model and the noise threshold bias the answer.

The user's three questions:

1. **Do I need to correct the data?** Is there an effect of instrument, day or operator, and is it large enough to matter?
2. **Should the unstimulated tube be subtracted from the stimulated one?**
3. **Where do the differences between samples come from?** One value per factor:
   - unwanted: instrument, day, operator;
   - wanted: patient, weeks since vaccination, prior peptide vaccination, peptide cocktail used or not (and which pool).

Requirements: generalisable to other cohorts; tables and plots aimed at the user; every method backed by clinical, regulatory or methods publications.

Framing: CLSI H62 (2021) and O'Hara et al. 2011 ask for intermediate precision across operators, days and instruments, and for a documented decision when an instrument is changed (Selliah et al. 2019, "transfer" validation). No bridge samples exist, so the precision is estimated retrospectively from variance components, as in the HIPC multicentre study (Finak et al. 2016). The report states this limitation.

## A. Model fixes (needed by all three questions)

| Fix | Why | Where |
|---|---|---|
| Week as a spline (brms `s(week_c, k = 4)`, lme4 `splines::ns(week_c, 3)`) plus the per-patient slope | Responses rise after each boost and then wane. A straight line pushes the curve into the day effects, because each visit is a new day | `design_formula()` in `src/harmonization/drift/drift_models.R` |
| `stimulated` indicator as a wanted factor; pool identity stays `(1 \| cocktail)` (patient-specific) | Makes "cocktail used or not" estimable | design block, frames |
| Replicate noise from the same sample (patient + day + week); lme4 tube-level residual as fallback; Maecker et al. 2005 CV curve as the external benchmark | Today's noise mixes pools and visits, so it is too high and hides drift | `replicate_sd()` in `drift_stats.R`, `09_decision.R` |
| Operator as an unwanted factor (sheet `tech_*`) | Built 2026-10-12, uncommitted | done |
| Every biological covariate generic | `prior_vacc` is still hard-coded in 4 places | `05_gated_frequencies.R`, `drift_models.R` (`frame_cfg`, `merge_tubes`, layer BG columns), `09_decision.R`, `10_report.Rmd`: use `design_columns(cfg)` and the `group_*` names |

## B. Q1: do I need to correct?

The existing decision step, rerun with operator and the fixes:

- **Rule 1:** P(|effect| > same-sample noise) ≥ 0.95.
- **Rule 2:** responder calls flipped when the effect is removed.
- Applied per unwanted factor level, for:
  - neg-tube technical markers;
  - neg background;
  - endpoints, both raw stimulated and stimulated together with neg.

Outcome table per factor:

| Rule 1 | Rule 2 | Action |
|---|---|---|
| fails | — | no correction |
| passes | 0 flips | keep the factor as a covariate downstream |
| passes | flips > 0 | correct |

When a factor reaches "correct":

- For frequencies: a per-level logit offset, written to `gated_frequencies_corrected.csv`. The originals are never overwritten.
- For the instrument alone: cyCombine or CytoNorm 2.0 with the bridge patients (16, 17, 61, 74) as a sensitivity check, confirming that the stimulated vs unstimulated difference survives (Schuyler et al. 2019).
- Per-day correction tools are not used. With no anchors and day ≈ visit, they would remove week biology.

Supporting evidence in the report:

- precision by operator, day and instrument as CVs, against the CLSI H62 / O'Hara targets (≤25 %; 30–35 % for rare events);
- gate-position drift (gating by one person; Welters et al. 2012 name gating as the main analytical source of variability).

## C. Q2: should the neg be subtracted?

Representations of the same endpoints, each fitted with all seven factors (part D):

| Code | Representation | Backed by |
|---|---|---|
| R3 | raw stimulated frequency | — |
| R4 | stimulated % minus neg % | HVTN practice (Horton et al. 2007); MIATA asks for raw and net values to be reported |
| R5 | stimulated and neg counts in one beta-binomial model, plus technical terms | MIMOSA (Finak et al. 2014) |
| R6 | responder call under the lab rule (≥ 0.05 %, ≥ 3×) | Moodie et al. 2010 (pre-specified rule) |
| R6b | call by one-sided Fisher's exact test, stim vs neg | HVTN (Horton et al. 2007) |
| R6c | call by MIMOSA posterior probability | Finak et al. 2014 |
| MFI | arcsinh median against stim minus neg, per cytokine marker | Kalina et al. 2012 (MFI comparability) |

Verdict rule:

- Subtraction (or the joint model) is recommended if it lowers the unwanted variance share against R3 without lowering the wanted share beyond its interval. This is the batch-removal vs bio-conservation trade-off scored by scIB (Luecken et al. 2022).
- The neg background is partly biology, so its own variance shares show what subtraction removes. Halvorson et al. 2025 (preprint; a second literature search on 2026-10-13 could not find it, so check the reference before citing) is the only paper comparing subtraction, ratio and Box-Cox; a ratio is added as a sensitivity representation.

## D. Q3: where do the differences come from?

**Readouts:**
- MFI of every marker within CD4, CD8 and live cells, for all tubes (new step 11, CytoML);
- the 8 cytokine frequencies.

**Model, one per readout and representation:**

```
y ~ instrument + operator + (1 | acq_day)                               # unwanted
  + (1 | patient) + s(week) + prior_vacc + stimulated + (1 | cocktail)  # wanted
```

**Variance shares per posterior draw**, following variancePartition (Hoffman & Schadt 2016) for the method and HIPC (Finak et al. 2016) for the cytometry precedent:
- fixed factors: the variance of their contribution across tubes;
- random factors: their variance;
- residual: latent-scale for beta-binomial and Bernoulli models (Nakagawa, Johnson & Schielzeth 2017).

**Confounding guard:** each model is refitted without operator (lme4). The unwanted total is the robust number; the split among instrument, day and operator is approximate, because operator and instrument are nested in calendar period.

## E. Plots and tables (each modelled on a published one)

| # | Output | Answers | Modelled on |
|---|---|---|---|
| 1 | Answer box per question (one line: correct yes/no per factor; subtract yes/no; top 3 sources of variation) | all | MIATA reporting |
| 2 | Acquisition map: patient × day, coloured by instrument, shaped by operator | design | existing report |
| 3 | Technical-factor overlap table (`technical_factors.csv`) | design | new (2026-10-12) |
| 4 | Stacked bars of variance share per readout, 7 factors + residual, unwanted in warm and wanted in cool colours | Q3 | Finak 2016 Fig 3; Brodin et al. 2015 |
| 5 | Violin of variance share per factor across all readouts, one panel per representation | Q3 | Hoffman & Schadt 2016 |
| 6 | Factor × readout heatmap of variance share | Q3 | Patin et al. 2018 |
| 7 | Scatter: unwanted share (x) against wanted share (y), one point per representation, arrows from R3 | Q2 | Luecken et al. 2022 (scIB) |
| 8 | Precision table: CV by operator, day and instrument per readout, against the CLSI H62 / O'Hara targets | Q1 | O'Hara et al. 2011; CLSI H62 |
| 9 | Noise against response size: our between-tube SD, Maecker 2005 CV values overlaid, 0.05 % line | Q1, Q2 | Maecker et al. 2005 Table 2 |
| 10 | Forest plot of every technical contrast in replicate-SD units; day heatmap with flip counts | Q1 | existing report |
| 11 | Per-day MFI trend of technical markers in negs, Levey-Jennings style (mean ± 2 SD bands) | Q1 | Perfetto et al. 2012; Green et al. 2011 |
| 12 | Stim against neg frequency scatter, coloured by call, one panel per calling rule; agreement table between rules | Q2 | Finak et al. 2014 (MIMOSA) |
| 13 | Week trajectory per patient and readout, with and without technical terms | biology | existing `06_biology.R` |
| 14 | MDS of per-tube medians, coloured by day/instrument/operator and separately by patient | Q1, Q3 | Trussart et al. 2020 (CytofRUV) |
| 15 | Calls heatmap patient-visit × pool, annotated by day/instrument/operator (do calls follow batch?) | Q1, Q2 | Lin et al. 2015 (COMPASS) |

All tables are also written as CSV next to `report.html`.

## F. Generalisable

- **Sample sheet** declares the categories:
  - `group_*` columns are wanted (protected);
  - `tech_*` columns are unwanted;
  - instrument and acquisition date come from the FCS header;
  - an FCS keyword such as `$OP` can be named in the config instead of a sheet column.
- **Readouts** come from the config: `parent_populations`, `neg_population`, markers and the positivity rule. The representations to run are listed in a new `representations:` block (default all).
- **No cohort name in code.** Fix the four `prior_vacc` hard-codes. Report titles take `cohort`.
- **`drift.R`** registers the new steps (11 MFI table, 12 variance partition, 13 calling-rule comparison) and keeps the resume and status machinery.
- **Template:** `config/drift_template.yaml` gains the new blocks with comments.

## G. Lab actions (outside the code, recommended by the clinical sources)

- **Export the FACSDiva CS&T bead history** (Levey-Jennings) for all 23 days and both instruments, as instrument-stability evidence (Perfetto et al. 2012; Green et al. 2011).
- **Run a bridge sample from now on:** one frozen PBMC aliquot on every acquisition day, or at least on any instrument change (Selliah et al. 2019; Le Lann et al. 2020). The pipeline already plans `stimulation = bridge` rows.

## Implementation

- **New:**
  - `src/harmonization/drift/drift_partition.R` (frames R3–R6c, `variance_shares()`);
  - `scripts/drift/11_mfi_table.R` (reuses `find_population()`, `asinh_cof()`, `marker_summary()` from step 4);
  - `scripts/drift/12_variance_partition.R` (brms or `--engine lme4`, resumable);
  - `scripts/drift/13_calling_rules.R` (Fisher and MIMOSA; MIMOSA installed into `mtb_drift` from Bioconductor).
- **Changed:** `drift_models.R` (spline, stimulated, generic covariates), `drift_stats.R` (`replicate_sd` grouping), `09_decision.R`, `10_report.Rmd` (three answer sections first, the catalogue above), `drift.R`, `config/drift_template.yaml`, `config/drift_mtb.yaml`.
- **Docs:**
  - `docs/plans/drift_three_questions_plan.md` (this plan, with the reference list and DOIs);
  - `docs/scripts/drift.md` (sheet categories, representations, reading the report).
- **Optional later:** COMPASS (needs Boolean combination counts); CytoBatchFlagR cross-check.

## References (DOIs checked by search; figure details partly from abstracts)

CLSI H62-Ed1 (2021) · O'Hara 2011 10.1016/j.jim.2010.09.036 · Green 2011 10.1016/j.jim.2010.07.004 · Wood 2013 10.1002/cyto.b.21108 · Selliah 2019 10.1002/cpcy.53 · Welters 2012 10.1007/s00262-012-1282-9 · Horton 2007 10.1016/j.jim.2007.03.002 · Maecker 2005 10.1186/1471-2172-6-13 · Finak 2016 10.1038/srep20686 · Finak 2014 (MIMOSA) 10.1093/biostatistics/kxt024 · Lin 2015 (COMPASS) 10.1038/nbt.3187 · Moodie 2010 10.1007/s00262-010-0875-4 · Janetzki 2009 (MIATA) 10.1016/j.immuni.2009.09.007 · Perfetto 2012 10.1038/nprot.2012.126 · Kalina 2012 10.1038/leu.2012.122 · Le Lann 2020 10.1038/s41598-020-68468-3 · Hoffman & Schadt 2016 10.1186/s12859-016-1323-z · Nakagawa 2017 10.1098/rsif.2017.0213 · Luecken 2022 10.1038/s41592-021-01336-8 · Brodin 2015 10.1016/j.cell.2014.12.020 · Patin 2018 10.1038/s41590-018-0049-7 · Trussart 2020 10.7554/eLife.59630 · Schuyler 2019 10.3389/fimmu.2019.02367 · Pedersen 2022 (cyCombine) 10.1038/s41467-022-29383-5 · Quintelier 2025 (CytoNorm 2.0) 10.1002/cyto.a.24910 · Halvorson 2025 (preprint, unverified on 2026-10-13) 10.1101/2025.07.15.664847

## Verification

1. **Unit tests** (`test_data/test_drift_partition.R`, base R):
   - `variance_shares()` recovers known components within 10 % and sums to 1;
   - R4 equals stim minus matched neg;
   - R6 equals `positivity_call()`;
   - the spline formula builds for brms and lme4;
   - a cohort config with a different covariate name runs (generalisation).
2. **Synthetic data** (`drift_synthetic.R` with operator effect and curved week):
   - shares recovered;
   - a pure week curve gives no detected day effect;
   - an injected operator shift that flips calls yields "correct", and the corrected table removes it.
3. **Real data:** lme4 on everything first (minutes), then brms in a screen on the dev box into `results_mtb_drift_2026-10-12_operator/`. Check:
   - R-hat ≤ 1.01 and no divergences;
   - brms and lme4 agree;
   - a verdict-by-verdict comparison with 2026-09-24.
4. **Regression:** existing 182 drift tests green.

## Status on 2026-10-12 (implementation day)

Implemented and tested (206 base-R tests, 24 of them new in `test_data/test_drift_partition.R`):

- Sheet `tech_*` columns, operator as a fixed technical factor, `design.reference_levels`.
- **Time model.** `design.time_df` (default 5): natural spline of log(1 + week) next to the linear slope. A 3-df spline on the week scale was tried first and did not help. Synthetic cohort with the MTB design, a response peak around week 12 and no day drift, 9 glmer fits per model:

  | Time model | Day variance (truth 0) | Operator effect (truth 0.8) |
  |---|---|---|
  | straight line | 0.13 to 0.94 | 1.16 to 1.85 |
  | log-time spline | 0.002 to 0.018 | 0.42 to 0.85 |

  Without a curve, both models give a day variance of about 0. Through `06_lme4.R` the layer B day variance fell from about 0.27–0.30 (line) to 0.004–0.015 (spline). The real-data design check with operator and spline reads `separable`: 0 % false day detection (13 % with operator and the straight line).
- **Rule 1 noise.**
  - Technical markers: the SD between tubes of the same sample.
  - Counts: the lme4 tube-level residual.
  - The source used is recorded in `decisions.csv` (`noise_source`).
- **Step 11** `mfi.csv`: 327 tubes × 13 readouts. flowWorkspace needs `TMPDIR` on the data volume; it crashed with an HDF5 read error when `/tmp` was used.
- **Step 12** variance partition, lme4 and brms. Counts are binomial with a per-tube random effect in both engines. A beta-binomial was tried first; its exact logit-scale tube variance explodes for rare events (0.62 residual share against 0.18 in lme4 on CD4 IFNg). With the per-tube term, brms and lme4 agree on CD4 IFNg (unwanted 0.37 against 0.33, wanted 0.49 against 0.49).
- **Step 13** calling rules. MIMOSA is not in Bioconductor 3.22 and was installed from GitHub (RGLab/MIMOSA). First counts:

  | Response | lab rule | Fisher | MIMOSA |
  |---|---|---|---|
  | CD4 IFNg | 127 | 151 | 164 |
  | CD4 TNFa | 72 | 147 | 147 |

  Fisher and MIMOSA call far more CD4 responses than the lab rule.
- **Report:** three answer sections, the precision table against CLSI H62 / O'Hara, the noise-vs-response-size figure with the Maecker 2005 values, and the Levey-Jennings view.

Not done:
- COMPASS and CytoBatchFlagR (optional);
- the instrument-only correction with cyCombine / CytoNorm 2.0, which is built only if part B returns "correct" for the instrument.
