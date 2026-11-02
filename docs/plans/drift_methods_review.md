[< Back to Documentation Index](../README.md)

# Drift analysis: critical review against the cytometry literature

Status: review written 2026-10-13. It covers the approach in `mtb_batch_drift_plan.md`, `drift_three_questions_plan.md` and `drift_technical_factors_plan.md`, the code under `scripts/drift/` and `src/harmonization/drift/`, and the results in `results_mtb_drift_2026-10-12_operator/`. The checks in section 2 were run on the development box on 2026-10-13 against that results folder, with the `mtb_drift` micromamba environment. Section 6 records the implementation of the recommendations and the rerun the same day.

## 1. Verdict in brief

The overall strategy fits the literature for a cohort with no anchors:

- keep the technical terms in the final model instead of correcting first (Nygaard et al. 2016);
- protect the biology;
- avoid anchor-based or balance-assuming correction tools (CytoNorm, CytoNorm 2.0 without controls, cyCombine, BatchAdjust, CytofRUV, iMUBAC).

Five problems weaken the conclusion "no correction is needed". All five are fixable without new data.

1. The endpoint and variance-partition models have no visit (patient x week) term. Visit-level variation therefore lands in the day term and in the "replicate noise" that rule 1 uses as its yardstick.
2. The verdict for frequencies rests on gates that the analyst moved by hand per acquisition. This applies to the cytokine gates as well as the CD4/CD8 gates. The frequency-level test therefore mostly measures how well the gating compensated, and that compensation cannot be audited.
3. Rule 1 treats "not detected" as "absent". With 13 of 23 days carrying one patient, most day effects are inconclusive rather than negligible.
4. Rule 2 counts flipped calls without a noise baseline. Counting noise alone flips about 130 of 2176 calls per resample, far more than any drift effect flips.
5. The positivity rule and the low parent event counts are a larger source of call instability than any technical factor found so far.

## 2. Checks run for this review

### 2.1 A visit term absorbs much of the "day" variance

The layer B lme4 model from `design_formula(cfg, "endpoint", "lme4")` (time spline, per-patient slope, instrument, operator, prior vaccination, pool, tube) was refitted with an added `(1 | patient:week)`. There are 56 visits with 3 to 5 stimulated tubes each, so this term is estimable.

| Response | Day variance, current | Day variance, with visit | Visit variance | Tube variance, current | Tube variance, with visit | AIC change |
|---|---|---|---|---|---|---|
| CD4 IFNg | 1.10 | 1.04 | 0.53 | 0.53 | 0.15 | -64 |
| CD4 IL2 | 0.73 | 0.63 | 0.22 | 0.22 | 0.11 | -38 |
| CD4 TNFa | 0.63 | 0.43 | 0.49 | 0.41 | 0.08 | -107 |
| CD4 CD107a | 0.68 | 0.34 | 0.40 | 0.36 | 0.11 | -90 |
| CD8 IFNg | 1.28 | 0.86 | 0.48 | 0.27 | 0.07 | -86 |
| CD8 IL2 | 0.47 | 0.26 | 0.49 | 0.51 | 0.11 | -63 |
| CD8 TNFa | 0.31 | 0.00 | 0.34 | 0.31 | 0.10 | -93 |
| CD8 CD107a | 0.70 | 0.44 | 0.41 | 0.28 | 0.06 | -197 |

All variances are on the logit scale.

The visit term is strongly preferred for every response. It does three things:

- It reduces the day variance by up to 100 percent (median about 35 percent).
- It cuts the tube-level residual by a factor of 2 to 5.
- It shows that two numbers the decision step relies on are misattributed:
  - `09_decision.R` uses the tube residual as the rule 1 threshold for counts;
  - the variance partition (`partition_terms()`, random intercepts only) books visit-level biology as acquisition day.

The "batch" share of 43 percent for logit frequencies in `partition_grid_summary.csv` is therefore an overestimate.

Why it happens: 13 of the 23 acquisition days hold a single patient. On those days, "day" and "visit" are the same set of tubes. Without a visit term, anything that is specific to that sample lands in the day effect. That includes deviation from the smooth time curve, thaw recovery and handling. The separability simulation (`simulate_endpoint()`) has no visit-level term and uses a tube noise SD of 0.2 against a real value of about 0.6. That is why it reported 0 percent false day detection.

### 2.2 The cytokine gates move with day and instrument

The cytokine gate bounds of all 328 tubes were read from the 20 workspaces. All gates are rectangles, and the bounds were converted to arcsinh with cofactor 150.

| Gate | SD of gate position across days | SD within a day | Gate shift on Symphony (p) | Neg-tube intensity shift on Symphony (layer BG) |
|---|---|---|---|---|
| CD4 TNFa | 0.43 | 0.07 | +0.69 (1e-35) | TNFa +0.80 |
| CD8 TNFa | 0.48 | 0.07 | +0.63 (2e-24) | |
| CD4 IL2 | 0.69 | 0.14 | -0.78 (7e-19) | IL2 -0.41 |
| CD8 IL2 | 0.54 | 0.10 | -0.73 (4e-37) | |
| CD4 CD107a | 0.41 | 0.16 | -0.65 (2e-36) | CD107a -0.47 |
| CD8 CD107a | 0.55 | 0.21 | -0.67 (5e-18) | |
| CD4 IFNg | 0.21 | 0.06 | -0.16 (2e-9) | IFNg -0.28 (not detected) |
| CD8 IFNg | 0.52 | 0.16 | +0.05 (0.47) | |

The gate shifts are adjusted for patient. Within a visit, the neg tube and every stimulated tube share one gate, which is good practice and makes background subtraction coherent.

Between acquisitions, the gates move 3 to 7 times more than within a day. On the Symphony they follow the measured intensity shift in direction and roughly in size. The existing gate-drift diagnostic only looks at the CD4/CD8 parent gates. Its two `gate_drift.csv` rows per marker are duplicates, because both parents are read from the same CD4 vs CD8 gate.

Consequences:

- "No instrument effect on frequencies" means "the analyst compensated the instrument effect at gating". That is a valid outcome, but it is a manual correction.
- If gates were set so that the neg tube shows a low background, the layer BG background frequencies partly measure the gate as well. The 2023-07-05 CD8 IFNg background event should be reviewed on the plot before it is attributed to the instrument.

Maecker et al. 2005, Jaimes et al. 2011 and Welters et al. 2012 all name gating as the largest analytical source of ICS variability. Finak et al. 2016 showed that central automated gating removes much of it.

### 2.3 Counting noise alone flips many calls

Each stimulated tube and its matched neg were resampled 500 times from binomials at their observed proportions. This is pure counting noise, with no overdispersion. The lab rule (at least 0.05 percent stim minus neg, at least 3x neg) was then reapplied.

- 2176 stimulated calls, 426 positive.
- About 133 calls flip per resample.
- 418 calls have a flip probability above 10 percent, and 236 above 25 percent.
- Parent counts are low:

  | Parent | Median events | 10th percentile |
  |---|---|---|
  | CD4 | 6004 | 762 |
  | CD8 | 2394 | 264 |

  At 264 parent events, one cytokine-positive cell is 0.38 percent, which is eight times the 0.05 percent threshold.
- 22 percent of matched negs hold 5 events or fewer, and 109 hold none. The 3x fold rule on such a neg is close to arbitrary (Moodie et al. 2010; Roederer 2008).

Against this baseline, "removing the day effect flips 2 calls" (CD4 IFNg, 2024-03-07) and "15 calls" (CD8 IFNg background, 2023-07-05) cannot be read without context. A flip count is only informative relative to the flips caused by noise of the same size.

### 2.4 A drift-free estimate of the week effect exists and disagrees

On 15 patient-days, two or three visits of one patient were acquired together (1304 stimulated tubes). Week contrasts within those patient-days are immune to day, instrument and operator drift by construction.

The test fitted two models to the same tubes:
- patient-day as a fixed effect (the within-day estimator);
- patient and day as random intercepts (the current approach).

The within-day slope was larger for all 8 responses. For example, CD4 IFNg was 5.8 against 1.2 logit per year, and CD4 IL2 was 3.3 against 1.7. These patient-days mostly cover weeks 0 to 17, so the slopes describe the early rise and are not comparable with the whole-study slope. The discrepancy is still a sign that the random-day model partly trades week signal against day effects.

The within-patient-day estimator is the standard sensitivity analysis for this design, because each patient-day acts as its own batch. It should be reported next to the model-based week effect.

## 3. Further issues, by part of the pipeline

### 3.1 Decision rule

| Issue | Why it matters | Improvement |
|---|---|---|
| Rule 1 has two outcomes: detected or not | Days with one or two patients have wide posteriors, so "not detected" mostly means "not identifiable". The current tables do not distinguish the two | Use three states with a region of practical equivalence: "detected" if P(\|effect\| > margin) ≥ 0.95, "negligible" if P(\|effect\| < margin) ≥ 0.95, otherwise "inconclusive" (Kruschke 2018; Lakens 2017 for the frequentist analogue). Report the inconclusive count prominently |
| The margin is the tube residual SD | The tube residual is inflated by visit variance (2.1), and it is a noise scale, not a relevance scale | Set the margin on an external relevance scale, for example the proficiency-panel CVs (Maecker 2005; Jaimes 2011: about 35 percent above 0.2 percent responses, worse below), converted to logit units |
| Rule 2 uses the posterior mean and an absolute flip count | Ignores effect uncertainty and the noise baseline (2.3) | Propagate posterior draws to get an expected flip count with an interval, and compare it with the counting-noise flips of 2.3. Call the effect relevant only if it adds flips beyond that baseline |
| 423 layer A day x statistic effects are judged at 0.95 | Some detections are expected by chance | Partial pooling already reduces this (Gelman, Hill and Yajima 2012). State the expected false-detection count alongside the result |

### 3.2 Layer A (technical markers in negs)

- Summaries are taken inside "Live cells", which is itself a hand-set gate on scatter and live/dead. The live/dead and scatter statistics are therefore conditioned on the gating. Summaries on ungated events, or on a fixed automated gate (openCyto, Finak et al. 2014 PLoS Comput Biol), would separate instrument from gater.
- Only the 56 neg tubes enter the layer A models. Step 4 with `--all_tubes` already summarises all 328 tubes, which is about six times more data on lineage and scatter markers. Six hours of peptide stimulation barely moves these markers; CD4 down-modulation on activated cells can be checked through the `stimulated` term. Fitting layer A on all tubes with a stimulation indicator would make day effects far better identified.
- One arcsinh cofactor (150) is used for both instruments. The Canto II and Symphony A1 have different detector noise, so equal arcsinh shifts are not equal signal shifts. Le Lann et al. 2020 and Solly et al. 2013 harmonise in bead (MESF or target value) units. Without beads, a per-channel landmark scale is the next best option, using the CD4- and CD8-negative peak positions and widths in negs (Hahne et al. 2010, `gaussNorm`).
- EMD uses two different references, depending on whether the patient has other Canto negs (`ref_kind`). The meaning of the EMD response therefore differs between tubes. An EMD-plus-MDS view of all neg tubes against each other, coloured by day, instrument and patient (Orlova et al. 2016; Hauchamps et al. 2025, CytoMDS), would show the same thing without a reference choice.

### 3.3 Model structure

- Add `(1 | patient:week)` to the endpoint, background and partition models (2.1).
- Day effects are exchangeable random intercepts. Drift in cytometry is often gradual or stepwise: laser ageing, tandem-dye degradation (CD4 is APC-Cy7, CD8 is PE-Cy7), reagent lot changes, service visits. A smooth term over calendar date (`s(acq_date)` in brms, or a Gaussian process), plus an indicator per lot or service date if the lab has the records, would borrow strength between neighbouring days. It would also separate trend from day-specific noise.
- The partition model has no per-patient time slope, so patient-specific kinetics also go to day. Use the same random-effect structure as the endpoint model.
- Operator, calendar period and instrument are close to aliased (`drift_technical_factors_plan.md` section 3). Report their variance as one technical block, and run the canonical-correlation check of variancePartition (Hoffman and Schadt 2016) before splitting it.
- The 2026-10-12 brms layer B fit for CD4 IFNg had 57 divergent transitions. Its day effects should not be quoted until it is refitted.

### 3.4 Positivity calling and background

- Lab rule, Fisher and MIMOSA agree poorly. Cohen's kappa between the lab rule and MIMOSA is 0.22 to 0.64, and for CD4 TNFa the counts are 72, 147 and 147 positive. The choice of rule changes more calls than any technical factor. The rule should be pre-specified (Moodie et al. 2010) and its sensitivity reported.
- Horton et al. 2007 pair the empirical rule with a count-based test. MIMOSA (Finak et al. 2014) models stim and neg counts jointly with overdispersion and shares strength across subjects. COMPASS (Lin et al. 2015) does the same for Boolean combinations. Both handle zero-event negs, which the fold rule cannot.
- MIMOSA and COMPASS have no batch term. A day effect that hits stim and neg alike cancels, but a stimulation-specific one does not. A beta-binomial model with a stimulated indicator and the technical terms is the natural extension, and it is what the `freq_raw` representation already fits.
- Set a minimum parent event count, below which a tube is reported as not evaluable. Roederer 2008 gives the Poisson argument. Jaimes et al. 2011 name event number as a main source of inter-laboratory variation.

### 3.5 Correction methods, if ever needed

The literature search confirms the current position:

- CytoNorm, BatchAdjust, CytofRUV and iMUBAC need anchors or control donors in every batch.
- CytoNorm 2.0 without controls and cyCombine assume comparable composition in every batch. On days that hold a single patient-week, that composition is the biology.

The two options that fit this cohort:

- For the instrument: landmark (peak) alignment of stable populations in negs, learned from the four bridge patients (Hahne et al. 2010). Finak et al. 2014 (Cytometry A) applied this inside gated subsets of HVTN ICS data and reduced false-positive vaccine calls. It is the closest published analogue to this cohort.
- For days: keep the terms in the model and do not pre-correct (Nygaard et al. 2016).

### 3.6 BayesFlow

Both variants failed their gates, and the plan already says that simulator quality is the limit. With the visit term missing from every simulator, further work there will not change the decision. Freeze it until anchor or bead data exist.

## 4. Recommended changes, in order

| # | Change | Effort | Verify |
|---|---|---|---|
| 1 | Add `(1 \| patient:week)` to `design_formula()` (background and endpoint classes) and `partition_terms()`. Add a visit SD to `simulate_endpoint()` and use a realistic tube SD | half a day | Synthetic cohort with visit noise and no day drift gives no detected day effect; real-data AIC table of 2.1 reproduced |
| 2 | Rule 1 with three states and an external margin; rule 2 as expected flips with interval, minus the counting-noise baseline | one day | Unit tests on `drift_decision()`; report shows inconclusive counts |
| 3 | Read every cytokine gate into `gate_bounds.csv` (step 4 already loops over the hierarchy); report gate shift against intensity shift per channel; fix the duplicate CD4/CD8 rows | half a day | Table of 2.2 reproduced by the pipeline |
| 4 | Within-patient-day week estimator (patient-day fixed effect) as a standard sensitivity output of `06_biology.R` | two hours | Table of 2.4 |
| 5 | Layer A on all tubes with a `stimulated` term; summaries on ungated or automatically gated events | one day | Day posterior SDs shrink; instrument effects agree with the neg-only fit |
| 6 | Re-gate cytokine gates with one fixed or landmark-normalised template per channel (openCyto), as the sensitivity analysis the plan already names (rank 12) | two days | Frequencies, calls and drift verdicts with manual against template gates |
| 7 | Smooth calendar-time term next to the day intercept; lot and service dates if available | one day | Compare LOO (brms `loo`) with and without |
| 8 | Minimum parent count; MIMOSA or the joint beta-binomial as the primary call, the lab rule as secondary | half a day plus a decision by the lab | Kappa table and counting-noise flip rate per rule |
| 9 | Refit CD4 IFNg (57 divergences) | part of 1 | R-hat ≤ 1.01, no divergences |

## 5. For future acquisitions (lab side)

The plan's lab actions are confirmed by every design guideline found (Cossarizza et al. 2019 and 2021; Kleinsteuber et al. 2016; Perfetto et al. 2012; Leek et al. 2010):

- one reference PBMC aliquot from a single large draw on every acquisition day, ideally two;
- CS&T and 8-peak bead records kept and exported;
- a written gating template with documented rules for moving a gate;
- a minimum number of acquired events per tube, set so that the parent populations clear the threshold of 3.4;
- all visits of a patient thawed and acquired together where the sample schedule allows it. This design made 2.4 possible.

## 6. Implementation and rerun (2026-10-13)

Code changes, all behind config fields so older configs behave as before:

| Recommendation | Implementation |
|---|---|
| 1 Visit term | `design.visit: [patient, week]` builds a `visit` factor in `design_frame()`. `(1 \| visit)` enters the endpoint class, the technical class on all tubes and the variance partition. The background classes are left without it: with one neg per visit it would be the tube residual under another name. In the partition the visit share is its own group (`sample`) and its own colour, neither biology nor batch. Both Stan programs take an optional visit intercept. `simulate_endpoint()` has a visit SD (0.6) and a realistic tube SD (0.3) |
| 2 Decision rule | `drift_decision()` has three states against `relevance_margin()`: log(1.35) on the logit scale for counts (`relevance.count_cv`), the same-sample SD otherwise. `flip_summary()` gives rule 2 over 200 posterior draws, and `counting_noise_flips()` gives the binomial baseline on the same tubes. "Detected, correct" now needs more flips than the 95th percentile of counting noise; flips below that are "detected, within counting noise" |
| 3 Gate audit | Step 4 records the cytokine gates under each parent with a `gate_marker` column. Step 9 measures each gate on its own dimension, which also removes the duplicated CD4/CD8 rows |
| 4 Within-day slope | `06_biology.R` writes `week_slope_within_day.csv` |
| 5 Layer A on all tubes | `technical_tubes: all`; a `stimulated` fixed term and the visit intercept. Ungated summaries: step 4 with `neg_population: root` into a subfolder, compared in step 15 |
| 6 Fixed gates | Uses the fixed-threshold frequencies of step 11 (99.5th percentile of the reference-instrument negs, one threshold per channel) in step 15. An openCyto re-gating was not built |
| 7 Calendar time | Step 15 refits with a natural spline of the acquisition date (4 df) next to the day intercept. lme4 AIC, not brms LOO |
| 8 Calling | Step 15 writes the counting-noise flip probability per call and the calls kept at minimum parent counts. The primary rule is unchanged (a decision for the lab) |
| 9 CD4 IFNg | Refitted with the visit term |

New step `15_sensitivity.R`; report sections for the three-state verdicts, the gate audit and the sensitivity analyses; 244 base-R tests pass (22 new).

Rerun: `results_mtb_drift_2026-10-13_review/` on the development box. Steps 1 to 3, 5, 11 and 13 and the layer BG fits were copied from `results_mtb_drift_2026-10-12_operator/`; steps 4, 1b, 6a and 6b (layers A and B), 6c, 6d, 12a, 12b, 15, 9, 14 and 10 were rerun. The command lines are in `run_review.sh` and `run_finish.sh` in that folder.

### 6.1 Convergence

- Layer B, brms with the visit term: R-hat ≤ 1.005 and no divergent transitions for all 8 endpoints, CD4 IFNg included (57 divergences before).
- Layer A on all tubes: R-hat ≤ 1.01. Divergences remain in three scatter and viability spreads (FSC 34, LiveDead 22, SSC 4), which are not correction candidates.
- brms and Stan agree on 885 of 885 day, 39 of 39 instrument and 78 of 78 operator effects.
- Variance partition: max R-hat 1.028, no divergences (as before).

### 6.2 Answer 1: correction

| Layer | Detected | Negligible | Inconclusive | Correct |
|---|---|---|---|---|
| A technical (19 statistics) | 70 | 0 | 416 | not applicable |
| BG background | 30 | 0 | 278 | 1 |
| B endpoints | 7 | 0 | 201 | 0 |

- No contrast is negligible. Every technical effect that was "not detected" before is inconclusive: the data cannot show that it is smaller than the margin. "No correction needed" is therefore "no correction shown to be needed".
- Layer B detects 5 day effects, the Symphony on CD8 TNFa (+1.09 logit, 95 percent interval 0.25 to 1.94) and operator Anna on CD8 TNFa. None changes more calls than counting noise. The day that flipped one call on 2026-10-12 (CD4 TNFa on 2024-07-12) is no longer detected.
- One correction candidate remains: the CD8 IFNg background on 2023-07-05 (+2.6 logit). Removing it changes 12.9 calls (95 percent interval 6 to 15), against a counting-noise 95th percentile of 4 on the same tubes. Flag or correct the CD8 IFNg calls of that day.

### 6.3 The hand gating carries the instrument correction

| Readout | Manual gates | One fixed threshold per channel |
|---|---|---|
| Symphony effect, CD8 CD107a | +0.47 (SE 0.51) | -5.41 (1.27) |
| CD8 IL2 | -0.12 (0.45) | -2.71 (0.71) |
| CD4 IL2 | -0.39 (0.53) | -1.96 (0.60) |
| CD4 CD107a | +0.08 (0.50) | -1.80 (0.61) |
| CD8 TNFa | +1.41 (0.37) | +1.99 (0.45) |
| Day variance, CD8 IFNg | 0.86 | 2.11 |
| Day variance, CD4 TNFa | 0.43 | 0.99 |

lme4, logit scale. In the lme4 partition, the instrument share of the logit frequencies is about 1 percent with the manual gates and 15 percent with fixed thresholds. The direction of the fixed-threshold effects matches the measured gate shifts (IL2 and CD107a gates 0.7 lower on the Symphony). The fixed threshold is set on the Canto negs, so it overstates the Symphony effect where the Symphony scale differs. The analyst's gates removed most of the instrument shift for six of eight readouts. The exception is CD8 TNFa, which keeps a Symphony effect in both versions.

The technical markers tell the same story one gate earlier. Without the live-cell gate, the Symphony shift of the CD8 median grows from -0.76 to -1.16 arcsinh, and that of the CD8 EMD from 0.45 to 0.65.

### 6.4 Visit, day and design

- With the visit term the brms day SD of the endpoints falls by up to 1.2 (CD8 TNFa 1.56 to 0.33), and the visit SD is 0.56 to 0.78. In the partition, the visit takes 10 to 32 percent of the frequency variance. The batch share of logit frequencies falls from 43 to 30 percent, and that of subtracted percentages from 18 to 15 percent.
- The design check with realistic visit and tube noise changes its verdict from "separable" to "separable only through the unstimulated layer". The week slope is still recovered (1 percent leakage), but the day SD is underestimated by 32 percent. A pure week effect produces a day SD above 0.1 in 33 percent of simulations. Day effects estimated from the endpoints alone are weakly identified.
- A smooth of calendar time does not explain the endpoint day effects (AIC worse for all 8). Side scatter drifts gradually (AIC -12 and -13 for the SSC median and EMD); no other technical marker improves by more than 4 AIC units.

### 6.5 Biology

- The whole-study week slope stays within 1.2 SE when the technical terms are removed, for every endpoint.
- On the 15 subject-days with several visits, the drift-free early-visit slope of the CD4 responses is 3.6 to 5.9 logit per year. On the same tubes, the random-day model gives -0.19 to 0.76 logit per year for IFNg, IL2 and TNFa, 3 to 4 SE lower. In this window the random-day model gives part of the rise to the day effects. The early response should be reported from the within-day contrasts.

### 6.6 Calls

- Counting noise alone changes 133 of 2176 lab-rule calls on average. 419 calls have a flip probability above 10 percent.
- Requiring 5000 parent events keeps 41 percent of the calls and cuts the expected flips to 28. Requiring 1000 keeps 76 percent with 78 expected flips.

### 6.7 Open

- The primary calling rule and a minimum parent count (lab decision).
- The CD8 IFNg calls of 2023-07-05, and the CD8 TNFa instrument effect (keep `instrument` in any CD8 TNFa analysis).
- Re-gating with a landmark-normalised template (Finak et al. 2014, Cytometry A), which would separate the instrument shift from the analyst's judgement better than one fixed threshold.

## 7. References

DOIs were resolved against Crossref on 2026-10-13.

- Finak G et al. 2014. High-throughput flow cytometry data normalization for clinical trials. Cytometry A 85:277. 10.1002/cyto.a.22433
- Finak G et al. 2014. Mixture models for single-cell assays with applications to vaccine studies (MIMOSA). Biostatistics 15:87. 10.1093/biostatistics/kxt024
- Finak G et al. 2014. OpenCyto. PLoS Comput Biol 10:e1003806. 10.1371/journal.pcbi.1003806
- Finak G et al. 2016. Standardizing flow cytometry immunophenotyping analysis from the Human ImmunoPhenotyping Consortium. Sci Rep 6:20686. 10.1038/srep20686
- Gelman A, Hill J, Yajima M 2012. Why we (usually) don't have to worry about multiple comparisons. J Res Educ Eff 5:189. 10.1080/19345747.2011.618213
- Hahne F et al. 2010. Per-channel basis normalization methods for flow cytometry data. Cytometry A 77:121. 10.1002/cyto.a.20823
- Hauchamps P et al. 2025. CytoMDS. Cytometry A 107:177. 10.1002/cyto.a.24921
- Hoffman GE, Schadt EE 2016. variancePartition. BMC Bioinformatics 17:483. 10.1186/s12859-016-1323-z
- Horton H et al. 2007. Optimization and validation of an 8-color ICS assay. J Immunol Methods 323:39. 10.1016/j.jim.2007.03.002
- Jaimes MC et al. 2011. Quality assurance of intracellular cytokine staining assays. J Immunol Methods 363:143. 10.1016/j.jim.2010.08.004
- Kleinsteuber K et al. 2016. Standardization and quality control for high-dimensional mass cytometry studies of human samples. Cytometry A 89:903. 10.1002/cyto.a.22935
- Kruschke JK 2018. Rejecting or accepting parameter values in Bayesian estimation. Adv Methods Pract Psychol Sci 1:270. 10.1177/2515245918771304
- Lakens D 2017. Equivalence tests. Soc Psychol Personal Sci 8:355. 10.1177/1948550617697177
- Le Lann L et al. 2020. Standardization procedure for flow cytometry data harmonization in prospective multicenter studies. Sci Rep 10:11567. 10.1038/s41598-020-68468-3
- Leek JT et al. 2010. Tackling the widespread and critical impact of batch effects. Nat Rev Genet 11:733. 10.1038/nrg2825
- Lin L et al. 2015. COMPASS. Nat Biotechnol 33:610. 10.1038/nbt.3187
- Maecker HT et al. 2005. Standardization of cytokine flow cytometry assays. BMC Immunol 6:13. 10.1186/1471-2172-6-13
- Moodie Z et al. 2010. Response definition criteria for ELISPOT assays revisited. Cancer Immunol Immunother 59:1489. 10.1007/s00262-010-0875-4
- Nygaard V, Rødland EA, Hovig E 2016. Methods that remove batch effects while retaining group differences may lead to exaggerated confidence. Biostatistics 17:29. 10.1093/biostatistics/kxv027
- Orlova DY et al. 2016. Earth Mover's Distance (EMD). PLoS One 11:e0151859. 10.1371/journal.pone.0151859
- Pedersen CB et al. 2022. cyCombine. Nat Commun 13:1698. 10.1038/s41467-022-29383-5
- Perfetto SP et al. 2012. Quality assurance for polychromatic flow cytometry using a suite of calibration beads. Nat Protoc 7:2067. 10.1038/nprot.2012.126
- Quintelier KLA et al. 2025. CytoNorm 2.0. Cytometry A 107:69. 10.1002/cyto.a.24910
- Roederer M 2008. How many events is enough? Are you positive? Cytometry A 73:384. 10.1002/cyto.a.20549
- Solly F et al. 2013. Comparable flow cytometry data can be obtained with two types of instruments, Canto II and Navios. Cytometry A 83:1066. 10.1002/cyto.a.22404
- Van Gassen S et al. 2020. CytoNorm. Cytometry A 97:268. 10.1002/cyto.a.23904
- Welters MJP et al. 2012. Harmonization of the intracellular cytokine staining assay. Cancer Immunol Immunother 61:967. 10.1007/s00262-012-1282-9
- Cossarizza A et al. 2019. Guidelines for the use of flow cytometry and cell sorting in immunological studies (second edition). Eur J Immunol 49:1457. 10.1002/eji.201970107
- Cossarizza A et al. 2021. Third edition. Eur J Immunol 51:2708. 10.1002/eji.202170126

Not found, despite being cited in earlier plans or searched for: Halvorson et al. 2025 (already flagged in `drift_three_questions_plan.md`). No published head-to-head comparison of the FACSCanto II and FACSymphony A1 for ICS was found.
