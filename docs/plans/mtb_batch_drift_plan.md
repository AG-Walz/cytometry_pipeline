[< Back to Documentation Index](../README.md)

# MTB ICS cohort: detecting acquisition drift and instrument shift

Status: agreed on 2026-09-21 after a structured interview; scripts implemented the same day (section 9). Everything runs CPU-only in a screen on highmem1.

## 1. Question

Over 2.5 years of acquisition, do the MTB intracellular cytokine staining (ICS) tubes carry non-biological variation attributable to acquisition day or to the cytometer, is it large enough to change a positivity call, and if so which correction is justified?

The filename field called `batch` is prior BCG vaccination status (1 = first vaccination in the study, 2 = previously vaccinated). It is a biological covariate. It must be protected in every model and must never be passed to a correction method as a batch label.

## 2. What the data look like

Source: `fcs_overview.xlsx` (from `scripts/fcs_overview.py`) joined with the FCS TEXT segments (`$DATE`, `$CYT`, `CREATOR`, `$PAR`). The folder is `/mnt/volume/data/MTB_FACS/` on the development box and `/mnt/volume/patrick/UMAP/` on highmem1, so every script takes the data root as a parameter and defaults to nothing.

| Item | Value |
|---|---|
| Tubes | 328 FCS files, 20 patients, 20 FlowJo workspaces |
| Neg tubes | 56, exactly one per patient-week (37 patient-days, 2 to 21 tubes each) |
| Peptide cocktail tubes | 272, mostly 3 to 5 per patient-week |
| Acquisition days | 23, from 2022-11-24 to 2025-06-04 |
| Instruments | FACSCanto II, Diva 6.1.2: 274 tubes on 19 days. FACSymphony A1, Diva 9.0.2: 54 tubes on 4 days in 2024 |
| Parameters | Canto 12 or 14, Symphony 12 or 13. Panel versions differ across runs and must be reconciled to the common marker set first |
| Markers | FSC, SSC, live/dead (Aqua), CD4 (APC-Cy7), CD8 (PE-Cy7), IFNg (PE), IL-2 (APC), TNFa (Pacific Blue), CD107a (FITC) |

Identifiability comes from repeated patients:

| Role | Patients |
|---|---|
| Acquired on more than one day (day drift is estimable within patient) | 12, 16, 17, 23, 34, 39, 48, 61, 69, 74 (10 patients, 2 to 5 days each) |
| Acquired on both instruments (instrument shift is estimable within patient) | 16, 17, 61, 74 |
| Symphony only (instrument fully confounded with patient) | 79 |
| Single day only (contribute to between-patient variance only) | 11, 57, 59, 65, 73, 79, 82, 85, 87, 88 |

Two caveats to state in every result:

- Day and instrument are partially confounded with calendar time and with patient. The four Symphony days fall in Feb to Aug 2024. A shift on those days is "instrument" only if the Canto days around them do not show the same shift. The Canto ran on 2024-07-12 between Symphony days, which is the one direct contemporaneous comparison.
- Week-in-study is confounded with day within a patient by design (each visit is a new day). Week is the biological effect of interest, so day drift is only estimable from the assumption that unstimulated (neg) marker distributions for lineage and viability markers do not change with vaccination. That assumption is the load-bearing one and is tested in step 4.2.

## 3. Decision rule

Correction is triggered per marker and per instrument or day contrast when both hold:

1. The technical effect on background-subtracted frequencies is larger than the same-day tube-to-tube replicate noise of that patient (posterior probability above 0.95 that the effect exceeds the replicate standard deviation).
2. Applying or removing the effect flips positivity calls under the fixed rule: background-subtracted frequency at or above 0.05 percent of the parent population and at least 3 times the neg tube frequency.

Both thresholds live in one config block so the lab can change them. Effects that pass rule 1 but not rule 2 are reported as "detected, not clinically relevant", not corrected.

## 4. Two layers of evidence

### 4.1 Layer A: neg tubes, marker distributions (pure technical signal)

Neg tubes have no antigen, so CD4, CD8, live/dead, FSC and SSC distributions should depend on patient and on nothing else. Cytokine channels in negs carry background only and are analysed as background level, not as biology.

Per neg tube, after compensation from the file's SPILL matrix and arcsinh transform (cofactor 150 for the Canto and Symphony PMT scales, checked on the data), extract for each marker on live singlets:

- median, MAD, 10th and 90th percentile
- for CD4 and CD8: position of the bimodal valley (density minimum between the two peaks) and fraction positive using that valley
- Earth Mover's Distance (EMD) to a pooled reference built from all Canto negs of the same patient where available, else to the pooled Canto neg distribution

Also record from the workspace, via CytoML: the CD4 and CD8 gate boundary coordinates per tube. Analysts adjust gates per run, so moving gates can mask an instrument shift at the frequency level. Gate coordinates regressed on day and instrument tell you whether that happened.

### 4.2 Layer B: gated frequencies (endpoint-level signal)

From the FlowJo gates read through CytoML/flowWorkspace: for every cocktail tube, the frequency of IFNg+, IL-2+, TNFa+, CD107a+ and their Boolean combinations within CD4+ and within CD8+, minus the matched neg tube of the same patient-week. Keep raw and subtracted values, and the parent counts, because the model works on counts.

## 5. Yardstick: Bayesian hierarchical models

This is the reference against which everything else, including BayesFlow, is judged. Two fits of the same model to check that the answer does not depend on the software.

Model for layer A, one per marker summary $y$:

```
y_tube ~ Normal(mu, sigma)
mu = alpha + a_patient + b_instrument + c_day
a_patient ~ Normal(0, tau_p)      # 20 levels
c_day     ~ Normal(0, tau_d)      # 23 levels, partially pooled
b_instrument: fixed effect, weakly informative prior Normal(0, 1) on the arcsinh scale
```

Model for layer B, per cytokine and parent population, on counts:

```
k_tube ~ BetaBinomial(n_parent, p, phi)
logit(p) = alpha + a_patient + w_week[patient] + v_vaccstatus + s_cocktail
           + b_instrument + c_day
```

Prior vaccination status enters as `v_vaccstatus` so nothing removes it. Week enters as a patient-specific slope or spline. `c_day` is the quantity of interest for drift and `b_instrument` for the shift. Report posterior means, 95 percent credible intervals, and `tau_d` relative to `tau_p` as the drift-to-biology ratio.

Fit once with `brms` (which writes Stan) and once with hand-written Stan through `cmdstanr`, both CPU, 4 chains. Agreement between them is a software check, not extra evidence.

Frequentist cross-check for reviewers who prefer it: `lme4::glmer` with the same structure, and the intraclass correlation of day within patient. Same conclusions expected, no new information.

## 6. BayesFlow: contained experiment in two variants

BayesFlow is amortised simulation-based inference. It does not detect batch effects itself. It learns a posterior estimator from simulated data, so its output is only meaningful if the simulator captures the real generative process. Both variants below are trained on CPU and produce posteriors over the same `b_instrument` and `c_day` parameters as the yardstick, so the comparison is direct.

### 6.1 Variant 1: tube-summary simulator

Simulator draws `a_patient`, `c_day`, `b_instrument`, `tau_p`, `tau_d`, `sigma` from the yardstick priors and generates the layer A summary table for the real design (same patients, days and instruments as the real 56 negs). Summary network: a set transformer over tubes, with day and instrument indicators as inputs. Inference network: coupling flow. Around 50k simulations, a few CPU hours.

Purpose: if this variant does not recover the brms posteriors on the real data and on held-out simulations, the pipeline is wrong and variant 2 is not worth running.

### 6.2 Variant 2: single-cell distribution simulator

Per tube, simulate 2000 cells from a mixture model per marker: two log-normal components for bimodal lineage markers (CD4, CD8), one for live/dead and scatter, with a low-frequency positive component for cytokine channels. Parameters:

- patient-level component means and weights (biology)
- day-level additive shift and multiplicative gain per fluorescence channel (drift, e.g. laser power and PMT voltage changes)
- instrument-level shift, gain and spread change per channel (Canto vs Symphony have different detectors and Diva versions)
- an event-rate and debris fraction term so the summary network sees realistic QC variation

Summary network: a DeepSet over cells producing a tube embedding, then a set transformer over tubes. Inference on `c_day`, `b_instrument` and the gain parameters. Training data around 100k simulated experiments, subsampled cells, CPU training overnight on highmem. Validation by simulation-based calibration and posterior predictive checks on the real negs (do simulated negs with the inferred parameters reproduce the real EMDs and valley positions).

Purpose: this variant can represent distribution-shape changes (spread, gain) that tube medians miss, and once trained it scores any new FCS file in seconds, which is a real QC asset for future runs. It drives the correction decision only if it passes calibration and agrees with the yardstick on the summaries both can see.

### 6.3 How BayesFlow is judged

| Check | Pass criterion |
|---|---|
| Simulation-based calibration | Rank histogram not rejected (chi-square p above 0.01) for at least 90 percent of parameters |
| Coverage on held-out simulations | 95 percent intervals contain the truth for 90 to 99 percent of test sets, for at least 90 percent of parameters |
| Agreement with yardstick | Posterior mean inside the brms 95 percent interval for at least 95 percent of the day and instrument effects |
| Posterior predictive on real negs (variant 2) | At least 85 percent of real per-tube medians inside the 90 percent predictive band |
| Recovery R squared | Reported, not gated. With one neg tube per patient-day, most day effects are weakly identified by any method, brms included, so a wide but calibrated posterior is the right answer and scores a low R squared |
| Cost | Reported honestly: wall-clock, simulations, tuning rounds |

If a variant fails any gated check, the decision falls back to the yardstick and the report says so.

## 7. Alternatives considered, and why each is or is not the primary path

Ranked from most to least suited to this dataset. Several are cheap and run as diagnostics regardless.

| Rank | Method | Verdict | Reason |
|---|---|---|---|
| 1 | Bayesian hierarchical model on neg summaries and on gated counts (section 5) | Primary | Handles the unbalanced repeated-measures design, partial pooling of 23 days, protects vaccination status and week as covariates, gives effect sizes with uncertainty. Reviewer-standard. |
| 2 | EMD between neg tubes, per marker, regressed on day and instrument | Diagnostic, run always | Assumption-free view of whole-distribution shift. Cannot separate patient from day without the model, but is the input for it and the visual that convinces. |
| 3 | Gate-position drift from workspaces | Diagnostic, run always | Only way to know whether analysts hand-corrected instrument shifts. Zero extra data needed. |
| 4 | Mixed-model ICC and variance components (`lme4`) | Cross-check | Same model as rank 1 without priors. Runs in seconds, useful as a frequentist appendix. |
| 5 | BayesFlow, both variants (section 6) | Contained experiment | Potentially the richest description of shift, reusable for future QC, but simulator-dependent and expensive. Must earn its place against rank 1. |
| 6 | kBET, LISI or silhouette on the neg-tube cell embedding | Secondary diagnostic | Measures mixing of cells across days in a UMAP or PCA. Sensitive to cell composition differences between patients, so a low mixing score does not prove a technical effect here. Report, do not decide on. |
| 7 | Cluster-abundance shifts (FlowSOM on pooled negs, then per-cluster mixed model) | Optional | Finds shifts localised to a subpopulation that global summaries miss. Adds a clustering step with its own instability. Run only if ranks 1 to 3 find something and you need to localise it. |
| 8 | CytoNorm quantile alignment | Not for detection, candidate for correction | Requires an anchor per batch. There is none. Could be run with pooled negs per day as pseudo-anchors, but negs differ by patient, so it would partially normalise away patient biology. Only if rule 3 triggers and only on lineage and viability channels. |
| 9 | Harmony on cell embeddings | Not for detection, weak for correction here | Corrects an embedding, not marker values, so gated frequencies cannot be recomputed from it. Existing wrapper in `src/harmonization/batch_correction.R` stays for UMAP use only. |
| 10 | ComBat on tube-level frequencies | Candidate for correction, not detection | Location-scale adjustment of the 23 day means would shrink real week effects because week and day are confounded within patient. Only defensible on neg-derived technical covariates, and only if the model says the day term is large. |
| 11 | Time-in-file drift checks (flowAI, PeacoQC) | Pre-processing QC, run first | Catches within-tube flow instability, clogging and fluidics drift. Not a batch method, but a tube that fails it should be excluded before any of the above. |
| 12 | Re-gating with openCyto | Rejected as primary | Would change the numbers the lab has reported. Kept as a sensitivity analysis only if gate-position drift (rank 3) is large. |

## 8. Correction strategy, conditional on findings

| Finding | Action |
|---|---|
| Instrument shift on lineage or viability markers, frequency calls flip | Per-channel quantile mapping of Symphony tubes to the Canto scale, learned from the four bridge patients' negs, applied before re-reading gates. Patient 79 corrected with the pooled mapping and flagged. |
| Instrument shift, calls do not flip | No correction. Include `instrument` as a covariate in the downstream statistical model of the vaccine response. |
| Day drift with `tau_d` above replicate noise, calls flip | Per-day additive offset from the neg model on the affected channels, applied to all tubes of that day. Report before and after. |
| Day drift, calls do not flip | No correction. Include `day` as a random effect downstream. |
| Cytokine background level differs by day or instrument | Never corrected on the marker scale. Handled by the existing background subtraction plus a day-level background term in the response model. |
| Gate-position drift large | Rerun layers A and B with openCyto gates as a sensitivity analysis and report both. |

## 9. Run order on highmem (CPU, screen session)

Implemented under `scripts/drift/` (thin CLIs) and `src/harmonization/drift/` (modules). Base-R modules are covered by `test_data/test_drift.R`; `test_data/drift_synthetic.R` generates step 4/5 tables with known effects so steps 6 to 9 can be exercised without CytoML.

```bash
bash scripts/drift/00_setup_env.sh cytometry_env bayesflow_env      # once
screen -S mtb_drift
bash scripts/drift/run_all.sh /mnt/volume/patrick/UMAP results/mtb_drift config/drift_mtb.yaml cytometry_env bayesflow_env
```

`run_all.sh` takes the cohort config as its third argument and passes it to every step as `--config`; the config holds the file-name regex, the neg rule, the instrument labels, the marker patterns and the gate names, so another cohort is a new YAML (plus `SAMPLE_SHEET=<csv>` when file names carry no metadata). It expects `fcs/` and `wsp/` under the data root (override with `FCS_DIR` and `WSP_DIR`), writes `run_all.log` into the results folder, and honours `SKIP_BF=1`, `CORES`, `WORKERS`, `BF_EPOCHS`, `BF_SIMS`, `BF_CELLS`.

| Step | Script | Output | Needs | Approx. time |
|---|---|---|---|---|
| 0 | `00_setup_env.sh` | brms, cmdstanr, yaml into the conda env; CPU BayesFlow conda env (Python 3.12) | conda, internet | 1 h |
| 1 | `01_extract_metadata.R` | `tubes.csv`, `bridge_report.csv` | base R | seconds (verified here on the 328 files) |
| 2 | `02_panel_reconcile.R` | `panel_map.csv`, `panel_signatures.csv`, `panel_gaps.csv` | base R, yaml | seconds (verified here) |
| 3 | `03_tube_qc.R` | `qc.csv`, `excluded_tubes.csv` | flowCore, PeacoQC or flowAI | 30 min |
| 4 | `04_neg_summaries.R --all_tubes` | `neg_summaries.csv`, `neg_emd.csv`, `gate_bounds.csv`, `neg_cells/*.csv.gz` | CytoML, flowWorkspace | 30 min. Run once with `--list_populations` first and set `neg_population` in the config |
| 5 | `05_gated_frequencies.R` | `all_gate_counts.csv`, `gated_frequencies.csv` | CytoML, flowWorkspace | 30 min |
| 6a | `06_yardstick_brms.R` | `effects_brms.csv`, `fits_brms_summary.csv`, `fits_brms/` | brms | 1 to 2 h |
| 6b | `06_yardstick_stan.R` | `effects_stan.csv` | cmdstanr, `stan/*.stan` | 1 h |
| 6c | `06_lme4.R` | `lme4_variance_components.csv` | lme4 | minutes |
| 7 | `07_bayesflow_v1.py` | `effects_bayesflow_v1.csv`, `bf_v1_recovery.csv`, `bf_v1_sbc.csv`, `bf_v1_timing.json` | BayesFlow venv | hours (see timing note) |
| 8 | `08_bayesflow_v2.py` | `effects_bayesflow_v2.csv`, `bf_v2_recovery.csv`, `bf_v2_sbc.csv`, `bf_v2_ppc.csv`, `bf_v2_timing.json` | BayesFlow venv, `neg_cells/` | overnight |
| 9 | `09_decision.R` | `decisions.csv`, `gate_drift.csv`, `fitter_agreement.csv`, BayesFlow pass/fail printed | base R | minutes |
| 10 | `10_report.Rmd` | `report.html` | rmarkdown, ggplot2 | minutes |

Steps 1 to 6 and 9 to 10 give a complete, defensible answer on their own. Steps 7 and 8 run afterwards and step 9 is repeated so the report includes them.

Verification status on 2026-09-21: every step ran end to end on the real data on the development box, in a micromamba environment (`/mnt/volume/micromamba`, env `mtb_drift`) built from the same package list as `00_setup_env.sh`. Fixes made during that run: CytoML needs `do.call` to pass a `subset` value; FlowJo-derived cytoframes name compensated channels `Comp-<channel>`; rectangle gates have one finite side per dimension; the QC exclusion rule was made lenient because PeacoQC flags 30 to 55 percent of events on eleven gateable tubes. The highmem1 run is a reproduction of this, not a first attempt.

### Results of the development-box run (2026-09-21)

Full outputs: `/mnt/volume/data/MTB_FACS/results_mtb_drift_2026-09-21/` (`report.html`, `decisions.csv`, `fitter_agreement.csv`, `gate_drift.csv`, `lme4_variance_components.csv`, BayesFlow acceptance tables, logs).

Yardstick. brms and hand-written Stan agree on 623 of 630 day and instrument effects (Stan's posterior mean inside the brms 95 percent interval); the seven exceptions are layer B day effects on days with few tubes. All brms fits have R-hat at or below 1.015 after the layer B refit with adapt_delta 0.99.

Instrument (FACSymphony A1 relative to FACSCanto II):

| Where | Finding | Rule 1 | Rule 2 |
|---|---|---|---|
| Neg tubes, scatter | FSC-A and SSC-A medians differ by tens of thousands of raw units | passes | not applicable: a detector scale difference, reported and not corrected; scatter is excluded from correction candidates |
| Neg tubes, fluorescence | CD8 median -0.62 (95 percent interval -1.07 to -0.11), CD4 valley -0.46 (-0.73 to -0.17), CD8 and live/dead EMD +0.4, CD8 spread and live/dead median borderline; all in arcsinh units | passes for 6 of 13 statistics | not applicable at this layer |
| Gate positions (rank 3 diagnostic) | analysts placed the CD4 gate 0.26 and the CD8 gate 0.70 arcsinh units lower on Symphony days (both p < 0.001 after adjusting for patient), i.e. they compensated the shift by hand at gating time | | |
| Gated frequencies, all 8 responses | instrument effect not detected for any cytokine in CD4 or CD8; intervals span zero; removing the posterior mean effect flips 0 positivity calls | fails | 0 flips |

Decision: no instrument correction of the frequencies. `instrument` enters the downstream response model as a covariate. Patient 79 (Symphony only) is flagged in the bridge report.

Acquisition day:

| Layer | Detected above replicate noise | Positivity flips |
|---|---|---|
| A, neg tubes (423 day-by-statistic effects) | 60, mostly scatter and live/dead medians, and CD8 spread | not applicable |
| B, frequencies (180 day-by-response effects) | 4, none clinically relevant | at most 4 calls per day out of about 270; 12 in total |

Decision: day drift is real in the raw signal and almost entirely absorbed by gating and background subtraction. No per-day correction; `acq_day` enters the response model as a random effect. The day variance on the logit scale is of the same order as the between-patient variance (brms tau_day 0.7 to 1.7 against tau_patient 0.3 to 3.5), which is the argument for keeping the day term rather than ignoring it.

BayesFlow, both variants fail the acceptance gate and the decision rests on the yardstick:

| Variant | Calibration | Agreement with brms | Other | Verdict |
|---|---|---|---|---|
| v1, tube-summary simulator, 30k simulations, 20 epochs, 2 h 20 min CPU | SBC pass on 96 percent, coverage pass on 100 percent | 0.94 (gate 0.95); the 27 misses are days with one or two neg tubes, where v1 shrinks towards zero more than brms | median recovery R squared 0.56 | fails narrowly; a larger budget on highmem may pass |
| v2, single-cell simulator, 15 epochs of 100 batches, 3 h 15 min CPU | SBC pass on 96 percent, coverage pass on 92 percent | 0.76 | recovers the three instrument shifts well (R squared 0.91 to 0.96) but no day effect (R squared 0), and the posterior predictive covers only 75 percent of real neg medians (gate 85 percent): the mixture simulator is too simple for the real distributions | fails |

What was learnt about BayesFlow here: it delivers calibrated posteriors and, in v2, a real single-cell view of the instrument shift, but with one neg tube per patient-day no method identifies day effects well, and the simulator quality, not the inference engine, is the limiting factor. It is not the right primary tool for this cohort. It would earn its place as a per-file QC scorer for future acquisitions once the v2 simulator reproduces the real neg distributions.

## 10. Open items and assumptions

- Cofactor and compensation: use the SPILL matrix in each file. Confirm the Symphony files carry a valid SPILL and that the Diva 9 export did not already apply it.
- Patient 16 carries both labels because the patient was treated twice (label 1 before the first treatment, 2 after relapse and retreatment). The per-tube label is correct and stays as is.
- The 13-parameter Symphony tubes (2 files) carry an extra `SP SSC-A` channel; it maps to nothing and is ignored. 172 files have no `$PnS` marker names at all, so channels are mapped by fluorochrome name, which is constant across the study (config `markers`).
- Positivity thresholds 0.05 percent and 3 times background are the working assumption. Adjust in the config if the lab uses another rule.
- No GPU anywhere. Variant 2 is sized for CPU. If training exceeds one night, reduce cells per tube to 1000 and simulations to 50k before changing the model.
