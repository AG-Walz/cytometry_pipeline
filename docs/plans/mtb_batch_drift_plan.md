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
| Simulation-based calibration | Rank histograms uniform for all inferred parameters |
| Recovery on held-out simulations | R squared above 0.9 for `b_instrument`, `c_day` |
| Agreement with yardstick | Posterior means within the brms 95 percent interval for all 23 days and the instrument contrast |
| Posterior predictive on real negs | Real summary statistics inside the 90 percent predictive band |
| Cost | Reported honestly: wall-clock, simulations, tuning rounds |

If variant 2 fails any check, the decision falls back to the yardstick and the report says so.

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

Implemented under `scripts/mtb_drift/` (thin CLIs) and `src/harmonization/drift/` (modules). Base-R modules are covered by `test_data/test_mtb_drift.R`; `test_data/mtb_drift_synthetic.R` generates step 4/5 tables with known effects so steps 6 to 9 can be exercised without CytoML.

```bash
bash scripts/mtb_drift/00_setup_env.sh cytometry_env ~/bfvenv      # once
screen -S mtb_drift
bash scripts/mtb_drift/run_all.sh /mnt/volume/patrick/UMAP results/mtb_drift cytometry_env ~/bfvenv
```

`run_all.sh` expects `fcs/` and `wsp/` under the data root (override with `FCS_DIR` and `WSP_DIR`), writes `run_all.log` into the results folder, and honours `SKIP_BF=1`, `CORES`, `WORKERS`, `BF_EPOCHS`, `BF_SIMS`, `BF_CELLS`.

| Step | Script | Output | Needs | Approx. time |
|---|---|---|---|---|
| 0 | `00_setup_env.sh` | brms, cmdstanr, yaml into the conda env; CPU BayesFlow venv | conda, internet | 1 h |
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

Verification status on 2026-09-21: steps 1 and 2 ran on the real data here; steps 7, 8 and 9 ran end to end on synthetic tables in smoke mode; steps 3 to 6 and 10 are syntax-checked only because this box has no CytoML, brms or Stan, so their first real run is on highmem.

## 10. Open items and assumptions

- Cofactor and compensation: use the SPILL matrix in each file. Confirm the Symphony files carry a valid SPILL and that the Diva 9 export did not already apply it.
- Patient 16 carries both vaccination labels across tubes. Confirm the correct label before layer B fitting.
- The 13-parameter Symphony tubes (2 files) carry an extra `SP SSC-A` channel; it maps to nothing and is ignored. 172 files have no `$PnS` marker names at all, so channels are mapped by fluorochrome name, which is constant across the study (config `markers`).
- Positivity thresholds 0.05 percent and 3 times background are the working assumption. Adjust in the config if the lab uses another rule.
- No GPU anywhere. Variant 2 is sized for CPU. If training exceeds one night, reduce cells per tube to 1000 and simulations to 50k before changing the model.
