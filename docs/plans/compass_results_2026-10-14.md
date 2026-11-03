# Polyfunctional T-cell responses in the MTB ICS cohort: results of step 16

Run on 2026-10-14 on the dev box (micromamba env `mtb_drift`, R 4.5, COMPASS 1.48.0, brms with CmdStan 2.39),
following `compass_plan_2026-10-14.md`. All outputs are in
`/mnt/volume/data/MTB_FACS/results_mtb_compass_2026-10-14/`. Code: `scripts/drift/16_compass.R`,
`src/harmonization/drift/drift_compass.R`, tests in `test_data/test_drift_compass.R` (40 pass).

## 1. Answers in brief

1. **Which marker combinations respond to the peptides.** In CD4 T cells, every combination of two or more
   markers that contains IFNg, IL-2 or TNFa is increased by peptide stimulation, at every visit after
   vaccination. Single-marker CD4 cells barely respond (IFNg-only cells never do). The strongest and most
   frequent responders are IFNg+IL-2+TNFa cells with or without CD107a. In CD8 T cells the response is
   small and confined to CD107a-containing polyfunctional subsets (mainly IFNg+TNFa+CD107a), seen at the
   first two or three visits; IFNg-only CD8 cells decrease at some visits.
2. **CD4 and CD8 both activated at later timepoints (H1).** Not supported. CD4 respond clearly at the late
   visits (net 1.7 % of CD4 T cells induced, P > 0.99). CD8 show no net response at the late visits; the
   probability that both respond late is 0.12 (visit model) and 0.50 (weeks model, i.e. a coin toss).
3. **More markers per responding cell later (H2, pre-specified as the share of induced cells with >= 2
   markers).** Not shown. In CD4 the share is already about 80 % at V1 and stays there (change +0.01,
   95 % interval -0.16 to +0.21). The weeks model leans towards an increase (+0.15, -0.05 to +0.31,
   P = 0.93) without excluding zero.
4. **More subsets respond later (H3).** The pre-specified measure turned out to be uninformative (its
   interval spans about -6 to +4 subsets even on simulated data with a planted effect), so it can neither
   support nor refute the hypothesis. Counting subsets that are clearly up (P >= 0.95), of 14 per parent:
   CD4 10 (V1), 9 (V2), 11 (V3), 9 (V4); CD8 1, 3, 1, 0. No rise in CD4 beyond V3, and none in CD8.
5. **Exploratory (chosen after seeing 2 to 4): the polyfunctional CD4 response grows.** The stimulation
   effect of CD4 subsets with >= 2 markers is larger at V3/V4 than at V1 (+0.24 log odds ratio, P = 0.986,
   visit model; +0.46, P > 0.999, weeks model), and grows more than that of single-marker subsets (+0.19,
   P = 0.988; +0.56, P > 0.999). So late CD4 responses are stronger and more polyfunctional relative to the
   unstimulated control. CD8 responses do not grow; in the weeks model they shrink (P(increase) = 0.004).
6. **Prior peptide vaccination** makes no difference to the CD4 response; for CD8 there is a weak hint of
   a larger response in previously vaccinated patients (section 7).

The two biggest caveats (section 9): the late visits are mostly different patients from the early ones,
and storage time is confounded with weeks. A refit in which every patient has their own stimulation
effect leaves all of the above unchanged (section 10); storage time remains open until freezing dates are
available.

## 2. What was done

| Stage | What | Output |
|---|---|---|
| 16a Boolean counts | For all 327 tubes and both parents, per-cell membership of the four FlowJo hand gates (IFNg, IL-2, TNFa, CD107a) coded into 16 Boolean subsets | `boolean_counts.csv` |
| 16b COMPASS | One COMPASS fit per parent, 272 stimulated tubes each against their matched neg, 14 shared subsets, 3 seeds, 10,000 iterations x 8 replications. The response map and the scores are the **mean model**, `mean_gamma` averaged over the three seeds, with the seeds' draws pooled for the intervals (section 3.1) | `compass_gamma.csv`, `compass_scores.csv`, `compass_mcmc.csv`, heatmaps |
| 16c subset model | Beta-binomial brms model of stimulated and neg counts of the 14 subsets x 2 parents (9,156 rows); visit or weeks as time; 4 chains x 2,000 iterations; about 5 h per fit | `subset_effects_*.csv`, `hypotheses_*.csv`, `compass_model/` |
| contrasts | Exploratory late vs V1 contrasts from the saved fits | `exploratory_late_vs_ref_*.csv` |
| summary | Patient x pool table, per-visit table, figures | `compass_pools.csv`, `compass_by_visit.csv`, `figures/` |

Subset model (visit version):

```
count | trials(parent_count) ~ parent * stim * visit + parent * stim * prior_vacc + instrument + operator
    + (1 + stim + stim:visit | ps) + (1 | ps:patient) + (1 | ps:visit_id) + (0 + stim | ps:cocktail) + (1 | acq_day)
phi ~ 1 + (1 | ps)                       # beta-binomial overdispersion per subset
```

`ps` is parent x subset. The stimulation effect of each subset is its log odds ratio stimulated vs neg
within the same sample (stimulated and neg tubes of one visit share `(1 | ps:visit_id)`), varying by
subset and visit with partial pooling towards the mean of its parent. The weeks version replaces visit by a
3-df natural spline of log(1 + week). Predictions are for a typical sample (patient, visit, pool and day
effects at zero, Canto, operator Patrick, not previously vaccinated).

Changes from the plan, all made before seeing real results:

- Beta-binomial instead of binomial with an observation-level effect: same purpose, sampled faster and
  mixed better in the smoke test.
- Fixed effects per parent (`parent * stim * ...`): with one shared mean, the smoke test pulled null CD8
  subsets towards the strong CD4 effects.
- One subset list for both parents (14 of 15 non-null subsets pass the COMPASS filter; only
  IFNg+IL-2 without TNFa and CD107a is too rare). Reading the COMPASS source showed that its functionality
  score divides by all 2^4 - 1 = 15 subsets, so CD4 and CD8 scores are on one scale.

## 3. Checks

| Check | Result |
|---|---|
| Boolean subsets sum to the parent | exact, all 654 tube x parent pairs |
| Marker marginals vs `gated_frequencies.csv` | exact, all 2,616 tube x parent x marker counts |
| Unit tests | 40 of 40 pass |
| COMPASS, 3 seeds | FS correlation 0.99992 (CD4), 0.99961 (CD8); per-tube response probabilities correlate 0.9998 / 0.9994, largest difference 0.07 / 0.12 |
| COMPASS hyperparameter R-hat | **Did not converge.** Gelman-Rubin on the Dirichlet hyperparameters across the three seeds reaches 6.3 (CD4) and 5.4 (CD8), against a target near 1. COMPASS gives no threshold of its own: `?COMPASSMCMCDiagnosis` says only that R-hat "should be near 1 but rarely are in practice. Very large values may be a concern". Six is a large value. The largest values sit on subsets in which almost no tube responds, whose stimulated-side parameter is barely informed, but that is an explanation and not a remedy. The agreement between seeds in the row above is supporting evidence that the quantity we read off, the per-tube response probability, is stable; it is **not** a convergence diagnostic, because chains can agree on one marginal while the hyperparameters they condition on have not mixed. Treat the COMPASS response map as descriptive, and take the inference from the subset model, whose own diagnostics are clean |
| COMPASS vs MIMOSA (any IFNg-containing subset vs MIMOSA IFNg call) | kappa 0.74 (CD4), 0.82 (CD8); vs the lab rule 0.48 / 0.64 |
| COMPASS without the low-information tubes | Refitted dropping 34 CD4 and 107 CD8 tubes below 1,000 parent events: response probabilities correlate 0.999 / 0.994 with the full run, under 1 % of calls flip, no subset moves by more than 2.9 percentage points (section 3.2) |
| brms convergence, 6 fits | max R-hat 1.009 to 1.013, no divergent transitions, min effective-sample ratio 0.05 to 0.14 |
| Synthetic data, planted effect (CD4 IFNg+IL-2+TNFa rising 1 -> 1.5 -> 2.5 -> 3 log OR over V1 to V4, CD8 none) | recovered as 0.74 -> 1.24 -> 1.99 -> 2.19 (large effects shrink somewhat); 95 % intervals cover the truth in 98 % of subset x visit cells; 3 of 132 null cells called (borderline); H2 detects the planted shift (+0.30, P = 0.976) |
| Synthetic null (constant responses, day drift) | 99 % coverage; 6 of 132 null cells called, 4 of them at V0; H2 near 0 for both parents (P 0.39, 0.26) |
| Lesson from the synthetic runs | H3 as defined is uninformative (interval about -6 to +4 subsets in every scenario); H1's net response has little power because it sums noise from large background subsets; V0 effects are pulled towards V1 (few V0 tubes), so V0 calls are not reliable |

### 3.1 Mean model over the seeds, and a bug in the package's own averaging

COMPASS's documentation recommends running several seeds and reports that the averaged model "should be
better then any of the individual models". The response map and the scores in section 4.1 are therefore no
longer taken from the first chain: `mean_gamma` is averaged over the three seeds, and the three seeds'
gamma draws are pooled for the FS and PFS intervals. Earlier versions of this document used seed 100 alone.

It changes almost nothing, which is the reassuring outcome: the mean absolute change in a per-tube response
probability is 0.001, the largest is 0.051, and 3 of 7,616 calls at the 0.5 cut flip. The number of tubes
with any response is identical for both parents (CD4 207, CD8 60), as are the agreement statistics with
MIMOSA and the lab rule. The seed-1 value is kept beside the mean in `compass_gamma.csv` as
`mean_gamma_seed1`, and `compass_mcmc.csv` records the largest difference and the number of flipped calls.

`COMPASS::COMPASSMCMCDiagnosis()` also advertises a `mean_model`, but **it does not average**. It calls
`abind()` through `Map()`, so each fit is bound into its own three-dimensional array and `[[1]]` keeps only
the first; the returned `mean_model$fit$mean_gamma` is identical to seed 1. Verified against COMPASS 1.48.0
on these fits: `identical(pkg_mean_model, seed1_mean_gamma)` is `TRUE`. The averaging here is done directly
in `mean_gamma_across_seeds()` instead. The R-hat values this function returns are unaffected and are still
used. Worth reporting upstream.

### 3.2 Sensitivity: COMPASS without the low-information tubes

The subset model was already refitted without tubes below 1,000 parent events (section 5), but COMPASS
itself had not been, and its CD8 map in particular rests on many small denominators. COMPASS was therefore
refitted with `--variant no_low_info`, dropping a tube when its own or its matched control's parent count
is below 1,000: 34 of 272 tubes for CD4 and **107 of 272 for CD8**, so this is a real test rather than a
formality. Three seeds again, the same 14 subsets as the full run so the maps are directly comparable
(`run_16b_no_low_info.sh`, `compare_16b_variants.R`).

The response map barely moves. On the tubes the two runs share:

| | CD4 | CD8 |
|---|---|---|
| tubes kept of 272 | 238 | 165 |
| per-tube response probabilities, correlation | 0.9989 | 0.9941 |
| calls flipped at the 0.5 cut | 30 of 3,332 (0.9 %) | 8 of 2,310 (0.35 %) |
| largest change in a single probability | 0.12 | 0.19 |
| tubes with any response, full to no_low_info | 193 to 196 | 53 to 58 |
| functionality score, correlation | 0.9995 | 0.9960 |
| polyfunctionality score, correlation | 0.9988 | 0.9994 |

Per subset, the share of tubes called responding shifts by at most 2.9 percentage points (CD4 IFNg+IL-2+TNFa
and IFNg+TNFa+CD107a), and every shift is upward: dropping the noisy tubes sharpens the shared prior slightly
rather than changing what it says. The ordering of subsets is unchanged in both parents, the same subsets sit
at zero, and CD8 remains confined to CD107a-containing combinations led by IFNg+TNFa+CD107a. Nothing in
section 4.1 changes.

One diagnostic does improve. The hyperparameter R-hat of section 3 falls from 6.3 to 3.5 for CD4 once the
34 small tubes are gone, which supports the reading that the poor mixing comes from subsets and tubes that
carry almost no information. It does not improve for CD8 (5.4), where 107 tubes were dropped and the
remaining counts are still small, so the convergence problem is not solved by this refit either.

Full tables in `compass_low_info_sensitivity.csv` and `compass_low_info_sensitivity_by_subset.csv`; the
variant's own figures are the `*_no_low_info.png` set, including its heatmaps.

## 4. Which subsets respond

### 4.1 Per tube (COMPASS)

Share of the 272 stimulated tubes in which a subset responds (posterior probability > 0.5):

| CD4 subset | tubes | CD8 subset | tubes |
|---|---|---|---|
| IFNg+IL-2+TNFa+CD107a | 54 % | IFNg+TNFa+CD107a | 12 % |
| IFNg+IL-2+TNFa | 44 % | IFNg+CD107a | 6 % |
| IL-2+TNFa+CD107a | 34 % | TNFa | 5 % |
| IFNg+TNFa+CD107a | 32 % | TNFa+CD107a | 4 % |
| IFNg+TNFa | 30 % | IFNg+IL-2+TNFa+CD107a | 4 % |
| TNFa+CD107a | 24 % | CD107a | 2 % |
| TNFa | 22 % | IFNg | 1 % |
| IL-2+TNFa | 21 % | all others | 0 % |
| IL-2 | 10 % | | |
| CD107a, IFNg+CD107a | 5 % each | | |
| IFNg, IL-2+CD107a, IFNg+IL-2+CD107a | 0 % | | |

- At least one CD4 subset responds in 87 % (V1), 78 % (V2), 81 % (V3) and 42 % (V4) of tubes; for CD8,
  32 %, 24 %, 15 % and 3 %.
- Among responding CD4 subset calls, 84 to 91 % involve >= 2 markers at every post-vaccination visit.
- Per patient x pool (`compass_pools.csv`): 119 of 135 pools have at least one responding CD4 tube, and
  every one of the 20 patients has at least one responding pool. For CD8, 36 pools in 13 patients.
- Figures: `figures/16_compass_heatmap_CD4.png`, `16_compass_heatmap_CD8.png`, `16_pfs_by_visit.png`,
  `16_responding_degree_by_visit.png`.

**How to read the heatmaps, and what COMPASS assumes.** COMPASS pairs every individual with its own
unstimulated control and treats individuals as independent and exchangeable; the shared hyperparameters
are what let it borrow strength. Here an individual is one stimulated tube, which is the unit the package
itself documents, and it is the right choice: had the patient been the individual, COMPASS would have
silently *summed* that patient's visits and pools into one pseudo-sample ("There are multiple samples per
individual; these will be aggregated"). But neither assumption holds exactly in this cohort, and the
heatmap rows are not 272 independent observations:

- **One control serves several stimulated tubes.** The 272 stimulated tubes of each parent draw on only
  55 distinct negative tubes. 265 of 272 share their control with at least one other stimulated tube, the
  median tube shares with five others, and one control is the reference for nine tubes. Tubes sharing a
  control rise and fall together: a negative that happened to be high makes every pool of that visit look
  unresponsive, and a low one makes them all look like responders. The model has no term for this, so the
  rows are correlated in blocks that the figure cannot show.
- **One patient appears many times.** 20 patients produced the 272 tubes, a median of 13.5 and up to 34
  tubes each. The shared prior, which encodes how often a subset responds across the cohort, is therefore
  estimated partly from pseudo-replicates, and patients with many pools and visits pull it towards their
  own behaviour. Consistency across rows overstates the cohort-level evidence.
- **COMPASS is one-sided.** Its hypothesis is that the stimulated proportion exceeds the unstimulated one.
  A probability near zero means "no evidence of an increase", never "a decrease". Every decrease reported
  in this document comes from the subset model, not from these figures.
- **The 0.5 cut is a convention.** The paper fixes no threshold; 0.7 and 0.9 are equally common in the
  literature, and the per-tube counts above shrink under a stricter cut.

Both heatmap figures now carry a short version of this note in their caption, so the structure travels with
the picture. The per-tube control sharing is in `compass_scores.csv` as `n_stim_sharing_neg`.

### 4.2 Cohort level (subset model, visit)

Stimulation effect as log odds ratio stimulated vs neg (exp gives the fold change in odds); + / -: posterior
probability of an increase >= 0.95 / <= 0.05. Figure: `figures/16_subset_effects_visit.png`.

| CD4 subset | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| IFNg+IL-2+TNFa+CD107a | 0.67+ | 1.00+ | 0.97+ | 1.46+ | 1.03+ |
| IFNg+IL-2+TNFa | 0.63+ | 0.99+ | 0.87+ | 1.54+ | 1.08+ |
| IFNg+IL-2+CD107a | 0.56+ | 0.71+ | 0.93+ | 1.12+ | 0.75+ |
| IFNg+TNFa+CD107a | 0.51+ | 0.82+ | 0.85+ | 1.32+ | 0.81+ |
| IL-2+TNFa+CD107a | 0.47+ | 0.78+ | 0.59+ | 1.31+ | 0.82+ |
| IL-2+TNFa | 0.48+ | 0.79+ | 0.61+ | 1.46+ | 0.86+ |
| IFNg+TNFa | 0.43+ | 0.71+ | 0.68+ | 1.16+ | 0.74+ |
| TNFa+CD107a | 0.24 | 0.55+ | 0.40+ | 1.04+ | 0.55+ |
| IFNg+CD107a | 0.34 | 0.50+ | 0.62+ | 0.86+ | 0.50+ |
| IL-2+CD107a | -0.04 | 0.00 | 0.07 | 0.14 | -0.01 |
| TNFa | 0.05 | 0.28+ | 0.11 | 0.54+ | 0.25 |
| IL-2 | -0.09 | 0.10 | -0.11 | 0.33+ | 0.07 |
| IFNg | -0.04 | 0.03 | 0.17 | 0.06 | -0.00 |
| CD107a | -0.15 | -0.05 | -0.24- | -0.00 | -0.09 |

A log odds ratio of 1.0 is a 2.7-fold, 1.46 a 4.3-fold increase over the matched neg.

| CD8 subset | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| IFNg+TNFa+CD107a | 0.30 | 0.34+ | 0.40+ | 0.44+ | 0.09 |
| IFNg+CD107a | 0.17 | 0.18 | 0.52+ | 0.32 | -0.09 |
| IL-2+CD107a | 0.03 | -0.05 | 0.45+ | 0.04 | -0.36 |
| IFNg+IL-2+CD107a | 0.21 | 0.12 | 0.57 | 0.30 | -0.16 |
| IFNg+IL-2+TNFa+CD107a | 0.16 | 0.17 | 0.27 | 0.27 | -0.10 |
| IFNg | -0.03 | -0.02 | -0.43- | 0.14 | -0.26 |
| CD107a | -0.03 | -0.01 | -0.08 | -0.06 | -0.27- |
| all other CD8 subsets | -0.3 to 0.4, none called | | | | |

The weeks model gives the same CD4 picture (all >= 2-marker subsets up at every time after vaccination,
rising with time) and for CD8 an early response of IFNg+IL-2+TNFa+CD107a (1.27+ at V1, 0.92+ at V2) and
IFNg+TNFa+CD107a (0.73+, 0.51+) that is gone at later times, plus a decrease of IFNg-only CD8 cells late
(-0.74-).

### 4.3 Size of the response

Net induced cells (sum over subsets of stimulated minus neg probability, i.e. the drop of the
all-negative subset), percent of the parent, typical sample:

| | V0 | V1 | V2 | V3 | V4 |
|---|---|---|---|---|---|
| CD4, visit model | 0.20 (-0.63, 1.07) | 0.88 (0.35, 1.62) | 0.49 (-0.01, 1.07) | 2.83 (1.57, 4.61) | 0.64 (0.08, 1.42) |
| CD8, visit model | -0.05 (-1.64, 1.58) | 0.09 (-0.76, 0.96) | -0.35 (-1.26, 0.45) | 0.10 (-0.75, 1.01) | -0.93 (-2.32, 0.05) |
| CD4, weeks model | 0.00 | 1.23 (0.52, 2.26) | 0.90 (0.22, 1.81) | 1.75 (0.90, 3.02) | 2.94 (1.52, 4.95) |

The weeks model predicts at the median week of each visit's tubes (V1 6, V2 15, V3 148, V4 228 weeks).
The two time axes disagree about V4: the visit model puts V4 below V3, the weeks model above. V4 rests on
three patients (12, 16, 17), and weeks and visit number are only loosely related (V2 spans 11 to 188
weeks), so neither trajectory beyond V3 is well determined.

## 5. Pre-specified hypotheses

| | Model | CD4 | CD8 |
|---|---|---|---|
| H1 net response at V3/V4 (fraction of parent) | visit | 0.0174 (0.0097, 0.0282), P = 1.00 | -0.0041 (-0.0128, 0.0027), P = 0.12 |
| | weeks | 0.0202 (0.0107, 0.0343), P = 1.00 | 0.0001 (-0.0076, 0.0081), P = 0.50 |
| H1 P(both respond late) | visit / weeks | 0.12 / 0.50 | |
| H2 share of induced cells with >= 2 markers, late - V1 | visit | +0.01 (-0.16, 0.21), P = 0.50 | +0.17 (-0.33, 0.66), P = 0.72 |
| | weeks | +0.15 (-0.05, 0.31), P = 0.93 | +0.01 (-0.49, 0.67), P = 0.42 |
| H3 subsets up, late - V1 | visit | +0.1 (-2.5, 2.5) | -1.8 (-7.0, 3.5) |
| | weeks | +0.3 (-2.0, 3.0) | -3.8 (-8.0, 1.0), P = 0.03 |

Share of induced CD4 cells with >= 2 markers by visit (visit model): 0.79 (V1), 0.91 (V2), 0.78 (V3),
0.82 (V4); weeks model 0.69 (V1) to 0.84 (late). The share is near its ceiling from the first visit,
which is why H2 has little room to move: CD4 responses are polyfunctional from the start.

Sensitivity runs, visit model: without tubes below 1,000 parent events (6,930 rows) and without the four
fallback-neg tubes (9,044 rows), the CD4 results stay within 0.01 of the full model (H1 0.0179 and 0.0171;
H2 +0.01 and +0.00). The CD8 values move inside their wide intervals (H2 +0.11 and +0.21; P(both respond
late) 0.14 and 0.06). No conclusion changes.

## 6. Exploratory: does the response grow from V1 to late visits?

Not pre-specified; added after H2 and H3 turned out insensitive. Late = V3/V4 (visit model) or the median
week of the V3/V4 tubes, 148 weeks (weeks model).

| Contrast, late - V1 | Visit model | Weeks model |
|---|---|---|
| CD4 mean log OR of >= 2-marker subsets | +0.24 (0.03, 0.46), P = 0.986 | +0.46 (0.24, 0.69), P = 0.9998 |
| CD4 mean log OR of single-marker subsets | +0.05 (-0.14, 0.24), P = 0.71 | -0.10 (-0.29, 0.08), P = 0.14 |
| CD4 polyfunctional minus single | +0.19 (0.02, 0.40), P = 0.988 | +0.56 (0.32, 0.82), P = 1.000 |
| CD4 net response (fraction of CD4) | +0.0086 (0.0018, 0.0174), P = 0.993 | +0.0079 (-0.0005, 0.0192), P = 0.97 |
| CD8 mean log OR of >= 2-marker subsets | -0.06 (-0.28, 0.17), P = 0.30 | -0.37 (-0.63, -0.10), P = 0.004 |
| CD8 mean log OR of single-marker subsets | -0.10 (-0.28, 0.10), P = 0.17 | -0.30 (-0.50, -0.11), P = 0.001 |
| CD8 polyfunctional minus single | +0.04 (-0.12, 0.20), P = 0.71 | -0.07 (-0.35, 0.22), P = 0.32 |

Read together: the CD4 response becomes stronger relative to background in the polyfunctional subsets,
and only there. This is the closest the data come to "more cells expressing multiple cytokines later",
but measured as fold over the unstimulated control rather than as the share of induced cells. CD8 shows
no late gain. Per-subset contrasts are in `exploratory_late_vs_ref_*.csv`; the largest CD4 gains are
IL-2+TNFa (+0.37, visit; +0.52, weeks), IFNg+IL-2+TNFa (+0.32; +0.53) and IL-2+TNFa+CD107a (+0.29; +0.81).

## 7. Prior vaccination, technical factors, pools and patients

- **Prior peptide vaccination** (40 tubes from patients 12, 16, 48): stimulation x prior_vacc -0.01
  (-0.19, 0.17) for CD4, i.e. no effect, and +0.17 (-0.02, 0.36, P = 0.96) for CD8 on the log OR scale, a
  weak hint of a larger CD8 response in previously vaccinated patients that rests on three patients.
- **Instrument and operator** shift the overall subset frequencies by amounts indistinguishable from zero
  (Symphony +0.08, -0.72 to 0.86; operator Anna -0.31, -0.92 to 0.27). They act on stimulated and neg
  tubes alike and largely cancel in the stimulation effect, because both tubes of a visit were acquired
  together. Acquisition day SD 0.70 on the logit scale.
- **Variation of the stimulation effect between pools** (SD 0.31 on the log OR scale per subset) is smaller
  than the variation of subset frequencies between patients (0.68) and visits (0.61).
- **Patient x pool table**: `compass_pools.csv` lists, for every patient, pool and parent, the tubes,
  visits, mean and maximum functionality score, polyfunctionality score and the subsets that respond in any
  tube.

## 8. Baseline (V0)

V0 has 8 stimulated tubes from 4 patients (59, 74, 79, 87), all acquired on the same day as their V1 and V2.
COMPASS finds small but consistent responses in two of them: patient 74 in CD4 IFNg+IL-2+TNFa+CD107a
(22 to 28 cells per 10,000 CD4 in the stimulated tubes against 8 in the neg) and patient 59 in CD8
IFNg+TNFa+CD107a (19 to 52 per 10,000 against 3). Patients 79 and 87 show nothing. Either these patients
had pre-existing T cells against the peptides, or V0 is not a clean baseline for them. The subset model's
V0 calls (CD4 polyfunctional subsets up) are not reliable: on simulated data the few V0 tubes were pulled
towards the V1 effect, and 4 of 6 false calls in the synthetic null were at V0.

## 9. Limitations

1. **Late visits are mostly other patients.** V3 has 7 patients, of whom 11, 39 and 57 have no earlier
   visit; V4 has 3, of whom 12 and 17 have no earlier visit. The model gives each pool its own response, but
   in the main fits not each patient, so a strongly responding patient seen only late could make "late"
   look stronger. Within the five patients with an early and a late visit, CD4 polyfunctionality goes up in
   two (16, 23), down in two (61, 74) and stays flat in one (69). The refit with a patient-specific
   stimulation effect (section 10) gives the same results, so this does not drive the findings; with 5
   patients seen early and late, the time contrasts still rest largely on comparisons between patients.
2. **Storage time.** Visits of one patient were often thawed and acquired together, so a later visit was
   stored for a shorter time. Better-preserved cells could respond more strongly, which would mimic the
   growth in section 6. Freezing dates are needed to separate the two.
3. **Exploratory contrasts.** Section 6 was chosen after seeing sections 4 and 5. Treat it as a hypothesis
   for the next cohort, not as a confirmed result.
4. **Pre-specified measures.** H3 cannot discriminate, and H1's net response is noisy (section 3). H2 is
   sound but has little room because CD4 responses are polyfunctional from V1.
5. **Killing.** CD107a marks degranulation, not killing; there is no perforin or granzyme in the panel.
   The CD107a single-positive gate has a high background (1.2 % of CD4 and 3.0 % of CD8 in negs) and does
   not respond.
6. **Hand gates.** All subsets inherit the hand-set gates, placed per visit by one person.
7. **COMPASS treats tubes as independent.** 265 of the 272 stimulated tubes share their unstimulated
   control with another stimulated tube (55 distinct controls, up to 9 tubes each), and the 272 tubes come
   from 20 patients. COMPASS has no term for either, so its response map overstates how much independent
   evidence there is. The subset model carries the dependence through its sample, patient and pool effects;
   COMPASS output is descriptive. Full statement of what this does and does not allow: section 4.1.

Data note found while plotting: patient 16's file names disagree on prior vaccination (`MTB016_2_*` at
V1/V2, `MTB16_1_*` at V4), so the model treats that patient's V4 tubes as not previously vaccinated. Please
check which is right.

## 10. Sensitivity with a patient-specific stimulation effect

The visit and weeks models refitted with `(1 + stim | ps:patient)` in place of `(1 | ps:patient)`, so
each patient has their own stimulation effect per subset on top of the pool effect (finished 2026-10-15
00:53 UTC; max R-hat 1.011, no divergent transitions). Patients differ in responsiveness by SD 0.27
(0.19 to 0.36) on the log OR scale, about as much as pools (0.32).

| Quantity | Visit, main | Visit, patient effect | Weeks, main | Weeks, patient effect |
|---|---|---|---|---|
| H1 CD4 net response late | 0.0174 (P 1.00) | 0.0169 (P 1.00) | 0.0202 (P 1.00) | 0.0215 (P 1.00) |
| H1 CD8 net response late | -0.0041 (P 0.12) | -0.0058 (P 0.07) | 0.0001 (P 0.50) | -0.0007 (P 0.42) |
| H1 P(both respond late) | 0.12 | 0.07 | 0.50 | 0.42 |
| H2 CD4 share change | +0.01 (P 0.50) | +0.01 (P 0.50) | +0.15 (P 0.93) | +0.10 (P 0.80) |
| CD4 >= 2-marker log OR, late - V1 | +0.24 (0.03, 0.46), P 0.986 | +0.23 (0.01, 0.46), P 0.977 | +0.46 (0.24, 0.69), P 0.9998 | +0.48 (0.25, 0.72), P 1.000 |
| CD4 polyfunctional minus single, late - V1 | +0.19, P 0.988 | +0.18, P 0.978 | +0.56, P 1.000 | +0.61, P 1.000 |
| CD8 >= 2-marker log OR, late - V1 | -0.06, P 0.30 | -0.06, P 0.31 | -0.37, P 0.004 | -0.38, P 0.006 |
| CD4 subsets clearly up V1 / V2 / V3 / V4 | 10 / 9 / 11 / 9 | 10 / 9 / 11 / 9 | 10 / 10 / 10 / 11 | 10 / 10 / 9 / 11 |

Every conclusion of sections 4 to 6 holds. The H2 lean in the weeks model weakens (P 0.93 to 0.80).

## 11. Files

| File | Content |
|---|---|
| `boolean_counts.csv` | filename, parent, subset, degree, count, parent_count (327 tubes x 2 parents x 16 subsets) |
| `boolean_marginal_check.csv` | Boolean marginals next to the step 5 counts |
| `compass_gamma.csv` | per tube and subset, COMPASS posterior probability of response |
| `compass_scores.csv` | per tube: FS and PFS with 95 % intervals, responding subsets, low-information flag, neg match |
| `compass_mcmc.csv`, `compass_fits/` | COMPASS diagnostics and fits, 3 seeds per parent |
| `compass_calls_agreement.csv` | COMPASS vs MIMOSA vs lab rule |
| `compass_pools.csv`, `compass_by_visit.csv` | per patient x pool and per visit summaries |
| `subset_effects_<model>.csv` | per subset and time: log OR, difference in percent, P(up), call |
| `hypotheses_<model>.csv`, `hypotheses_all_models.csv` | H1 to H3 and per-time levels |
| `exploratory_late_vs_ref_<model>.csv` | section 6 |
| `compass_model/` | brms fits (.rds), diagnostics, synthetic truths |
| `figures/` | heatmaps, PFS by visit (cohort mean with patient-resampled interval), PFS per patient and pool (`16_pfs_by_patient_CD4.png`, `_CD8.png`), responding subsets by degree, subset effects grid |
| `compass_low_info_sensitivity*.csv` | section 3.2, the COMPASS refit without low-information tubes |
| `compass_gamma_no_low_info.csv`, `compass_scores_no_low_info.csv`, `compass_mcmc_no_low_info.csv`, `figures/*_no_low_info.png` | that refit's own outputs |
| `pre_mean_model/` | the seed-1 tables the mean model replaced (section 3.1) |
| `run_16c.sh`, `run_16c_pstim.sh`, `run_16b_no_low_info.sh`, `compare_16b_variants.R`, `16*.log` | exact commands and logs |

## References

- Lin L, Finak G, Ushey K, et al. COMPASS identifies T-cell subsets correlated with clinical outcomes.
  Nat Biotechnol 2015;33:610-616.
- Finak G, McDavid A, Chattopadhyay P, et al. Mixture models for single-cell assays with applications to
  vaccine studies. Biostatistics 2014;15:87-101.
- Gelman A, Hill J, Yajima M. Why we (usually) don't have to worry about multiple comparisons.
  J Res Educ Eff 2012;5:189-211.
- Bürkner P-C. brms: an R package for Bayesian multilevel models using Stan. J Stat Softw 2017;80(1):1-28.
