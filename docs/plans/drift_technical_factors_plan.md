[< Back to Documentation Index](../README.md)

# Drift analysis: user-declared technical factors (operator first), how correction works, and what is tested

Status: written 2026-10-12. Parts 1 and 2 are implemented and uncommitted; parts 3 to 6 are open.

## 1. What the user asked for

Users should be able to add any metadata they have to the sample sheet and say whether each column is technical (drift to detect and, if needed, remove) or biological (to keep and separate from drift). The first such column is the person who measured each tube.

## 2. What exists now (uncommitted)

| Piece | Behaviour |
|---|---|
| Sample sheet | `group_*` columns stay protected biology. New `tech_*` columns are technical factors (prefix stripped). A blank `tech_*` cell becomes the level `unknown`, with a warning. A `tech_*` column may not reuse a built-in name (`instrument`, `acq_date`, ...). |
| FCS header | A `tech_*` column replaces a header field of the same name, so `tech_operator` overrides `$OP`. A technical term listed in the config without a sheet column is read from `tubes.csv`, so `operator` works from `$OP` alone. |
| Config | `design.technical`: the first term other than `instrument` is the day (random intercept); every other term, instrument included, is a fixed factor contrasted against its reference. `design.reference_levels` sets the reference of a non-instrument factor (default: first level in sort order). `config/drift_mtb.yaml`: `[acq_day, instrument, operator]`, reference Patrick. |
| Models | brms, Stan (through the design matrix X, no Stan change) and lme4 take the extra fixed factors from `design_formula()`. Effects are written with `term = operator`. |
| Decision | Rules 1 and 2 run for every technical contrast. Rule 2 for an operator effect shifts the tubes of that operator (layer B) or the matched negs of that operator (layer BG, found through `neg_filename`). |
| Checks and report | `check_design.R` writes `technical_factors.csv` (tubes, days, days shared with another level, subjects per level). The report's forest plot and tables show every technical contrast. |
| Step 5 | Carries every design column. Older `gated_frequencies.csv` files are back-filled from `tubes.csv`. |
| Tests | 182 base-R tests pass (15 new). |
| Behaviour change | An extra technical term used to become a random intercept in brms and lme4, but Stan ignored it and step 9 never judged it. It is now a fixed factor in every fitter. |

## 3. MTB facts established on 2026-10-12

- The operator list matches the FCS `$OP` keyword on all 328 tubes (the files spell Anna as `AnnaD`). It confirms the header; it adds no new information.
- Operator is almost nested in day: one operator per day except 2023-06-07 (Anna and Yacine). It also follows calendar time: Yacine 2022-11 to 2023-06 and one Symphony day in 2024-08; Anna 2023-06 to 2024-03; Patrick 2024-05 to 2025-06. An operator effect cannot be told apart from anything else that changed between those periods (reagent lots, instrument service). Gating was done by one person, so gating style is not part of it.
- Five patients were measured by more than one operator (12, 16, 23, 61, 74). These within-patient contrasts carry the operator estimate.
- Separability simulation, 30 repetitions per scenario, on the real design:

| | without operator | with operator |
|---|---|---|
| week slope recovered (true 1.04) | 1.02 | 1.02 |
| day-SD leakage | 8 % | 12 % |
| false day detection, pure week effect | 0 % | 13 % |
| week bias, pure day shift (logit per year) | 0.002 | 0.011 |

  Adding operator costs some precision on the day term, but leaves the biology (week slope) intact. The simulation does not yet inject an operator effect (part 4.1).

## 4. How "correction" works today

No value in the pipeline has ever been corrected, and no step writes corrected frequencies or intensities. The analysis decides whether a correction would be needed.

1. The model estimates each technical effect: instrument, operator, each day.
2. Rule 1: is the effect larger than the same-day tube-to-tube noise, with posterior probability at least 0.95?
3. Rule 2, a counterfactual computed in memory: the posterior mean effect is subtracted on the logit scale from the affected tubes. For a background effect, only the matched neg is shifted. Then the positivity calls (≥0.05 % and ≥3× neg) are recounted. The number of changed calls is the "flips".
4. The verdict is a label in `decisions.csv`: `not detected`, `detected, not clinically relevant` (0 flips) or `detected, correct` (flips > 0).

Where a shift was removed on the MTB data, it was removed by people, not by the pipeline. The analysts placed the CD4 and CD8 gates lower on Symphony days, and the gate-position diagnostic shows this. Day drift in the raw signal disappears through gating and background subtraction.

The correction methods in section 8 of `mtb_batch_drift_plan.md` (quantile mapping of Symphony tubes, per-day offsets) are written down but not implemented, because no verdict on the 2026-09-21 and 2026-09-24 runs required them. The one relevant background day (2023-07-05, CD8 IFNg) was handled by flagging, not by changing values.

## 5. What has been tested, and what has not

| Comparison | Status | Result |
|---|---|---|
| Unit tests of every module (formulas, frames, control matching, flips, sheet) | done | 182 pass |
| Synthetic data with known instrument, day and background effects, steps 6 to 9 (`drift_synthetic.R`) | done for instrument and day | effects recovered (2026-09-24) |
| brms against hand-written Stan, real data | done | 423/423 A, 272/272 BG, 182/184 B |
| lme4 against brms | done | same conclusions |
| Week slope with and without technical terms, leave-one-day-out | done | CD8 IFNg moves 3.8 SE, others ≤1.3 SE |
| Separability of week and day by simulation | done, with and without operator | table in section 3 |
| BayesFlow calibration and agreement | done | both variants fail the gate |
| Operator list against FCS `$OP` | done | 328/328 agree |
| Step 1 with operator sheet against the 2026-09-24 tube table | done | identical apart from the new columns |
| Synthetic recovery of an operator effect | open | 6.1 |
| Real-data fits with operator in the model | open | 6.3 |
| Decisions with operator against decisions without (same data) | open | 6.3 |
| Positivity calls before and after an applied correction | not possible yet: nothing is ever applied | 6.5 |
| External check with bridging controls (same PBMC aliquot every day) | not possible: none acquired | lab action, plan of 2026-09-24 section 7 |

## 6. Plan

| Step | Work | Verify | Effort |
|---|---|---|---|
| 6.1 | Inject an operator shift into `simulate_endpoint()` and `drift_synthetic.R` (truth rows with `term = operator`) | Operator effect recovered inside its 95 % interval; leakage into week slope below 5 %; decision step labels it detected when injected above replicate noise and not detected when zero | half a day |
| 6.2 | Re-fit the two unconverged layer B models of 2026-09-24 (CD4 TNFa R-hat 1.10, CD8 IFNg 19 divergences) with more iterations | R-hat ≤ 1.01, 0 divergences | part of 6.3 |
| 6.3 | Real-data rerun in `results_mtb_drift_2026-10-12_operator/`: 6c lme4, 6d biology, 6a brms (about 1.5 h), 6b Stan, 9, 10. Then a side-by-side table: every instrument and day verdict with and without operator, plus the operator contrasts (Anna and Yacine against Patrick) per response | brms/Stan agreement as before; every changed verdict listed and explained | half a day, mostly waiting |
| 6.4 | Interpretation rule for operator: a detected operator effect is reported as "operator or period" because of the nesting in section 3. It is corrected only if rules 1 and 2 both pass and the five multi-operator patients show the same direction | Written into the report text | one hour |
| 6.5 | Only if 6.3 yields `detected, correct` for any technical term: implement the correction (per-level offset on the logit scale for frequencies, written to `gated_frequencies_corrected.csv`, never overwriting the originals), recount calls and report before and after | Recounted calls equal the rule 2 flip count; synthetic correction removes an injected shift | one day |
| 6.6 | Documentation: sample sheet section of `docs/scripts/drift.md` (`tech_*` columns, `reference_levels`, header override), the README of `scripts/drift`, and section 9 of `mtb_batch_drift_plan.md` (the inconsistent count of layer B day effects, 9 against 13, gets resolved from `decisions.csv`) | Walkthrough runs as written | two hours |

Decisions for the user before 6.3:

1. Keep operator in the default MTB config, or run it as a sensitivity analysis next to the 2026-09-24 model? The simulation favours keeping the primary model without operator and reporting operator as a sensitivity analysis, because adding it raises false day detection from 0 to 13 percent while the operator contrast is confounded with calendar period.
2. Should 6.5 be built in any case, so that the pipeline can produce corrected values for future cohorts, or only when a verdict requires it?
