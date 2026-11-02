# Handover: COMPASS workflow for the MTB ICS cohort

Written 2026-10-14. Nothing in this plan is implemented yet. Read sections 1, 6 and 9 before writing code.

## 1. Question and why COMPASS

The user wants to know **which T-cell population drives the biological effect**: CD4 or CD8, and which cytokine combination, and whether it changes with weeks since first vaccination (primary) and prior peptide vaccination (`prior_vacc`, secondary).

Today the pipeline cannot answer this:

- Step 6d (`06_biology.R`) fits each of the 8 single-cytokine responses (CD4/CD8 × IFNg, IL-2, TNFa, CD107a) separately. Ranking separate fits by significance is not a comparison between populations (Gelman & Stern 2006).
- A cell producing IFNg and TNFa is counted in both single gates, so the 8 responses overlap.

COMPASS (Lin et al. 2015, Nat Biotechnol 33:610) works on the Boolean combinations of the cytokine gates, so each cell is counted once. Per stimulated sample it gives:

- the posterior probability that each subset responds above its unstimulated control;
- a functionality score (FS) and a polyfunctionality score (PFS).

It was built for this question in ICS vaccine studies. It is listed as "not done" in `drift_three_questions_plan.md`.

## 2. What the current outputs already say

From `results_mtb_drift_2026-10-13_review/week_slope_within_day.csv` (logit per year):

| Response | Within patient-day | Between days |
|---|---|---|
| CD4 IFNg | 5.9 ± 1.3 | 0.3 ± 0.5 |
| CD4 TNFa | 5.0 ± 1.2 | −0.2 ± 0.5 |
| CD4 IL-2 | 3.6 ± 0.8 | 0.8 ± 0.3 |
| CD8, all four | 0.4 to 1.0 ± 1.1 to 1.6 | 0.3 to 1.1 ± 0.3 to 0.6 |

- In the full model (`week_slope_sensitivity.csv`) no week slope is more than about 1.7 SE from zero.
- The within-day signal points to CD4 Th1 cytokines. It rests on only 15 patient-days with more than one visit.
- COMPASS should show whether this is one polyfunctional CD4 subset or several single-cytokine ones.

## 3. Data facts (checked 2026-10-14)

- **Tubes.** 328 FCS files, 327 after exclusion: 272 stimulated, 55 negs.
- **Samples.** 56 patient × week samples carry stimulated tubes, with 2 to 10 tubes each.
- **Peptide pools.** The pools in `cocktail` are patient-specific, about 130 distinct IDs, mostly used 1 to 5 times. COMPASS therefore cannot be run per pool. Treat all peptide pools as one stimulation, as step 13 already does for MIMOSA.
- **Gating tree.** `/Lymphocytes/Single Cells/Live cells/{CD4,CD8}/{CD107a, IFNy|INFy, IL-2, TNF}`.
  - All four cytokine gates exist under both parents in every tube.
  - FlowJo spells IFNg as `IFNy` or `INFy`. `canon_gate()` in `05_gated_frequencies.R` maps the names onto `cfg$cytokine_markers`. Reuse that mapping.
- **Parent sizes.** CD4 counts range from 4 to 142,475 (median 6,004). CD8 counts range from 0 to 40,927 (median 2,394). Tiny parents must be dropped. Step 15 in the review tree has `min_parent_sensitivity.csv` to choose the cut.
- **Neg matching.** `match_controls()` gives 268 exact matches, 3 `same_patient_day` and 1 `same_patient_reference_instrument`. One neg serves every stimulated tube of its sample.
- **Data locations.**
  - Raw data: `/mnt/volume/data/MTB_FACS/{fcs,wsp}`.
  - Latest full results: `/mnt/volume/data/MTB_FACS/results_mtb_drift_2026-10-13_review/`.
- **COMPASS package.** Version 1.48.0 is in Bioconductor 3.22 (the `mtb_drift` env) but **not installed**. MIMOSA is installed (from GitHub).

## 4. COMPASS facts verified in the 1.48.0 source

- **Entry point.** `SimpleCOMPASS(n_s, n_u, meta, individual_id, iterations = 1e4, replications = 8, seed = 100)` fits from count matrices. The `COMPASS()` path needs a GatingSet container and is not needed here.
- **Column names and count.** `n_s` and `n_u` need one column per Boolean subset, named `IFNg&!IL2&TNFa&!CD107a`, i.e. all `2^M` combinations.
- **Null category last.** The all-negative category must be the **last column**: `COMPASS-discrete.R` treats column K as the null. `COMPASS()` builds it as `total − sum(kept subsets)`.
- **Filtering.** `COMPASS()` drops sparse subsets with `category_filter = function(x) colSums(x > 5) > 2` on the stimulated counts, and their cells go into the null. `SimpleCOMPASS` does no filtering, so apply the same filter before calling it and rebuild the null column.
- **Row order.** Rows of `n_s` and `n_u` and the rows of `meta` must share the individual IDs. If they are not identical, `SimpleCOMPASS` silently reorders by name, so pass them already sorted.
- **Outputs.**
  - `fit$fit$mean_gamma` is the I × K posterior response probability.
  - `fit$fit$gamma` holds the gamma draws (I × K × N) of the last replication, so FS and PFS can be computed per draw.
  - `FunctionalityScore()`, `PolyfunctionalityScore()`, `COMPASSMCMCDiagnosis()` (R-hat on the alpha hyperparameters).
- Source tarball: `https://bioconductor.org/packages/3.22/bioc/src/contrib/COMPASS_1.48.0.tar.gz`.

## 5. Proposed workflow

New files only, so nothing collides with the unmerged review patch:

- `scripts/drift/16_compass.R`
- `src/harmonization/drift/drift_compass.R`
- `test_data/test_drift_compass.R`

Step number 16, because the review tree already uses 15.

### 5.1 Boolean counts (needs workspaces)

For every tube and each parent in `cfg$parent_populations`:

1. Open the GatingSets with `open_workspace()`, the same loop as steps 5 and 11.
2. Take the parent's event indices and the four cytokine child indices with `flowWorkspace::gh_pop_get_indices()`.
3. Map the gate names with the step 5 mapping, encode each parent event as a 4-bit code, and tabulate all 16 codes.

Write `boolean_counts.csv` with columns `filename, parent, subset, count, parent_count`, using COMPASS subset names and the null named explicitly.

The step stops with an error unless both checks pass, since they cost nothing:

- the 16 counts sum to `parent_count`;
- each cytokine's marginal equals `count` in `gated_frequencies.csv` exactly.

Put the counting in a pure function (indices in, counts out) so it can be tested without flowWorkspace.

### 5.2 COMPASS fit (minutes)

1. Run one `SimpleCOMPASS` fit per parent, with `individual_id = filename` of the stimulated tube.
2. Each `n_u` row is the Boolean counts of that tube's matched neg, taken from `gated_frequencies.csv` columns `neg_filename` and `neg_match`.
3. Drop tubes whose stimulated or neg parent count is below the chosen minimum, and record them.
4. Apply the category filter and put the null column last.
5. Run 3 seeds, and record R-hat and the correlation of FS across seeds (target > 0.99).

Outputs:

- `compass_gamma.csv`: `filename, parent, subset, mean_gamma`.
- `compass_scores.csv`: `filename, parent, FS, PFS`, plus FS and PFS 2.5 and 97.5 % from the gamma draws.
- `compass_fits/*.rds`.
- A COMPASS heatmap per parent, annotated with week bin, `prior_vacc` and instrument.

### 5.3 Biology (brms, the yardstick engine)

1. **Primary model.** Model FS jointly over both parents:
   ```
   FS ~ parent * (week_c + week_nl1..5 + prior_vacc) + instrument + operator
        + (1 | acq_day) + (0 + parent | patient) + (1 | patient:week) + (1 | cocktail)
   ```
   - FS lies in [0, 1]: use a Beta family on FS squeezed into (0, 1), or logit-normal.
   - Propagate COMPASS uncertainty by passing the per-tube FS draw SD with `se()`, or by refitting on several gamma draws.
   - Build the formula through `design_formula()`-style helpers and the config `design` block, as the other steps do.
2. **Answer to "which population".** Report posterior contrasts: the CD4 minus CD8 week slope with its interval and P(> 0), and the same for `prior_vacc`.
3. **Subset level.** Keep subsets with `mean_gamma > 0.5` in at least 10 % of tubes of a parent.
   - Fit their stimulated counts jointly: binomial on the parent, observation-level random effect, neg counts entered as an offset or covariate.
   - Fixed terms: `subset * (week terms + prior_vacc)`, plus the same technical and random terms.
   - This shows whether one polyfunctional subset or several carry the week effect.
4. **PFS.** Same as the primary model, as a secondary readout.

lme4 versions of the same models are a quick cross-check, as step 6c already provides.

## 6. Confounds to settle before drawing conclusions

1. **Storage time versus week.**
   - Several visits of one patient acquired on the same day means stored samples thawed together. Within such a day, a later visit was stored for less time, so week and storage time are perfectly collinear.
   - Between days, week mixes with day drift instead.
   - The within-day CD4 slopes in section 2 could be storage loss rather than vaccine biology.
   - Collection or freezing dates are needed to add `storage_years` as a covariate. Ask the user (section 9).
2. **Hand gating.** Cytokine gates were set by hand per workspace, and the 2026-10-13 review found this removed the Symphony shift for 6 of 8 readouts. Boolean subsets inherit the hand gates.
   - Sensitivity run: Boolean counts from one fixed threshold per channel. Reuse the thresholds that `fixed_gate_frequencies()` in the review tree derives from per-tube histograms, and apply them per cell in 5.1, since Boolean combinations need joint per-cell values.
3. **Dependence.** Tubes of one sample share their neg and their visit. COMPASS treats tubes as independent individuals. The downstream model must therefore carry `(1 | patient:week)`, the `design.visit` term in the review tree, and `(1 | cocktail)`.
4. **Fallback negs.** Four stimulated tubes use a fallback neg. Report them and rerun without them.

## 7. Validation (success criteria)

| Check | Pass |
|---|---|
| Unit tests in `test_data/test_drift_compass.R` (base R, existing harness): 4-bit coding, subset names, null last, filter and null rebuild, neg pairing | all pass |
| Boolean marginals versus `gated_frequencies.csv` on the real data | exact for every tube, parent and cytokine |
| Synthetic cohort (`test_data/drift_synthetic.R` design): one CD4 polyfunctional subset with a week effect, CD8 without | COMPASS flags that subset; the CD4 − CD8 week contrast excludes 0; CD8 contrast covers 0 |
| Synthetic null: day drift, no week effect | week contrast covers 0 |
| MCMC: 3 seeds | FS correlation > 0.99, R-hat reported |
| Agreement with step 13 | kappa between "any IFNg-containing subset responds" (COMPASS) and MIMOSA IFNg calls, reported in `calls_agreement.csv` style |

## 8. Practicalities and house rules

- **Environment.**
  - Install with `BiocManager::install("COMPASS")` in `mtb_drift` (`export MAMBA_ROOT_PREFIX=/mnt/volume/micromamba/root; /mnt/volume/micromamba/bin/micromamba run -n mtb_drift ...`).
  - Add COMPASS to `environment_drift.yml`, `00_setup_env.sh` and the `doctor` package list in `drift.R`.
  - On highmem1 the env is conda `cytometry_env` (rstan, not cmdstanr). COMPASS does not need Stan.
- **Temporary files.** flowWorkspace steps need `TMPDIR=/mnt/volume/tmp_drift`, because `/tmp` crashed step 11 with an HDF5 error and the root disk has about 1 GB free.
- **Long runs.** Detach them with `setsid nohup`; the user drives the session from a PC that may sleep.
- **Step registry.** Register step 16 in `STEPS` in `drift.R` with a small hunk. The unmerged review patch (`/mnt/volume/drift_review/work/review.patch`) also edits `drift.R`.
- **Config.** Add a `compass:` block to `config/drift_mtb.yaml` and `config/drift_template.yaml`: parents, cytokines (default `cytokine_markers`), category filter thresholds, minimum parent count, iterations, replications, seeds.
- **Style.** Match the scripts: `common.R`, `parse_args()`, `write_table()`, everything driven by the config, base-R tests, British spelling, no new dependencies besides COMPASS.
- **Git.** The agent never commits or pushes. Stage nothing unasked; give the user the commit commands as text.

## 9. Open questions for the user

1. Are sample collection or freezing dates available? Without them the week effect cannot be separated from storage time within day (section 6.1).
2. Unit of analysis: one COMPASS individual per stimulated tube (recommended, comparable with step 13), or pools summed per sample (sensitivity)?
3. Primary gates: the FlowJo hand gates (recommended, what the lab reports) with fixed thresholds as sensitivity, or the other way round?
4. Keep CD107a in the Boolean set? It is a degranulation marker, commonly included; dropping it halves the number of subsets.

## References

- Lin L, Finak G, Ushey K, et al. COMPASS identifies T-cell subsets correlated with clinical outcomes. Nat Biotechnol 2015;33:610–616.
- Finak G, McDavid A, Chattopadhyay P, et al. Mixture models for single-cell assays with applications to vaccine studies. Biostatistics 2014;15:87–101 (MIMOSA, step 13).
- Gelman A, Stern H. The difference between "significant" and "not significant" is not itself statistically significant. Am Stat 2006;60:328–331.
- Gelman A, Hill J, Yajima M. Why we (usually) don't have to worry about multiple comparisons. J Res Educ Eff 2012;5:189–211.
