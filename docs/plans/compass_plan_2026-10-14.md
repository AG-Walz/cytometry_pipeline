# Polyfunctional T-cell responses to peptide vaccination in the MTB ICS cohort: analysis plan

Version of 2026-10-14. Agreed with the study scientist in a plan review; nothing is implemented yet.
The command for `16_compass.R` in section 7 shows the intended interface of a script still to be written.
Extends `compass_workflow_handover.md` (same folder), whose goals it replaces where they differ.

## 1. Study and question

Patients received patient-specific peptide vaccines. PBMC from each visit were stimulated in vitro with
the patient's peptide pools and analysed by intracellular cytokine staining (ICS) on two cytometers
(FACSCanto II, FACSymphony A1). Each visit has one unstimulated control tube ("neg"). All ICS samples were
previously confirmed as responders by ELISpot, so a response is assumed to exist.

The question is which T-cell response the vaccine induces:

1. Which Boolean combinations of IFNg, IL-2, TNFa and CD107a are increased or decreased by peptide
   stimulation, in CD4 and in CD8 T cells.
2. How these responses change over time (visit number and weeks since first vaccination).
3. Whether polyfunctionality, i.e. cells co-expressing several markers, rises at later timepoints.

### Pre-specified hypotheses (fixed before fitting)

| ID | Hypothesis | Contrast |
|---|---|---|
| H1 | CD4 and CD8 both respond at late visits | P(stimulation effect > 0 in CD4 and in CD8 at V3/V4) |
| H2 | More markers per responding cell later | Share of stimulation-induced cells with >= 2 markers, late visits vs V1, per parent |
| H3 | More subsets respond later | Expected number of subsets with a positive stimulation effect, late vs V1, per parent |

"Not detectable" is an acceptable outcome for each; the contrasts are fixed so the result counts either way.

## 2. Why single-gate analysis is not enough

The existing pipeline (step 6) fits the 8 single-marker responses (CD4/CD8 x IFNg, IL-2, TNFa, CD107a)
separately. A cell producing IFNg and TNFa is counted in both gates, so the 8 responses overlap, and
ranking separate fits by significance does not compare populations (Gelman & Stern 2006). Boolean
combinations count each cell once and can show shifts between single and multiple producers that
single gates hide.

What the single-gate models found (2026-10-13 review, logit per year):

- Full model: no week slope more than about 1.7 SE from zero for any of the 8 responses.
- Within patient-day (15 patient-days with several visits acquired together): CD4 IFNg 5.9 +/- 1.3,
  CD4 TNFa 5.0 +/- 1.2, CD4 IL-2 3.6 +/- 0.8; CD8 0.4 to 1.0 +/- 1.1 to 1.6. This within-day signal is
  confounded with storage time (section 6).
- Lab-rule positivity (>= 0.05 % background-subtracted and >= 3-fold over neg), stimulated tubes after
  week 0: CD4 IFNg 47 %, IL-2 29 %, TNFa 27 %, CD107a 8 %; CD8 IFNg 15 %, IL-2 13 %, TNFa 14 %, CD107a 6 %.

## 3. Data

- 327 FCS files after exclusion: 272 stimulated, 55 neg. 20 patients.
- Gating (FlowJo, by hand per workspace): `/Lymphocytes/Single Cells/Live cells/{CD4,CD8}/{CD107a, IFNy, IL-2, TNF}`.
  All four marker gates exist under both parents in every tube.
- Peptide pools are patient-specific: 135 patient-pool combinations, 52 used once, 83 tested at more
  than one visit. All pools are treated as one stimulation, with pool identity as a random effect.
- Neg matching: 268 stimulated tubes have a neg from the same patient and visit; 4 use a fallback
  (3 same patient and day, 1 same patient on the reference instrument).
- Covariates: instrument, operator, acquisition day, prior peptide vaccination (`prior_vacc`, 2 levels).
- Visit V0 is a pre-vaccination baseline.

Stimulated tubes by visit (visit number = vaccinations received):

| Visit | Weeks (min, median, max) | Patients | Tubes prior_vacc 1 / 2 |
|---|---|---|---|
| V0 | 0, 0, 0 | 4 | 8 / 0 |
| V1 | 5, 6, 8 | 12 | 51 / 9 |
| V2 | 11, 15, 188 | 15 | 87 / 24 |
| V3 | 32, 148, 315 | 7 | 62 / 0 |
| V4 | 70, 228, 279 | 3 | 24 / 7 |

Only 7 patients have three or more visits, so late visits are mostly different patients from early ones.

Parent sizes, smaller of stimulated and matched neg count (stimulated tubes):

| Parent | 1-49 | 50-99 | 100-249 | 250-999 | >= 1,000 |
|---|---|---|---|---|---|
| CD4 | 2 | 8 | 2 | 22 | 238 |
| CD8 | 2 | 3 | 27 | 75 | 165 |

No tube has an empty parent. Small parents cluster in a few patients (patient 23 most visits; patients
85, 17, 16 several visits in CD8).

## 4. Decisions taken in the review

| Topic | Decision |
|---|---|
| Unit of analysis | One stimulated tube (patient x visit x pool), paired with its neg |
| Gates | FlowJo hand gates (what the lab reports); a fixed-threshold rerun is not part of this version |
| Markers | IFNg, IL-2, TNFa, CD107a: 16 Boolean subsets per parent |
| Small parents | No tubes dropped. Tubes < 1,000 parent events flagged; sensitivity rerun without them |
| Time axis | Visit (V0-V4, categorical) and weeks since first vaccination (spline), reported side by side |
| Direction | Two-sided: subsets may increase or decrease with stimulation |
| CD4 vs CD8 functionality score | Side table only |
| Storage dates | Not available yet; analysis built without them, week results flagged as confounded |
| ELISpot | Data exist; to be supplied later. Config hook for an agreement table, off by default |

## 5. Methods

### 5.1 Boolean counts

For every tube and parent, per-cell membership of the four marker gates is read from the FlowJo
workspace (`flowWorkspace::gh_pop_get_indices`), encoded as a 4-bit code, and tabulated into 16 counts.
Hard checks: the 16 counts sum to the parent count, and each marker's marginal equals the single-gate
count already in the pipeline.

### 5.2 COMPASS (per-tube response map)

COMPASS (Lin et al. 2015, Bioconductor 1.48.0), `SimpleCOMPASS` on count matrices, one fit per parent,
one individual per stimulated tube, `n_u` from the matched neg. Sparse subsets are merged into the
all-negative category with the package's default filter (subset kept if > 5 cells in > 2 tubes), and
the all-negative category is the last column. Three seeds; R-hat on the hyperparameters and the
across-seed correlation of the functionality score (target > 0.99).

Outputs per tube: posterior response probability per subset, functionality and polyfunctionality score
(PFS) with 95 % intervals from the posterior draws. Heatmaps ordered by patient and visit.

COMPASS only detects increases over the neg and pools information across tubes: a tube with few events
gets response probabilities close to the cohort average, with wide intervals. It is therefore used as a
descriptive map and cross-check, not as the test of the hypotheses.

### 5.3 Two-sided subset model (primary)

Stimulated and neg counts of each non-null subset, CD4 and CD8 together, in one binomial model
(brms). Each neg tube enters once; stimulated and neg tubes of one visit share the visit random
effect, so the stimulation effect is estimated within sample:

```
count | trials(parent_count) ~ stim * time * subset      # subset terms partially pooled
                             + instrument + operator + prior_vacc
                             + (1 | acq_day) + (1 | patient) + (1 | patient:visit)
                             + (0 + stim | cocktail) + (1 | obs)
```

- `subset` is the parent x Boolean subset (up to 30 levels); stimulation and stimulation x time effects
  vary by subset with partial pooling, which also addresses multiple comparisons
  (Gelman, Hill & Yajima 2012).
- `time` is either visit (V0-V4) or a weeks spline; two fits.
- `prior_vacc` enters as a main effect and with `stim`; its interaction with visit is not estimable
  (no prior_vacc 2 tubes at V0 or V3).
- `(1 | obs)` absorbs overdispersion.
- About 10,000 rows; runtime may be hours. An lme4 version of the same model is a quick cross-check.

Hypothesis contrasts are computed per posterior draw:

- H1: stimulation effect at V3/V4 per parent; on the weeks axis, at the median week of V3/V4 tubes.
- H2: stimulation-induced cells per subset, `p(stim) - p(neg)` on the probability scale, aggregated by
  number of positive markers (1-4); share with >= 2 markers among induced cells, late vs V1.
- H3: expected number of subsets with a positive stimulation effect, late vs V1.

### 5.4 Checks and sensitivity analyses

- V0 baseline tubes as an over-calling check.
- Rerun without tubes < 1,000 parent events.
- Rerun without the 4 fallback-neg tubes.
- Agreement (kappa) between COMPASS "any IFNg-containing subset responds" and the existing MIMOSA
  IFNg calls.
- Synthetic data (existing simulator): one polyfunctional CD4 subset rising over time with CD8 flat
  must be recovered; a null with day drift only must give contrasts covering zero.

## 6. Limitations

1. **Storage time.** Several visits of one patient were often thawed and acquired on the same day;
   within such a day the later visit was stored for less time, so weeks and storage time are collinear.
   Between days, weeks mix with day-to-day drift. Collection or freezing dates are requested.
2. **Few late patients.** V3 has 7 patients and V4 has 3; between-patient differences and visit effects
   are partly confounded.
3. **Hand gating.** Boolean subsets inherit the hand-set gates, which were placed per visit by one person.
4. **Killing.** The panel has no perforin or granzyme; cytotoxic potential can only be read as
   CD107a degranulation, which is rarely positive (6-8 % of tubes by the lab rule).
5. **Dependence.** Tubes of one visit share their neg; COMPASS treats tubes as independent, the
   subset model carries visit, patient and pool random effects.

## 7. Where the data and code are, and how to run it

There are two machines. Both hold the same raw data; the code and the latest results differ.

### Dev box (de.NBI VM, where the plan was written)

| What | Path |
|---|---|
| Pipeline checkout | `/mnt/volume/cytometry_pipeline` (branch `feat/pbiapplauncher-entry`, with uncommitted work incl. steps 11-15) |
| Raw FCS (328 files) | `/mnt/volume/data/MTB_FACS/fcs/` |
| FlowJo workspaces (20 `.wsp`) | `/mnt/volume/data/MTB_FACS/wsp/` |
| File overview | `/mnt/volume/data/MTB_FACS/fcs_overview.csv` (filename, patient, batch = prior_vacc, vaccine = visit, week, cocktail) |
| Latest full results (2026-10-13 review) | `/mnt/volume/data/MTB_FACS/results_mtb_drift_2026-10-13_review/` |
| Key inputs for step 16 | `gated_frequencies.csv` (single-gate counts, neg pairing), `tubes.csv`, `calls_compare.csv` (step 13 MIMOSA calls), `min_parent_sensitivity.csv` |
| Cohort config | `config/drift_mtb.yaml` (in the checkout) |
| R environment | micromamba env `mtb_drift`, R 4.5, Bioconductor 3.22: flowWorkspace, CytoML, brms + cmdstanr, lme4, MIMOSA. COMPASS not yet installed |
| Scratch for flowWorkspace | `/mnt/volume/tmp_drift` (root disk has under 3 GB free) |

Running R there:

```bash
export MAMBA_ROOT_PREFIX=/mnt/volume/micromamba/root
export TMPDIR=/mnt/volume/tmp_drift
/mnt/volume/micromamba/bin/micromamba run -n mtb_drift Rscript scripts/drift/drift.R status --config config/drift_mtb.yaml
```

### highmem1 (large-memory server for long runs)

| What | Value |
|---|---|
| Connection | SSH from the dev box with `ssh peptimine_highmem1` (user `ubuntu`, key-based; host entry in `~/.ssh/config` on the dev box) |
| Resources | 28 cores, 503 GB RAM, `/mnt/volume` 4.9 TB (about 414 GB free) |
| Pipeline checkout | `/mnt/volume/apps/cytometry_pipeline` (branch `feat/pbiapplauncher-entry` at commit `db19418`, one commit behind the dev box and without its uncommitted work) |
| Raw FCS and workspaces | `/mnt/volume/patrick/UMAP/fcs/` (328) and `/mnt/volume/patrick/UMAP/wsp/` (20); same files as on the dev box |
| Drift results | `results/mtb_drift/` inside the checkout (an earlier run; the 2026-10-13 review results exist only on the dev box) |
| R environment | conda env `cytometry_env` at `/mnt/volume/bin/miniconda/envs/cytometry_env`, R 4.3.3: brms, rstan, cmdstanr, flowWorkspace, CytoML, lme4. COMPASS and MIMOSA not installed (R 4.3 needs the Bioconductor 3.18 build of COMPASS) |

Running a long step there, detached so it survives a closed laptop:

```bash
ssh peptimine_highmem1
cd /mnt/volume/apps/cytometry_pipeline
conda activate cytometry_env
setsid nohup Rscript scripts/drift/16_compass.R --config config/drift_mtb.yaml \
  --fcs_dir /mnt/volume/patrick/UMAP/fcs --wsp_dir /mnt/volume/patrick/UMAP/wsp \
  --out results/mtb_drift > results/mtb_drift/compass.log 2>&1 &
```

`16_compass.R` does not exist yet. To run it on highmem1, the new code must first reach that checkout (commit and push on the dev box, then `git pull` on highmem1), Data paths are not in the config: every step takes `--fcs_dir`, `--wsp_dir` and `--out` on the command line (as `scripts/drift/run_all.sh` does), so the same config works on both machines. Agents working on this project only read from highmem1; the user runs commands that change it.

## 8. Points where a review would help most

1. Is the two-sided binomial subset model (5.3) a sound primary analysis, or is a multinomial /
   Dirichlet-multinomial model over all 16 subsets preferable, given that subsets share one parent?
2. Is the H2 definition (share of induced cells with >= 2 markers, from `p(stim) - p(neg)`)
   appropriate when some subsets decrease with stimulation?
3. With 7 and 3 patients at V3 and V4, is visit as a categorical term too ambitious, and would pooling
   V3/V4 into one "late" level be better?
4. Any better use of the ELISpot data once available (validation of calls versus covariate)?

## References

- Lin L, Finak G, Ushey K, et al. COMPASS identifies T-cell subsets correlated with clinical outcomes.
  Nat Biotechnol 2015;33:610-616.
- Finak G, McDavid A, Chattopadhyay P, et al. Mixture models for single-cell assays with applications to
  vaccine studies. Biostatistics 2014;15:87-101.
- Gelman A, Stern H. The difference between "significant" and "not significant" is not itself
  statistically significant. Am Stat 2006;60:328-331.
- Gelman A, Hill J, Yajima M. Why we (usually) don't have to worry about multiple comparisons.
  J Res Educ Eff 2012;5:189-211.
