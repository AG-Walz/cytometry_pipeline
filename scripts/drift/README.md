# scripts/drift: acquisition drift and instrument shift

Implements `docs/plans/mtb_batch_drift_plan.md`. Thin CLIs here, logic in `src/harmonization/drift/`.
`drift.R` is the entry point; the numbered step scripts are what it runs and can still be called on
their own. Everything cohort-specific lives in one YAML passed to every step with `--config`;
`config/drift_mtb.yaml` is the MTB ICS cohort and `config/drift_template.yaml` the starting point that
`drift.R init` copies for a new cohort. The walkthrough for a new cohort, from a folder of FCS and
.wsp files to the report, is `docs/scripts/drift.md`.

Quick start:

```bash
conda env create -f environment_drift.yml                       # R side; or 00_setup_env.sh for both envs
Rscript scripts/drift/drift.R doctor                            # packages, CmdStan, optional --bf_env
Rscript scripts/drift/drift.R init --fcs_dir <fcs> --out results/cohort --cohort NAME   # sample_sheet.csv + config
#   complete the sheet and the EDIT fields of the config, harmonise marker names (below)
Rscript scripts/drift/drift.R run --sample_sheet results/cohort/sample_sheet.csv --fcs_dir <fcs> --out results/cohort --config results/cohort/drift_NAME.yaml --only 1
Rscript scripts/drift/drift.R check --out results/cohort --config results/cohort/drift_NAME.yaml     # verdict: separable?
screen -S drift
Rscript scripts/drift/drift.R run --sample_sheet results/cohort/sample_sheet.csv --fcs_dir <fcs> --wsp_dir <wsp> \
  --out results/cohort --config results/cohort/drift_NAME.yaml [--with_bayesflow --bf_env bayesflow_env]
Rscript scripts/drift/drift.R status --out results/cohort --config results/cohort/drift_NAME.yaml
```

`run` executes, in order, 1 metadata, 1b check_design, 2 panel, 3 tube QC, 4 neg summaries, 5 frequencies,
6a brms, 6b Stan, 6c lme4, 6d biology (week slope with and without technical terms, leave-one-day-out),
7 and 8 BayesFlow (opt-in), 9 decision, 10 report. A step runs when one of its outputs is missing or
older than one of its inputs (the config counts as an input); `--only 6c,9,10`, `--from 6`, `--force`
override. Output is tee'd to `<out>/run_all.log` and every step leaves `<out>/.done/<step>.json` (time,
duration, git commit, config checksum, command line). `bash scripts/drift/run_all.sh DATA_ROOT OUT CFG
[r_env] [bf_env]` is the old positional form and calls `drift.R run` with the same meaning of
`FCS_DIR`, `WSP_DIR`, `SAMPLE_SHEET`, `SKIP_BF`, `CORES`, `WORKERS`. Outputs and timings per step are in
section 9 of the plan and in `docs/scripts/drift.md`.

Before the first full run, harmonise the marker names in the FCS files so the `markers` block only
needs the canonical names: `Rscript scripts/harmonize_marker_names.R --fcs_dir <fcs> --probe` shows every
`$PnS` spelling per channel, a run with `--canonical config/drift_mtb.yaml` previews the renames, and
`--apply --in_place` writes them with an undo ledger (see `docs/scripts/harmonize-marker-names.md`). Run
it before step 1; step 2 then sees one spelling per marker.

Before the first full run, list the FlowJo gate tree and set `neg_population`, `parent_populations`
and `cytokine_gate_regex` in the config to match it:

```bash
Rscript scripts/drift/drift.R run --sample_sheet results/cohort/sample_sheet.csv --fcs_dir <fcs> --out results/cohort --config results/cohort/drift_NAME.yaml --only 1,2
Rscript scripts/drift/04_neg_summaries.R --wsp_dir <wsp> --fcs_dir <fcs> --out results/cohort --config results/cohort/drift_NAME.yaml --list_populations
```

For the MTB cohort the file names carry the metadata (`tubes.filename_regex` in `config/drift_mtb.yaml`),
so step 1 also runs without `--sample_sheet`; `drift.R init --config config/drift_mtb.yaml` turns that
regex into a pre-filled sheet, which is the migration path.

Tests without cytometry packages: `Rscript test_data/test_drift.R`, `Rscript test_data/test_drift_controls.R`
(control matching and background flips) and `Rscript test_data/test_drift_design.R`. Synthetic end-to-end
check of steps 6 to 10: run steps 1 and 2 (`drift.R run ... --only 1,2`), then
`Rscript test_data/drift_synthetic.R --out results/syn --config config/drift_mtb.yaml`, then
`drift.R run --out results/syn --config config/drift_mtb.yaml --from 6`; recovered effects should match
`synthetic_truth.csv`, which lists the injected technical (A), background (BG) and endpoint (B) effects per
`layer` and `response` (NA response = shared day drift). Both BayesFlow scripts accept `--smoke` for a
one-minute end-to-end check.
