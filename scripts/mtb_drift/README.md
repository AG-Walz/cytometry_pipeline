# scripts/mtb_drift — acquisition drift and instrument shift, MTB ICS cohort

Implements `docs/plans/mtb_batch_drift_plan.md`. Thin CLIs here, logic in `src/harmonization/drift/`,
settings in `config/mtb_drift.yaml`. Run order, outputs and timings are in section 9 of the plan.

Quick start on highmem1:

```bash
bash scripts/mtb_drift/00_setup_env.sh cytometry_env bayesflow_env
screen -S mtb_drift
bash scripts/mtb_drift/run_all.sh /mnt/volume/patrick/UMAP results/mtb_drift cytometry_env bayesflow_env
```

Before the first full run, list the FlowJo gate tree and set `neg_population`, `parent_populations`
and `cytokine_gate_regex` in the config to match it:

```bash
Rscript scripts/mtb_drift/01_extract_metadata.R --fcs_dir /mnt/volume/patrick/UMAP/fcs --out results/mtb_drift
Rscript scripts/mtb_drift/02_panel_reconcile.R --out results/mtb_drift
Rscript scripts/mtb_drift/04_neg_summaries.R --wsp_dir /mnt/volume/patrick/UMAP/wsp --fcs_dir /mnt/volume/patrick/UMAP/fcs --out results/mtb_drift --list_populations
```

Tests without cytometry packages: `Rscript test_data/test_mtb_drift.R`. Synthetic end-to-end check of
steps 6 to 9: run steps 1 and 2, then `Rscript test_data/mtb_drift_synthetic.R --out results/syn`, then
the step 6 to 9 scripts on `results/syn`; recovered effects should match `synthetic_truth.csv`.
Both BayesFlow scripts accept `--smoke` for a one-minute end-to-end check.
