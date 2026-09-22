#!/usr/bin/env bash
# Run the whole drift analysis. Steps 1-6 and 9-10 give the complete yardstick answer; 7-8 are the
# BayesFlow experiment and run last. Intended for a screen session on highmem1:
#   screen -S mtb_drift
#   bash scripts/mtb_drift/run_all.sh /mnt/volume/patrick/UMAP results/mtb_drift [r_conda_env] [bayesflow_conda_env]
# DATA_ROOT must contain fcs/ and wsp/ (or set FCS_DIR and WSP_DIR). Set SKIP_BF=1 to skip 7-8.
# Without screen: nohup bash scripts/mtb_drift/run_all.sh ... > /dev/null 2>&1 &   (progress in $OUT/run_all.log)
set -euo pipefail
DATA_ROOT="${1:?data root}"; OUT="${2:?results dir}"; ENV="${3:-cytometry_env}"; BFENV="${4:-bayesflow_env}"
FCS_DIR="${FCS_DIR:-$DATA_ROOT/fcs}"; WSP_DIR="${WSP_DIR:-$DATA_ROOT/wsp}"
HERE="$(cd "$(dirname "$0")" && pwd)"; R="conda run --no-capture-output -n $ENV Rscript"; PY="conda run --no-capture-output -n $BFENV python"
mkdir -p "$OUT"; exec > >(tee -a "$OUT/run_all.log") 2>&1
step() { echo; echo "===== $(date '+%F %T') $* ====="; }

step 1 metadata;        $R "$HERE/01_extract_metadata.R" --fcs_dir "$FCS_DIR" --out "$OUT"
step 2 panel;           $R "$HERE/02_panel_reconcile.R" --out "$OUT"
step 3 tube QC;         $R "$HERE/03_tube_qc.R" --out "$OUT" --workers "${WORKERS:-8}"
step 4 neg summaries;   $R "$HERE/04_neg_summaries.R" --wsp_dir "$WSP_DIR" --fcs_dir "$FCS_DIR" --out "$OUT" --all_tubes
step 5 frequencies;     $R "$HERE/05_gated_frequencies.R" --wsp_dir "$WSP_DIR" --fcs_dir "$FCS_DIR" --out "$OUT"
step 6a brms;           $R "$HERE/06_yardstick_brms.R" --out "$OUT" --cores "${CORES:-4}"
step 6b stan;           $R "$HERE/06_yardstick_stan.R" --out "$OUT" --cores "${CORES:-4}"
step 6c lme4;           $R "$HERE/06_lme4.R" --out "$OUT"
step 9 decision;        $R "$HERE/09_decision.R" --out "$OUT"
step 10 report;         $R -e "rmarkdown::render('$HERE/10_report.Rmd', params = list(out = normalizePath('$OUT')), output_file = 'report.html', output_dir = normalizePath('$OUT'))"
if [ -z "${SKIP_BF:-}" ]; then
  # Budgets are the ones that ran on the dev box (v1: 2 h 20 min; v2 training: 3 h 15 min, both on 28 CPU cores
  # while sharing them with other jobs). Posterior sampling is chunked to keep v2 under a few GB of RAM.
  step 7 bayesflow v1;  $PY "$HERE/07_bayesflow_v1.py" --out "$OUT" --epochs "${BF_EPOCHS:-30}" --num_sims "${BF_SIMS:-30000}" --batch_size 128
  step 8 bayesflow v2;  $PY "$HERE/08_bayesflow_v2.py" --out "$OUT" --epochs "${BF2_EPOCHS:-15}" --batches_per_epoch 100 --batch_size 8 \
                          --n_cells "${BF_CELLS:-800}" --chunk 4 --n_sbc 100 --posterior_samples 1000
  step 9 decision again, with BayesFlow; $R "$HERE/09_decision.R" --out "$OUT"
  step 10 report again; $R -e "rmarkdown::render('$HERE/10_report.Rmd', params = list(out = normalizePath('$OUT')), output_file = 'report.html', output_dir = normalizePath('$OUT'))"
fi
step done
