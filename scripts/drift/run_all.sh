#!/usr/bin/env bash
# Run the whole drift analysis. Thin wrapper around `Rscript scripts/drift/drift.R run`, kept so the
# old positional call still works:
#   bash scripts/drift/run_all.sh DATA_ROOT OUT CFG [r_conda_env] [bayesflow_conda_env]
# maps to
#   conda run -n r_conda_env Rscript scripts/drift/drift.R run --fcs_dir DATA_ROOT/fcs --wsp_dir DATA_ROOT/wsp \
#     --out OUT --config CFG --cores $CORES --workers $WORKERS --bf_env bayesflow_conda_env [--with_bayesflow] [--sample_sheet $SAMPLE_SHEET]
# DATA_ROOT must contain fcs/ and wsp/ (or set FCS_DIR and WSP_DIR). Set SAMPLE_SHEET to a CSV when the
# file names carry no metadata. Steps 7-8 (BayesFlow) run unless SKIP_BF=1, as before. drift.R skips
# steps whose outputs are up to date, logs to $OUT/run_all.log and writes $OUT/.done/<step>.json.
# Intended for a screen session on highmem1:
#   screen -S mtb_drift
#   bash scripts/drift/run_all.sh /mnt/volume/patrick/UMAP results/mtb_drift config/drift_mtb.yaml [r_conda_env] [bayesflow_conda_env]
# Without screen: nohup bash scripts/drift/run_all.sh ... > /dev/null 2>&1 &   (progress in $OUT/run_all.log)
set -euo pipefail
DATA_ROOT="${1:?data root}"; OUT="${2:?results dir}"; CFG="${3:?cohort config yaml}"; ENV="${4:-cytometry_env}"; BFENV="${5:-bayesflow_env}"
FCS_DIR="${FCS_DIR:-$DATA_ROOT/fcs}"; WSP_DIR="${WSP_DIR:-$DATA_ROOT/wsp}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ARGS=(run --fcs_dir "$FCS_DIR" --wsp_dir "$WSP_DIR" --out "$OUT" --config "$CFG"
      --cores "${CORES:-4}" --workers "${WORKERS:-8}" --bf_env "$BFENV")
[ -n "${SAMPLE_SHEET:-}" ] && ARGS+=(--sample_sheet "$SAMPLE_SHEET")
[ -z "${SKIP_BF:-}" ] && ARGS+=(--with_bayesflow)
exec conda run --no-capture-output -n "$ENV" Rscript "$HERE/drift.R" "${ARGS[@]}"
