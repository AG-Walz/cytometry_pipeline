#!/usr/bin/env bash
# =============================================================================
# app/run.sh — launch the FlowJo -> UMAP app (PBIAppLauncher entry point)
# =============================================================================
# Usage: bash app/run.sh <port>
#
# The conda environment is ACTIVATED rather than called through its bin/
# directory. That matters here: the app runs `Rscript` as a subprocess, and
# setup_env.sh installs an activation hook that pins R_LIBS_USER to the
# environment's own library and neutralises a global ~/.Rprofile / ~/.Renviron.
# Invoking $CONDA_PREFIX/bin/gunicorn directly skips that hook, leaving the
# subprocess to find the system Rscript — which has no CytoML — or to load
# ABI-incompatible packages from ~/R. See docs/scripts/umap-from-flowjo.md.
#
# Override CONDA_BASE or UMAP_APP_ENV if the install lives elsewhere.
# =============================================================================
set -euo pipefail

PORT="${1:?usage: run.sh <port>}"
CONDA_BASE="${CONDA_BASE:-/mnt/volume/bin/miniconda}"
UMAP_APP_ENV="${UMAP_APP_ENV:-cytometry_env}"

if [ ! -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
    echo "conda not found at $CONDA_BASE — set CONDA_BASE to the install root." >&2
    exit 1
fi

# Activation hooks are not nounset-safe (cmdstan's reads an unset $CMDSTAN), so
# relax -u around them or the script dies before gunicorn starts.
# shellcheck disable=SC1091
. "$CONDA_BASE/etc/profile.d/conda.sh"
set +u
conda activate "$UMAP_APP_ENV"
set -u

# Repo root, so `app.app` is importable and the R script is found relative to it.
cd "$(dirname "$(readlink -f "$0")")/.."

if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript is not on PATH after activating '$UMAP_APP_ENV'." >&2
    exit 1
fi

# --workers 1 is REQUIRED, not a tuning choice: the active run lives in a
# module-level object, so a second worker would own a second, unrelated runner
# and the browser would poll whichever process its request happened to reach.
# Concurrency comes from threads instead.
exec gunicorn \
    --workers 1 \
    --threads 8 \
    --timeout 120 \
    --bind "0.0.0.0:${PORT}" \
    app.app:server
