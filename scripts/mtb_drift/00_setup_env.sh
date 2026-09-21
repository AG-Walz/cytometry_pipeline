#!/usr/bin/env bash
# Step 0: add what the drift analysis needs on top of the existing cytometry_env conda environment,
# and create a CPU-only Python venv for BayesFlow (kept separate so keras/torch never touch conda R).
#   bash scripts/mtb_drift/00_setup_env.sh [conda_env_name] [bf_venv_dir]
set -euo pipefail
ENV="${1:-cytometry_env}"
BFVENV="${2:-$HOME/bfvenv}"

echo "== R packages into conda env $ENV =="
conda install -y -n "$ENV" -c conda-forge -c bioconda \
  r-brms r-cmdstanr cmdstan r-yaml r-lme4 r-rmarkdown r-knitr r-ggplot2 \
  bioconductor-cytoml bioconductor-flowworkspace bioconductor-flowcore bioconductor-peacoqc bioconductor-flowai
conda run -n "$ENV" Rscript -e 'cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN", cmdstanr::cmdstan_path())); cat("cmdstan:", cmdstanr::cmdstan_version(), "\n")' \
  || echo "cmdstanr cannot see CmdStan; run: conda run -n $ENV Rscript -e 'cmdstanr::install_cmdstan(cores = 8)'"

echo "== BayesFlow venv at $BFVENV (CPU torch) =="
python3 -m venv "$BFVENV"
"$BFVENV/bin/pip" install -q --upgrade pip
"$BFVENV/bin/pip" install -q --index-url https://download.pytorch.org/whl/cpu torch
"$BFVENV/bin/pip" install -q "bayesflow>=2.0.14,<3" pandas pyyaml scipy matplotlib
"$BFVENV/bin/python" -c "import bayesflow as bf, keras; print('bayesflow', bf.__version__, 'keras', keras.__version__, keras.backend.backend())"

echo "== check =="
conda run -n "$ENV" Rscript -e 'for (p in c("CytoML","flowWorkspace","flowCore","PeacoQC","flowAI","brms","cmdstanr","lme4","yaml","rmarkdown")) cat(sprintf("%-14s %s\n", p, requireNamespace(p, quietly = TRUE)))'
