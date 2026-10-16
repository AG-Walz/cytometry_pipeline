#!/usr/bin/env bash
# Step 0: add what the drift analysis needs on top of the existing cytometry_env conda environment,
# and create a separate CPU-only conda env for BayesFlow (so keras/torch never touch the R env).
#   bash scripts/drift/00_setup_env.sh [r_conda_env] [bayesflow_conda_env]
# The R package list is also in environment_drift.yml at the repository root, for a fresh env:
#   conda env create -f environment_drift.yml   (then `Rscript scripts/drift/drift.R doctor` to check it)
set -euo pipefail
ENV="${1:-cytometry_env}"
BFENV="${2:-bayesflow_env}"

echo "== R packages into conda env $ENV =="
conda install -y -n "$ENV" -c conda-forge -c bioconda \
  r-brms r-cmdstanr cmdstan r-yaml r-lme4 r-rmarkdown r-knitr r-ggplot2 \
  bioconductor-cytoml bioconductor-flowworkspace bioconductor-flowcore bioconductor-peacoqc bioconductor-flowai
conda run -n "$ENV" Rscript -e 'cmdstanr::set_cmdstan_path(Sys.getenv("CMDSTAN", cmdstanr::cmdstan_path())); cat("cmdstan:", cmdstanr::cmdstan_version(), "\n")' \
  || echo "cmdstanr cannot see CmdStan; run: conda run -n $ENV Rscript -e 'cmdstanr::install_cmdstan(cores = 8)'"

echo "== BayesFlow conda env $BFENV (Python 3.12, CPU torch) =="
# BayesFlow 2 needs Python 3.12 or 3.13, which the system python on highmem does not provide,
# so the interpreter comes from conda-forge and the packages from pip inside that env.
conda create -y -n "$BFENV" -c conda-forge python=3.12 pip
conda run -n "$BFENV" pip install -q --upgrade pip
conda run -n "$BFENV" pip install -q --index-url https://download.pytorch.org/whl/cpu torch
conda run -n "$BFENV" pip install -q "bayesflow>=2.0.14,<3" pandas pyyaml scipy matplotlib
conda run -n "$BFENV" python -c "import bayesflow as bf, keras; print('bayesflow', bf.__version__, 'keras', keras.__version__, keras.backend.backend())"

echo "== check =="
conda run -n "$ENV" Rscript -e 'for (p in c("CytoML","flowWorkspace","flowCore","PeacoQC","flowAI","brms","cmdstanr","lme4","yaml","rmarkdown")) cat(sprintf("%-14s %s\n", p, requireNamespace(p, quietly = TRUE)))'
