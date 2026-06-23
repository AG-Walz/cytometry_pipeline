#!/bin/bash
# setup_env.sh
# Automates the creation/update of the cytometry conda environment.
#
# NOTE on R library isolation:
#   A global ~/.Rprofile or ~/.Renviron that hard-codes .libPaths()/R_LIBS_USER
#   to a specific R version's library (e.g. ~/R/x86_64-pc-linux-gnu-library/4.4)
#   will shadow the conda environment's own R packages. Loading packages built
#   for one R version into another causes "undefined symbol" / ABI errors
#   (e.g. SETLENGTH, NO_REFERENCES) and can break `conda env create` at a
#   Bioconductor post-link step. We therefore neutralise those global files for
#   every R invocation here (R_PROFILE_USER/R_ENVIRON_USER=/dev/null) and install
#   the same isolation as a conda activation hook, so the env stays self-contained.

ENV_NAME="cytometry_env"
YML_FILE="environment.yml"
EXTRA_R_SCRIPT="install_extra_packages.R"

echo "==================================================="
echo "Cytometry Pipeline Environment Setup"
echo "==================================================="

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Error: conda could not be found."
    echo "Please install Miniconda or Anaconda first."
    exit 1
fi

# --- Isolate R from any global ~/.Rprofile / ~/.Renviron during setup ---
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null
unset R_LIBS_USER R_LIBS_SITE 2>/dev/null || true

CONDA_BASE="$(conda info --base)"
ENV_PREFIX="$CONDA_BASE/envs/$ENV_NAME"

create_env() {
    echo "Creating environment '$ENV_NAME'..."
    conda env create -f "$YML_FILE"
}

# Decide create / recreate / update.
if conda env list | grep -qE "[/ ]$ENV_NAME$"; then
    if [ ! -x "$ENV_PREFIX/bin/R" ]; then
        echo "Environment '$ENV_NAME' exists but is incomplete (no R found)."
        echo "Removing the broken environment and recreating it..."
        conda env remove -n "$ENV_NAME" -y
        create_env
    else
        echo "Environment '$ENV_NAME' already exists."
        read -p "Do you want to update it? (y/n) " -n 1 -r
        echo ""
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            echo "Updating environment..."
            conda env update -f "$YML_FILE" --prune
        fi
    fi
else
    create_env
fi

# Fail fast: if R is still missing, the env did not build correctly.
if [ ! -x "$ENV_PREFIX/bin/R" ]; then
    echo "==================================================="
    echo "✗ Environment creation failed (no R in $ENV_PREFIX)."
    echo "  Scroll up for the conda error. Common cause: a global"
    echo "  ~/.Rprofile or ~/.Renviron pinning .libPaths()/R_LIBS_USER."
    echo "==================================================="
    exit 1
fi

# --- Install conda activation hooks so the env is self-contained at runtime ---
echo "Installing R-library isolation activation hook..."
ACT_DIR="$ENV_PREFIX/etc/conda/activate.d"
DEACT_DIR="$ENV_PREFIX/etc/conda/deactivate.d"
mkdir -p "$ACT_DIR" "$DEACT_DIR"

cat > "$ACT_DIR/zz-r-libs-isolation.sh" <<'EOF'
# Keep conda R self-contained: ignore any global ~/.Rprofile / ~/.Renviron that
# pin a different R version's package library (avoids ABI / undefined-symbol
# errors). Set by setup_env.sh.
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null
export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
EOF

cat > "$DEACT_DIR/zz-r-libs-isolation.sh" <<'EOF'
unset R_PROFILE_USER R_ENVIRON_USER R_LIBS_USER
EOF

echo "==================================================="
echo "Installing Extra R Packages (GitHub/Non-Conda)"
echo "==================================================="

# Activate the environment (activation hook applies the isolation above) and
# run the R installer.
eval "$(conda shell.bash hook)"
conda activate "$ENV_NAME"

if [ $? -eq 0 ]; then
    echo "✓ Environment activated"
    echo "Running R installer script..."
    Rscript "$EXTRA_R_SCRIPT"
else
    echo "✗ Failed to activate environment."
    exit 1
fi

echo "==================================================="
echo "Setup Complete!"
echo "To use the environment, run:"
echo "  conda activate $ENV_NAME"
echo "==================================================="
