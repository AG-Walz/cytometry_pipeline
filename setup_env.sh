#!/bin/bash
# setup_env.sh
# Automates the creation/update of the cytometry conda environment

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

# Check if environment exists
if conda info --envs | grep -q "^$ENV_NAME"; then
    echo "Environment '$ENV_NAME' already exists."
    read -p "Do you want to update it? (y/n) " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Updating environment..."
        conda env update -f $YML_FILE --prune
    fi
else
    echo "Creating environment '$ENV_NAME'..."
    conda env create -f $YML_FILE
fi

echo "==================================================="
echo "Installing Extra R Packages (GitHub/Non-Conda)"
echo "==================================================="

# Activate environment to run the R script
# Note: 'conda activate' might not work in scripts without eval
eval "$(conda shell.bash hook)"
conda activate $ENV_NAME

if [ $? -eq 0 ]; then
    echo "✓ Environment activated"
    echo "Running R installer script..."
    Rscript $EXTRA_R_SCRIPT
else
    echo "✗ Failed to activate environment."
    exit 1
fi

echo "==================================================="
echo "Setup Complete!"
echo "To use the environment, run:"
echo "  conda activate $ENV_NAME"
echo "==================================================="
