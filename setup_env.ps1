# =============================================================================
# setup_env.ps1 — Windows PowerShell setup script for cytometry environment
# =============================================================================
# Usage:  Open PowerShell, cd to project root, then run:
#           .\setup_env.ps1
#
# Prerequisites:
#   1. Miniconda/Anaconda installed and in PATH
#   2. Rtools 4.3+ installed from https://cran.r-project.org/bin/windows/Rtools/
# =============================================================================

$ENV_NAME = "cytometry_env"
$YML_FILE = "environment_windows.yml"
$EXTRA_R_SCRIPT = "install_extra_packages.R"

Write-Host "===================================================" -ForegroundColor Cyan
Write-Host "Cytometry Pipeline Environment Setup (Windows)"     -ForegroundColor Cyan
Write-Host "===================================================" -ForegroundColor Cyan

# -------------------------------------------------------------------
# 1. Configure R library path
# -------------------------------------------------------------------
Write-Host ""
Write-Host "--- R Library Configuration ---" -ForegroundColor Yellow
Write-Host ""
Write-Host "R needs a single folder where ALL packages are installed."
Write-Host "This avoids packages scattering across multiple locations."
Write-Host ""

# Check if R_LIBS_USER is already set persistently
$existingLib = [Environment]::GetEnvironmentVariable("R_LIBS_USER", "User")
if ($existingLib) {
    Write-Host "Current R_LIBS_USER: $existingLib" -ForegroundColor Green
    $changeLib = Read-Host "Keep this library path? (y/n)"
    if ($changeLib -match "^[Yy]") {
        $R_LIB_PATH = $existingLib
    } else {
        $R_LIB_PATH = $null
    }
} else {
    $R_LIB_PATH = $null
}

if (-not $R_LIB_PATH) {
    Write-Host ""
    Write-Host "Enter the FULL path to your R library folder."
    Write-Host "Examples:"
    Write-Host "  D:\R_libs"
    Write-Host "  C:\Users\$env:USERNAME\R\library"
    Write-Host "  C:\R\win-library\4.3"
    Write-Host ""
    $R_LIB_PATH = Read-Host "R library path"

    # Trim quotes if the user pasted a quoted path
    $R_LIB_PATH = $R_LIB_PATH.Trim('"').Trim("'").Trim()

    if (-not $R_LIB_PATH) {
        Write-Host "Error: No path provided. Exiting." -ForegroundColor Red
        exit 1
    }
}

# Create the directory if it doesn't exist
if (-not (Test-Path $R_LIB_PATH)) {
    Write-Host "Directory does not exist. Creating: $R_LIB_PATH"
    New-Item -ItemType Directory -Path $R_LIB_PATH -Force | Out-Null
}

# Normalize to full path
$R_LIB_PATH = (Resolve-Path $R_LIB_PATH).Path
Write-Host "[OK] R library path: $R_LIB_PATH" -ForegroundColor Green

# --- Persist R_LIBS_USER as a permanent User environment variable ---
[Environment]::SetEnvironmentVariable("R_LIBS_USER", $R_LIB_PATH, "User")
$env:R_LIBS_USER = $R_LIB_PATH
Write-Host "[OK] R_LIBS_USER set permanently for user $env:USERNAME" -ForegroundColor Green

# --- Write ~/.Renviron so R always picks it up (even outside conda) ---
$renvironPath = Join-Path $env:USERPROFILE ".Renviron"
# Use forward slashes in .Renviron (R convention on all platforms)
$R_LIB_PATH_FWDSLASH = $R_LIB_PATH -replace '\\', '/'
$renvironLine = "R_LIBS_USER=$R_LIB_PATH_FWDSLASH"

if (Test-Path $renvironPath) {
    # Read existing content, replace or append R_LIBS_USER line
    $content = Get-Content $renvironPath -Raw
    if ($content -match "(?m)^R_LIBS_USER=") {
        $content = $content -replace "(?m)^R_LIBS_USER=.*$", $renvironLine
    } else {
        $content = $content.TrimEnd() + "`n$renvironLine`n"
    }
    Set-Content -Path $renvironPath -Value $content -NoNewline
} else {
    Set-Content -Path $renvironPath -Value "$renvironLine`n"
}
Write-Host "[OK] Written to $renvironPath" -ForegroundColor Green
Write-Host ""

# -------------------------------------------------------------------
# 2. Check prerequisites
# -------------------------------------------------------------------
Write-Host "--- Checking Prerequisites ---" -ForegroundColor Yellow

# Check conda
$condaCmd = Get-Command conda -ErrorAction SilentlyContinue
if (-not $condaCmd) {
    Write-Host "Error: conda not found in PATH." -ForegroundColor Red
    Write-Host "Install Miniconda from https://docs.conda.io/en/latest/miniconda.html"
    exit 1
}
Write-Host "[OK] conda found: $($condaCmd.Source)" -ForegroundColor Green

# Check Rtools
$rtoolsPaths = @(
    "C:\rtools43",
    "C:\rtools44",
    "C:\rtools45",
    "$env:LOCALAPPDATA\Programs\R\rtools43",
    "$env:LOCALAPPDATA\Programs\R\rtools44"
)
$rtoolsFound = $false
foreach ($p in $rtoolsPaths) {
    if (Test-Path $p) {
        Write-Host "[OK] Rtools found: $p" -ForegroundColor Green
        $rtoolsFound = $true
        break
    }
}
if (-not $rtoolsFound) {
    Write-Host "[WARN] Rtools not found in standard locations." -ForegroundColor Yellow
    Write-Host "  Some R packages may fail to compile from source."
    Write-Host "  Download from: https://cran.r-project.org/bin/windows/Rtools/"
}

# Check Git
$gitCmd = Get-Command git -ErrorAction SilentlyContinue
if ($gitCmd) {
    Write-Host "[OK] git found: $($gitCmd.Source)" -ForegroundColor Green
} else {
    Write-Host "[WARN] git not found. GitHub package installs will fail." -ForegroundColor Yellow
    Write-Host "  Download from: https://git-scm.com/download/win"
}

# -------------------------------------------------------------------
# 3. Create or update conda environment
# -------------------------------------------------------------------
Write-Host ""
Write-Host "--- Conda Environment ---" -ForegroundColor Yellow

$envList = conda info --envs 2>&1 | Out-String
if ($envList -match $ENV_NAME) {
    Write-Host "Environment '$ENV_NAME' already exists." -ForegroundColor Yellow
    $reply = Read-Host "Update it? (y/n)"
    if ($reply -match "^[Yy]") {
        Write-Host "Updating environment..."
        conda env update -f $YML_FILE --prune
    }
} else {
    Write-Host "Creating environment '$ENV_NAME'..."
    conda env create -f $YML_FILE
}

# -------------------------------------------------------------------
# 4. Activate and install extra R packages
# -------------------------------------------------------------------
Write-Host ""
Write-Host "===================================================" -ForegroundColor Cyan
Write-Host "Installing Extra R Packages (GitHub/Bioconductor)"  -ForegroundColor Cyan
Write-Host "===================================================" -ForegroundColor Cyan

# Initialize conda for this shell session
$condaHook = conda shell.powershell hook 2>&1 | Out-String
Invoke-Expression $condaHook
conda activate $ENV_NAME

if ($LASTEXITCODE -eq 0) {
    Write-Host "[OK] Environment activated" -ForegroundColor Green

    # Verify R sees the library path
    Write-Host ""
    Write-Host "Verifying R library path..."
    Rscript -e "cat('R sees library:', .libPaths()[1], '\n'); cat('R_LIBS_USER:', Sys.getenv('R_LIBS_USER'), '\n')"

    Write-Host ""
    Write-Host "Running R installer script..."
    Rscript $EXTRA_R_SCRIPT
} else {
    Write-Host "[FAIL] Could not activate environment." -ForegroundColor Red
    Write-Host "Try manually: conda activate $ENV_NAME && Rscript $EXTRA_R_SCRIPT"
    exit 1
}

# -------------------------------------------------------------------
# 5. Try installing problematic Bioconductor packages
# -------------------------------------------------------------------
Write-Host ""
Write-Host "Attempting Bioconductor packages that often fail on Windows..."
$biocScript = @"
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
pkgs <- c('flowWorkspace', 'openCyto', 'CytoML', 'cytolib')
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat(sprintf('  Installing %s...\n', p))
    tryCatch({
      BiocManager::install(p, update = FALSE, ask = FALSE, type = 'binary')
      cat(sprintf('  [OK] %s installed.\n', p))
    }, error = function(e) {
      cat(sprintf('  [SKIP] %s failed: %s\n', p, e[['message']]))
      cat(sprintf('         This is expected on Windows. FlowJo import (Pipeline D) will be unavailable.\n'))
    })
  } else {
    cat(sprintf('  [OK] %s already installed.\n', p))
  }
}
"@
Rscript -e $biocScript

# -------------------------------------------------------------------
# 6. Summary
# -------------------------------------------------------------------
Write-Host ""
Write-Host "===================================================" -ForegroundColor Cyan
Write-Host "Setup Complete!"                                     -ForegroundColor Green
Write-Host "===================================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Configuration summary:" -ForegroundColor White
Write-Host "  R library path : $R_LIB_PATH"
Write-Host "  Conda env      : $ENV_NAME"
Write-Host "  .Renviron      : $renvironPath"
Write-Host ""
Write-Host "The R library path has been set in TWO places:" -ForegroundColor White
Write-Host "  1. Windows env variable R_LIBS_USER (persistent, all R sessions)"
Write-Host "  2. ~/.Renviron file (picked up by R/RStudio on startup)"
Write-Host ""
Write-Host "To use the environment:" -ForegroundColor White
Write-Host "  conda activate $ENV_NAME"
Write-Host ""
Write-Host "To run a pipeline:" -ForegroundColor White
Write-Host "  cd pipelines"
Write-Host '  Rscript pipeline_A_raw_fcs.R --config ..\config\pipeline_A_config.R'
Write-Host ""
Write-Host "To verify your library path from R:" -ForegroundColor White
Write-Host '  Rscript -e ".libPaths()"'
Write-Host ""
