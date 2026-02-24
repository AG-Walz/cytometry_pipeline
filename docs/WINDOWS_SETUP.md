# Windows Setup Guide — Cytometry Pipeline

A step-by-step guide to install and run the cytometry pipeline on a fresh Windows 10 or 11 computer. No prior programming experience is assumed.

**Time required**: approximately 45-90 minutes (mostly waiting for downloads and installs).

**What you will install** (all free):

| Software | What it is | Why you need it |
|----------|-----------|-----------------|
| Git | Version control tool | To download ("clone") the pipeline code |
| R | Statistical programming language | The pipeline is written in R |
| Rtools | C/C++ compiler for R | Some R packages need to be compiled from source code |
| Miniconda | Package and environment manager | Installs all the R packages the pipeline needs, in an isolated environment |

---

## Part 1: Install the four required programs

You must install these **in this exact order**. Do not skip any step.

---

### Step 1 — Install Git

Git is the tool that lets you download the pipeline code from the online repository.

1. Open your web browser (Edge, Chrome, Firefox — any is fine).
2. Go to: **https://git-scm.com/download/win**
3. The download should start automatically. If it does not, click the link that says **"Click here to download manually"**.
4. Once the file (something like `Git-2.xx.x-64-bit.exe`) has downloaded, **double-click it** to start the installer.
5. The installer will show several screens. Use these settings:
   - **Select Components**: leave everything at the defaults. Click **Next**.
   - **Choosing the default editor**: leave it as "Use Vim" (you won't need to use it). Click **Next**.
   - **Adjusting the name of the initial branch**: leave as "Let Git decide". Click **Next**.
   - **Adjusting your PATH environment**: select **"Git from the command line and also from 3rd-party software"** (this is usually already selected). Click **Next**.
   - For all remaining screens, leave defaults and keep clicking **Next**, then **Install**.
6. When the installer finishes, click **Finish**.

**How to verify it worked:**
- Press the `Windows` key on your keyboard, type `cmd`, and press Enter. A black window opens.
- Type `git --version` and press Enter.
- You should see something like: `git version 2.44.0.windows.1`
- Type `exit` and press Enter to close the window.

---

### Step 2 — Install R

R is the programming language the pipeline is written in.

1. Go to: **https://cran.r-project.org/bin/windows/base/**
2. Click the link at the top that says **"Download R-4.x.x for Windows"** (the version number may be different; any 4.3+ version is fine).
3. Once the file (something like `R-4.4.1-win.exe`) has downloaded, **double-click it**.
4. The installer will show several screens:
   - **Select language**: choose English. Click **OK**.
   - **Information**: Click **Next**.
   - **Select destination location**: leave the default (`C:\Program Files\R\R-4.x.x`). Click **Next**.
   - **Select components**: leave everything checked. Click **Next**.
   - **Startup options**: select **"No (accept defaults)"**. Click **Next**.
   - **Select Start Menu folder**: leave default. Click **Next**.
   - **Additional tasks**:
     - Check **"Save version number in registry"** (should already be checked).
     - Check **"Associate R with .RData files"** (should already be checked).
     - Click **Next**.
   - Click **Install**. Wait for it to finish. Click **Finish**.

**How to verify it worked:**
- Press the `Windows` key, type `cmd`, press Enter.
- Type `where Rscript` and press Enter.
- You should see a path like: `C:\Program Files\R\R-4.4.1\bin\Rscript.exe`
- If you see `INFO: Could not find files`, R was not added to your PATH. See Troubleshooting at the end.
- Type `exit` and press Enter.

---

### Step 3 — Install Rtools

Rtools provides the compilers that some R packages need to build from source code.

1. Go to: **https://cran.r-project.org/bin/windows/Rtools/**
2. Click the version of Rtools that **matches your R version**:
   - If you installed R 4.3.x → click **Rtools43**
   - If you installed R 4.4.x → click **Rtools44**
   - If you installed R 4.5.x → click **Rtools45**
3. On the next page, click the **"Rtools4x installer"** link to download the `.exe` file.
4. **Double-click** the downloaded file to start the installer.
5. Settings:
   - **Select destination location**: leave the default (e.g. `C:\rtools44`). Click **Next**.
   - **Select additional tasks**: make sure **"Add rtools to system PATH"** is **checked**. Click **Next**.
   - Click **Install**. Wait. Click **Finish**.

**How to verify it worked:**
- Press the `Windows` key, type `cmd`, press Enter.
- Type `gcc --version` and press Enter.
- You should see something like: `gcc.exe (GCC) 13.2.0`
- Type `exit` and press Enter.

---

### Step 4 — Install Miniconda

Miniconda is a lightweight package manager that will install all the R packages the pipeline needs.

1. Go to: **https://docs.conda.io/en/latest/miniconda.html**
2. Under "Latest Miniconda Installer Links", find the row for **Windows** and click **Miniconda3 Windows 64-bit**.
3. **Double-click** the downloaded file (`Miniconda3-latest-Windows-x86_64.exe`).
4. Settings:
   - **License agreement**: click **I Agree**.
   - **Install type**: select **"Just Me (recommended)"**. Click **Next**.
   - **Install location**: leave the default (something like `C:\Users\YourName\miniconda3`). Click **Next**.
   - **Advanced options**:
     - **IMPORTANT**: Check **"Add Miniconda3 to my PATH environment variable"** (the installer warns against this — ignore the warning and check it anyway. You need this).
     - Leave "Register Miniconda3 as my default Python" checked.
     - Click **Install**.
   - Wait for it to finish. Click **Finish**.

**How to verify it worked:**
- **Close and reopen** any command prompt or PowerShell windows you have open (this is important — the PATH changes only take effect in new windows).
- Press the `Windows` key, type `powershell`, press Enter.
- Type `conda --version` and press Enter.
- You should see something like: `conda 24.1.2`
- Type `exit` and press Enter.

---

## Part 2: Download the pipeline code

Now you will use Git to download the pipeline code to your computer.

1. **Decide where to put the project.** For example, you might want it at `C:\cytometry`. Avoid paths with spaces (do NOT put it in `C:\My Documents\cytometry`).

2. Open **PowerShell**:
   - Press the `Windows` key on your keyboard.
   - Type `powershell`
   - Click on **"Windows PowerShell"** (the blue icon). **Do not** click "Run as administrator".

3. In the PowerShell window, type the following command and press Enter. Replace `YOUR_REPO_URL` with the actual Git repository URL you were given:

   ```
   git clone YOUR_REPO_URL C:\cytometry
   ```

   For example, if the repository is on GitHub:
   ```
   git clone https://github.com/your-org/cytometry.git C:\cytometry
   ```

   You will see progress messages as it downloads. Wait until it says something like `done.` and you get a new prompt line.

4. Now navigate into the project folder by typing:
   ```
   cd C:\cytometry
   ```

5. Verify you are in the right place by typing:
   ```
   dir
   ```
   You should see folders like `src`, `pipelines`, `config`, `docs`, and files like `environment_windows.yml`, `setup_env.ps1`, etc.

---

## Part 3: Run the setup script

The setup script installs all R packages and configures your environment. It will ask you one question.

1. You should still be in PowerShell in the `C:\cytometry` folder from Part 2. If not, open PowerShell and type:
   ```
   cd C:\cytometry
   ```

2. Windows blocks unsigned scripts by default. You need to allow this script to run. Type:
   ```
   Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
   ```
   If it asks **"Do you want to change the execution policy?"**, type `Y` and press Enter.

3. Now run the setup script:
   ```
   .\setup_env.ps1
   ```

4. The script will print a header and then ask:

   ```
   --- R Library Configuration ---

   R needs a single folder where ALL packages are installed.
   This avoids packages scattering across multiple locations.

   Enter the FULL path to your R library folder.
   Examples:
     D:\R_libs
     C:\Users\YourName\R\library
     C:\R\win-library\4.3

   R library path: _
   ```

   **What to type here:** This is the folder where you want all R packages to live. If you already have one, type its full path. If you don't have one yet, type a path where you want to create one, for example:
   ```
   C:\R_libs
   ```
   Press Enter. The script will create the folder if it doesn't exist.

5. The script will now:
   - Check that conda, Rtools, and Git are installed (you should see green `[OK]` messages).
   - Create the conda environment (this downloads many packages — **it can take 10-30 minutes**).
   - Install extra R packages from CRAN, Bioconductor, and GitHub.
   - Try to install some packages that may fail on Windows (this is expected — it will print `[SKIP]` for those).

6. **What to expect during installation:**
   - Many lines of text will scroll by. This is normal.
   - You may see some warnings in yellow or red. Warnings about packages being "skipped" are expected.
   - If you see `[FAIL] Could not activate environment`, see the Troubleshooting section below.
   - The entire process may take **20-60 minutes** depending on your internet speed.

7. When it finishes, you will see:

   ```
   ===================================================
   Setup Complete!
   ===================================================

   Configuration summary:
     R library path : C:\R_libs
     Conda env      : cytometry_env
     .Renviron      : C:\Users\YourName\.Renviron
   ```

**Keep this PowerShell window open** — you will use it in the next parts.

---

## Part 4: Verify the installation

Still in the same PowerShell window:

1. **Activate the environment** (you must do this every time you open a new PowerShell window to work with the pipeline):
   ```
   conda activate cytometry_env
   ```
   Your prompt should change from `PS C:\cytometry>` to something like `(cytometry_env) PS C:\cytometry>`. The `(cytometry_env)` at the beginning means the environment is active.

2. **Check that R is working:**
   ```
   Rscript -e "cat('R version:', R.version$version.string, '\n')"
   ```
   You should see: `R version: R version 4.x.x (2024-xx-xx)` or similar.

3. **Check that your library path is set correctly:**
   ```
   Rscript -e ".libPaths()"
   ```
   The FIRST path listed should be your library folder (e.g. `C:/R_libs`). If it is, all packages will install and load from there.

4. **Check a key package loads:**
   ```
   Rscript -e "library(data.table); cat('data.table OK\n'); library(ggplot2); cat('ggplot2 OK\n'); library(flowCore); cat('flowCore OK\n')"
   ```
   You should see three `OK` messages. If any fail, see Troubleshooting.

---

## Part 5: How to run a pipeline

Every time you want to use the pipeline, you need to do two things first: **open PowerShell** and **activate the conda environment**. Then you can run any pipeline or script.

### Step-by-step: running a pipeline

1. **Open PowerShell:**
   - Press the `Windows` key.
   - Type `powershell`.
   - Click "Windows PowerShell".

2. **Activate the environment** (you must do this every time):
   ```
   conda activate cytometry_env
   ```
   Wait until you see `(cytometry_env)` at the start of the line.

3. **Navigate to the pipelines folder:**
   ```
   cd C:\cytometry\pipelines
   ```

4. **Run the pipeline you need.** Replace the config path with your actual config file:

   **Pipeline A** — Full analysis from raw .fcs files:
   ```
   Rscript pipeline_A_raw_fcs.R --config ..\config\pipeline_A_config.R
   ```

   **Pipeline B** — Quick statistics from CSV expression files:
   ```
   Rscript pipeline_B_csv_expression.R --config ..\config\pipeline_B_config.R
   ```

   **Pipeline C** — Cluster abundance from CSV files with gating columns:
   ```
   Rscript pipeline_C_csv_gating.R --config ..\config\pipeline_C_config.R
   ```

   **Pipeline E** — Fast UMAP/clustering from preprocessed .fcs files:
   ```
   Rscript pipeline_E_preprocessed_umap.R --config ..\config\pipeline_E_config.R
   ```

5. The pipeline will run and print progress messages. When it finishes, your results will be in the output directory specified in the config file.

### Step-by-step: running an individual script

Individual scripts (the building blocks inside `src/`) must be run from **their own folder**. This is because each script finds shared utilities relative to its working directory.

1. Open PowerShell and activate the environment (steps 1-2 above).

2. Navigate to the script's folder. For example, to run the Leiden clustering script:
   ```
   cd C:\cytometry\src\clustering
   ```

3. Run it:
   ```
   Rscript leiden.R
   ```

**Important rule:** always `cd` into the folder that contains the script before running it. Do NOT run it from the project root like this:
```
# WRONG — will fail with "file not found" error:
cd C:\cytometry
Rscript src\clustering\leiden.R
```

---

## Part 6: Editing a config file before running

Each pipeline reads a config file that tells it where your data is, where to save results, and which methods to use. You need to edit this file to point to YOUR data.

1. Open **File Explorer** (the folder icon in your taskbar, or press `Windows + E`).

2. Navigate to `C:\cytometry\config\`.

3. Right-click on the config file you want to edit (e.g. `pipeline_A_config.R`) and choose **"Open with" > "Notepad"**. If you don't see "Open with", choose "Open with" > "Choose another app" > "Notepad".

4. The config file looks like this (simplified example):
   ```r
   params <- list(
     input_dir      = "C:/path/to/your/fcs/files",
     output_dir     = "C:/path/to/your/output",
     markers        = c("CD3", "CD4", "CD8", "CD19", "CD56"),
     clustering_method = "leiden",
     dimred_method     = "umap",
     batch_method      = "none"
   )
   ```

5. **Change the paths to match your data:**
   - `input_dir`: the folder where your `.fcs` or `.csv` data files are.
   - `output_dir`: the folder where you want results saved (it will be created if it doesn't exist).
   - Use **forward slashes** (`/`) in paths, not backslashes. For example:
     - Correct: `"C:/Users/Maria/Data/fcs_files"`
     - Wrong: `"C:\Users\Maria\Data\fcs_files"`

6. Save the file (`Ctrl+S`) and close Notepad.

7. Now run the pipeline as described in Part 5.

---

## Part 7: What each pipeline expects as input

| Pipeline | Input files | Where to put them |
|----------|------------|-------------------|
| **A** (Raw FCS) | `.fcs` files from the cytometer | One folder with all `.fcs` files |
| **B** (CSV Expression) | `.csv` files (rows = cells, columns = markers, no gating columns) | One folder with all `.csv` files |
| **C** (CSV Gating) | `.csv` files that have expression columns AND binary gating/cluster columns | One folder with all `.csv` files |
| **D** (FlowJo) | `.wsp` workspace files + `.fcs` files | Not supported on Windows — see note below |
| **E** (Preprocessed) | `.fcs` files that were already preprocessed (compensated, transformed) | One folder with all `.fcs` files |

**Note on Pipeline D:** Pipeline D requires FlowJo workspace import, which depends on the `flowWorkspace` R package. This package has C++ dependencies that do not compile reliably on Windows. If you need to compare gating strategies, export your gating from FlowJo as CSV files and use Pipeline C instead.

---

## Part 8: Every-day use (cheat sheet)

After initial setup is complete, this is all you need to do each time:

```powershell
# 1. Open PowerShell (Windows key > type "powershell" > click it)

# 2. Activate the environment
conda activate cytometry_env

# 3. Go to the pipelines folder
cd C:\cytometry\pipelines

# 4. Run the pipeline you need (pick one):
Rscript pipeline_A_raw_fcs.R --config ..\config\pipeline_A_config.R
Rscript pipeline_B_csv_expression.R --config ..\config\pipeline_B_config.R
Rscript pipeline_C_csv_gating.R --config ..\config\pipeline_C_config.R
Rscript pipeline_E_preprocessed_umap.R --config ..\config\pipeline_E_config.R
```

To run an individual script instead:
```powershell
conda activate cytometry_env
cd C:\cytometry\src\clustering
Rscript leiden.R
```

---

## Part 9: Windows compatibility — what works and what doesn't

### Works perfectly

| Feature | Notes |
|---------|-------|
| Loading CSV / FCS files | All data loading works |
| UMAP, t-SNE visualization | Works (Rtools needed for compilation — you installed it in Step 3) |
| Leiden, FlowSOM, DBSCAN clustering | All work |
| Statistical analysis (GLMM, limma, diffcyt, survival) | All work |
| All ggplot2-based plots (heatmaps, barplots, violin, etc.) | All work |
| Memory monitoring | Cross-platform support added |
| ComBat batch correction | Works (`batch_method = "combat"` in config) |

### May need workarounds

| Feature | Issue | What to do |
|---------|-------|------------|
| PhenoGraph clustering | The `Rphenograph` package sometimes fails to compile | Use `clustering_method = "leiden"` or `"flowsom"` instead |
| Harmony batch correction | Has a Python dependency that can be tricky | Use `batch_method = "combat"` instead |
| CytoNorm batch correction | Depends on flowWorkspace | Use `batch_method = "combat"` instead |
| GUI pop-up dialogs | The `tcltk` GUI can be unstable on some Windows versions | Already configured to fall back to text prompts automatically |

### Does not work on Windows

| Feature | Why | Alternative |
|---------|-----|-------------|
| Diffusion maps | The `destiny` R package does not have a Windows version | Use `dimred_method = "umap"` or `dimred_method = "tsne"` in your config |
| FlowJo workspace import (Pipeline D) | The `flowWorkspace` package C++ build fails on Windows | Export data from FlowJo as CSV files, then use Pipeline C |

---

## Troubleshooting

### "Rscript is not recognized as an internal or external command"

R was not added to your system PATH during installation.

**Fix:**
1. Press `Windows + R`, type `sysdm.cpl`, press Enter.
2. Click the **Advanced** tab.
3. Click **Environment Variables** at the bottom.
4. In the bottom section ("System variables"), find the row called **Path**. Click it, then click **Edit**.
5. Click **New** and type the path to R's `bin` folder. It will be something like:
   ```
   C:\Program Files\R\R-4.4.1\bin
   ```
   (Check in File Explorer what version you have.)
6. Click **OK** on all three windows.
7. **Close and reopen** PowerShell. Try `Rscript --version` again.

### "conda is not recognized as an internal or external command"

Miniconda was not added to your PATH.

**Fix:**
1. Press `Windows + R`, type `sysdm.cpl`, press Enter.
2. Click **Advanced** > **Environment Variables**.
3. In the bottom section ("System variables"), find **Path**, click **Edit**.
4. Click **New** and add these three paths (adjust the username):
   ```
   C:\Users\YourName\miniconda3
   C:\Users\YourName\miniconda3\Scripts
   C:\Users\YourName\miniconda3\Library\bin
   ```
5. Click **OK** on all windows. Close and reopen PowerShell. Try `conda --version`.

### "gcc is not recognized"

Rtools was not added to your PATH, or you didn't check the "Add rtools to system PATH" box during installation.

**Fix:** Uninstall Rtools from "Add or Remove Programs", then re-install following Step 3 above, making sure to check the PATH checkbox.

### The setup script says "running scripts is disabled on this system"

You need to allow PowerShell to run scripts.

**Fix:** In PowerShell, type:
```
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```
Type `Y` when asked, then try `.\setup_env.ps1` again.

### "conda env create" hangs for a very long time

The default conda solver can be slow. Install the faster solver:

```
conda install -n base -c conda-forge conda-libmamba-solver
conda config --set solver libmamba
```

Then try the setup script again.

### "Error: package 'xxxxx' was installed under R version 4.x.x"

Your R version was upgraded but the old packages remain. Reinstall the package:
```
Rscript -e "install.packages('xxxxx')"
```

### "Cannot open file 'src/R_shared/io_utils.R': No such file or directory"

You ran the script from the wrong folder. Each script must be run from its own directory:

```powershell
# WRONG:
cd C:\cytometry
Rscript src\clustering\leiden.R

# CORRECT:
cd C:\cytometry\src\clustering
Rscript leiden.R
```

### A package fails with "loading failed for 'x86_64-w64-mingw32'"

The package could not compile. Try installing the pre-compiled binary:
```
Rscript -e "install.packages('package_name', type = 'binary')"
```

### R scripts show strange errors like "unexpected input" or "incomplete string"

Git may have changed the line endings in the files. Fix:
```
cd C:\cytometry
git config core.autocrlf true
git checkout -- .
```

### "Warning: Memory usage could not be determined"

This is not an error. It means the pipeline could not query how much RAM your computer has. It will continue running normally. This can happen on very old Windows versions.

### Your project path has spaces and things break

Move the project to a path without spaces. For example, instead of:
```
C:\My Documents\lab stuff\cytometry    <-- spaces will cause problems
```
Use:
```
C:\cytometry                            <-- no spaces, works perfectly
```

---

## Appendix A: How to check which R packages are installed

```powershell
conda activate cytometry_env
Rscript -e "pkgs <- installed.packages()[,'Package']; cat(sort(pkgs), sep='\n')"
```

## Appendix B: How to install a single missing R package

```powershell
conda activate cytometry_env

# From CRAN (most packages):
Rscript -e "install.packages('package_name')"

# From Bioconductor (bioinformatics packages):
Rscript -e "BiocManager::install('package_name')"

# From GitHub:
Rscript -e "devtools::install_github('username/package_name')"
```

## Appendix C: How to update the pipeline code

When a new version of the pipeline is released:

```powershell
cd C:\cytometry
git pull
```

If the environment.yml changed, update the environment too:
```powershell
conda activate cytometry_env
conda env update -f environment_windows.yml --prune
Rscript install_extra_packages.R
```

## Appendix D: How to completely start over

If something is badly broken and you want to reset everything:

```powershell
# Remove the conda environment
conda deactivate
conda env remove -n cytometry_env

# Re-run setup
cd C:\cytometry
.\setup_env.ps1
```

## Appendix E: Folder structure after setup

After installation, your computer will have:

```
C:\cytometry\                       <-- the pipeline code (from Git)
    src\                            <-- all R scripts (the building blocks)
    pipelines\                      <-- pipeline orchestrators (A, B, C, E)
    config\                         <-- config files (edit these for your data)
    docs\                           <-- documentation (you are reading this)
    test_data\                      <-- example test files
    environment_windows.yml         <-- conda package list
    setup_env.ps1                   <-- the setup script you ran
    install_extra_packages.R        <-- extra package installer

C:\R_libs\                          <-- your R library (or wherever you chose)
    data.table\                     <-- installed R packages live here
    ggplot2\
    flowCore\
    ... (many more)

C:\Users\YourName\miniconda3\       <-- Miniconda installation
    envs\cytometry_env\             <-- the isolated conda environment

C:\Program Files\R\R-4.x.x\        <-- R installation
C:\rtools44\                        <-- Rtools (compiler) installation
C:\Users\YourName\.Renviron         <-- R config file (points to your library)
```
