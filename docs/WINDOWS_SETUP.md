# Running the FlowJo → UMAP app on a Windows PC

This guide starts from the point where you have **downloaded the ZIP from GitHub,
unpacked it, and have VS Code installed**.

**Time required:** ~45–75 minutes, almost all of it waiting for the R packages.

At the end you run one command, `python app/app.py`, and use the tool in your
normal Windows browser.

---

## Read this first: why we use WSL

The app itself is Python, but it is a front-end — every button shells out to
`Rscript scripts/umap_from_flowjo.R`, which reads FlowJo `.wsp` workspaces using
the Bioconductor packages **CytoML** and **flowWorkspace**. Those have a C++
backend that does not build reliably on Windows (which is why
`environment_windows.yml` excludes them). Two smaller things break natively too:
the Cancel button uses POSIX process groups, and `environment_windows.yml`
contains no Python at all.

The fix is **WSL2** — a real Ubuntu running inside Windows, shipped by Microsoft,
free, no dual boot, no VM to configure. Your Windows files stay where they are,
and the app's web page opens in your ordinary Windows browser.

> You do **not** need Nextflow, Docker, Rtools, or a separate Windows R install.
> Nextflow only orchestrates the batch pipelines A–E; the app does not use it.
> `setup_env.ps1` and `environment_windows.yml` are obsolete — ignore them.

---

## Step 1 — Install WSL2 with Ubuntu

*Skip if you already have a WSL terminal open.*

1. `Windows` key → type `powershell` → **right-click** "Windows PowerShell" →
   **Run as administrator**.
2. Run:

   ```powershell
   wsl --install
   ```

3. **Reboot.** An Ubuntu window opens by itself and asks you to create a UNIX
   username and password (new, unrelated to your Windows login — the password
   shows no characters as you type, which is normal).

Reopen it any time with `Windows` key → `Ubuntu`.

---

## Step 2 — Copy the project into Ubuntu

Ubuntu sees your Windows drives under `/mnt/c/...`, but running from there is
slow and causes permission problems. Copy the project into the Linux home
directory instead. In the Ubuntu terminal, adjusting the source path:

```bash
cp -r "/mnt/c/Users/YOUR_WINDOWS_USER/Downloads/cytometry_pipeline-main" ~/cytometry_pipeline
cd ~/cytometry_pipeline
ls          # expect: app/  src/  scripts/  docs/  environment.yml  setup_env.sh
```

If `ls` shows no `app/` folder, you downloaded a branch other than `main` — get
the ZIP from the repository's default branch.

**VS Code:** install the **WSL** extension (`Ctrl+Shift+X`, search "WSL",
publisher Microsoft), then run `code .` from this folder. The bottom-left corner
should read `WSL: Ubuntu`, and VS Code's terminal (``Ctrl+` ``) becomes an Ubuntu
terminal you can use for everything below.

---

## Step 3 — Install Miniconda (inside Ubuntu)

```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Enter through the licence, type `yes` to accept, accept the default location,
and answer **`yes`** to *"update your shell profile to automatically initialize
conda?"*.

Then **close and reopen the terminal** and check:

```bash
conda --version
conda config --set solver libmamba     # optional, but much faster
```

---

## Step 4 — Build the environment

```bash
cd ~/cytometry_pipeline
bash setup_env.sh
```

This creates `cytometry_env` with R, Bioconductor, **and** the app's Python
dependencies (`environment.yml` includes `python`, `dash`,
`dash-mantine-components`, `dash-iconify`, `gunicorn`). It takes **20–60
minutes** and ends with `Setup Complete!`.

**Verify both halves:**

```bash
conda activate cytometry_env
Rscript -e "library(CytoML); library(flowWorkspace); library(uwot); cat('R OK\n')"
python -c "import dash, dash_mantine_components; print('Python OK')"
```

---

## Step 5 — Run the app

```bash
conda activate cytometry_env
cd ~/cytometry_pipeline
python app/app.py
```

It prints:

```
FlowJo → UMAP app on http://127.0.0.1:8050
```

Open **http://127.0.0.1:8050** in Edge or Chrome **on Windows** — WSL2 forwards
localhost automatically, so no extra configuration is needed.

Stop the app with `Ctrl+C` in the terminal.

| Flag | Default | Purpose |
|------|---------|---------|
| `--port` | `8050` | Use another port if 8050 is taken |
| `--debug` | off | Hot-reload and full error pages |

Leave `--host` alone. The file browser is rooted at `/` with no authentication,
so binding to anything but `127.0.0.1` hands a filesystem browser and a job
runner to anyone who can reach the port.

---

## Step 6 — Point it at your data

Put your `.wsp` workspaces and `.fcs` files **inside the Linux side**, so the app's
file browser and R both see plain Linux paths:

```bash
mkdir -p ~/data
cp -r "/mnt/c/Users/YOUR_WINDOWS_USER/Desktop/my_experiment" ~/data/
```

`/mnt/c/...` paths do work in the browser if you prefer not to copy — expect it
to be noticeably slower on large FCS files.

Then follow the four steps in the UI: **Data** (workspace folder + FCS folder →
Load workspaces) → **Population & markers** → **Parameters** → **Run**. The exact
command is shown before it runs, and output streams into the Log pane.

**Getting results back to Windows:** paste `\\wsl$\Ubuntu\home\YOUR_UNIX_USER\`
into File Explorer's address bar, or copy explicitly:

```bash
cp -r ~/cytometry_results "/mnt/c/Users/YOUR_WINDOWS_USER/Desktop/"
```

Full feature documentation: [`app/README.md`](../app/README.md).

---

## Troubleshooting

**`wsl --install` is not recognised**
Windows too old. Update Windows, or install "Ubuntu" from the Microsoft Store and
enable *Virtual Machine Platform* under "Turn Windows features on or off".

**`python: command not found` / `ModuleNotFoundError: dash`**
The environment is not active. `conda activate cytometry_env` — needed in every
new terminal window.

**`Cannot find .../scripts/umap_from_flowjo.R`**
You are not in the repository root, or you unpacked a branch without `scripts/`.
`cd ~/cytometry_pipeline` first.

**The page does not load in the Windows browser**
Check the terminal still shows the app running, and use `127.0.0.1:8050` rather
than `localhost:8050`. If the port is busy: `python app/app.py --port 8060`.

**"Loading workspaces" takes minutes**
Expected — it builds a full `GatingSet` per workspace. Watch the Log tab.

**Everything is very slow**
You are running from `/mnt/c/...`. Move the project (and ideally the data) into
`~` as in Step 2.

**WSL eats too much RAM**
Create `C:\Users\YOUR_WINDOWS_USER\.wslconfig` containing:
```
[wsl2]
memory=12GB
```
then `wsl --shutdown` in PowerShell and reopen Ubuntu.

---

## Everyday cheat sheet

```bash
# Windows key → "Ubuntu", then:
conda activate cytometry_env
cd ~/cytometry_pipeline
python app/app.py
# → open http://127.0.0.1:8050 in your Windows browser
```

---

## If you also need the batch pipelines (A–E)

Not required for the app. Those are Nextflow workflows; install it into the same
environment with `conda install -c bioconda -c conda-forge nextflow` and run e.g.

```bash
nextflow run nextflow/main.nf -c nextflow/conf/local.config \
    -params-file nextflow/conf/pipeline_A_params.yaml
```
