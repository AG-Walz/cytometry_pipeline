"""Paths, defaults and small shared constants for the FlowJo-UMAP Dash app.

The DEFAULTS dict mirrors the defaults of ``scripts/umap_from_flowjo.R`` so the
GUI never silently disagrees with the CLI. If an option's default changes there,
change it here too.
"""

from __future__ import annotations

import os
import shlex
import shutil
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
REPO_ROOT = APP_DIR.parent
CLI_SCRIPT = REPO_ROOT / "scripts" / "umap_from_flowjo.R"
SYNC_SCRIPT = REPO_ROOT / "scripts" / "sync_fcs_names.R"

# Where the last-used form state is remembered between launches.
CACHE_DIR = Path(
    os.environ.get("UMAP_APP_CACHE", Path.home() / ".cache" / "umap_flowjo_app")
)
SESSION_FILE = CACHE_DIR / "last_session.json"

# File-browser filters.
WSP_EXTENSIONS = (".wsp",)

# Plot files the results gallery knows about. PNG is preferred for inline
# preview: an SVG scatter of a large embedding carries one node per cell and can
# be tens of MB, which browsers handle badly.
PLOT_EXTENSIONS = (".png", ".svg", ".pdf")
PREVIEW_EXTENSIONS = (".png", ".svg")
EXTENSION_PREFERENCE = (".png", ".svg", ".pdf")
#: An SVG larger than this is offered for download instead of embedded: with
#: one node per cell, a cohort-sized embedding is a file the browser cannot
#: render inline, and one such <img> freezes the whole gallery.
MAX_INLINE_SVG_BYTES = 8 * 1024 * 1024

#: Non-plot outputs offered for download in the gallery: cells.csv.gz,
#: combination_frequencies_<col>.csv, separation_metrics.csv, run_info.json and
#: friends. ".gz" is the suffix of "cells.csv.gz".
DATA_EXTENSIONS = (".csv", ".tsv", ".json", ".gz")

SUBPOP_MODES = ("furthest", "cd4cd8", "direct", "deepest", "combinations",
                "none")
#: --layout values, in the order the form lists them.
LAYOUTS = (
    ("gates", "Gate-separated islands (recommended)"),
    ("markers", "Marker values"),
)
NORMALIZE_METHODS = ("none", "zscore", "minmax")
# uwot's supported distance metrics.
UMAP_METRICS = ("euclidean", "cosine", "manhattan", "hamming", "correlation")
OUTPUT_FORMATS = ("svg", "png", "pdf")

#: The R script's own --output_formats default. The form default below differs
#: on purpose (PNG previews; SVG of a large embedding does not), so run_argv
#: compares against THIS to decide whether the flag needs to be passed.
R_DEFAULT_OUTPUT_FORMATS = ("svg",)

#: Defaults for every form field, mirroring the R CLI.
DEFAULTS: dict = {
    "wsp": "",
    "fcs_dir": "",
    # Text, not a number: the R side takes an index OR a sample-group name.
    "group": "1",
    "allow_missing_samples": False,
    # How a workspace's samples are paired with FCS files. "content" (the R
    # default) tries the file name first, then lets CytoML pair the rest on the
    # $FIL keyword inside each file plus its event count. "strict" makes a name
    # mismatch final. See --match_mode in scripts/umap_from_flowjo.R.
    "match_mode": "content",
    # "INFy=IFN;TNFa=TNF" — spellings to treat as the same marker. Empty by
    # default: merging two markers is the user's call, never a guess.
    "marker_aliases": "",
    "recover_markers": True,
    "population": "",
    # {workspace label: gate}. Empty means "one population for every
    # workspace", the original behaviour; non-empty switches the run to
    # --population_map. A workspace absent from it is excluded from the run.
    "population_map": {},
    "sample_metadata": "",
    "sample_pattern": "",
    "annotate": [],
    "markers": [],
    # Plots split by a metadata column, the order of its levels ("col=a,b;
    # col2=c"; the first level is the reference of the fold-change maps), and
    # one "keep only samples where" clause, which run_params turns into
    # --filter "col=a,b" or "col!=a,b".
    "split_by": [],
    "level_order": "",
    "filter_column": "",
    "filter_values": [],
    "filter_mode": "include",
    "subpop_mode": "furthest",
    "explicit_subpop": False,
    "subpop_select": [],
    # R's default is "markers"; the form starts on "gates" (FORM_START below).
    # run_argv passes --layout always, so the omission rule cannot turn the
    # form's choice back into R's default.
    "layout": "markers",
    # Gate layout only. Keep these equal to the R defaults: they are what
    # run_argv compares against to decide whether a flag is needed.
    "gate_gap": 2,
    "gate_tier_ratio": 2,
    "gate_softness": 4,
    "supervise": 0,
    "cells_per_combination": 2000,
    "min_combination_cells": 50,
    # Size of the balanced fit set in the gate layout; 0 there means 100000.
    "umap_subsample": 0,
    # Markers layout only: what UMAP is fitted on before every cell is
    # projected (--fit_subset): "random" cells (random_fit_cells of them, sent
    # as --umap_subsample; 0 = all) or SuperCellCyto "supercell" means.
    "fit_subset": "random",
    "random_fit_cells": 0,
    "supercell_gamma": 20,
    "supercell_k_knn": 5,
    # At most this many supercells are fitted on; R raises gamma to match.
    "supercell_max": 100000,
    "normalize": "none",
    "downsample": 0,
    # Worker processes for the per-file plots, which are ~90% of a run.
    # 0 = decide automatically (all cores but one, capped by the cgroup quota).
    "n_cores": 0,
    "n_neighbors": 15,
    # R resolves an omitted --min_dist per layout (0.1 markers, 0.25 gates);
    # run_argv always passes the form's value, so this is the markers value
    # and the one the layout switch swaps against.
    "min_dist": 0.1,
    "metric": "euclidean",
    "seed": 42,
    "point_size": 0.4,
    "point_alpha": 0.6,
    "plot_width": 10.0,
    "plot_height": 8.0,
    "output_formats": ["png"],
    "output_parent": "",
    "output_name": "",
    # Reuse the steps an earlier run into the same folder saved (R's default);
    # off passes --no_resume.
    "resume": True,
    # Folder the Results tab shows. Set by the last finished run, or by hand.
    "results_dir": "",
    # FCS-names tab. The folder only. The write mode is deliberately NOT
    # remembered: a refresh must never bring the app back pre-armed to apply.
    "fcs_sync_dir": "",
}

#: Where the form deliberately starts somewhere other than the R default. Used
#: for a fresh session and by Reset form; a remembered choice still wins.
FORM_START: dict = {
    "layout": "gates",
    # The gate layout's tuned min_dist (rounder islands). The layout switch
    # swaps 0.1 and 0.25 while the user has not set a value of their own.
    "min_dist": 0.25,
}

#: min_dist each layout starts from; see FORM_START and steps.layout_min_dist.
LAYOUT_MIN_DIST: dict = {"markers": 0.1, "gates": 0.25}


def rscript_command() -> list[str]:
    """The argv prefix used to invoke R.

    Defaults to plain ``Rscript`` because the app is meant to run inside the
    same ``cytometry_env`` conda environment as the pipeline. Override with
    ``UMAP_RSCRIPT`` (e.g. ``conda run -n cytometry_env Rscript``) when Python
    and R live in different environments.
    """
    override = os.environ.get("UMAP_RSCRIPT", "").strip()
    return shlex.split(override) if override else ["Rscript"]


def line_buffer_prefix() -> list[str]:
    """``stdbuf`` prefix that forces line-buffered output, when available.

    R block-buffers stdout when it is a pipe rather than a terminal, which would
    make the live log arrive only once the run finished. ``stdbuf -oL -eL``
    restores line buffering; where it does not exist (macOS without coreutils)
    the log simply updates in chunks.
    """
    if shutil.which("stdbuf"):
        return ["stdbuf", "-oL", "-eL"]
    return []
