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
# Sample-metadata tables. The R side sniffs "," vs ";" itself.
METADATA_EXTENSIONS = (".csv", ".tsv", ".txt")

# Plot files the results gallery knows about. PNG is preferred for inline
# preview: an SVG scatter of a large embedding carries one node per cell and can
# be tens of MB, which browsers handle badly.
PLOT_EXTENSIONS = (".png", ".svg", ".pdf")
PREVIEW_EXTENSIONS = (".png", ".svg")
EXTENSION_PREFERENCE = (".png", ".svg", ".pdf")

SUBPOP_MODES = ("furthest", "cd4cd8", "direct", "deepest", "none")
NORMALIZE_METHODS = ("none", "zscore", "minmax")
# uwot's supported distance metrics.
UMAP_METRICS = ("euclidean", "cosine", "manhattan", "hamming", "correlation")
OUTPUT_FORMATS = ("svg", "png", "pdf")

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
    "subpop_mode": "furthest",
    "explicit_subpop": False,
    "subpop_select": [],
    "normalize": "none",
    "downsample": 0,
    "n_neighbors": 15,
    "min_dist": 0.1,
    "metric": "euclidean",
    "seed": 42,
    "point_size": 0.4,
    "point_alpha": 0.6,
    "plot_width": 10.0,
    "plot_height": 8.0,
    "output_formats": ["svg"],
    "output_parent": "",
    "output_name": "",
    # FCS-names tab. The folder only. The write mode is deliberately NOT
    # remembered: a refresh must never bring the app back pre-armed to apply.
    "fcs_sync_dir": "",
}


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
