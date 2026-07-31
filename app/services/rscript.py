"""Build and parse the ``scripts/umap_from_flowjo.R`` command line.

Everything the GUI does ultimately becomes one of two invocations:

* **inspect** — ``--list_populations``, which loads the workspace and prints the
  sample list, the gating tree and the common marker panel;
* **run** — the full analysis.

The parser below reads the inspect output. Its block headers are a contract
shared with ``src/visualization/umap_from_flowjo.R``; both sides carry a comment
saying so.
"""

from __future__ import annotations

import re
import shlex

from app import config


def _base_argv(script=None) -> list[str]:
    """The stdbuf + Rscript + script-path prefix.

    ``script`` defaults to the UMAP CLI so the two existing callers are
    unchanged; the FCS-names tab passes ``config.SYNC_SCRIPT``.
    """
    return (
        config.line_buffer_prefix()
        + config.rscript_command()
        + [str(script or config.CLI_SCRIPT)]
    )


def inspect_argv(wsp: str, fcs_dir: str, group: str,
                 allow_missing_samples: bool = False) -> list[str]:
    """Argv that loads the workspace(s) and prints samples/populations/markers.

    ``group`` is passed through as text: the R side takes either an index or a
    sample-group *name*, which is what lets one value cover a folder of
    workspaces that order their groups differently.

    ``allow_missing_samples`` has to be honoured HERE and not only on the run:
    the R side resolves every FCS file before it will print a population tree,
    so without it a single absent file leaves the form with no populations to
    pick from and the toggle further down unreachable.
    """
    argv = _base_argv() + [
        "--wsp", wsp,
        "--fcs_dir", fcs_dir,
        "--group", str(group).strip(),
        "--list_populations",
    ]
    if allow_missing_samples:
        argv.append("--allow_missing_samples")
    return argv


def run_argv(params: dict) -> list[str]:
    """Argv for a full run built from the form state.

    Only options the user actually set away from the R default are appended, so
    the printed command stays close to what someone would type by hand — with
    the exception of the paths and the population, which are always explicit.
    """
    argv = _base_argv() + [
        "--wsp", params["wsp"],
        "--fcs_dir", params["fcs_dir"],
    ]
    # --population is omitted entirely when a map is given: the R side ignores
    # it there, and printing a gate the run does not use would misrepresent the
    # command someone is about to read or paste.
    if not (params.get("population_map") or {}):
        argv.extend(["--population", params["population"]])
    argv.extend(["--output_dir", params["output_dir"]])

    def add(flag: str, value, default) -> None:
        if value is not None and value != default:
            argv.extend([flag, str(value)])

    group = params.get("group")
    add("--group", str(group).strip() if group is not None else None,
        config.DEFAULTS["group"])

    # A bare flag, so `add` (which appends a value) does not apply.
    if params.get("allow_missing_samples"):
        argv.append("--allow_missing_samples")

    # Per-workspace populations override --population entirely on the R side, so
    # the two are never both meaningful. ';' separates pairs because a FlowJo
    # gate name may contain a comma but not, in practice, a semicolon.
    pop_map = params.get("population_map") or {}
    if pop_map:
        argv.extend(["--population_map",
                     ";".join(f"{w}={p}" for w, p in pop_map.items() if p)])

    # Annotations. --sample_pattern is only meaningful without a CSV, and the
    # R side ignores it in that case, so it is not sent — the printed command
    # should not imply an option is doing something it is not.
    metadata = (params.get("sample_metadata") or "").strip()
    pattern = (params.get("sample_pattern") or "").strip()
    if metadata:
        argv.extend(["--sample_metadata", metadata])
    elif pattern:
        argv.extend(["--sample_pattern", pattern])
    annotate = params.get("annotate") or []
    if annotate and (metadata or pattern):
        argv.extend(["--annotate", ",".join(annotate)])

    markers = params.get("markers") or []
    if markers:
        argv.extend(["--markers", ",".join(markers)])

    # --subpop_select overrides --subpop_mode on the R side, so only one of the
    # two is ever passed.
    subpop_select = params.get("subpop_select") or []
    if subpop_select:
        argv.extend(["--subpop_select", ",".join(subpop_select)])
    else:
        add("--subpop_mode", params.get("subpop_mode"),
            config.DEFAULTS["subpop_mode"])

    add("--normalize", params.get("normalize"), config.DEFAULTS["normalize"])
    add("--downsample", params.get("downsample"), config.DEFAULTS["downsample"])
    add("--n_neighbors", params.get("n_neighbors"), config.DEFAULTS["n_neighbors"])
    add("--min_dist", params.get("min_dist"), config.DEFAULTS["min_dist"])
    add("--metric", params.get("metric"), config.DEFAULTS["metric"])
    add("--seed", params.get("seed"), config.DEFAULTS["seed"])
    add("--point_size", params.get("point_size"), config.DEFAULTS["point_size"])
    add("--point_alpha", params.get("point_alpha"), config.DEFAULTS["point_alpha"])
    add("--plot_width", params.get("plot_width"), config.DEFAULTS["plot_width"])
    add("--plot_height", params.get("plot_height"), config.DEFAULTS["plot_height"])

    formats = params.get("output_formats") or []
    if sorted(formats) != sorted(config.DEFAULTS["output_formats"]):
        argv.extend(["--output_formats", ",".join(formats)])

    return argv


# --- sync_fcs_names.R --------------------------------------------------------
#
# Five modes, ordered read-only first. That ordering is the safety design: the
# two harmless modes are the ones you land on, and reaching a write means moving
# past two labels that say "writes nothing".
#
# Only 4 of the script's 16 flags get a UI. The rest (--set, --unset, --guid,
# --strip_extension, --allow, --ledger, --no_recursive, --no_original) are
# deliberately CLI-only: they are expert options, and --allow in particular
# downgrades a refusal that exists to stop the tool writing to a file it does
# not fully understand. That must stay a typed act, not a checkbox.
SYNC_MODES: dict[str, dict] = {
    "probe": {
        "label": "Census — read every file and report. Writes nothing.",
        "flags": ["--probe"], "button": "Run census", "color": "blue",
    },
    "dry": {
        "label": "Dry run — show exactly what would change. Writes nothing.",
        "flags": [], "button": "Run dry run", "color": "blue",
    },
    "inplace": {
        "label": "APPLY — patch the original files where they sit.",
        "flags": ["--apply", "--in_place"], "button": "Apply in place",
        "color": "red", "writes": True,
    },
    "copy": {
        "label": "APPLY — write patched copies to another folder.",
        "flags": ["--apply"], "button": "Apply to copies", "color": "red",
        "writes": True, "needs_output": True,
    },
    "revert": {
        "label": "Undo — revert the most recent run from its ledger.",
        "flags": ["--revert"], "button": "Revert last run", "color": "orange",
        "writes": True,
    },
}


def sync_argv(fcs_dir: str, mode: str, also: str = "",
              output_dir: str = "") -> list[str]:
    """Argv for one FCS-names run."""
    spec = SYNC_MODES[mode]
    argv = _base_argv(config.SYNC_SCRIPT) + ["--fcs_dir", fcs_dir]
    argv += list(spec["flags"])
    if spec.get("needs_output") and output_dir.strip():
        argv += ["--output_dir", output_dir.strip()]
    # Revert replays a ledger; every other option would be noise on that command.
    if mode != "revert" and also.strip():
        argv += ["--also", also.strip()]
    return argv


def keyword_table_argv(fcs_dir: str, csv_path: str) -> list[str]:
    """Argv that dumps every keyword of every file to a CSV. Read-only."""
    return _base_argv(config.SYNC_SCRIPT) + [
        "--fcs_dir", fcs_dir, "--probe", "--keyword_table", csv_path,
    ]


def command_string(argv: list[str]) -> str:
    """The command as you would paste it into a shell."""
    return shlex.join(argv)


# --- inspect output parsing --------------------------------------------------

_SAMPLES_HEADER = re.compile(r"^Samples in loaded group \(\s*(\d+)\s*\):\s*$")
_POPS_HEADER = re.compile(r"^Population tree \(\s*(\d+)\s* nodes\):\s*$")
_MARKERS_HEADER = re.compile(r"^Common markers \(\s*(\d+)\s*\):\s*$")
_WORKSPACES_HEADER = re.compile(r"^Workspaces resolved \(\s*(\d+)\s*\):\s*$")

_SAMPLE_ITEM = re.compile(r"^  - (.+?)\s*$")
# Not anchored on "/": gs_get_pop_paths() lists the bare "root" node first, and
# requiring a leading slash made _collect stop on that very first line — so the
# population list came back EMPTY for every real workspace.
_POP_ITEM = re.compile(r"^  (\S.*?)\s*$")
_MARKER_ITEM = re.compile(r"^  (\S.*?)\s*$")
# Four-space indent, deliberately: it keeps the two-space item patterns above
# from swallowing this block, which is printed after them.
_WORKSPACE_ITEM = re.compile(r"^    (\S.*?)\s*$")

_POPS_BY_WSP_HEADER = re.compile(r"^Populations by workspace \(\s*(\d+)\s*\):\s*$")
# Inside that block the workspace rows keep the four-space indent and its gates
# take six, which is what lets the two levels be told apart on the way back in.
_PBW_WORKSPACE = re.compile(r"^    (\S.*?)\s*\|\s*\d+ node\(s\)\s*$")
_PBW_GATE = re.compile(r"^      (\S.*?)\s*$")


def _collect(lines: list[str], start: int, item_re: re.Pattern) -> list[str]:
    """Consume consecutive ``item_re`` matches starting at ``start``."""
    out: list[str] = []
    for line in lines[start:]:
        m = item_re.match(line.rstrip("\n"))
        if not m:
            break
        out.append(m.group(1))
    return out


def parse_inspect(stdout: str) -> dict:
    """Extract samples, populations, markers and workspaces from inspect output.

    Missing blocks come back as empty lists rather than raising: the caller
    decides what is fatal (no populations is; an empty marker list is only a
    warning, since the R side prints an explanatory line in that case).

    With several workspaces the samples are the run-wide display labels and the
    populations are the gates *shared by every workspace* — the only ones a run
    can use, since the R side refuses a population it cannot resolve everywhere.
    """
    lines = stdout.splitlines()
    result: dict = {
        "samples": [], "populations": [], "markers": [], "workspaces": [],
        # {workspace label: its own gates}. Each workspace's FULL tree, unlike
        # "populations", which is the intersection and collapses to a few nodes
        # when the workspaces were gated independently.
        "populations_by_workspace": {},
    }

    for i, raw in enumerate(lines):
        line = raw.rstrip("\n")
        if _SAMPLES_HEADER.match(line):
            result["samples"] = _collect(lines, i + 1, _SAMPLE_ITEM)
        elif _POPS_HEADER.match(line):
            # The R side prints a parenthetical when no gate is shared; prose,
            # not a gate name (same convention as the marker block).
            result["populations"] = [
                p for p in _collect(lines, i + 1, _POP_ITEM)
                if not p.startswith("(")
            ]
        elif _MARKERS_HEADER.match(line):
            markers = _collect(lines, i + 1, _MARKER_ITEM)
            # The R side prints an explanatory parenthetical when the set is
            # empty; that is prose, not a marker name.
            result["markers"] = [m for m in markers if not m.startswith("(")]
        elif _POPS_BY_WSP_HEADER.match(line):
            result["populations_by_workspace"] = _collect_pops_by_workspace(
                lines, i + 1)
        elif _WORKSPACES_HEADER.match(line):
            result["workspaces"] = [
                _parse_workspace_row(row)
                for row in _collect(lines, i + 1, _WORKSPACE_ITEM)
            ]
    return result


def _collect_pops_by_workspace(lines: list[str], start: int) -> dict:
    """Two-level block: a workspace row, then its gates, until neither matches."""
    out: dict[str, list[str]] = {}
    current: str | None = None
    for line in lines[start:]:
        line = line.rstrip("\n")
        ws = _PBW_WORKSPACE.match(line)
        if ws:
            current = ws.group(1)
            out[current] = []
            continue
        gate = _PBW_GATE.match(line)
        if gate and current is not None:
            out[current].append(gate.group(1))
            continue
        break
    return out


def _parse_workspace_row(row: str) -> dict:
    """One ``label | path | group N "name" | K sample(s)`` row.

    Anything unexpected degrades to a label-only entry rather than raising —
    this block is for display, and a formatting drift should not break a load
    that otherwise succeeded.
    """
    parts = [p.strip() for p in row.split("|")]
    entry = {"label": parts[0] if parts else row, "wsp": "", "group": "",
             "samples": ""}
    if len(parts) >= 2:
        entry["wsp"] = parts[1]
    if len(parts) >= 3:
        entry["group"] = parts[2]
    if len(parts) >= 4:
        entry["samples"] = parts[3]
    return entry
