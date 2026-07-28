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


def _base_argv() -> list[str]:
    return (
        config.line_buffer_prefix()
        + config.rscript_command()
        + [str(config.CLI_SCRIPT)]
    )


def inspect_argv(wsp: str, fcs_dir: str, group: int) -> list[str]:
    """Argv that loads the workspace and prints samples/populations/markers."""
    return _base_argv() + [
        "--wsp", wsp,
        "--fcs_dir", fcs_dir,
        "--group", str(int(group)),
        "--list_populations",
    ]


def run_argv(params: dict) -> list[str]:
    """Argv for a full run built from the form state.

    Only options the user actually set away from the R default are appended, so
    the printed command stays close to what someone would type by hand — with
    the exception of the paths and the population, which are always explicit.
    """
    argv = _base_argv() + [
        "--wsp", params["wsp"],
        "--fcs_dir", params["fcs_dir"],
        "--population", params["population"],
        "--output_dir", params["output_dir"],
    ]

    def add(flag: str, value, default) -> None:
        if value is not None and value != default:
            argv.extend([flag, str(value)])

    add("--group", params.get("group"), config.DEFAULTS["group"])
    add("--expected_samples", params.get("expected_samples"),
        config.DEFAULTS["expected_samples"])

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


def command_string(argv: list[str]) -> str:
    """The command as you would paste it into a shell."""
    return shlex.join(argv)


# --- inspect output parsing --------------------------------------------------

_SAMPLES_HEADER = re.compile(r"^Samples in loaded group \(\s*(\d+)\s*\):\s*$")
_POPS_HEADER = re.compile(r"^Population tree \(\s*(\d+)\s* nodes\):\s*$")
_MARKERS_HEADER = re.compile(r"^Common markers \(\s*(\d+)\s*\):\s*$")

_SAMPLE_ITEM = re.compile(r"^  - (.+?)\s*$")
_POP_ITEM = re.compile(r"^  (/.*?)\s*$")
_MARKER_ITEM = re.compile(r"^  (\S.*?)\s*$")


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
    """Extract samples, populations and markers from the inspect output.

    Missing blocks come back as empty lists rather than raising: the caller
    decides what is fatal (no populations is; an empty marker list is only a
    warning, since the R side prints an explanatory line in that case).
    """
    lines = stdout.splitlines()
    result: dict[str, list[str]] = {"samples": [], "populations": [], "markers": []}

    for i, raw in enumerate(lines):
        line = raw.rstrip("\n")
        if _SAMPLES_HEADER.match(line):
            result["samples"] = _collect(lines, i + 1, _SAMPLE_ITEM)
        elif _POPS_HEADER.match(line):
            result["populations"] = _collect(lines, i + 1, _POP_ITEM)
        elif _MARKERS_HEADER.match(line):
            markers = _collect(lines, i + 1, _MARKER_ITEM)
            # The R side prints an explanatory parenthetical when the set is
            # empty; that is prose, not a marker name.
            result["markers"] = [m for m in markers if not m.startswith("(")]
    return result
