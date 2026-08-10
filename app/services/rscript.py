"""Build the ``scripts/umap_from_flowjo.R`` command line, and read what it wrote.

Everything the GUI does ultimately becomes one of two invocations:

* **inspect** — ``--list_populations``, which loads the workspaces and, with
  ``--report_json``, writes the whole load out as data;
* **run** — the full analysis.

This module builds those command lines and reads the JSON back. It does not
compute anything: every number, gate, marker and diagnosis in the report is R's,
because R is the only side that can see a GatingSet.

There used to be a ~120-line regular-expression parser here that recovered all of
that from the run's stdout, against four block headers and three significant
indentation levels held in sync by a comment on each side. It is gone. The R
side writes ``--report_json`` instead — see
``src/visualization/umap_load_report.R`` for why.
"""

from __future__ import annotations

import json
import shlex
from pathlib import Path

from app import config

#: Major schema version this reader understands. R writes ``schema`` into every
#: report; a mismatch here means the two have drifted and is reported rather
#: than papered over with ``.get()`` calls that quietly return nothing.
REPORT_SCHEMA_MAJOR = 1


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
                 allow_missing_samples: bool = False,
                 match_mode: str = "content",
                 marker_aliases: str = "",
                 recover_markers: bool = True,
                 report_json: str = "") -> list[str]:
    """Argv that loads the workspace(s) and writes the load report.

    Deliberately has NO ``--workspaces``: a load reads every workspace it finds.
    The health panel's on/off switches decide what goes into the RUN, and they
    can only decide it for workspaces the report actually contains — filtering
    the load too would delete a switched-off workspace's row from the panel it
    is switched off in, leaving nothing to switch back on. See run_argv, which
    is where the selection belongs and where it is converted from the exclusion
    set the UI stores to the inclusion list the R flag takes.

    ``group`` is passed through as text: the R side takes either an index or a
    sample-group *name*, which is what lets one value cover a folder of
    workspaces that order their groups differently.

    ``allow_missing_samples`` has to be honoured HERE and not only on the run:
    the R side resolves every FCS file before it will print a population tree,
    so without it a single absent file leaves the form with no populations to
    pick from and the toggle further down unreachable.

    ``match_mode`` and ``marker_aliases`` are passed here too, and for the same
    reason: both change WHICH workspaces load and WHICH markers they agree on,
    so an inspect that omitted them would populate the form from a different
    dataset than the run it is configuring.
    """
    argv = _base_argv() + [
        "--wsp", wsp,
        "--fcs_dir", fcs_dir,
        "--group", str(group).strip(),
        "--list_populations",
    ]
    if allow_missing_samples:
        argv.append("--allow_missing_samples")
    if match_mode and match_mode != config.DEFAULTS["match_mode"]:
        argv.extend(["--match_mode", match_mode])
    if (marker_aliases or "").strip():
        argv.extend(["--marker_aliases", marker_aliases.strip()])
    if not recover_markers:
        argv.append("--no_recover_markers")
    if report_json:
        argv.extend(["--report_json", report_json])
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

    # Both of these decide WHICH data the run sees, not how it is drawn, so they
    # are passed whenever they differ from the R default — never inferred.
    add("--match_mode", params.get("match_mode"), config.DEFAULTS["match_mode"])
    aliases = (params.get("marker_aliases") or "").strip()
    if aliases:
        argv.extend(["--marker_aliases", aliases])
    if not params.get("recover_markers", True):
        argv.append("--no_recover_markers")
    # Only sent when something is actually switched off: the flag's absence
    # means "every workspace", which is both the R default and what the printed
    # command should say when nothing was excluded.
    selected = params.get("workspaces") or []
    if selected:
        argv.extend(["--workspaces", ",".join(selected)])

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
    add("--n_cores", params.get("n_cores"), config.DEFAULTS["n_cores"])
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
# Three modes, read-only first, exactly one of which writes to your data. The
# CLI has more (a census that skips the per-file preview, copies into another
# folder, --set/--unset/--guid/--strip_extension/--allow); they are expert
# options and the printed command is copy-pasteable, so adding one by hand is
# easy. --allow especially: it downgrades a refusal that exists to stop the
# tool writing to a file it does not fully understand, which should stay a
# deliberate typed act rather than a checkbox.
SYNC_MODES: dict[str, dict] = {
    "dry": {
        "label": "Preview — list what would change. Writes nothing.",
        "flags": [], "button": "Preview", "color": "blue",
    },
    "inplace": {
        "label": "Apply — write $FIL into the original files.",
        "flags": ["--apply", "--in_place"], "button": "Apply",
        "color": "red", "writes": True,
    },
    "revert": {
        "label": "Undo — restore the names from the last run's ledger.",
        "flags": ["--revert"], "button": "Undo", "color": "orange",
        "writes": True,
    },
}


def sync_argv(fcs_dir: str, mode: str, also: str = "") -> list[str]:
    """Argv for one FCS-names run."""
    spec = SYNC_MODES[mode]
    argv = _base_argv(config.SYNC_SCRIPT) + ["--fcs_dir", fcs_dir]
    argv += list(spec["flags"])
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


# --- reading the load report -------------------------------------------------


class ReportError(Exception):
    """The report is missing, unreadable, or of a schema we do not know."""


def read_report(path: str | Path) -> dict:
    """Load one ``--report_json`` document.

    Raises ``ReportError`` with something a user can act on. The caller treats
    that as "the load did not produce a usable result", which is the same thing
    it did when the stdout parser came back with no populations — except that
    now the reason is specific.
    """
    p = Path(path)
    if not p.is_file():
        raise ReportError(
            "R did not write its load report. The workspaces may have failed "
            "before finishing — see the Log tab."
        )
    try:
        report = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ReportError(f"The load report could not be read: {exc}") from exc

    schema = str(report.get("schema") or "")
    major = schema.split(".", 1)[0]
    if major != str(REPORT_SCHEMA_MAJOR):
        raise ReportError(
            f"The load report is schema {schema or '(none)'}, but this app "
            f"reads schema {REPORT_SCHEMA_MAJOR}.x. The R side and the app are "
            f"from different versions of the repo."
        )
    return report


def inspect_from_report(report: dict) -> dict:
    """The report, reshaped into what the form's callbacks already expect.

    A translation layer, deliberately: the rest of the app was written against
    the keys the old parser produced, and changing R's vocabulary and every
    callback in one step would make both harder to review. The report itself is
    kept under ``"report"`` so the health panel can render the parts that have
    no equivalent in the old shape.
    """
    loaded = [w for w in report.get("workspaces") or []
              if w.get("status") == "loaded"]

    return {
        # Gates every workspace shares — the only ones --population can use.
        "populations": [g["path"] for g in report.get("populations") or []],
        # {path: label}, so a gate FlowJo named "  " is visible in the tree.
        # R decides the label; see gate_display_label().
        "population_labels": {g["path"]: g["label"]
                              for g in report.get("populations") or []},
        "markers": list((report.get("markers") or {}).get("common") or []),
        "samples": [s["display"] for s in report.get("samples") or []],
        "workspaces": [
            {"label": w["label"], "wsp": w.get("wsp") or "",
             "group": w.get("group") or "",
             "samples": f"{w.get('n_loaded', 0)} sample(s)"}
            for w in loaded
        ],
        "populations_by_workspace": {
            w["label"]: [g["path"] for g in w.get("gates") or []]
            for w in loaded
        },
        "gate_labels_by_workspace": {
            w["label"]: {g["path"]: g["label"] for g in w.get("gates") or []}
            for w in loaded
        },
        "markers_by_workspace": {w["label"]: list(w.get("markers") or [])
                                 for w in loaded},
        "annotations": report.get("annotations") or {},
        "report": report,
    }


def alias_string(pairs: dict) -> str:
    """``{written: used}`` as the ``--marker_aliases`` CLI form."""
    return ";".join(f"{k}={v}" for k, v in pairs.items() if k and v)


# --- annotation preview ------------------------------------------------------


def annotation_preview_argv(samples_path: str, out_path: str, metadata: str,
                            pattern: str, id_column: str = "") -> list[str]:
    """Argv for ``--annotation_preview``: no workspace, no FCS, no GatingSet."""
    argv = _base_argv() + [
        "--annotation_preview", out_path,
        "--preview_samples", samples_path,
    ]
    if (metadata or "").strip():
        argv.extend(["--sample_metadata", metadata.strip()])
    elif (pattern or "").strip():
        argv.extend(["--sample_pattern", pattern.strip()])
    if (id_column or "").strip():
        argv.extend(["--sample_id_column", id_column.strip()])
    return argv


def describe_annotations(metadata: str, pattern: str,
                         samples: list[str]) -> tuple[list[str], dict]:
    """``(columns, preview)`` for the current metadata/pattern settings.

    Runs R synchronously — the same choice, and for the same reason, as the FCS
    keyword table: it opens no FCS file and no workspace, so it answers in
    milliseconds and the caller is waiting for it.

    This used to be a Python reimplementation of the matching rules in
    ``src/loading/sample_annotations.R``, which its own docstring described as a
    mirror that R would win against if the two disagreed. Calling the real thing
    removes the question.
    """
    if not (metadata or "").strip() and not (pattern or "").strip():
        return [], {}

    import subprocess
    import tempfile

    with tempfile.TemporaryDirectory(prefix="umap_annot_") as tmp:
        samples_path = str(Path(tmp) / "samples.txt")
        out_path = str(Path(tmp) / "preview.json")
        Path(samples_path).write_text("\n".join(samples), encoding="utf-8")
        argv = annotation_preview_argv(samples_path, out_path, metadata, pattern)
        try:
            proc = subprocess.run(argv, cwd=str(config.REPO_ROOT),
                                  capture_output=True, text=True, timeout=120)
        except (OSError, subprocess.SubprocessError) as exc:
            return [], {"error": f"Could not run R: {exc}"}
        if not Path(out_path).is_file():
            detail = (proc.stderr or proc.stdout or "").strip().splitlines()
            return [], {"error": "R could not resolve the annotations: "
                                 + (detail[-1] if detail else "no output")}
        try:
            preview = json.loads(Path(out_path).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            return [], {"error": f"Could not read the annotation preview: {exc}"}

    if preview.get("error"):
        return [], {"error": preview["error"]}
    columns = list(preview.get("columns") or [])
    # The form wants a level COUNT per column; R sends the levels themselves,
    # which is the more useful thing to have in the document.
    preview["levels"] = {c: len(v or [])
                         for c, v in (preview.get("levels") or {}).items()}
    preview["unmatched"] = list(preview.get("unmatched") or [])
    if preview.get("source"):
        preview["source"] = f"From {preview['source']}"
    return columns, preview
