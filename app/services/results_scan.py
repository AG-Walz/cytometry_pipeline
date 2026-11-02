"""Turn an output directory into the structure the results gallery renders.

File names come from ``create_safe_filename()`` on the R side, which composes
them as::

    umap_<pop>_by_file.<ext>
    umap_<pop>_by_workspace.<ext>          (only when >1 workspace was loaded)
    umap_<pop>_by_subpopulation.<ext>
    umap_<pop>_marker_<M>_ALL.<ext>
    umap_<pop>_marker_<M>_file_<sample>.<ext>
    umap_<pop>_by_combination.<ext>        (gate layout)
    umap_<pop>_split_<col>_<level>.<ext>   (--split_by)
    umap_<pop>_split_<col>_panel.<ext>
    umap_<pop>_diff_<col>_<level>.<ext>

and then replaces every character outside ``[A-Za-z0-9_.-]`` with ``_``. That
sanitising means the marker and sample names recovered here are the *sanitised*
ones — good enough to label and filter by, but not identical to the originals.
Anything that does not match falls into "other" rather than being dropped, so
the gallery can never quietly hide a file the run produced.

The split plots are the exception to reading names: ``<col>_<level>`` cannot be
split back apart when a column name has an underscore, and a long name is
truncated to 100 characters. So the run writes ``plot_index.json`` (file stem,
kind, column, level, reference per plot), and when it is there it decides those
plots' kind, column, level and order. The name rules stay as the fallback for a
folder without one.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

from app.config import (
    DATA_EXTENSIONS,
    EXTENSION_PREFERENCE,
    MAX_INLINE_SVG_BYTES,
    PLOT_EXTENSIONS,
    PREVIEW_EXTENSIONS,
)

PLOT_INDEX = "plot_index.json"
#: Folders of a run that hold no results: the resume checkpoints and the
#: app's own run log.
SKIP_DIRS = {"checkpoints", "_run"}
#: Index kinds the gallery files by the index; any other kind (marker plots,
#: the overview plots) is filed by its name as before.
INDEX_KINDS = ("combination", "split", "split_panel", "diff")


@dataclass
class Plot:
    """One logical plot, possibly present in several file formats."""

    key: str
    title: str
    section: str = "other"
    marker: str = ""
    sample: str = ""
    variants: dict[str, Path] = field(default_factory=dict)
    #: Position in plot_index.json, which is the order R drew the levels in.
    #: None when the plot is not in an index.
    rank: int | None = None

    @property
    def preview(self) -> Path | None:
        """The best file to show inline, or None when nothing embeddable exists.

        An SVG over ``MAX_INLINE_SVG_BYTES`` does not count: it carries one
        node per cell and would stall the page, so it is download-only.
        """
        for ext in EXTENSION_PREFERENCE:
            if ext in PREVIEW_EXTENSIONS and ext in self.variants:
                path = self.variants[ext]
                if ext == ".svg" and self._size_of(path) > MAX_INLINE_SVG_BYTES:
                    continue
                return path
        return None

    @property
    def oversized_svg(self) -> bool:
        """True when the only embeddable file is an SVG too large to embed."""
        svg = self.variants.get(".svg")
        return (self.preview is None and svg is not None
                and self._size_of(svg) > MAX_INLINE_SVG_BYTES)

    @staticmethod
    def _size_of(path: Path) -> int:
        try:
            return path.stat().st_size
        except OSError:
            return 0

    @property
    def downloads(self) -> list[tuple[str, Path]]:
        return [
            (ext, self.variants[ext])
            for ext in EXTENSION_PREFERENCE
            if ext in self.variants
        ]

    @property
    def size(self) -> int:
        best = self.preview or next(iter(self.variants.values()), None)
        try:
            return best.stat().st_size if best else 0
        except OSError:
            return 0


@dataclass
class Results:
    directory: Path
    overview: list[Plot] = field(default_factory=list)
    marker_all: dict[str, Plot] = field(default_factory=dict)
    marker_file: dict[str, dict[str, Plot]] = field(default_factory=dict)
    #: Sample-annotation plots. These reuse ``Plot.marker`` for the annotation
    #: COLUMN (patient, timepoint) and ``Plot.sample`` for the LEVEL, since the
    #: file-name scheme is the same shape as the marker one.
    annot_all: dict[str, Plot] = field(default_factory=dict)
    annot_level: dict[str, dict[str, Plot]] = field(default_factory=dict)
    annotations: list[str] = field(default_factory=list)
    #: Metadata splits, per column (``Plot.marker`` = column, ``Plot.sample``
    #: = level, as for annotations): the side-by-side panel, the per-level
    #: plots and the fold-change maps, each list in the order R drew them.
    split_panel: dict[str, Plot] = field(default_factory=dict)
    split_level: dict[str, list[Plot]] = field(default_factory=dict)
    split_diff: dict[str, list[Plot]] = field(default_factory=dict)
    split_columns: list[str] = field(default_factory=list)
    #: CSV/JSON outputs (cells.csv.gz, run_info.json, ...), for download.
    data_files: list[Path] = field(default_factory=list)
    other: list[Plot] = field(default_factory=list)
    markers: list[str] = field(default_factory=list)
    samples: list[str] = field(default_factory=list)
    formats: set[str] = field(default_factory=set)
    n_files: int = 0

    @property
    def is_empty(self) -> bool:
        return self.n_files == 0

    @property
    def has_preview(self) -> bool:
        """False when the run produced only PDFs (nothing embeddable)."""
        return bool(self.formats & set(PREVIEW_EXTENSIONS))


def natural_key(value: str) -> list:
    """Sort key that compares digit runs as numbers ("W5" < "W12" < "W110").

    Same order as natural_order() in src/visualization/umap_gate_layout.R,
    which is the order R draws split levels in unless told otherwise.
    """
    return [("", float(t)) if t.isdigit() else (t.lower(), float("-inf"))
            for t in re.findall(r"[0-9]+|[^0-9]+", str(value))]


def _text(value) -> str:
    """One index field as text; jsonlite may write a scalar as [x], NA as null."""
    if isinstance(value, list):
        value = value[0] if value else None
    return "" if value is None else str(value)


def read_plot_index(directory: str | os.PathLike) -> dict[str, dict]:
    """``{stem: {kind, column, level, reference, rank}}`` from plot_index.json.

    {} when there is no index or it cannot be read; the gallery then files
    everything by name. Accepts what jsonlite writes for a data.frame: an array
    of row objects (NA fields left out, or null), or an object of columns,
    either of them optionally under a "plots" key.
    """
    try:
        doc = json.loads((Path(directory) / PLOT_INDEX).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    if isinstance(doc, dict) and isinstance(doc.get("plots"), (list, dict)):
        doc = doc["plots"]
    if isinstance(doc, dict):
        cols = {k: v for k, v in doc.items() if isinstance(v, list)}
        n = max((len(v) for v in cols.values()), default=0)
        rows = [{k: (v[i] if i < len(v) else None) for k, v in cols.items()}
                for i in range(n)]
    elif isinstance(doc, list):
        rows = [r for r in doc if isinstance(r, dict)]
    else:
        return {}

    index: dict[str, dict] = {}
    for rank, row in enumerate(rows):
        stem = Path(_text(row.get("file"))).name
        if Path(stem).suffix.lower() in PLOT_EXTENSIONS:
            stem = Path(stem).stem
        if not stem or stem in index:
            continue
        entry = {k: _text(row.get(k))
                 for k in ("kind", "column", "level", "reference")}
        entry["rank"] = rank
        index[stem] = entry
    return index


def _split_token(stem: str) -> str:
    """"_split_" or "_diff_" when the name says this is a split plot, else "".

    Whichever of the four family tokens comes first after the population
    wins, so a marker or annotation plot whose own marker or level happens to
    contain "_split_" stays with the markers or annotations.
    """
    hits = {t: stem.find(t) for t in ("_split_", "_diff_", "_marker_", "_annot_")}
    hits = {t: i for t, i in hits.items() if i >= 0}
    if not hits:
        return ""
    first = min(hits, key=hits.get)
    return first if first in ("_split_", "_diff_") else ""


def _column_and_level(tail: str, columns: list[str]) -> tuple[str, str]:
    """Split "<col>_<level>" using the columns known from the folder.

    The longest known column that prefixes the tail wins (so
    "prior_peptide_vaccination_yes" is not read as column "prior"); without
    one, the first "_" separates them, which is right for any column name
    without an underscore.
    """
    for column in sorted(columns, key=len, reverse=True):
        if tail.startswith(column + "_"):
            return column, tail[len(column) + 1:]
    column, _, level = tail.partition("_")
    return column, level


def _split_columns(stems: list[str], data_names: list[str]) -> list[str]:
    """Sanitised split column names the folder itself spells out in full."""
    columns = set()
    for stem in stems:
        if _split_token(stem) == "_split_" and stem.endswith("_panel"):
            columns.add(stem.split("_split_", 1)[1][: -len("_panel")])
    for name in data_names:
        m = re.match(r"^combination_frequencies_(.+)\.csv$", name)
        if m:
            columns.add(m.group(1))
    return sorted(c for c in columns if c)


def _classify_indexed(entry: dict) -> tuple[str, str, str, str] | None:
    """``(section, column, level, title)`` from a plot_index.json entry."""
    kind, column = entry.get("kind"), entry.get("column") or ""
    level, reference = entry.get("level") or "", entry.get("reference") or ""
    if kind == "combination":
        return "combination", "", "", "Coloured by gate combination"
    if not column:
        return None
    if kind == "split_panel":
        return "split_panel", column, "", f"{column}: all levels"
    if kind == "split":
        return "split_level", column, level, f"{column}: {level}"
    if kind == "diff":
        against = f" vs {reference}" if reference else ", change"
        return "split_diff", column, level, f"{column}: {level}{against}"
    return None


def _classify(stem: str, entry: dict | None = None,
              split_columns: list[str] | None = None) -> tuple[str, str, str, str]:
    """Map a file stem to ``(section, marker, sample, title)``.

    ``entry`` is the stem's plot_index.json row, which decides when it names
    one of INDEX_KINDS. The split and combination names are tested before the
    older suffix rules: a level called "by_file" must not make a split plot an
    overview.
    """
    if entry and entry.get("kind") in INDEX_KINDS:
        indexed = _classify_indexed(entry)
        if indexed:
            return indexed
    if stem.endswith("_by_combination"):
        return "combination", "", "", "Coloured by gate combination"
    token = _split_token(stem)
    if token:
        tail = stem.split(token, 1)[1]
        if token == "_split_" and tail.endswith("_panel"):
            column = tail[: -len("_panel")]
            return "split_panel", column, "", f"{column}: all levels"
        column, level = _column_and_level(tail, split_columns or [])
        if token == "_split_":
            return "split_level", column, level, f"{column}: {level}"
        return ("split_diff", column, level,
                f"{column}: {level}, change against the reference")

    if stem.endswith("_by_file"):
        return "overview", "", "", "Coloured by source file"
    if stem.endswith("_by_workspace"):
        return "overview", "", "", "Coloured by FlowJo workspace"
    if stem.endswith("_by_subpopulation"):
        return "overview", "", "", "Coloured by FlowJo subpopulation"

    if "_marker_" in stem:
        tail = stem.split("_marker_", 1)[1]
        if "_file_" in tail:
            # rsplit: the sample name is the last segment, so a marker
            # containing "_file_" still parses correctly.
            marker, sample = tail.rsplit("_file_", 1)
            return "marker_file", marker, sample, f"{marker} — {sample}"
        if tail.endswith("_ALL"):
            marker = tail[: -len("_ALL")]
            return "marker_all", marker, "", f"{marker} — all samples"

    # Sample annotations (patient, timepoint, …) mirror the marker scheme, so
    # the same rsplit reasoning applies: the level is the last segment.
    if "_annot_" in stem:
        tail = stem.split("_annot_", 1)[1]
        if "_level_" in tail:
            column, level = tail.rsplit("_level_", 1)
            return "annot_level", column, level, f"{column}: {level}"
        if tail.endswith("_ALL"):
            column = tail[: -len("_ALL")]
            return "annot_all", column, "", f"Coloured by {column}"

    return "other", "", "", stem


def _run_files(directory: Path) -> list[Path]:
    """Every file of a run, in its subfolders too (overview/, splits/<col>/,
    data/, ...; runs before 2026-10-15 were flat), minus SKIP_DIRS."""
    found: list[Path] = []
    for root, dirs, files in os.walk(directory):
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS)
        found.extend(Path(root) / f for f in files)
    return sorted(found, key=lambda p: p.name.lower())


def scan(output_dir: str | os.PathLike) -> Results:
    """Group every plot file of the run in ``output_dir`` for display, and
    list its data files. Plots are identified by file name wherever they sit:
    the names are unique within a run."""
    directory = Path(output_dir)
    results = Results(directory=directory)
    if not directory.is_dir():
        return results

    # stem -> Plot, so the three formats of one figure collapse into one card.
    plots: dict[str, Plot] = {}
    order: list[str] = []

    try:
        entries = _run_files(directory)
    except OSError:
        return results

    plot_files: list[Path] = []
    for path in entries:
        if not path.is_file():
            continue
        ext = path.suffix.lower()
        if ext in PLOT_EXTENSIONS:
            plot_files.append(path)
        elif ext in DATA_EXTENSIONS:
            results.data_files.append(path)

    index = read_plot_index(directory)
    columns = _split_columns([p.stem for p in plot_files],
                             [p.name for p in results.data_files])
    for path in plot_files:
        ext = path.suffix.lower()
        results.n_files += 1
        results.formats.add(ext)
        stem = path.stem
        if stem not in plots:
            entry = index.get(stem)
            section, marker, sample, title = _classify(stem, entry, columns)
            plots[stem] = Plot(
                key=stem, title=title, section=section, marker=marker, sample=sample,
                rank=entry["rank"] if entry else None,
            )
            order.append(stem)
        plots[stem].variants[ext] = path

    for stem in order:
        plot = plots[stem]
        section = plot.section
        if section in ("overview", "combination"):
            results.overview.append(plot)
        elif section == "split_panel":
            results.split_panel[plot.marker] = plot
        elif section == "split_level":
            results.split_level.setdefault(plot.marker, []).append(plot)
        elif section == "split_diff":
            results.split_diff.setdefault(plot.marker, []).append(plot)
        elif section == "marker_all":
            results.marker_all[plot.marker] = plot
        elif section == "marker_file":
            results.marker_file.setdefault(plot.marker, {})[plot.sample] = plot
        elif section == "annot_all":
            results.annot_all[plot.marker] = plot
        elif section == "annot_level":
            results.annot_level.setdefault(plot.marker, {})[plot.sample] = plot
        else:
            results.other.append(plot)

    results.annotations = sorted(set(results.annot_all) | set(results.annot_level))
    results.markers = sorted(set(results.marker_all) | set(results.marker_file))
    samples: set[str] = set()
    for per_sample in results.marker_file.values():
        samples.update(per_sample)
    results.samples = sorted(samples)

    # Split plots in the order R drew them (the index order, which follows
    # --level_order); without an index, natural order of the level.
    def _split_rank(plot: Plot) -> tuple:
        return (plot.rank is None, plot.rank or 0, natural_key(plot.sample))

    for per_column in (results.split_level, results.split_diff):
        for plots_ in per_column.values():
            plots_.sort(key=_split_rank)
    split_plots = list(results.split_panel.values()) + [
        p for d in (results.split_level, results.split_diff)
        for ps in d.values() for p in ps]
    first_rank: dict[str, tuple] = {}
    for plot in split_plots:
        rank = (plot.rank is None, plot.rank or 0)
        first_rank[plot.marker] = min(first_rank.get(plot.marker, rank), rank)
    results.split_columns = sorted(first_rank,
                                   key=lambda c: (first_rank[c], natural_key(c)))

    # The gate combinations first: in the gate layout that is the figure the
    # run is for. Then sample, workspace, subpopulation, which reads more
    # naturally than alphabetical: from the finest grouping to the coarsest.
    _overview_order = ("_by_file", "_by_workspace", "_by_subpopulation")

    def _overview_rank(plot: Plot) -> int:
        if plot.section == "combination":
            return -1
        for i, suffix in enumerate(_overview_order):
            if plot.key.endswith(suffix):
                return i
        return len(_overview_order)

    results.overview.sort(key=_overview_rank)
    return results


def human_size(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{n:.1f} GB"
