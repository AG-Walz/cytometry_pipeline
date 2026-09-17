"""Turn an output directory into the structure the results gallery renders.

File names come from ``create_safe_filename()`` on the R side, which composes
them as::

    umap_<pop>_by_file.<ext>
    umap_<pop>_by_workspace.<ext>          (only when >1 workspace was loaded)
    umap_<pop>_by_subpopulation.<ext>
    umap_<pop>_marker_<M>_ALL.<ext>
    umap_<pop>_marker_<M>_file_<sample>.<ext>

and then replaces every character outside ``[A-Za-z0-9_.-]`` with ``_``. That
sanitising means the marker and sample names recovered here are the *sanitised*
ones — good enough to label and filter by, but not identical to the originals.
Anything that does not match falls into "other" rather than being dropped, so
the gallery can never quietly hide a file the run produced.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

from app.config import (
    EXTENSION_PREFERENCE,
    MAX_INLINE_SVG_BYTES,
    PLOT_EXTENSIONS,
    PREVIEW_EXTENSIONS,
)


@dataclass
class Plot:
    """One logical plot, possibly present in several file formats."""

    key: str
    title: str
    section: str = "other"
    marker: str = ""
    sample: str = ""
    variants: dict[str, Path] = field(default_factory=dict)

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


def _classify(stem: str) -> tuple[str, str, str, str]:
    """Map a file stem to ``(section, marker, sample, title)``."""
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


def scan(output_dir: str | os.PathLike) -> Results:
    """Group every plot file in ``output_dir`` for display.

    Non-recursive: the R script writes flat into the output directory.
    """
    directory = Path(output_dir)
    results = Results(directory=directory)
    if not directory.is_dir():
        return results

    # stem -> Plot, so the three formats of one figure collapse into one card.
    plots: dict[str, Plot] = {}
    order: list[str] = []

    try:
        entries = sorted(directory.iterdir(), key=lambda p: p.name.lower())
    except OSError:
        return results

    for path in entries:
        if not path.is_file():
            continue
        ext = path.suffix.lower()
        if ext not in PLOT_EXTENSIONS:
            continue
        results.n_files += 1
        results.formats.add(ext)
        stem = path.stem
        if stem not in plots:
            section, marker, sample, title = _classify(stem)
            plots[stem] = Plot(
                key=stem, title=title, section=section, marker=marker, sample=sample
            )
            order.append(stem)
        plots[stem].variants[ext] = path

    for stem in order:
        plot = plots[stem]
        section = plot.section
        if section == "overview":
            results.overview.append(plot)
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

    # Sample, then workspace, then subpopulation reads more naturally than
    # alphabetical — it goes from the finest grouping to the coarsest.
    _overview_order = ("_by_file", "_by_workspace", "_by_subpopulation")

    def _overview_rank(plot: Plot) -> int:
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
