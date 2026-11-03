"""The results gallery: whatever the run actually wrote to the output directory.

A run over a full panel with several samples produces hundreds of files, and
Dash's html.Img exposes no lazy-loading attribute, so the gallery controls the
load itself: per-marker plots live in an accordion whose panels are only built
once opened (Mantine keeps closed panels mounted, so merely collapsing them
would not help), and each section is capped.
"""

from __future__ import annotations

from pathlib import Path
from urllib.parse import quote

import dash_mantine_components as dmc
from dash import ALL, Input, Output, callback, ctx, html, no_update

from app.services import results_scan, serving

#: Cards rendered eagerly per marker section before the browser would suffer.
MAX_CARDS_PER_SECTION = 60


def _download_links(plot: results_scan.Plot):
    return dmc.Group(
        [
            dmc.Anchor(
                ext.lstrip("."),
                href=f"/plot/download?path={quote(str(path))}",
                size="xs",
                target="_blank",
            )
            for ext, path in plot.downloads
        ],
        gap="xs",
    )


def _card(plot: results_scan.Plot):
    preview = plot.preview
    if preview is not None:
        # Dash's html.Img exposes no `loading` attribute, so images cannot be
        # marked lazy here; the collapsed per-marker accordions and
        # MAX_CARDS_PER_SECTION are what keep the number rendered at once sane.
        body = html.Img(
            src=serving.url_for_plot(preview),
            id={"comp": "plot-img", "path": str(preview)},
            n_clicks=0,
            style={
                "width": "100%",
                "height": "auto",
                "cursor": "zoom-in",
                "background": "#fff",
                "borderRadius": "4px",
            },
        )
    elif plot.oversized_svg:
        body = dmc.Alert(
            "SVG too large to preview — download it, or tick PNG next run.",
            color="gray", variant="light",
        )
    else:
        body = dmc.Alert(
            "PDF only — no inline preview.", color="gray", variant="light"
        )

    return dmc.Card(
        [
            body,
            dmc.Space(h=6),
            dmc.Text(plot.title, size="xs", fw=600,
                     style={"wordBreak": "break-word"}),
            dmc.Space(h=4),
            dmc.Group(
                [
                    dmc.Text(results_scan.human_size(plot.size), size="xs",
                             c="dimmed"),
                    _download_links(plot),
                ],
                justify="space-between",
            ),
        ],
        withBorder=True,
        radius="sm",
        p="xs",
    )


def _grid(plots: list[results_scan.Plot]):
    shown = plots[:MAX_CARDS_PER_SECTION]
    children = [dmc.SimpleGrid([_card(p) for p in shown],
                               cols={"base": 1, "sm": 2, "xl": 3}, spacing="sm")]
    if len(plots) > len(shown):
        children.append(
            dmc.Text(
                f"{len(plots) - len(shown)} more file(s) in this section are "
                "not shown — narrow the filters to see them.",
                size="xs",
                c="orange",
                mt=8,
            )
        )
    return html.Div(children)


def _split_panel(results: results_scan.Results, column: str):
    """One split column: all levels side by side, each level, the fold-change maps."""
    levels = results.split_level.get(column) or []
    diffs = results.split_diff.get(column) or []
    panel = results.split_panel.get(column)
    children = []
    if panel is not None:
        children += [_card(panel), dmc.Space(h=10)]
    if levels:
        children += [dmc.Text("Each level", size="xs", c="dimmed", mb=4),
                     _grid(levels), dmc.Space(h=10)]
    if diffs:
        children += [dmc.Text("Change in frequency against the reference level", size="xs",
                              c="dimmed", mb=4),
                     _grid(diffs)]
    return html.Div(children)


def _data_files(results: results_scan.Results):
    """The run's tables, as plain download links."""
    return dmc.Group(
        [
            dmc.Anchor(
                f"{p.name} ({results_scan.human_size(_size(p))})",
                href=f"/plot/download?path={quote(str(p))}",
                size="xs",
                target="_blank",
            )
            for p in results.data_files
        ],
        gap="md",
    )


def _settings_section(results: results_scan.Results) -> list:
    """"Settings used": how the run ended, its code version, and the command
    that repeats it with every option spelt out. [] for runs without a record.
    """
    if not results.command and not results.settings:
        return []
    s = results.settings
    git = s.get("git") or {}
    status = s.get("status", "")
    if status == "running":
        # R writes the outcome when it ends; "running" on a finished folder
        # means the process was killed before it could.
        status = "running (or stopped without finishing)"
    facts = [
        f"Started {s.get('started', '?')} on {(s.get('host') or {}).get('name', '?')}",
        f"status: {status or '?'}"
        + (f", {s['minutes']} min" if s.get("minutes") is not None else ""),
        f"code: {git.get('describe') or git.get('commit') or '?'}"
        + (" (with local changes)" if git.get("dirty") else ""),
    ]
    children = [
        dmc.Divider(label="Settings used", labelPosition="left"),
        dmc.Space(h=8),
        dmc.Text(" · ".join(facts), size="xs", c="dimmed"),
    ]
    if s.get("error"):
        children.append(dmc.Text(f"Error: {s['error']}", size="xs", c="red"))
    if results.command:
        children += [
            dmc.Space(h=6),
            dmc.Text("Repeat this run: bash _run/command.sh (every option is "
                     "written out, defaults included).", size="xs"),
            dmc.Code(results.command, block=True,
                     style={"fontSize": "11px", "maxHeight": "260px",
                            "overflow": "auto"}),
        ]
    children += [
        dmc.Space(h=6),
        dmc.Group(
            [dmc.Anchor(p.name, href=f"/plot/download?path={quote(str(p))}",
                        size="xs", target="_blank")
             for p in results.settings_files],
            gap="md"),
        dmc.Space(h=16),
    ]
    return children


def _data_section(results: results_scan.Results) -> list:
    """"Settings used", then the "Data files" heading and links; [] when there
    is neither."""
    settings = _settings_section(results)
    if not results.data_files:
        return settings
    return settings + [
        dmc.Divider(label="Data files", labelPosition="left"),
        dmc.Space(h=8),
        _data_files(results),
        dmc.Space(h=16),
    ]


def _size(path: Path) -> int:
    try:
        return path.stat().st_size
    except OSError:
        return 0


def register_callbacks() -> None:
    @callback(
        Output("gallery", "children"),
        Output("split-header", "children"),
        Output("split-accordion", "children"),
        Output("gallery-rest", "children"),
        Output("per-marker-header", "children"),
        Output("marker-accordion", "children"),
        Output("annot-accordion", "children"),
        Output("results-path", "children"),
        Output("filter-marker", "data"),
        Output("filter-sample", "data"),
        Input({"comp": "fb-input", "bid": "results"}, "value"),
        Input("store-run-generation", "data"),
        Input("btn-refresh-results", "n_clicks"),
        Input("filter-marker", "value"),
        Input("filter-sample", "value"),
        Input("marker-accordion", "value"),
        Input("annot-accordion", "value"),
        Input("split-accordion", "value"),
    )
    def render_gallery(directory, _generation, _refresh, marker_filter,
                       sample_filter, open_markers, open_annotations,
                       open_splits):
        """Render the gallery for ``directory``.

        The accordions' ``value`` props are inputs so that only the panels the
        user has actually opened are built. Writing a component's children
        while reading its value is not a loop: setting children never changes
        value.
        """
        open_markers = set(open_markers or [])
        open_annotations = set(open_annotations or [])
        open_splits = set(open_splits or [])
        directory = (directory or "").strip()
        if not directory:
            return (
                dmc.Alert(
                    "Run something, or pick a results folder above.",
                    color="gray",
                    variant="light",
                ),
                None,
                [],
                None,
                None,
                [],
                [],
                "No results yet.",
                [],
                [],
            )

        serving.allow(directory)
        results = results_scan.scan(directory)

        path_text = dmc.Text(f"{results.n_files} file(s)", size="sm",
                             c="dimmed")

        marker_options = [{"value": m, "label": m} for m in results.markers]
        sample_options = [{"value": s, "label": s} for s in results.samples]

        if results.is_empty:
            # The gate layout writes its tables before any plot, so a run
            # whose plotting failed leaves a folder of data files only.
            if not Path(directory).is_dir():
                note = "This folder does not exist."
            elif results.data_files:
                note = "No plot files in this folder; its data files are below."
            else:
                note = "No plot files in this folder."
            return (
                dmc.Alert(note, color="gray", variant="light"),
                None,
                [],
                html.Div(_data_section(results)),
                None,
                [],
                [],
                path_text,
                marker_options,
                sample_options,
            )

        sections: list = []

        # A map that failed is said above the plots, where it is seen.
        if results.embedding.get("suspect"):
            reasons = results.embedding.get("reasons") or []
            if isinstance(reasons, str):
                reasons = [reasons]
            sections.append(dmc.Alert(
                dmc.Stack([dmc.Text(r, size="sm") for r in reasons], gap=4),
                title="This embedding looks wrong",
                color="orange", variant="light"))
            sections.append(dmc.Space(h=10))

        if not results.has_preview:
            sections.append(
                dmc.Alert(
                    "This run produced PDFs only, which cannot be embedded. "
                    "Use the download links, or tick PNG in step 3 next time.",
                    color="blue",
                    variant="light",
                )
            )
            sections.append(dmc.Space(h=10))

        if results.overview:
            sections.append(dmc.Divider(label="Overview", labelPosition="left"))
            sections.append(dmc.Space(h=8))
            sections.append(_grid(results.overview))
            sections.append(dmc.Space(h=16))

        # Split by metadata, straight after the overview: one collapsed panel
        # per column, built only once opened (a column can have dozens of
        # levels, and each level has a plot and a fold-change map). Hidden by a
        # sample filter, like the annotation plots: neither is per sample.
        split_items = []
        if not sample_filter:
            for column in results.split_columns:
                n_levels = len(results.split_level.get(column) or [])
                n_diffs = len(results.split_diff.get(column) or [])
                label = f"{column} — {n_levels} level(s)" + (
                    f", {n_diffs} fold-change map(s)" if n_diffs else "")
                split_items.append(
                    dmc.AccordionItem(
                        [
                            dmc.AccordionControl(label),
                            dmc.AccordionPanel(
                                _split_panel(results, column)
                                if column in open_splits
                                else dmc.Text("Opening…", size="xs", c="dimmed")
                            ),
                        ],
                        value=column,
                    )
                )
        split_header = (
            html.Div([dmc.Divider(label="Split by metadata",
                                  labelPosition="left"),
                      dmc.Space(h=8)])
            if split_items else None
        )
        gallery_top, sections = sections, []

        # Sample annotations (patient, timepoint, …): the "coloured by" plot for
        # every column, then each column's levels in its own collapsed panel.
        # A cohort can have dozens of patients, so the levels are not rendered
        # until you open them — same reasoning as the per-marker section.
        if results.annotations and not sample_filter:
            annot_overviews = [results.annot_all[c] for c in results.annotations
                               if c in results.annot_all]
            if annot_overviews:
                sections.append(
                    dmc.Divider(label="Sample annotations",
                                labelPosition="left")
                )
                sections.append(dmc.Space(h=8))
                sections.append(_grid(annot_overviews))
                sections.append(dmc.Space(h=12))


        # One collapsed panel per annotation column holding its level plots.
        # Built lazily for the same reason as the per-marker panels: a cohort
        # of 33 patients is 33 full-size images, and Mantine keeps closed
        # panels mounted, so filling them all would load every one at once.
        annot_items = []
        if not sample_filter:
            for column in results.annotations:
                levels = results.annot_level.get(column) or {}
                if not levels:
                    continue
                panel = (
                    _grid([levels[k] for k in sorted(levels)])
                    if column in open_annotations
                    else dmc.Text("Opening…", size="xs", c="dimmed")
                )
                annot_items.append(
                    dmc.AccordionItem(
                        [
                            dmc.AccordionControl(
                                f"{column} — {len(levels)} level(s), each "
                                "highlighted against the whole embedding"
                            ),
                            dmc.AccordionPanel(panel),
                        ],
                        value=column,
                    )
                )

        markers = [m for m in results.markers
                   if not marker_filter or m == marker_filter]

        all_sample_plots = [results.marker_all[m] for m in markers
                            if m in results.marker_all]
        if all_sample_plots and not sample_filter:
            sections.append(
                dmc.Divider(label="Markers — all samples", labelPosition="left")
            )
            sections.append(dmc.Space(h=8))
            sections.append(_grid(all_sample_plots))
            sections.append(dmc.Space(h=16))

        per_marker_items = []
        for marker in markers:
            per_sample = results.marker_file.get(marker) or {}
            plots = [
                plot
                for sample, plot in sorted(per_sample.items())
                if not sample_filter or sample == sample_filter
            ]
            if not plots:
                continue
            # Build the images only for panels the user has opened. Mantine
            # keeps closed panels mounted, so filling them all would load every
            # marker x sample plot at once — hundreds of images for a full
            # panel across several samples.
            panel = (
                _grid(plots)
                if marker in open_markers
                else dmc.Text("Opening…", size="xs", c="dimmed")
            )
            per_marker_items.append(
                dmc.AccordionItem(
                    [
                        dmc.AccordionControl(
                            f"{marker} — {len(plots)} sample plot(s)"
                        ),
                        dmc.AccordionPanel(panel),
                    ],
                    value=marker,
                )
            )

        if results.other:
            sections.append(dmc.Divider(label="Other files", labelPosition="left"))
            sections.append(dmc.Space(h=8))
            sections.append(_grid(results.other))
            sections.append(dmc.Space(h=16))

        sections += _data_section(results)

        # The per-marker accordion lives outside #gallery so its open set can
        # drive which panels are built; only its header travels with the rest.
        header = (
            html.Div(
                [
                    dmc.Divider(label="Markers — per sample", labelPosition="left"),
                    dmc.Space(h=8),
                ]
            )
            if per_marker_items
            else None
        )

        return (
            html.Div(gallery_top),
            split_header,
            split_items,
            html.Div(sections),
            header,
            per_marker_items,
            annot_items,
            path_text,
            marker_options,
            sample_options,
        )

    @callback(
        Output("lightbox", "opened"),
        Output("lightbox-body", "children"),
        Output("lightbox", "title"),
        Input({"comp": "plot-img", "path": ALL}, "n_clicks"),
        prevent_initial_call=True,
    )
    def open_lightbox(clicks):
        # Dash fires this whenever the set of images changes, with all-zero
        # click counts; only a real click should open the modal.
        if not clicks or not any(clicks):
            return no_update, no_update, no_update
        triggered = ctx.triggered_id
        if not isinstance(triggered, dict):
            return no_update, no_update, no_update
        path = triggered.get("path", "")
        return (
            True,
            html.Img(
                src=serving.url_for_plot(path),
                style={"width": "100%", "height": "auto", "background": "#fff"},
            ),
            Path(path).name,
        )
