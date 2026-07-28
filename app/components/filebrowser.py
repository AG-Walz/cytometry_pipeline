"""A reusable path field: a text input plus a dmc.Tree browser in a modal.

Used four times (workspace file, FCS directory, output parent, results
directory), so every part is addressed by a pattern-matching id keyed on
``bid``. The tree is rooted at ``/`` and loads each directory the first time it
is expanded.

All the modal's state transitions — open, expand, jump-to-path, confirm,
cancel — are handled by one callback. That is deliberate: several of them write
the same tree props, and routing them through a single ``ctx.triggered_id``
switch keeps that honest instead of scattering ``allow_duplicate`` outputs.
"""

from __future__ import annotations

from pathlib import Path

import dash_mantine_components as dmc
from dash import (
    MATCH,
    Input,
    Output,
    State,
    callback,
    clientside_callback,
    ctx,
    dcc,
    html,
    no_update,
)

from app.services import fs_tree


def _id(comp: str, bid: str | dict) -> dict:
    return {"comp": comp, "bid": bid}


def browser_field(
    bid: str,
    label: str,
    description: str,
    mode: str,
    extensions: tuple[str, ...] | None = None,
    placeholder: str = "",
    value: str = "",
) -> html.Div:
    """A labelled path input with a Browse button and its modal.

    ``mode`` is ``"file"`` or ``"dir"``; ``extensions`` filters the files shown
    in file mode (directories are always listed so you can navigate through
    them).
    """
    dirs_only = mode == "dir"
    what = "directory" if dirs_only else "file"

    return html.Div(
        [
            dcc.Store(
                id=_id("fb-cfg", bid),
                data={"dirs_only": dirs_only, "extensions": list(extensions or [])},
            ),
            # Holds the node the tree should be scrolled to after opening or
            # jumping; a clientside callback does the scrolling.
            dcc.Store(id=_id("fb-scroll", bid), data=""),
            html.Div(id=_id("fb-scroll-sink", bid), style={"display": "none"}),
            dmc.Group(
                [
                    dmc.TextInput(
                        id=_id("fb-input", bid),
                        label=label,
                        description=description,
                        placeholder=placeholder,
                        value=value,
                        style={"flex": 1},
                    ),
                    dmc.Button(
                        "Browse",
                        id=_id("fb-open", bid),
                        variant="default",
                        mt=28 if description else 24,
                    ),
                ],
                align="flex-start",
                gap="xs",
                wrap="nowrap",
            ),
            dmc.Modal(
                id=_id("fb-modal", bid),
                title=f"Select a {what} — {label}",
                size="xl",
                opened=False,
                zIndex=2000,
                children=[
                    dmc.Group(
                        [
                            dmc.TextInput(
                                id=_id("fb-goto", bid),
                                placeholder="Jump to a path, e.g. /mnt/volume/data",
                                style={"flex": 1},
                            ),
                            dmc.Button("Go", id=_id("fb-goto-btn", bid),
                                       variant="light"),
                            dmc.Switch(
                                id=_id("fb-hidden", bid),
                                label="Hidden",
                                checked=False,
                                size="sm",
                            ),
                        ],
                        align="center",
                        gap="xs",
                        wrap="nowrap",
                    ),
                    dmc.Space(h=10),
                    dmc.ScrollArea(
                        dmc.Tree(
                            id=_id("fb-tree", bid),
                            data=fs_tree.root_data(),
                            expanded=[],
                            selected=[],
                            selectOnClick=True,
                            expandOnClick=True,
                            allowRangeSelection=False,
                            renderNode={"function": "fsNode"},
                        ),
                        h=420,
                        type="auto",
                        offsetScrollbars=True,
                        style={"border": "1px solid var(--mantine-color-default-border)",
                               "borderRadius": "var(--mantine-radius-sm)",
                               "padding": "8px"},
                    ),
                    dmc.Space(h=10),
                    dmc.Text(
                        f"No {what} selected.",
                        id=_id("fb-selected", bid),
                        size="sm",
                        c="dimmed",
                        style={"wordBreak": "break-all"},
                    ),
                    dmc.Space(h=10),
                    dmc.Group(
                        [
                            dmc.Button("Cancel", id=_id("fb-cancel", bid),
                                       variant="subtle", color="gray"),
                            dmc.Button("Use this", id=_id("fb-confirm", bid),
                                       disabled=True),
                        ],
                        justify="flex-end",
                    ),
                ],
            ),
        ]
    )


def register_callbacks() -> None:
    """Wire the browser. Call once at import time of the app module."""

    @callback(
        Output(_id("fb-modal", MATCH), "opened"),
        Output(_id("fb-tree", MATCH), "data"),
        Output(_id("fb-tree", MATCH), "expanded"),
        Output(_id("fb-input", MATCH), "value"),
        Output(_id("fb-scroll", MATCH), "data"),
        Output(_id("fb-tree", MATCH), "selected"),
        Input(_id("fb-open", MATCH), "n_clicks"),
        Input(_id("fb-tree", MATCH), "expanded"),
        Input(_id("fb-goto-btn", MATCH), "n_clicks"),
        Input(_id("fb-goto", MATCH), "n_submit"),
        Input(_id("fb-hidden", MATCH), "checked"),
        Input(_id("fb-confirm", MATCH), "n_clicks"),
        Input(_id("fb-cancel", MATCH), "n_clicks"),
        State(_id("fb-tree", MATCH), "data"),
        State(_id("fb-tree", MATCH), "selected"),
        State(_id("fb-input", MATCH), "value"),
        State(_id("fb-goto", MATCH), "value"),
        State(_id("fb-cfg", MATCH), "data"),
        prevent_initial_call=True,
    )
    def browser_state(
        _open_clicks,
        expanded,
        _goto_clicks,
        _goto_submit,
        show_hidden,
        _confirm_clicks,
        _cancel_clicks,
        data,
        selected,
        current_value,
        goto_value,
        cfg,
    ):
        comp = (ctx.triggered_id or {}).get("comp")
        dirs_only = bool(cfg.get("dirs_only"))
        extensions = tuple(cfg.get("extensions") or ()) or None
        show_hidden = bool(show_hidden)
        keep = (no_update,) * 6

        if comp == "fb-cancel":
            return (False, *keep[1:])

        if comp == "fb-confirm":
            chosen = (selected or [None])[0]
            ok, message = fs_tree.validate_selection(chosen, dirs_only, extensions)
            if not ok:
                return keep
            return False, no_update, no_update, message, no_update, no_update

        if comp in ("fb-open", "fb-goto-btn", "fb-goto"):
            if comp == "fb-open":
                # Open where the field already points, so re-browsing is one
                # click; otherwise start from the working directory.
                target = current_value or str(Path.cwd())
            else:
                target = (goto_value or "").strip()
                if not target:
                    return keep

            tree, expand = fs_tree.expand_to(target, dirs_only, extensions,
                                             show_hidden)
            landed = expand[-1] if expand else "/"
            # In directory mode the place you navigated to is almost always the
            # place you want, so preselect it and let "Use this" be one click.
            preselect = [landed] if dirs_only else []
            opened = True if comp == "fb-open" else no_update
            return opened, tree, expand, no_update, landed, preselect

        if comp == "fb-hidden":
            # Re-list everything currently open under the new hidden-file rule.
            tree = fs_tree.root_data()
            for value in sorted(expanded or [], key=len):
                fs_tree.load_into(tree, value, dirs_only, extensions, show_hidden)
            return no_update, tree, *keep[2:]

        if comp == "fb-tree":
            # Lazy load: fill in whichever newly-expanded nodes are still
            # placeholders. `expanded` is not echoed back, which would loop.
            tree = data or fs_tree.root_data()
            changed = False
            for value in expanded or []:
                if fs_tree.load_into(tree, value, dirs_only, extensions, show_hidden):
                    changed = True
            return no_update, (tree if changed else no_update), *keep[2:]

        return keep

    # The tree renders the whole path from "/" downwards, so after opening or
    # jumping the node of interest is usually well below the fold. Scroll it
    # into view clientside — there is no server-side way to move a scrollbar.
    clientside_callback(
        """
        function(target) {
            if (!target) { return window.dash_clientside.no_update; }
            setTimeout(function () {
                var el = document.querySelector(
                    '[data-value="' + target.replace(/"/g, '\\\\"') + '"]'
                );
                if (el && el.scrollIntoView) {
                    el.scrollIntoView({block: 'center', behavior: 'smooth'});
                }
            }, 200);
            return window.dash_clientside.no_update;
        }
        """,
        Output(_id("fb-scroll-sink", MATCH), "children"),
        Input(_id("fb-scroll", MATCH), "data"),
        prevent_initial_call=True,
    )

    @callback(
        Output(_id("fb-selected", MATCH), "children"),
        Output(_id("fb-selected", MATCH), "c"),
        Output(_id("fb-confirm", MATCH), "disabled"),
        Input(_id("fb-tree", MATCH), "selected"),
        State(_id("fb-cfg", MATCH), "data"),
        prevent_initial_call=True,
    )
    def selection_preview(selected, cfg):
        dirs_only = bool(cfg.get("dirs_only"))
        extensions = tuple(cfg.get("extensions") or ()) or None
        chosen = (selected or [None])[0]
        ok, message = fs_tree.validate_selection(chosen, dirs_only, extensions)
        return message, ("dimmed" if not ok else None), not ok
