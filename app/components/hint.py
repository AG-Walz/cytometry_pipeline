"""Explanations on hover, not on the page.

Every setting used to carry its explanation as a grey paragraph under its
label, so the form was mostly text and the settings were hard to find. The
explanation now sits behind a small info icon next to the label and opens on
hover, keyboard focus or touch -- the Mantine pattern for field help (a
Tooltip on an icon inside the label; Input labels take components).
"""

from __future__ import annotations

import dash_mantine_components as dmc
from dash import html
from dash_iconify import DashIconify

#: Width of the explanation bubble, in px.
HINT_WIDTH = 320


def hint_icon(text, width: int = HINT_WIDTH):
    """The info icon with ``text`` in a tooltip (None when there is no text)."""
    if not text:
        return None
    return dmc.Tooltip(
        # A focusable span, so the explanation is reachable by keyboard too.
        html.Span(
            DashIconify(icon="tabler:info-circle", width=15,
                        color="var(--mantine-color-dimmed)"),
            tabIndex=0,
            style={"display": "inline-flex", "cursor": "help",
                   "verticalAlign": "middle"},
        ),
        label=text,
        multiline=True,
        w=width,
        withArrow=True,
        openDelay=150,
        position="top-start",
        events={"hover": True, "focus": True, "touch": True},
        zIndex=4000,
    )


def hint(label, text, width: int = HINT_WIDTH):
    """``label`` followed by an info icon that explains it on hover.

    Returns ``label`` unchanged when there is no text, so it can wrap any
    label unconditionally.
    """
    icon = hint_icon(text, width)
    if icon is None:
        return label
    return html.Span(
        [label, icon],
        style={"display": "inline-flex", "alignItems": "center", "gap": "4px"},
    )
