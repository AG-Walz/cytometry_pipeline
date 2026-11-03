"""Remember the form state between launches.

Everything the user typed is written to a small JSON file so reopening the app
(or an accidental browser refresh) does not mean retyping paths. Only known keys
are restored, and a corrupt or unreadable file is silently ignored — a broken
cache must never stop the app from starting.
"""

from __future__ import annotations

import json

from app.config import CACHE_DIR, DEFAULTS, FORM_START, SESSION_FILE


def load() -> dict:
    """The remembered state, filled in with defaults for anything missing.

    A key never remembered starts at ``FORM_START`` where that has one, else at
    the R default.

    A session saved before the gate layout existed has no "layout" key. It
    opens on ``FORM_START``'s layout, so it also takes that layout's min_dist:
    the remembered one was chosen for the marker layout and would otherwise
    pack the gate islands with it. Every other remembered value is kept.
    """
    state = {**DEFAULTS, **FORM_START}
    try:
        with open(SESSION_FILE, encoding="utf-8") as fh:
            saved = json.load(fh)
    except (OSError, ValueError):
        return state
    if not isinstance(saved, dict):
        return state
    for key, value in saved.items():
        if key in DEFAULTS:
            state[key] = value
    if "layout" not in saved:
        state.update(FORM_START)
    return state


def save(state: dict) -> None:
    """Merge ``state`` into the saved file. Failures are ignored.

    Merging rather than replacing matters: the form is saved continuously, and
    at page load — before the workspace has been read — fields like the chosen
    population are necessarily empty. Omitting them from ``state`` then leaves
    the remembered value intact instead of erasing it.
    """
    payload = load()
    payload.update({k: v for k, v in state.items() if k in DEFAULTS})
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        tmp = SESSION_FILE.with_suffix(".json.tmp")
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)
        tmp.replace(SESSION_FILE)
    except (OSError, TypeError, ValueError):
        pass
