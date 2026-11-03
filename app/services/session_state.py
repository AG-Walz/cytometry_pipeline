"""Remember the form state between launches.

Everything the user typed is written to a small JSON file so reopening the app
(or an accidental browser refresh) does not mean retyping paths. Only known keys
are restored, and a corrupt or unreadable file is silently ignored — a broken
cache must never stop the app from starting.
"""

from __future__ import annotations

import json

from app.config import (CACHE_DIR, DEFAULTS, FORM_START, FORM_START_VERSION,
                        SESSION_FILE)


def load() -> dict:
    """The remembered state, filled in with defaults for anything missing.

    A key never remembered starts at ``FORM_START`` where that has one, else at
    the R default.

    A session with no "layout" key, or saved under an older FORM_START
    (``form_start_version``), opens on ``FORM_START``'s layout and that
    layout's min_dist. Every other remembered value is kept.
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
    # A session saved under an older FORM_START remembers that start, not a
    # choice: until version 2 the form opened on the gate layout, and nearly
    # every remembered "gates" is that default. Those keys start over once;
    # everything else remembered is kept.
    version = saved.get("form_start_version") or 1
    if "layout" not in saved or version < 2:
        state.update(FORM_START)
    if version < 3:
        for key, old_default in (("point_size", 0.4), ("point_alpha", 0.6)):
            if state.get(key) == old_default:
                state[key] = None
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
    payload["form_start_version"] = FORM_START_VERSION
    try:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        tmp = SESSION_FILE.with_suffix(".json.tmp")
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)
        tmp.replace(SESSION_FILE)
    except (OSError, TypeError, ValueError):
        pass
