#!/usr/bin/env python
"""Dash front-end for scripts/umap_from_flowjo.R.

Run from the repository root with the pipeline environment active::

    conda activate cytometry_env
    python app/app.py

The app builds and executes the same command line you would type by hand — it
shows that command before running it — and reads back whatever the script wrote.
It binds to 127.0.0.1 by default: the file browser is rooted at ``/`` and there
is no authentication, so exposing it on a network would hand anyone who can
reach it a filesystem browser and a job runner.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow `python app/app.py` from anywhere by making the repo root importable.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import dash_mantine_components as dmc  # noqa: E402
from dash import Dash  # noqa: E402

from app import config, layout  # noqa: E402
from app.callbacks import execution, fcsnames, results, steps  # noqa: E402
from app.components import filebrowser  # noqa: E402
from app.services import serving, session_state  # noqa: E402


def create_app() -> Dash:
    app = Dash(
        __name__,
        title="FlowJo → UMAP",
        assets_folder=str(config.APP_DIR / "assets"),
        # Components inside the stepper's inactive panes are not in the initial
        # DOM, and the pattern-matching browser ids are generated per field.
        suppress_callback_exceptions=True,
        external_stylesheets=dmc.styles.ALL,
    )
    # A callable layout is re-evaluated on every page load, so the remembered
    # form state is restored on a browser refresh too — not only when the
    # server happens to start.
    app.layout = lambda: layout.build_layout(session_state.load())

    filebrowser.register_callbacks()
    steps.register_callbacks()
    fcsnames.register_callbacks()
    execution.register_callbacks()
    results.register_callbacks()
    serving.register_route(app.server)
    return app


#: Built at import so a WSGI server can find it, and reused by main().
APP = create_app()

#: WSGI entry point, e.g. ``gunicorn --workers 1 --threads 8 app.app:server``.
#:
#: ONE worker only. The active run lives in a module-level object
#: (services.runner.RUNNER), so a second worker would own a second, unrelated
#: runner: the browser would poll whichever process the request landed in and
#: see an idle app while a run was in progress elsewhere. Threads, not workers,
#: are what this app scales with — and it is a single-user tool regardless.
server = APP.server


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1",
                        help="interface to bind (default: %(default)s)")
    parser.add_argument("--port", type=int, default=8050,
                        help="port to bind (default: %(default)s)")
    parser.add_argument("--debug", action="store_true",
                        help="enable Dash hot-reload and error pages")
    args = parser.parse_args()

    if not config.CLI_SCRIPT.is_file():
        parser.error(
            f"Cannot find {config.CLI_SCRIPT}. Run the app from a checkout of "
            "the cytometry_pipeline repository."
        )

    print(f"FlowJo → UMAP app on http://{args.host}:{args.port}")
    APP.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
