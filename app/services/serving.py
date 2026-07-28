"""Serve plot files that live outside Dash's ``assets`` directory.

Output directories are chosen at runtime, so their files cannot be static
assets. A small Flask route streams them instead, restricted to directories the
app has actually displayed. The app binds to localhost, but a route that will
hand back any path on request is still worth not writing.
"""

from __future__ import annotations

import mimetypes
from pathlib import Path
from urllib.parse import quote

from flask import Response, abort, request, send_file

from app.config import PLOT_EXTENSIONS

#: Directories whose plot files may be served (populated as galleries render).
_ALLOWED: set[Path] = set()


def allow(directory: str | Path) -> None:
    try:
        _ALLOWED.add(Path(directory).resolve())
    except OSError:
        pass


def _is_allowed(path: Path) -> bool:
    try:
        resolved = path.resolve()
    except OSError:
        return False
    if resolved.suffix.lower() not in PLOT_EXTENSIONS:
        return False
    return any(resolved.parent == allowed for allowed in _ALLOWED)


def url_for_plot(path: str | Path) -> str:
    """A cache-busting URL for ``path`` (its mtime is the version)."""
    p = Path(path)
    try:
        version = int(p.stat().st_mtime)
    except OSError:
        version = 0
    return f"/plot?path={quote(str(p))}&v={version}"


def register_route(server) -> None:
    @server.route("/plot")
    def serve_plot():  # noqa: D401 - Flask view
        raw = request.args.get("path", "")
        if not raw:
            abort(400)
        path = Path(raw)
        if not _is_allowed(path) or not path.is_file():
            abort(404)
        mime, _ = mimetypes.guess_type(path.name)
        # SVGs are served inline so <img> can render them rather than
        # triggering a download.
        return send_file(path, mimetype=mime or "application/octet-stream")

    @server.route("/plot/download")
    def download_plot():  # noqa: D401 - Flask view
        raw = request.args.get("path", "")
        path = Path(raw)
        if not raw or not _is_allowed(path) or not path.is_file():
            abort(404)
        return send_file(path, as_attachment=True, download_name=path.name)


__all__ = ["allow", "url_for_plot", "register_route", "Response"]
