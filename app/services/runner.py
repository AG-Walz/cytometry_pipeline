"""A single-slot subprocess runner with a live, pollable log.

The app is deliberately one-run-at-a-time, so one module-level ``Runner``
instance owns whichever R process is active — an inspect or a full run. Both
stream into the same log, and neither can start while the other is alive, which
is exactly the guarantee the UI promises.

Output is read on a background thread into a bounded buffer; the UI polls
``snapshot()`` from a ``dcc.Interval``. Nothing here is safe for multi-user
serving, and it is not meant to be: the app binds to localhost.
"""

from __future__ import annotations

import os
import re
import signal
import subprocess
import threading
import time

#: Hard cap on retained log lines, so a pathological run cannot exhaust memory.
MAX_LINES = 20000

#: Progress lines emitted by ``emit_progress()`` in src/R_shared/io_utils.R.
#: They are intercepted here rather than shown: they are for the progress bar,
#: and one per sample would bury the log they share a stream with. R owns the
#: counting; this only reads it. Keep in step with PROGRESS_PREFIX there.
_PROGRESS_RE = re.compile(
    r"^::progress:: phase=(\S+) done=(\d+) total=(\d+) label=(.*)$"
)

#: How each phase reads in the UI, in the order they occur.
PROGRESS_PHASES = {
    "load": "Building gating sets",
    "extract": "Extracting cells",
    "embed": "Computing the embedding",
    "plot": "Drawing plots",
}

IDLE = "idle"
RUNNING = "running"
DONE = "done"
FAILED = "failed"
CANCELLED = "cancelled"


class Runner:
    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._proc: subprocess.Popen | None = None
        self._thread: threading.Thread | None = None
        self._lines: list[str] = []
        self._truncated = False
        self._status = IDLE
        self._kind = ""
        self._argv: list[str] = []
        self._meta: dict = {}
        self._returncode: int | None = None
        self._started_at: float | None = None
        self._finished_at: float | None = None
        self._cancelling = False
        #: Bumped on every completion so the UI can act once per finished run.
        self._generation = 0
        #: Latest ::progress:: line, and when its phase began (for the ETA).
        self._progress: dict = {}
        self._phase_started_at: float | None = None

    # -- lifecycle ---------------------------------------------------------

    def start(self, argv: list[str], cwd: str, kind: str, meta: dict | None = None):
        """Launch ``argv``. Returns ``(ok, message)``; refuses while busy."""
        with self._lock:
            if self._status == RUNNING:
                return False, "A run is already in progress."
            self._lines = []
            self._truncated = False
            self._status = RUNNING
            self._kind = kind
            self._argv = list(argv)
            self._meta = dict(meta or {})
            self._returncode = None
            self._progress = {}
            self._phase_started_at = None
            self._started_at = time.time()
            self._finished_at = None
            self._cancelling = False

            try:
                # start_new_session puts the child in its own process group so
                # cancelling kills R and anything it spawned, not just the
                # stdbuf wrapper.
                self._proc = subprocess.Popen(
                    argv,
                    cwd=cwd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                    start_new_session=True,
                )
            except (OSError, ValueError) as exc:
                self._status = FAILED
                self._finished_at = time.time()
                self._generation += 1
                self._lines = [f"Failed to launch: {exc}", "", "Command was:"]
                self._lines += ["  " + " ".join(argv)]
                return False, str(exc)

            self._thread = threading.Thread(
                target=self._pump, args=(self._proc,), daemon=True
            )
            self._thread.start()
            return True, "started"

    def _pump(self, proc: subprocess.Popen) -> None:
        """Drain the child's output, then record its exit status."""
        assert proc.stdout is not None
        try:
            for line in proc.stdout:
                line = line.rstrip("\n")
                hit = _PROGRESS_RE.match(line)
                if hit:
                    phase, done, total, label = hit.groups()
                    with self._lock:
                        # Each phase restarts the clock: they are different work
                        # at different speeds, so extrapolating "Drawing plots"
                        # from how long "Building gating sets" took would give a
                        # confidently wrong estimate.
                        if phase != self._progress.get("phase"):
                            self._phase_started_at = time.time()
                        self._progress = {
                            "phase": phase,
                            "done": int(done),
                            "total": int(total),
                            "label": label,
                        }
                    continue
                with self._lock:
                    if len(self._lines) < MAX_LINES:
                        self._lines.append(line)
                    elif not self._truncated:
                        self._truncated = True
                        self._lines.append(
                            f"… output truncated at {MAX_LINES} lines; "
                            "the run continues and files are still being written."
                        )
        except (ValueError, OSError):
            # Pipe closed underneath us (e.g. during a cancel).
            pass
        finally:
            returncode = proc.wait()
            with self._lock:
                self._returncode = returncode
                self._finished_at = time.time()
                if self._cancelling:
                    self._status = CANCELLED
                elif returncode == 0:
                    self._status = DONE
                else:
                    self._status = FAILED
                self._generation += 1
                self._proc = None

    def cancel(self) -> bool:
        """Terminate the active process group. No-op when idle."""
        with self._lock:
            proc = self._proc
            if proc is None or self._status != RUNNING:
                return False
            self._cancelling = True
            self._lines.append("")
            self._lines.append("── cancelled by user ──")
            pid = proc.pid

        try:
            os.killpg(os.getpgid(pid), signal.SIGTERM)
        except (ProcessLookupError, PermissionError, OSError):
            try:
                proc.terminate()
            except OSError:
                return False

        def _hard_kill() -> None:
            time.sleep(5)
            with self._lock:
                still = self._proc
            if still is not None and still.pid == pid:
                try:
                    os.killpg(os.getpgid(pid), signal.SIGKILL)
                except OSError:
                    pass

        threading.Thread(target=_hard_kill, daemon=True).start()
        return True

    # -- inspection --------------------------------------------------------

    def snapshot(self) -> dict:
        """A consistent view of the run for the UI."""
        with self._lock:
            elapsed = None
            if self._started_at is not None:
                end = self._finished_at or time.time()
                elapsed = end - self._started_at
            return {
                "status": self._status,
                "kind": self._kind,
                "log": "\n".join(self._lines),
                "n_lines": len(self._lines),
                "returncode": self._returncode,
                "elapsed": elapsed,
                "argv": list(self._argv),
                "meta": dict(self._meta),
                "generation": self._generation,
                "progress": self._progress_view(),
            }

    def _progress_view(self) -> dict:
        """The current phase as a fraction, with an estimate of what is left.

        Caller holds the lock. Returns {} when nothing has reported yet, which
        is the honest state for the first minutes of a load: better a bar that
        appears late than one that sits at 0% claiming to know something.

        The estimate is linear in the CURRENT phase only — see _pump — and is
        withheld until a few items are done, because extrapolating a whole run
        from one sample produces a number that swings by tens of minutes
        between polls and is worse than no number at all.
        """
        p = self._progress
        if not p or not p.get("total"):
            return {}
        done, total = p["done"], p["total"]
        fraction = max(0.0, min(1.0, done / total))
        remaining = None
        if (self._phase_started_at is not None and done >= 3 and done < total
                and self._status == RUNNING):
            per_item = (time.time() - self._phase_started_at) / done
            remaining = per_item * (total - done)
        return {
            "phase": p["phase"],
            "phase_label": PROGRESS_PHASES.get(p["phase"], p["phase"]),
            "done": done,
            "total": total,
            "label": p.get("label", ""),
            "percent": round(fraction * 100),
            "remaining": remaining,
        }

    @property
    def busy(self) -> bool:
        with self._lock:
            return self._status == RUNNING


#: The app's single runner.
RUNNER = Runner()


def format_elapsed(seconds: float | None) -> str:
    if seconds is None:
        return "—"
    seconds = int(seconds)
    return f"{seconds // 60:02d}:{seconds % 60:02d}"
