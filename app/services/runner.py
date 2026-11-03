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

from app.services import detached

#: Hard cap on retained log lines, so a pathological run cannot exhaust memory.
MAX_LINES = 20000

#: Progress lines emitted by ``emit_progress()`` in src/R_shared/io_utils.R.
#: They are intercepted here rather than shown: they are for the progress bar,
#: and one per sample would bury the log they share a stream with. R owns the
#: counting; this only reads it. Keep in step with PROGRESS_PREFIX there.
#: ``elapsed`` is optional so a log written before R started emitting it still
#: parses — the bar then loses its timings, not its position. It sits before
#: ``label``, which stays last because a label may contain anything at all.
_PROGRESS_RE = re.compile(
    r"^::progress:: phase=(\S+) done=(\d+) total=(\d+)"
    r"(?: elapsed=([0-9.]+))? label=(.*)$"
)

#: Every phase a run can report, IN ORDER, with how it reads in the UI and how
#: much of the run it is expected to take.
#:
#: The weights are what make the bar honest about the shape of a run rather than
#: just its position: the four embedding phases are not four equal quarters, and
#: a bar that draws them as such tells someone watching the 20-minute neighbour
#: search that they are 5% through when they are nearer half. They are only
#: PRIORS — as each phase finishes, its measured duration replaces its weight
#: (see Runner._phase_weights), so the estimate converges on the machine and the
#: data actually in front of it. Derived from a 9.8M-cell, 328-sample run and
#: the 500k-cell measurements recorded in src/R_shared/stats_utils.R.
PROGRESS_PHASES = {
    "load": "Building gating sets",
    "extract": "Extracting cells",
    "embed_nn": "Neighbour search",
    "embed_init": "Initialising the layout",
    "embed": "Optimising the embedding",
    "plot": "Drawing plots",
}

#: Prior share of a run, per phase. Relative, not absolute — only the ratios are
#: used. Kept beside PROGRESS_PHASES so a new phase cannot be added without one.
PHASE_WEIGHTS = {
    "load": 8.0,
    "extract": 10.0,
    "embed_nn": 45.0,
    "embed_init": 7.0,
    "embed": 20.0,
    "plot": 10.0,
}

#: The order phases occur in, which is the order they are drawn in. Taken from
#: PROGRESS_PHASES rather than repeated, so the two cannot disagree.
PHASE_ORDER = list(PROGRESS_PHASES)

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
        #: {phase: [first_elapsed, last_elapsed]} in the R process's own
        #: seconds. What turns the segmented bar from a guess into a
        #: measurement as the run proceeds; see _phase_weights.
        self._phase_times: dict[str, list[float]] = {}
        #: Set when following a run that is not our child (see attach()).
        self._detached_dir: str | None = None
        self._detached_pid: int = 0

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
            self._phase_times = {}
            self._started_at = time.time()
            self._finished_at = None
            self._cancelling = False
            # This run's output comes down a pipe. Any detached run we were
            # following is no longer what the UI is showing, and leaving the
            # pointer set made snapshot() re-read THAT run's log on every poll —
            # overwriting this one's lines and re-deciding this one's status from
            # a finished process's pid. That is what surfaced as an inspect
            # reporting "failed (exit code None)" while printing a complete log.
            self._detached_dir = None
            self._detached_pid = 0

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
                    phase, done, total, elapsed, label = hit.groups()
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
                        self._note_phase_time(phase, elapsed)
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
        if self._detached_dir:
            stopped = detached.stop(self._detached_dir)
            if stopped:
                with self._lock:
                    self._cancelling = True
                    self._status = CANCELLED
                    self._finished_at = time.time()
                    self._generation += 1
            return stopped
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

    # -- detached runs -----------------------------------------------------

    def start_detached(self, argv, cwd, output_dir, kind, meta=None):
        """Launch a run that outlives the app, and follow its log.

        Same single-slot rule as start(): one thing at a time. The difference is
        where the output goes — a file in the run's own output directory rather
        than a pipe — which is what makes it survivable and re-readable. See
        app/services/detached.py.
        """
        with self._lock:
            if self._status == RUNNING:
                return False, "A run is already in progress."
        try:
            status = detached.launch(argv, cwd, output_dir, kind, meta)
        except OSError as exc:
            with self._lock:
                self._status = FAILED
                self._lines = [f"Failed to launch: {exc}", "", "Command was:",
                               "  " + " ".join(argv)]
                self._finished_at = time.time()
                self._generation += 1
            return False, str(exc)
        self.attach(output_dir)
        return True, "started"

    def attach(self, output_dir) -> bool:
        """Follow a detached run — ours, or one a previous app instance began.

        Returns False when there is nothing there to follow. Refuses while
        something else is running, so attaching cannot silently displace the run
        the user is currently watching.
        """
        status = detached.read_status(output_dir)
        if not status:
            return False
        with self._lock:
            if self._status == RUNNING and self._detached_dir != str(output_dir):
                return False
            self._detached_dir = str(output_dir)
            self._argv = list(status.get("argv") or [])
            self._kind = status.get("kind") or "run"
            self._meta = dict(status.get("meta") or {})
            self._started_at = status.get("started_at")
            self._detached_pid = int(status.get("pid") or 0)
            self._returncode = status.get("exit_code")
            self._status = RUNNING if status.get("running") else DONE
            self._finished_at = None if status.get("running") else time.time()
            self._progress = {}
            self._phase_started_at = None
            self._phase_times = {}
            self._lines = []
        self._refresh_detached()
        return True

    def _refresh_detached(self) -> None:
        """Re-read the detached run's log and re-derive status and progress.

        The log is the only channel — there is no pipe to a process we may not
        even have started — so progress is recovered by scanning it. Cheap
        enough at poll rate: the scan is a regex over the tail, and only the
        LAST progress line of each phase matters.
        """
        directory = self._detached_dir
        if not directory:
            return
        # Aliveness FIRST, then the log. The other order loses a race that
        # happens on every fast run: read a log that is still being written,
        # then observe a process that has since exited, and conclude from an
        # incomplete log that it stopped early. Asking about the process first
        # means a log read after "not alive" is a log nobody can still append
        # to.
        alive = detached.is_alive(self._detached_pid)
        text = detached.read_log(directory)
        code = None if alive else detached.read_exit_code(directory)
        lines = text.splitlines()

        keep: list[str] = []
        progress: dict = {}
        seen: list[tuple[str, str | None]] = []
        for line in lines:
            hit = _PROGRESS_RE.match(line)
            if hit:
                phase, done, total, elapsed, label = hit.groups()
                progress = {"phase": phase, "done": int(done),
                            "total": int(total), "label": label}
                # Rebuilt from the log every poll, because for a detached run
                # the log is the only memory there is: a browser reopened an
                # hour into a run has to recover how long each earlier phase
                # took from what R wrote, not from anything held here.
                seen.append((phase, elapsed))
                continue
            keep.append(line)

        with self._lock:
            if len(keep) > MAX_LINES:
                self._truncated = True
                keep = keep[-MAX_LINES:]
            self._lines = keep
            self._phase_times = {}
            for phase, elapsed in seen:
                self._note_phase_time(phase, elapsed)
            if progress.get("phase") != self._progress.get("phase"):
                # Re-derived after an attach, so the phase clock starts now
                # rather than pretending to know when the phase began. The ETA
                # is therefore absent for a moment after attaching, which is
                # better than an ETA computed from a start time we invented.
                self._phase_started_at = time.time()
            self._progress = progress
            if alive:
                return
            finished = code == 0 if code is not None else self._looks_finished(keep)
            verdict = DONE if finished else FAILED
            if self._status == RUNNING:
                self._returncode = code
                self._status = verdict
                self._finished_at = time.time()
                self._generation += 1
            elif self._status == FAILED and verdict == DONE:
                # A verdict taken from a log that turned out to be incomplete.
                # Correcting it matters more than never changing our mind: the
                # alternative is an app that says a finished run failed and has
                # no way to stop saying it.
                self._returncode = code
                self._status = DONE
                self._generation += 1

    def _note_phase_time(self, phase: str, elapsed: str | None) -> None:
        """Record when a phase was first and last heard from.

        Caller holds the lock. ``elapsed`` is R's own seconds-since-start (see
        emit_progress in src/R_shared/io_utils.R), so it needs no shared clock
        and survives the app being restarted mid-run. A line without it — an
        older R, or another emitter — is simply not counted; the phase keeps its
        prior weight rather than getting a wrong one.
        """
        if elapsed is None:
            return
        try:
            value = float(elapsed)
        except ValueError:
            return
        span = self._phase_times.get(phase)
        if span is None:
            self._phase_times[phase] = [value, value]
        else:
            span[1] = value

    def _phase_weights(self) -> dict[str, float]:
        """Relative width of each phase's segment.

        Caller holds the lock. A phase that has FINISHED is worth what it
        actually took; one that has not is worth its prior, rescaled so the
        priors and the measurements are in the same units. That rescaling is the
        whole trick: without it, a run whose measured seconds are large next to
        the priors would draw every unfinished phase as a sliver, and the bar
        would lurch every time a phase ended.

        With nothing measured yet this returns the priors unchanged, which is
        the right answer for the first seconds of a run.
        """
        measured = {p: max(0.0, s[1] - s[0]) for p, s in self._phase_times.items()}
        current = self._progress.get("phase")
        # The current phase is still growing, so its measurement is a lower
        # bound, not a duration. It keeps its prior until it is behind us.
        done_phases = {p: v for p, v in measured.items()
                       if p != current and v > 0}
        if not done_phases:
            return dict(PHASE_WEIGHTS)
        prior_of_done = sum(PHASE_WEIGHTS.get(p, 1.0) for p in done_phases)
        scale = (sum(done_phases.values()) / prior_of_done) if prior_of_done else 1.0
        return {
            p: done_phases.get(p, PHASE_WEIGHTS.get(p, 1.0) * scale)
            for p in PHASE_ORDER
        }

    def _phase_view(self) -> list[dict]:
        """One entry per phase, in order: its width, its state, its fill.

        Caller holds the lock. This is the whole segmented bar — the callback
        that draws it works nothing out.
        """
        weights = self._phase_weights()
        total_weight = sum(weights.values()) or 1.0
        current = self._progress.get("phase")
        seen_current = False
        out = []
        for phase in PHASE_ORDER:
            if phase == current:
                seen_current = True
                done, total = self._progress.get("done", 0), self._progress.get("total", 1)
                fill = max(0.0, min(1.0, done / total)) if total else 0.0
                state = "running"
            elif phase in self._phase_times and not seen_current:
                fill, state = 1.0, "done"
            elif seen_current:
                fill, state = 0.0, "pending"
            else:
                # Never reported and already behind us: a phase this run does
                # not have (no plots, or an embedding fed pre-computed
                # coordinates). Drawn as skipped rather than silently dropped,
                # so the shape of the bar still matches the run.
                fill, state = 0.0, "skipped"
            out.append({
                "phase": phase,
                "label": PROGRESS_PHASES.get(phase, phase),
                "weight": round(100.0 * weights.get(phase, 0.0) / total_weight, 2),
                "fill": round(fill, 4),
                "state": state,
            })
        return out

    @staticmethod
    def _looks_finished(lines: list[str]) -> bool:
        """Did the run reach its end, or did it stop?

        Only for runs started before the launcher recorded an exit code — see
        detached.read_exit_code(), which is the real answer when it exists.
        Completion is judged from what the run printed; the R side prints one
        of these last and only on success.
        """
        tail = "\n".join(lines[-40:])
        return "Done. Plots written to:" in tail or "Load report:" in tail

    # -- inspection --------------------------------------------------------

    def snapshot(self) -> dict:
        """A consistent view of the run for the UI."""
        # A detached run has no pipe feeding the buffer, so its log and status
        # are re-read here, at poll rate, before the view is taken.
        if self._detached_dir:
            self._refresh_detached()
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

        phases = self._phase_view()
        # How much of the WHOLE run is behind us: every phase's share of the bar
        # times how full it is. This is the number the estimate comes from, and
        # it is why the estimate no longer collapses when a phase changes — the
        # old one extrapolated the current phase alone, so "Drawing plots" at
        # 1/300 predicted hours a second after a 20-minute search finished.
        overall = sum(s["weight"] * s["fill"] for s in phases) / 100.0

        remaining = None
        run_elapsed = max((s[1] for s in self._phase_times.values()), default=0.0)
        if self._status == RUNNING and overall > 0.02 and run_elapsed > 5:
            remaining = run_elapsed * (1.0 - overall) / overall
        elif (self._phase_started_at is not None and done >= 3 and done < total
                and self._status == RUNNING):
            # No usable elapsed in the log (an older R). Fall back to the
            # phase-local estimate rather than showing nothing.
            per_item = (time.time() - self._phase_started_at) / done
            remaining = per_item * (total - done)

        return {
            "phase": p["phase"],
            "phase_label": PROGRESS_PHASES.get(p["phase"], p["phase"]),
            "done": done,
            "total": total,
            "label": p.get("label", ""),
            "percent": round(fraction * 100),
            "overall_percent": round(overall * 100),
            "phases": phases,
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
