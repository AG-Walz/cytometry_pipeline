"""How long a run has left, from what it has done and what the last run took.

Pure functions: the Runner collects the facts (the current run's progress
lines, the reference run, the clock) and this works out the bar and the
estimate. Kept apart from the Runner so it can be tested without a process.

The old estimate extrapolated ``elapsed * (1 - done) / done`` over prior phase
weights. Two things made that wrong by a factor of ten on real runs:

* the run log of an output folder holds EVERY run made into it, and all of
  them were read as one run (a 1 h 30 run showed "16 h 31 min left");
* a phase with one step ("Initialising the layout", 0 / 1) adds elapsed time
  but no progress, so the estimate went UP for as long as the phase ran. On
  the J445 data that phase is 52 of 93 minutes, not the 7 % the priors said.

The estimate now comes from a REFERENCE: the phase durations of the last
complete run with the same settings (see phase_history), or else of the last
complete run in the same log. Without one there is no number, because the
priors cannot predict a phase whose cost depends on the data this much.
"""

from __future__ import annotations

import time

#: Every phase a run can report, IN ORDER, with how it reads in the UI.
#: Keep in step with emit_progress() calls in src/.
PROGRESS_PHASES = {
    "load": "Building gating sets",
    "extract": "Extracting cells",
    "embed_nn": "Neighbour search",
    "embed_init": "Initialising the layout",
    "embed": "Optimising the embedding",
    "plot": "Drawing plots",
}

#: The word under each part of the bar. Short, because a part can be narrow.
SHORT_LABELS = {
    "load": "Gating",
    "extract": "Cells",
    "embed_nn": "Neighbours",
    "embed_init": "Layout start",
    "embed": "Optimise",
    "plot": "Plots",
}

#: Prior share of a run, per phase. Used ONLY to draw the bar when there is no
#: reference run; it is never turned into a time.
PHASE_WEIGHTS = {
    "load": 8.0,
    "extract": 10.0,
    "embed_nn": 45.0,
    "embed_init": 7.0,
    "embed": 20.0,
    "plot": 10.0,
}

PHASE_ORDER = list(PROGRESS_PHASES)

#: "Resuming: '<stage>' reused from checkpoints/..." (ckpt_load in
#: src/visualization/umap_from_flowjo.R) -> the phases that then do not run.
#: A reused "model" still projects every cell, so "embed" still runs.
RESUMED_PHASES = {
    "extract": ("load", "extract"),
    "model": ("embed_nn", "embed_init"),
    "embedding": ("embed_nn", "embed_init", "embed"),
}

#: The spread shown around the estimate. Not a statistical interval: the same
#: settings on the same host differ by about this much between runs.
SPREAD = 0.15

#: The speed of this run relative to the reference, from the phases already
#: done, is clamped to this range. One slow file read must not triple the rest.
SPEED_LIMITS = (0.67, 1.5)

#: Phases with no more than this many steps are timed by the clock, not by
#: their steps: "0 / 1" says nothing about how far the phase is.
FEW_STEPS = 3


def durations(phase_times: dict[str, list[float]]) -> dict[str, float]:
    """Seconds per phase, from the start of the phase to the start of the next.

    ``phase_times`` is {phase: [first_elapsed, last_elapsed]} for ONE run. The
    last line of a phase is sent when its last item STARTS, so first-to-last
    leaves out the last item (for "load", a whole workspace). Start-to-start
    does not. The last phase runs to its last line.
    """
    ordered = sorted(phase_times.items(), key=lambda kv: kv[1][0])
    out = {}
    for i, (phase, (first, last)) in enumerate(ordered):
        end = ordered[i + 1][1][0] if i + 1 < len(ordered) else last
        out[phase] = max(0.0, end - first)
    return out


def _speed(done: dict[str, float], reference: dict[str, float]) -> float:
    """This run's time per phase over the reference's, for the phases done."""
    common = [p for p in done if reference.get(p)]
    ref = sum(reference[p] for p in common)
    took = sum(done[p] for p in common)
    if ref < 60 or took <= 0:
        return 1.0
    lo, hi = SPEED_LIMITS
    return max(lo, min(hi, took / ref))


def estimate(progress: dict, phase_times: dict[str, list[float]],
             resumed: set[str], reference: dict[str, float] | None,
             run_now: float | None, last_output_at: float | None,
             now: float | None = None) -> dict:
    """The progress view: one entry per phase, and the time left.

    ``progress`` is the latest progress line of the current run, ``run_now``
    the current time in the R process's seconds, ``reference`` the phase
    durations of the reference run (or None).
    """
    now = time.time() if now is None else now
    current = progress.get("phase")
    done_n, total_n = progress.get("done", 0), progress.get("total", 0)
    item_fill = max(0.0, min(1.0, done_n / total_n)) if total_n else 0.0
    measured = durations(phase_times)
    last_elapsed = max((s[1] for s in phase_times.values()), default=0.0)
    run_now = max(run_now or 0.0, last_elapsed)
    t_in = run_now - phase_times[current][0] if current in phase_times else 0.0
    t_in = max(0.0, t_in)

    # The state of each phase, before any width or time is given to it.
    states, seen_current = {}, False
    for phase in PHASE_ORDER:
        if phase == current:
            seen_current = True
            states[phase] = "running"
        elif phase in resumed:
            states[phase] = "reused"
        elif seen_current:
            states[phase] = "pending"
        elif phase in phase_times:
            states[phase] = "done"
        else:
            # Never reported and already behind us: a phase this run does not
            # have. Drawn as skipped, so the bar still matches the run.
            states[phase] = "skipped"
    done = {p: measured[p] for p, s in states.items() if s == "done"}

    ref = {p: v for p, v in (reference or {}).items() if v > 0}
    speed = _speed(done, ref) if ref else 1.0
    expected_cur = ref[current] * speed if current in ref else None
    by_clock = total_n <= FEW_STEPS or item_fill < 0.1

    # Fill of the running phase. A phase with few steps fills with the clock
    # against the reference, so the bar moves during a 50-minute silent step.
    fill = item_fill
    if by_clock and expected_cur:
        fill = max(fill, min(0.97, t_in / expected_cur))

    # Time left in the running phase.
    overrun = False
    if not by_clock and item_fill > 0:
        cur_left = t_in * (1.0 - item_fill) / item_fill
    elif expected_cur is not None:
        cur_left = max(0.0, expected_cur - t_in)
        overrun = t_in > expected_cur
    else:
        cur_left = None

    # Widths. With a reference: what each phase took or is expected to take.
    # Without one: the priors, rescaled to what the done phases measured.
    weights = {}
    if ref:
        for phase, state in states.items():
            if state == "done":
                weights[phase] = done[phase]
            elif state == "running":
                weights[phase] = max(expected_cur or 0.0, t_in)
            elif state == "pending":
                weights[phase] = ref.get(phase, 0.0) * speed
            else:
                weights[phase] = 0.0
    else:
        prior_done = sum(PHASE_WEIGHTS[p] for p in done)
        scale = (sum(done.values()) / prior_done) if prior_done and \
            sum(done.values()) > 0 else 1.0
        for phase, state in states.items():
            if state == "done" and done[phase] > 0:
                weights[phase] = done[phase]
            elif state in ("reused", "skipped"):
                weights[phase] = 0.0
            else:
                weights[phase] = PHASE_WEIGHTS[phase] * scale
    total_w = sum(weights.values()) or 1.0

    phases, steps = [], []
    for phase in PHASE_ORDER:
        state = states[phase]
        f = {"done": 1.0, "reused": 1.0, "running": fill}.get(state, 0.0)
        phases.append({
            "phase": phase,
            "label": PROGRESS_PHASES[phase],
            "short": SHORT_LABELS[phase],
            "weight": round(100.0 * weights[phase] / total_w, 2),
            "fill": round(f, 4),
            "state": state,
        })
        # The steps this run has: not reused, not skipped, and (with a
        # reference) not a pending phase the reference run never had.
        if state in ("done", "running") or (
                state == "pending" and (not ref or phase in ref)):
            steps.append(phase)
    overall = sum(s["weight"] * s["fill"] for s in phases) / 100.0

    remaining = None
    if ref and cur_left is not None:
        later = sum(ref.get(p, 0.0) * speed for p, s in states.items()
                    if s == "pending")
        remaining = cur_left + later

    return {
        "phase": current,
        "phase_label": PROGRESS_PHASES.get(current, current),
        "done": done_n,
        "total": total_n,
        "label": progress.get("label", ""),
        "percent": round(item_fill * 100),
        "overall_percent": round(overall * 100),
        "phases": phases,
        "step": steps.index(current) + 1 if current in steps else None,
        "steps": len(steps),
        "few_steps": total_n <= FEW_STEPS,
        "run_elapsed": run_now,
        "phase_elapsed": t_in,
        "phase_expected": expected_cur,
        "overrun": overrun,
        "has_reference": bool(ref),
        "remaining": remaining,
        "remaining_low": remaining * (1 - SPREAD) if remaining is not None else None,
        "remaining_high": remaining * (1 + SPREAD) if remaining is not None else None,
        "finish_at": now + remaining if remaining is not None else None,
        "silent_for": (now - last_output_at) if last_output_at else None,
    }


def resumed_phases(lines: list[str]) -> set[str]:
    """The phases a run skips because it reused a saved step."""
    out: set[str] = set()
    for line in lines:
        if line.startswith("Resuming: '"):
            stage = line[len("Resuming: '"):].split("'", 1)[0]
            out.update(RESUMED_PHASES.get(stage, ()))
    return out
