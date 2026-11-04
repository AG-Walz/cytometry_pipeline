"""Phase durations predicted from the size of the job, for unseen settings.

A run with settings the app has seen before is estimated from that run (see
phase_history). This module is for the first run of a settings combination:
the time of each phase is predicted from what the run is about to do, which
R states in its ``::runinfo::`` line (emit_run_info in src/R_shared/io_utils.R)
and in the totals of its progress lines.

The formulas follow how each phase scales, and the constants are measured:

* gating sets and cell extraction: per sample;
* neighbour search: (fitted cells x neighbours) ** 0.9 / threads ** 0.5
  (the search is limited by memory, so threads help less than linearly);
* layout initialisation (spectral, one thread): (fitted cells x neighbours)
  ** 1.3;
* optimisation: fitted cells x neighbours x epochs / optimiser threads (1 by
  default, see --embed_sgd_cores), plus the projection of the other cells:
  cells x neighbours / threads;
* plots: per plot and per cell.

Constants from the J445 runs on highmem2 of 2026-11-03 (27 search threads,
1 optimiser thread, 15 neighbours): one fitted on 300,000 of 2.86M cells,
one fitted on all of them. The powers are from these two runs: 9.5 times the
cells took 7.4 times as long to search and 19 times as long to initialise.

These are rough on other data and other hosts. On a 16-core laptop, a
114,000-cell example ran about 3 times faster than predicted (53 s against
about 150 s). So each host corrects the constants with its own complete runs
(calibrate): the median of measured / predicted per phase.
"""

from __future__ import annotations

import statistics

#: Seconds per unit of work, per phase (units as in predict()).
SECONDS = {
    "load": 0.25,          # per sample
    "extract": 0.25,       # per sample
    "embed_nn": 95.7,      # per (M fitted cells x neighbours) ** NN_POWER / threads ** 0.5
    "embed_init": 23.5,    # per (M fitted cells x neighbours) ** INIT_POWER
    "embed": 0.180,        # per (M fitted cells x neighbour x epoch / optimiser thread)
    "plot": 1.8,           # per (plot x M cells)
}
NN_POWER = 0.9
INIT_POWER = 1.3
#: Seconds per (M projected cells x neighbour / thread).
PROJECT_SECONDS = 165.0
#: Before the plots start their number is unknown: seconds per M cells.
PLOTS_PER_MILLION_CELLS = 60.0
#: uwot's epochs for a large fit, used until the optimisation reports its own.
DEFAULT_EPOCHS = 200
#: A host's correction of a constant is kept within this range.
CALIBRATION_LIMITS = (0.33, 3.0)


def predict(info: dict, totals: dict) -> dict[str, float] | None:
    """Predicted seconds per phase, or None before R has said the job's size.

    ``info`` is the ``::runinfo::`` line, ``totals`` the latest ``total`` of
    each phase's progress lines (workspaces, samples, epochs, plots).
    """
    if not info or "cells" not in info:
        return None
    cells = info["cells"] / 1e6
    fit = info.get("fit_cells", info["cells"]) / 1e6
    k = info.get("neighbors", 15)
    threads = max(1.0, info.get("threads", 1))
    sgd = max(1.0, info.get("sgd_threads", 1))
    epochs = totals.get("embed", DEFAULT_EPOCHS)
    out = {
        "embed_nn": SECONDS["embed_nn"] * (fit * k) ** NN_POWER / threads ** 0.5,
        "embed_init": SECONDS["embed_init"] * (fit * k) ** INIT_POWER,
        "embed": (SECONDS["embed"] * fit * k * epochs / sgd
                  + PROJECT_SECONDS * max(0.0, cells - fit) * k / threads),
        "plot": (SECONDS["plot"] * totals["plot"] * max(cells, 0.05)
                 if "plot" in totals else PLOTS_PER_MILLION_CELLS * cells),
    }
    if "extract" in totals:
        out["load"] = SECONDS["load"] * totals["extract"]
        out["extract"] = SECONDS["extract"] * totals["extract"]
    return out


def calibrate(prediction: dict[str, float], entries: list[dict]) -> dict[str, float]:
    """Correct a prediction with what this host measured on complete runs.

    ``entries`` are phase_history entries; the ones saved with their run's
    size are compared with what the constants would have predicted for them.
    """
    ratios: dict[str, list[float]] = {}
    for entry in entries:
        guess = predict(entry.get("info") or {}, entry.get("totals") or {})
        for phase, took in (entry.get("phases") or {}).items():
            if guess and guess.get(phase, 0) > 1 and took > 0:
                ratios.setdefault(phase, []).append(took / guess[phase])
    lo, hi = CALIBRATION_LIMITS
    return {phase: seconds * max(lo, min(hi, statistics.median(ratios[phase])))
            if phase in ratios else seconds
            for phase, seconds in prediction.items()}
