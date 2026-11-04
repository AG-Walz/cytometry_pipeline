"""Phase durations predicted from the size of the job, for unseen settings.

A run with settings the app has seen before is estimated from that run (see
phase_history). This module is for the first run of a settings combination:
the time of each phase is predicted from what the run is about to do, which
R states in its ``::runinfo::`` line (emit_run_info in src/R_shared/io_utils.R)
and in the totals of its progress lines.

The formulas and constants are fitted on the 15 complete runs logged on
highmem1 and highmem2 up to 2026-11-04 (J445 and the September "live cells"
runs; fit sets of 8,000 to 2.86M cells, 15 to 80 neighbours, 1 to 27
threads). M = million cells; k = neighbours; the edges of the neighbour
graph are 1.19-1.28 x fit x k in every run, so 1.23 x fit x k here.

=========================  ===========================================  ======
phase                      seconds                                      error (*)
=========================  ===========================================  ======
gating sets                0.33 x samples                               ~10 %
cell extraction            0.2 x samples                                x0.3-2
neighbour search           212 x fit^1.16 x k^0.40 / threads^0.43       ~10 %
initialisation, spectral   1.26 x (fit x k)^2.13                        ~40 %
initialisation, PCA        20 x fit
optimisation               118 x edges x epochs / 1000 x fit^0.2        ~15-35 %
                           / optimiser threads
projection                 20.3 x cells^0.92 x k^0.63 x (27/threads)^0.43  ~7 %
plots                      (14 + 1.9 x plots) x cells
=========================  ===========================================  ======

(*) typical error when each run is predicted from the other runs only. The
spectral initialisation is the weak point: it grows with the square of the
graph, and only one logged run is large (2.86M cells, 52 min); predicted
from the others alone it came out 2.7 times too long. The same run on the
same host also differs between days (highmem1, 2.59M cells: neighbour search
1196 s once and 2704 s once), so no formula is better than about x2 on a
busy host. Each host therefore corrects the constants with its own complete
runs (calibrate): the median of measured / predicted per phase.
"""

from __future__ import annotations

import statistics

EDGES_PER_FIT_K = 1.23
#: Before the plots start their number is unknown.
DEFAULT_PLOTS = 20
#: uwot's epochs for a large fit (serial optimiser), and with --embed_sgd_cores
#: above 1 (batch optimiser; see run_umap in src/R_shared/stats_utils.R).
DEFAULT_EPOCHS, BATCH_EPOCHS = 200, 500
#: A host's correction of a constant is kept within this range.
CALIBRATION_LIMITS = (0.33, 3.0)


def predict(info: dict, totals: dict) -> dict[str, float] | None:
    """Predicted seconds per phase, or None before R has said the job's size.

    ``info`` is the ``::runinfo::`` line, ``totals`` the latest ``total`` of
    each phase's progress lines (samples, epochs, plots).
    """
    if not info or "cells" not in info:
        return None
    cells = info["cells"] / 1e6
    fit = info.get("fit_cells", info["cells"]) / 1e6
    k = max(1.0, info.get("neighbors", 15))
    threads = max(1.0, info.get("threads", 1))
    sgd = max(1.0, info.get("sgd_threads", 1))
    epochs = totals.get("embed", BATCH_EPOCHS if sgd > 1 else DEFAULT_EPOCHS)
    spectral = info.get("spectral", 1) >= 1
    projected = max(0.0, cells - fit)
    edges = EDGES_PER_FIT_K * fit * k
    plots = totals.get("plot", DEFAULT_PLOTS)
    out = {
        "embed_nn": 212 * fit ** 1.16 * k ** 0.40 / threads ** 0.43 if fit else 0.0,
        "embed_init": (1.26 * (fit * k) ** 2.13 if spectral else 20 * fit)
                      if fit else 0.0,
        "embed": (118 * edges * epochs / 1000 * fit ** 0.2 / sgd if fit else 0.0)
                 + (20.3 * projected ** 0.92 * k ** 0.63 * (27 / threads) ** 0.43
                    if projected else 0.0),
        "plot": (14 + 1.9 * plots) * max(cells, 0.05),
    }
    if "extract" in totals:
        out["load"] = 0.33 * totals["extract"]
        out["extract"] = 0.2 * totals["extract"]
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
