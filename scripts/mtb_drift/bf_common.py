"""Shared pieces for the BayesFlow experiments (steps 7 and 8).

Both variants infer the SAME parameters the brms/Stan yardstick estimates - an
instrument shift `b_instr` and one shift per acquisition day `c_day[d]` - so the
posteriors are directly comparable in step 9 (effects_bayesflow_v1/v2.csv, same
long format as effects_brms.csv). CPU only; keras backend is torch.
"""
import os
os.environ.setdefault("KERAS_BACKEND", "torch")
import argparse, json, time
import numpy as np
import pandas as pd


def read_design(results_dir, layer="A"):
    """The real design: one row per neg tube with integer codes for patient, day, instrument."""
    tubes = pd.read_csv(os.path.join(results_dir, "tubes.csv"))
    negs = tubes[tubes.is_neg == True].copy()  # noqa: E712
    if os.path.exists(os.path.join(results_dir, "excluded_tubes.csv")):
        ex = pd.read_csv(os.path.join(results_dir, "excluded_tubes.csv")).filename
        negs = negs[~negs.filename.isin(ex)]
    negs = negs.sort_values(["patient", "acq_day"]).reset_index(drop=True)
    negs["patient_i"] = pd.factorize(negs.patient)[0]
    negs["day_i"] = pd.factorize(negs.acq_day)[0]
    ref = read_config()["reference_instrument"]
    negs["instr"] = (negs.instrument != ref).astype(float)
    return negs


def read_config():
    import yaml
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "..", "..", "config", "mtb_drift.yaml")) as fh:
        return yaml.safe_load(fh)


def day_levels(negs):
    return list(pd.unique(negs.acq_day))


def effects_long(model, layer, response, instr_draws, day_draws, day_names, instr_level):
    """Same columns as R effects_long(): model, layer, response, term, level, draw, value."""
    rows = [pd.DataFrame(dict(model=model, layer=layer, response=response, term="instrument",
                              level=instr_level, draw=np.arange(1, len(instr_draws) + 1),
                              value=instr_draws))]
    for j, name in enumerate(day_names):
        rows.append(pd.DataFrame(dict(model=model, layer=layer, response=response, term="day",
                                      level=name, draw=np.arange(1, day_draws.shape[0] + 1),
                                      value=day_draws[:, j])))
    return pd.concat(rows, ignore_index=True)


def base_parser(desc):
    p = argparse.ArgumentParser(description=desc)
    p.add_argument("--out", required=True, help="results dir written by steps 1-5")
    p.add_argument("--epochs", type=int, default=50)
    p.add_argument("--num_sims", type=int, default=50000, help="offline training simulations")
    p.add_argument("--batch_size", type=int, default=64)
    p.add_argument("--posterior_samples", type=int, default=2000)
    p.add_argument("--n_sbc", type=int, default=300, help="held-out simulations for calibration/recovery")
    p.add_argument("--seed", type=int, default=20260921)
    p.add_argument("--smoke", action="store_true", help="tiny run to check the pipeline end to end")
    p.add_argument("--restore", action="store_true", help="skip training, load the saved approximator from the checkpoint")
    p.add_argument("--chunk", type=int, default=25, help="test sets per posterior-sampling call (memory bound)")
    return p


def sample_chunked(workflow, conditions, num_samples, chunk):
    """workflow.sample over a dict of arrays in slices of `chunk` along axis 0, to bound memory.
    Pushing hundreds of test sets through a cell-level summary network at once needs tens of GB."""
    n = next(iter(conditions.values())).shape[0]
    parts = []
    for i in range(0, n, chunk):
        parts.append(workflow.sample(conditions={k: v[i:i + chunk] for k, v in conditions.items()},
                                     num_samples=num_samples))
    return {k: np.concatenate([p[k] for p in parts], axis=0) for k in parts[0]}


def maybe_restore(workflow, restore, path, custom_objects=None):
    """Load the saved approximator instead of training. Returns True when restored.
    `custom_objects` names any custom keras classes the saved model uses."""
    if not restore:
        return False
    import keras
    try:
        workflow.approximator = keras.saving.load_model(path, custom_objects=custom_objects)
    except Exception:
        workflow.load_approximator()
    print("restored approximator from", path)
    return True


def recovery_table(posterior, truth, names):
    """posterior: (n_test, n_samples, n_par); truth: (n_test, n_par)."""
    mean = posterior.mean(axis=1)
    rows = []
    for j, n in enumerate(names):
        t, m = truth[:, j], mean[:, j]
        r2 = 1 - ((t - m) ** 2).sum() / max(((t - t.mean()) ** 2).sum(), 1e-12)
        lo, hi = np.quantile(posterior[:, :, j], [0.025, 0.975], axis=1)
        rows.append(dict(parameter=n, r2=r2, coverage95=float(((t >= lo) & (t <= hi)).mean()),
                         rmse=float(np.sqrt(((t - m) ** 2).mean()))))
    return pd.DataFrame(rows)


def sbc_ranks(posterior, truth):
    """Rank of the truth among posterior samples, per test set and parameter (uniform if calibrated)."""
    return (posterior < truth[:, None, :]).sum(axis=1)


def sbc_uniformity(ranks, n_samples, n_bins=20):
    """Chi-square p-value of the rank histogram per parameter; > 0.05 = no evidence of miscalibration."""
    from scipy.stats import chisquare
    out = []
    for j in range(ranks.shape[1]):
        h, _ = np.histogram(ranks[:, j], bins=n_bins, range=(0, n_samples + 1))
        out.append(float(chisquare(h).pvalue))
    return out


class Timer:
    def __init__(self): self.t0 = time.time(); self.marks = {}
    def mark(self, k): self.marks[k] = round(time.time() - self.t0, 1)
    def dump(self, path):
        with open(path, "w") as fh: json.dump(self.marks, fh, indent=1)
