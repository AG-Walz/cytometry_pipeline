"""Shared pieces for the BayesFlow experiments (steps 7 and 8).

Both variants infer the SAME parameters the brms/Stan yardstick estimates - an
instrument shift `b_instr` and one shift per acquisition day `c_day[d]` - so the
posteriors are directly comparable in step 9 (effects_bayesflow_v1/v2.csv, same
long format as effects_brms.csv). CPU only; keras backend is torch.
"""
import os
os.environ.setdefault("KERAS_BACKEND", "torch")
import argparse, json, time
from types import SimpleNamespace
import numpy as np
import pandas as pd


def base_parser(desc):
    p = argparse.ArgumentParser(description=desc)
    p.add_argument("--out", required=True, help="results dir written by steps 1-5")
    p.add_argument("--config", required=True, help="cohort config (e.g. config/drift_mtb.yaml)")
    p.add_argument("--epochs", type=int, default=50)
    p.add_argument("--batch_size", type=int, default=64)
    p.add_argument("--posterior_samples", type=int, default=2000)
    p.add_argument("--n_sbc", type=int, default=300, help="held-out simulations for calibration/recovery")
    p.add_argument("--seed", type=int, default=20260921)
    p.add_argument("--smoke", action="store_true", help="tiny run to check the pipeline end to end")
    p.add_argument("--restore", action="store_true", help="skip training, load the saved approximator from the checkpoint")
    p.add_argument("--chunk", type=int, default=25, help="test sets per posterior-sampling call (memory bound)")
    return p


def setup(a, tag):
    """Seeds, the cohort config and the real design (one row per neg tube with integer codes for
    patient and day, and K instrument dummies against cfg reference_instrument, one per
    non-reference instrument in order of appearance), plus the output folder bf_<tag>/ and a timer."""
    import yaml
    rng = np.random.default_rng(a.seed)
    try:  # keras is only needed for seeding here; lets the design be built without it
        import keras; keras.utils.set_random_seed(a.seed)
    except ImportError:
        pass
    with open(a.config) as fh:
        cfg = yaml.safe_load(fh)
    tubes = pd.read_csv(os.path.join(a.out, "tubes.csv"))
    negs = tubes[tubes.is_neg == True].copy()  # noqa: E712
    if os.path.exists(os.path.join(a.out, "excluded_tubes.csv")):
        negs = negs[~negs.filename.isin(pd.read_csv(os.path.join(a.out, "excluded_tubes.csv")).filename)]
    negs = negs.sort_values(["patient", "acq_day"]).reset_index(drop=True)
    negs["patient_i"] = pd.factorize(negs.patient)[0]
    negs["day_i"] = pd.factorize(negs.acq_day)[0]
    other = [i for i in negs.instrument.unique() if i != cfg["reference_instrument"]]
    instr = np.stack([(negs.instrument == lv).values.astype(float) for lv in other], axis=1) \
        if other else np.zeros((len(negs), 0))  # (N, K); K = 0 for a single-instrument cohort
    outdir = os.path.join(a.out, f"bf_{tag}"); os.makedirs(outdir, exist_ok=True)
    return SimpleNamespace(rng=rng, cfg=cfg, negs=negs, N=len(negs), P=negs.patient_i.max() + 1,
                           D=negs.day_i.max() + 1, pat=negs.patient_i.values, day=negs.day_i.values,
                           instr=instr, K_instr=len(other), instr_levels=other, days=list(pd.unique(negs.acq_day)),
                           instr_level=other[0] if other else "none", outdir=outdir, T=Timer())


def instr_names(instr_levels, prefix="b_instr"):
    """Parameter names for the recovery/sbc tables: plain 'b_instr' when there is one non-reference
    instrument (the existing file format), 'b_instr[<name>]' per level otherwise."""
    return [prefix] if len(instr_levels) == 1 else [f"{prefix}[{lv}]" for lv in instr_levels]


def effects_long(model, layer, response, instr_draws, day_draws, day_names, instr_levels):
    """Same columns as R effects_long(): model, layer, response, term, level, draw, value.
    instr_draws is (draws,) or (draws, K) with one block of instrument rows per entry of
    instr_levels (a name or a list of names); K = 0 gives day rows only."""
    if isinstance(instr_levels, str):
        instr_levels = [instr_levels]
    instr_draws = np.asarray(instr_draws).reshape(len(instr_draws), -1)
    rows = []
    for j, lv in enumerate(instr_levels):
        rows.append(pd.DataFrame(dict(model=model, layer=layer, response=response, term="instrument",
                                      level=lv, draw=np.arange(1, instr_draws.shape[0] + 1),
                                      value=instr_draws[:, j])))
    for j, name in enumerate(day_names):
        rows.append(pd.DataFrame(dict(model=model, layer=layer, response=response, term="day",
                                      level=name, draw=np.arange(1, day_draws.shape[0] + 1),
                                      value=day_draws[:, j])))
    return pd.concat(rows, ignore_index=True)


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


def validate(workflow, simulator, a, names, k_instr, D, tag):
    """Recovery and simulation-based calibration on held-out simulations. The checked parameters
    are the reported ones: b_instr (k_instr columns; none when k_instr = 0) and the day effects tau_d * c_d.
    Writes bf_<tag>_recovery.csv and bf_<tag>_sbc.csv; returns the recovery table and the test set."""
    from scipy.stats import chisquare
    test = simulator.sample(a.n_sbc)
    post = sample_chunked(workflow, test, a.posterior_samples, a.chunk)
    S = (a.n_sbc, a.posterior_samples)
    posterior = np.concatenate(([post["b_instr"].reshape(*S, k_instr)] if k_instr else []) +
                               [post["tau_d"].reshape(*S, 1) * post["c_d"].reshape(*S, D)], axis=-1)
    truth = np.concatenate(([test["b_instr"].reshape(-1, k_instr)] if k_instr else []) +
                           [test["tau_d"].reshape(-1, 1) * test["c_d"].reshape(-1, D)], axis=1)
    mean = posterior.mean(axis=1)
    ranks = (posterior < truth[:, None, :]).sum(axis=1)  # uniform over test sets if calibrated
    rec, sbc = [], []
    for j, n in enumerate(names):
        t, m = truth[:, j], mean[:, j]
        lo, hi = np.quantile(posterior[:, :, j], [0.025, 0.975], axis=1)
        rec.append(dict(parameter=n, r2=1 - ((t - m) ** 2).sum() / max(((t - t.mean()) ** 2).sum(), 1e-12),
                        coverage95=float(((t >= lo) & (t <= hi)).mean()), rmse=float(np.sqrt(((t - m) ** 2).mean()))))
        h, _ = np.histogram(ranks[:, j], bins=20, range=(0, a.posterior_samples + 1))
        sbc.append(dict(parameter=n, sbc_p=float(chisquare(h).pvalue)))  # > 0.05: no evidence of miscalibration
    rec = pd.DataFrame(rec); rec.to_csv(os.path.join(a.out, f"bf_{tag}_recovery.csv"), index=False)
    pd.DataFrame(sbc).to_csv(os.path.join(a.out, f"bf_{tag}_sbc.csv"), index=False)
    return rec, test


class Timer:
    def __init__(self): self.t0 = time.time(); self.marks = {}
    def mark(self, k): self.marks[k] = round(time.time() - self.t0, 1)
    def dump(self, path):
        with open(path, "w") as fh: json.dump(self.marks, fh, indent=1)
