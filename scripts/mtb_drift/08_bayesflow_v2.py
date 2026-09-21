#!/usr/bin/env python
"""Step 8: BayesFlow variant 2 - single-cell distribution simulator, hierarchical summary network.

Generative model per simulated experiment, on the REAL neg-tube design (patients, days,
instruments) and the technical markers of config/mtb_drift.yaml, all on the arcsinh scale:

  per marker k and patient p: mixture location(s) mu_pk (bimodal markers: two components with
  patient-specific positive fraction), component SDs s_k
  day drift:        c_d ~ N(0, tau_d) shared across markers, marker loading lam_k (lam_1 = 1)
  instrument shift: b_instr_k ~ N(0, 1) per marker, additive; spread gain g_instr ~ logN(0, 0.2)
  cells:            n_cells per tube, drawn from the shifted mixture, plus 2 percent uniform debris

Observable: (n_tubes, n_cells, K + 3) - cell values and compact tube codes (patient, day, instr).
Summary network: DeepSet over cells -> tube embedding, concatenated with tube codes,
SetTransformer over tubes -> summary. Trained ONLINE (simulate on the fly), CPU.

Inferred and reported: b_instr_k per marker (comparable to the layer A 'median' instrument
effect for marker k) and c_d per day (comparable to the layer A 'median' day effects, scaled by
the marker loading). Writes effects_bayesflow_v2.csv, bf_v2_recovery.csv, bf_v2_sbc.csv,
bf_v2_ppc.csv (posterior predictive check on the real neg medians) and bf_v2_timing.json.

    python scripts/mtb_drift/08_bayesflow_v2.py --out results [--n_cells 1000 --epochs 40]
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bf_common import *  # noqa
import keras
from keras import ops
import bayesflow as bf


class CellsInTubes(bf.networks.SummaryNetwork):
    """x: (batch, n_tubes, n_cells, K + n_codes). Cells are pooled per tube by a DeepSet,
    the tube embedding is joined with the tube codes, and a SetTransformer pools the tubes."""

    def __init__(self, n_markers, summary_dim=32, cell_dim=16, **kwargs):
        super().__init__(**kwargs)
        self.n_markers = n_markers
        self.cells = bf.networks.DeepSet(summary_dim=cell_dim, embed_dim=32, mlp_widths=(64,), depth=1)
        self.tubes = bf.networks.SetTransformer(summary_dim=summary_dim, embed_dims=(64, 64), num_seeds=4)

    def build(self, input_shape):
        b, n, c, f = input_shape
        self.cells.build((None, c, self.n_markers))
        self.tubes.build((None, n, self.cells.summary_dim + f - self.n_markers))
        super().build(input_shape)

    def call(self, x, training=False, **kwargs):
        shp = ops.shape(x)
        b, n, c, f = shp[0], shp[1], shp[2], shp[3]
        cells = ops.reshape(x[..., :self.n_markers], (b * n, c, self.n_markers))
        e = self.cells(cells, training=training)
        e = ops.reshape(e, (b, n, self.cells.summary_dim))
        codes = x[:, :, 0, self.n_markers:]
        return self.tubes(ops.concatenate([e, codes], axis=-1), training=training)


def main():
    p = base_parser(__doc__)
    p.add_argument("--n_cells", type=int, default=1000)
    p.add_argument("--batches_per_epoch", type=int, default=200)
    a = p.parse_args()
    if a.smoke:
        a.epochs, a.batches_per_epoch, a.n_sbc, a.posterior_samples, a.n_cells, a.batch_size = 1, 4, 8, 100, 200, 4
    rng = np.random.default_rng(a.seed); keras.utils.set_random_seed(a.seed)
    cfg = read_config()
    negs = read_design(a.out)
    markers = [m for m in cfg["technical_markers"] if m not in ("FSC_A", "SSC_A")]  # fluorescence only
    bimodal = [m in cfg["bimodal_markers"] for m in markers]
    K, N = len(markers), len(negs)
    P, D = negs.patient_i.max() + 1, negs.day_i.max() + 1
    pat, day, instr = negs.patient_i.values, negs.day_i.values, negs.instr.values
    days = day_levels(negs)
    other = [i for i in negs.instrument.unique() if i != cfg["reference_instrument"]]
    instr_level = other[0] if other else "none"
    codes = np.stack([pat / max(P - 1, 1), day / max(D - 1, 1), instr], axis=1).astype("float32")  # (N, 3)
    outdir = os.path.join(a.out, "bf_v2"); os.makedirs(outdir, exist_ok=True)
    T = Timer()
    print(f"design: {N} tubes, {P} patients, {D} days, K={K} markers {markers}, {a.n_cells} cells/tube")

    # --- simulator -------------------------------------------------------------------------
    def prior():
        return dict(b_instr=rng.normal(0, 1, K), tau_d=abs(rng.normal(0, 0.5)), c_d=rng.normal(0, 1, D),
                    lam=np.concatenate([[1.0], rng.normal(1, 0.3, K - 1)]),
                    log_g_instr=rng.normal(0, 0.2))

    def likelihood(b_instr, tau_d, c_d, lam, log_g_instr):
        # patient biology: component means, positive fractions, component SDs
        mu_hi = rng.normal(3.5, 0.4, (P, K)); mu_lo = mu_hi - rng.uniform(2.0, 3.0, (P, K))
        frac = rng.beta(4, 5, (P, K)); sd = rng.uniform(0.3, 0.5, K)
        shift = (tau_d * c_d[day])[:, None] * lam[None, :] + instr[:, None] * b_instr[None, :]   # (N, K)
        gain = np.where(instr > 0, np.exp(log_g_instr), 1.0)[:, None]                            # (N, 1)
        cells = np.empty((N, a.n_cells, K), dtype="float32")
        for k in range(K):
            if bimodal[k]:
                pos = rng.random((N, a.n_cells)) < frac[pat, k][:, None]
                loc = np.where(pos, mu_hi[pat, k][:, None], mu_lo[pat, k][:, None])
            else:
                loc = mu_hi[pat, k][:, None] - 1.5
            cells[:, :, k] = loc + shift[:, k][:, None] + rng.normal(0, 1, (N, a.n_cells)) * sd[k] * gain
        debris = rng.random((N, a.n_cells)) < 0.02
        cells[debris] = rng.uniform(-1, 6, (debris.sum(), K)).astype("float32")
        obs = np.concatenate([cells, np.broadcast_to(codes[:, None, :], (N, a.n_cells, 3))], axis=-1)
        return dict(obs=obs.astype("float32"))

    simulator = bf.make_simulator([prior, likelihood])
    adapter = (bf.Adapter().convert_dtype("float64", "float32")
               .constrain("tau_d", lower=0)
               .concatenate(["b_instr", "tau_d", "c_d", "lam", "log_g_instr"], into="inference_variables")
               .rename("obs", "summary_variables"))
    workflow = bf.BasicWorkflow(
        simulator=simulator, adapter=adapter,
        inference_network=bf.networks.CouplingFlow(depth=6),
        summary_network=CellsInTubes(n_markers=K, summary_dim=48),
        initial_learning_rate=5e-4, checkpoint_filepath=outdir, checkpoint_name="v2")
    workflow.fit_online(epochs=a.epochs, num_batches_per_epoch=a.batches_per_epoch, batch_size=a.batch_size)
    T.mark("train")

    # --- calibration and recovery -----------------------------------------------------------
    test = simulator.sample(a.n_sbc)
    post = workflow.sample(conditions=test, num_samples=a.posterior_samples)
    names = [f"b_instr[{m}]" for m in markers] + [f"c_d[{d}]" for d in days]
    post_arr = np.concatenate([post["b_instr"].reshape(a.n_sbc, a.posterior_samples, K),
                               post["tau_d"].reshape(a.n_sbc, a.posterior_samples, 1) *
                               post["c_d"].reshape(a.n_sbc, a.posterior_samples, D)], axis=-1)
    truth = np.concatenate([test["b_instr"].reshape(-1, K), test["tau_d"].reshape(-1, 1) * test["c_d"].reshape(-1, D)], axis=1)
    rec = recovery_table(post_arr, truth, names); rec.to_csv(os.path.join(a.out, "bf_v2_recovery.csv"), index=False)
    pd.DataFrame(dict(parameter=names, sbc_p=sbc_uniformity(sbc_ranks(post_arr, truth), a.posterior_samples))).to_csv(
        os.path.join(a.out, "bf_v2_sbc.csv"), index=False)
    T.mark("validate")

    # --- real neg cells ------------------------------------------------------------------------
    real = np.empty((N, a.n_cells, K), dtype="float32")
    for i, fn in enumerate(negs.filename):
        df = pd.read_csv(os.path.join(a.out, "neg_cells", fn + ".csv.gz"))
        m = df.reindex(columns=markers).values.astype("float32")
        m = m[np.isfinite(m).all(axis=1)]
        idx = rng.choice(len(m), a.n_cells, replace=len(m) < a.n_cells)
        real[i] = m[idx]
    obs = np.concatenate([real, np.broadcast_to(codes[:, None, :], (N, a.n_cells, 3))], axis=-1)[None]
    ps = workflow.sample(conditions=dict(obs=obs), num_samples=a.posterior_samples)
    b = ps["b_instr"].reshape(-1, K); cd = ps["tau_d"].reshape(-1, 1) * ps["c_d"].reshape(-1, D)
    lam = ps["lam"].reshape(-1, K)
    eff = []
    for k, mk in enumerate(markers):   # day effect on marker k = lam_k * c_d, so it matches the per-marker yardstick
        eff.append(effects_long("bayesflow_v2", "A", f"{mk}:median", b[:, k], cd * lam[:, [k]], days, instr_level))
    eff = pd.concat(eff, ignore_index=True); eff.to_csv(os.path.join(a.out, "effects_bayesflow_v2.csv"), index=False)

    # posterior predictive: simulate from 100 posterior draws, compare real per-tube medians
    real_med = np.median(real, axis=1)                                      # (N, K)
    sims = []
    for j in rng.choice(a.posterior_samples, min(100, a.posterior_samples), replace=False):
        s = likelihood(b[j], float(ps["tau_d"].reshape(-1)[j]), ps["c_d"].reshape(-1, D)[j], lam[j],
                       float(ps["log_g_instr"].reshape(-1)[j]))["obs"][..., :K]
        sims.append(np.median(s, axis=1))
    sims = np.stack(sims)                                                    # (S, N, K)
    lo, hi = np.quantile(sims, [0.05, 0.95], axis=0)
    ppc = pd.DataFrame(dict(filename=np.repeat(negs.filename.values, K), marker=np.tile(markers, N),
                            real_median=real_med.ravel(), pred_lo=lo.ravel(), pred_hi=hi.ravel()))
    ppc["inside_90"] = (ppc.real_median >= ppc.pred_lo) & (ppc.real_median <= ppc.pred_hi)
    ppc.to_csv(os.path.join(a.out, "bf_v2_ppc.csv"), index=False)
    T.mark("real_data"); T.dump(os.path.join(a.out, "bf_v2_timing.json"))
    print(rec.to_string()); print(f"PPC: {ppc.inside_90.mean():.2f} of real medians inside the 90% band"); print(T.marks)


if __name__ == "__main__":
    main()
