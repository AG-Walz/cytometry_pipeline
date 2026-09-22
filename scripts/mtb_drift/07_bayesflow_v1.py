#!/usr/bin/env python
"""Step 7: BayesFlow variant 1 - tube-summary simulator on the REAL design.

Generative model (matches scripts/mtb_drift/stan/layer_a_normal.stan exactly):
    a_p ~ N(0, tau_p),  c_d ~ N(0, tau_d),  b_instr ~ N(0, 1),  alpha ~ N(0, 1)
    y_tube = alpha + a_p[patient] + c_d[day] + b_instr * instr + N(0, sigma)
for ONE standardised neg-tube statistic at a time (as in 06_yardstick_brms.R, layer A).
The observable is the vector of y over the real neg tubes plus the design codes, so the
summary network sees a set of (y, patient one-hot, day one-hot, instr) rows.

Inferred: b_instr and c_d[1..D] (plus alpha, tau_p, tau_d, sigma as nuisance).
Outputs effects_bayesflow_v1.csv (same long format as effects_brms.csv), bf_v1_recovery.csv,
bf_v1_sbc.csv, bf_v1_timing.json, and diagnostic PNGs under bf_v1/.

    python scripts/mtb_drift/07_bayesflow_v1.py --out results [--epochs 50 --num_sims 50000]
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bf_common import *  # noqa
import keras
import bayesflow as bf


def main():
    a = base_parser(__doc__).parse_args()
    if a.smoke:
        a.epochs, a.num_sims, a.n_sbc, a.posterior_samples = 2, 512, 20, 200
    rng = np.random.default_rng(a.seed); keras.utils.set_random_seed(a.seed)
    cfg = read_config()
    negs = read_design(a.out)
    P, D, N = negs.patient_i.max() + 1, negs.day_i.max() + 1, len(negs)
    pat, day, instr = negs.patient_i.values, negs.day_i.values, negs.instr.values
    days = day_levels(negs)
    other = [i for i in negs.instrument.unique() if i != cfg["reference_instrument"]]
    instr_level = other[0] if other else "none"
    outdir = os.path.join(a.out, "bf_v1"); os.makedirs(outdir, exist_ok=True)
    T = Timer()
    print(f"design: {N} neg tubes, {P} patients, {D} days, {int(instr.sum())} on {instr_level}")

    # design matrix is fixed, so it is a constant condition: one-hot patient & day, instr
    design = np.concatenate([np.eye(P)[pat], np.eye(D)[day], instr[:, None]], axis=1).astype("float32")

    def prior():
        return dict(alpha=rng.normal(0, 1), b_instr=rng.normal(0, 1),
                    tau_p=abs(rng.normal(0, 1)), tau_d=abs(rng.normal(0, 1)), sigma=abs(rng.normal(0, 1)),
                    c_d=rng.normal(0, 1, D))   # standard normal; scaled by tau_d in the likelihood

    def likelihood(alpha, b_instr, tau_p, tau_d, sigma, c_d):
        a_p = rng.normal(0, 1, P) * tau_p
        y = alpha + a_p[pat] + tau_d * c_d[day] + b_instr * instr + rng.normal(0, sigma, N)
        return dict(obs=np.concatenate([y[:, None].astype("float32"), design], axis=1))

    simulator = bf.make_simulator([prior, likelihood])

    # c_d is inferred in its standardised form; the reported day effect is tau_d * c_d
    adapter = (bf.Adapter()
               .convert_dtype("float64", "float32")
               .constrain(["tau_p", "tau_d", "sigma"], lower=0)
               .concatenate(["alpha", "b_instr", "tau_p", "tau_d", "sigma", "c_d"], into="inference_variables")
               .rename("obs", "summary_variables"))

    workflow = bf.BasicWorkflow(
        simulator=simulator, adapter=adapter,
        inference_network=bf.networks.CouplingFlow(depth=6),
        summary_network=bf.networks.SetTransformer(summary_dim=32, embed_dims=(64, 64), num_seeds=4),
        initial_learning_rate=5e-4, checkpoint_filepath=outdir, checkpoint_name="v1")
    if not maybe_restore(workflow, a.restore, os.path.join(outdir, "v1.keras")):
        train = simulator.sample(a.num_sims); T.mark("simulate")
        workflow.fit_offline(train, epochs=a.epochs, batch_size=a.batch_size,
                             validation_data=max(256, a.num_sims // 20))
    T.mark("train")

    # --- calibration and recovery on held-out simulations ---------------------------------
    test = simulator.sample(a.n_sbc)
    post = sample_chunked(workflow, test, a.posterior_samples, a.chunk)
    names = ["b_instr"] + [f"c_d[{d}]" for d in days]
    post_arr = np.concatenate([post["b_instr"].reshape(a.n_sbc, a.posterior_samples, 1),
                               post["tau_d"].reshape(a.n_sbc, a.posterior_samples, 1) *
                               post["c_d"].reshape(a.n_sbc, a.posterior_samples, D)], axis=-1)
    truth = np.concatenate([test["b_instr"].reshape(-1, 1), test["tau_d"].reshape(-1, 1) * test["c_d"].reshape(-1, D)], axis=1)
    rec = recovery_table(post_arr, truth, names); rec.to_csv(os.path.join(a.out, "bf_v1_recovery.csv"), index=False)
    ranks = sbc_ranks(post_arr, truth)
    pd.DataFrame(dict(parameter=names, sbc_p=sbc_uniformity(ranks, a.posterior_samples))).to_csv(
        os.path.join(a.out, "bf_v1_sbc.csv"), index=False)
    try:
        figs = workflow.plot_default_diagnostics(test_data=test, num_samples=min(500, a.posterior_samples))
        for k, f in figs.items(): f.savefig(os.path.join(outdir, f"{k}.png"), dpi=120)
    except Exception as e:  # diagnostics are optional
        print("diagnostic plots skipped:", e)
    T.mark("validate")

    # --- real data: every layer A statistic the yardstick fitted -----------------------------
    summ = pd.read_csv(os.path.join(a.out, "neg_summaries.csv"))
    emd = pd.read_csv(os.path.join(a.out, "neg_emd.csv"))
    long = []
    for mk in cfg["technical_markers"]:
        s = summ[summ.marker == mk]
        long += [("median", s.set_index("filename")["median"]), ("spread", s.set_index("filename")["spread"])]
        if "valley" in s and s.valley.notna().any():
            long += [("valley", s.set_index("filename")["valley"]), ("frac_pos", s.set_index("filename")["frac_pos"])]
        e = emd[emd.marker == mk].set_index("filename")["emd"]
        long.append(("emd", e))
        for stat, ser in long:
            y = ser.reindex(negs.filename).values.astype(float)
            ok = np.isfinite(y)
            if ok.sum() < N:  # a missing statistic breaks the fixed design; impute with the mean and say so
                print(f"  {mk}:{stat}: {N - ok.sum()} missing values imputed with the mean")
                y[~ok] = np.nanmean(y)
            sc, mu = y.std(), y.mean()
            obs = np.concatenate([((y - mu) / sc)[:, None].astype("float32"), design], axis=1)[None]
            ps = workflow.sample(conditions=dict(obs=obs), num_samples=a.posterior_samples)
            b = ps["b_instr"].reshape(-1) * sc
            cd = (ps["tau_d"].reshape(-1, 1) * ps["c_d"].reshape(-1, D)) * sc
            df = effects_long("bayesflow_v1", "A", f"{mk}:{stat}", b, cd, days, instr_level)
            df.to_csv(os.path.join(outdir, f"effects_{mk}_{stat}.csv"), index=False)
        long = []
    eff = pd.concat([pd.read_csv(os.path.join(outdir, f)) for f in os.listdir(outdir) if f.startswith("effects_")])
    eff.to_csv(os.path.join(a.out, "effects_bayesflow_v1.csv"), index=False)
    T.mark("real_data"); T.dump(os.path.join(a.out, "bf_v1_timing.json"))
    print(rec.to_string()); print("done", T.marks)


if __name__ == "__main__":
    main()
