#!/usr/bin/env python
"""Step 7: BayesFlow variant 1 - tube-summary simulator on the REAL design.

Generative model (matches scripts/drift/stan/layer_a_normal.stan exactly):
    a_p ~ N(0, tau_p),  c_d ~ N(0, tau_d),  b_instr ~ N(0, 1),  alpha ~ N(0, 1)
    y_tube = alpha + a_p[patient] + c_d[day] + b_instr * instr + N(0, sigma)
for ONE standardised neg-tube statistic at a time (as in 06_yardstick_brms.R, layer A).
The observable is the vector of y over the real neg tubes plus the design codes, so the
summary network sees a set of (y, patient one-hot, day one-hot, instr) rows.

Inferred: b_instr and c_d[1..D] (plus alpha, tau_p, tau_d, sigma as nuisance).
Outputs effects_bayesflow_v1.csv (same long format as effects_brms.csv), bf_v1_recovery.csv,
bf_v1_sbc.csv, bf_v1_timing.json, and diagnostic PNGs under bf_v1/.

    python scripts/drift/07_bayesflow_v1.py --out results --config config/drift_mtb.yaml [--epochs 50 --num_sims 50000]
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bf_common import *  # noqa
import bayesflow as bf


def main():
    p = base_parser(__doc__)
    p.add_argument("--num_sims", type=int, default=50000, help="offline training simulations")
    a = p.parse_args()
    if a.smoke:
        a.epochs, a.num_sims, a.n_sbc, a.posterior_samples = 2, 512, 20, 200
    s = setup(a, "v1")
    rng, P, D, N, pat, day, instr = s.rng, s.P, s.D, s.N, s.pat, s.day, s.instr
    print(f"design: {N} neg tubes, {P} patients, {D} days, {int(instr.sum())} on {s.instr_level}")
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
        initial_learning_rate=5e-4, checkpoint_filepath=s.outdir, checkpoint_name="v1")
    if not maybe_restore(workflow, a.restore, os.path.join(s.outdir, "v1.keras")):
        train = simulator.sample(a.num_sims); s.T.mark("simulate")
        workflow.fit_offline(train, epochs=a.epochs, batch_size=a.batch_size,
                             validation_data=max(256, a.num_sims // 20))
    s.T.mark("train")

    names = ["b_instr"] + [f"c_d[{d}]" for d in s.days]
    rec, test = validate(workflow, simulator, a, names, 1, D, "v1")
    try:
        figs = workflow.plot_default_diagnostics(test_data=test, num_samples=min(500, a.posterior_samples))
        for k, f in figs.items(): f.savefig(os.path.join(s.outdir, f"{k}.png"), dpi=120)
    except Exception as e:  # diagnostics are optional
        print("diagnostic plots skipped:", e)
    s.T.mark("validate")

    # --- real data: every layer A statistic the yardstick fitted -----------------------------
    summ = pd.read_csv(os.path.join(a.out, "neg_summaries.csv"))
    emd = pd.read_csv(os.path.join(a.out, "neg_emd.csv"))
    eff = []
    for mk in s.cfg["technical_markers"]:
        st = summ[summ.marker == mk].set_index("filename")
        stats = ["median", "spread"] + (["valley", "frac_pos"] if st.valley.notna().any() else [])
        series = [(stat, st[stat]) for stat in stats] + [("emd", emd[emd.marker == mk].set_index("filename")["emd"])]
        for stat, ser in series:
            y = ser.reindex(s.negs.filename).values.astype(float)
            ok = np.isfinite(y)
            if ok.sum() < N:  # a missing statistic breaks the fixed design; skip it (step 6 drops those tubes)
                print(f"  {mk}:{stat}: {N - ok.sum()} missing values, statistic skipped"); continue
            sc, mu = y.std(), y.mean()
            obs = np.concatenate([((y - mu) / sc)[:, None].astype("float32"), design], axis=1)[None]
            ps = workflow.sample(conditions=dict(obs=obs), num_samples=a.posterior_samples)
            b = ps["b_instr"].reshape(-1) * sc
            cd = (ps["tau_d"].reshape(-1, 1) * ps["c_d"].reshape(-1, D)) * sc
            eff.append(effects_long("bayesflow_v1", "A", f"{mk}:{stat}", b, cd, s.days, s.instr_level))
    pd.concat(eff, ignore_index=True).to_csv(os.path.join(a.out, "effects_bayesflow_v1.csv"), index=False)
    s.T.mark("real_data"); s.T.dump(os.path.join(a.out, "bf_v1_timing.json"))
    print(rec.to_string()); print("done", s.T.marks)


if __name__ == "__main__":
    main()
