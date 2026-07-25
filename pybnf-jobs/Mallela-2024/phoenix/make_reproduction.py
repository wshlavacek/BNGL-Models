#!/usr/bin/env python
"""Reproduction figure for the phoenix job (Mallela et al. 2024, Phoenix MSA).

Simulates phoenix.bngl over the paper's fit window (model day t = 0..648; day 0 = 2020-01-21,
day 648 = 2021-10-30) and overlays the model's DAILY new detected cases on the fit data
(phoenix.exp) -- panel A of Fig. 5 in the paper.

Two curves are drawn:

  * **published MAP** -- the .bngl nominals, i.e. the authors' Table 1 estimate
    (Phoenix/Output/adaptive_files/MLE_params.txt), which the ENTIRE 19-22-parameter
    adaptive-MCMC run produced;
  * **L-BFGS-B fit** -- this job's result, which re-estimates only beta and fD from the
    published MAP and leaves every other parameter (including every switch time) pinned.

The fit observable fDCs_Cum is the true cumulative detected-case count fD*C_S(t) (Eqs. 39-40);
PyBNF's per-observable `cumulative` flag differences it row-to-row into daily incidence before
scoring, which this script reproduces as C(t) - C(t-1) on the daily grid. The headline metric is
the job's own objective -- the negative-binomial NLL at the published dispersion r = 3.116234553368695
-- reported for both parameter sets, alongside relative-error / peak / cumulative summaries.

The neg-bin dispersion does not affect the mean trajectory, so this reproduction is
deterministic. Simulation goes through BNG2.pl (not bngsim) so the figure is reproducible
without the gradient toolchain.

Requires BNGPATH set (BNG2.pl) and matplotlib/numpy/scipy.
Usage: BNGPATH=... python make_reproduction.py [--params "beta fD"]
"""
import argparse
import glob
import os
import subprocess
import tempfile

import numpy as np
from scipy.special import gammaln
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BNG = os.path.join(os.environ["BNGPATH"], "BNG2.pl")
MODEL = os.path.join(HERE, "phoenix.bngl")
OBS = "fDCs_Cum"                 # the fit observable (a Molecules observable: cumulative counter)
T_END = 648                      # 2021-10-30; day 0 = 2020-01-21
FREE = ["beta", "fD"]

R_DISP = 3.116234553368695          # published negative-binomial dispersion (Table 1)
MAP = dict(beta=0.285816836838281, fD=0.538318135897164)   # published MAP (== .bngl nominals)
FIT = dict(beta=0.2931033097285776, fD=0.4590247441033879)   # this job's L-BFGS-B result
T0 = 13.95001678247903             # start of local transmission (pinned); the model predicts 0 cases before it


def _cols(gdat_path):
    with open(gdat_path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(gdat_path)
    return data[:, 0], {name: data[:, i] for i, name in enumerate(header)}


def simulate(params=None):
    """One BNG2.pl run: generate the finite network (42 species, 88 reactions) and integrate the
    ODE from t=0 to 648 at every integer day. `params` overrides the free-parameter nominals."""
    with open(MODEL) as fh:
        src = fh.read().split("end model")[0] + "end model\n"
    overrides = "".join('setParameter("%s",%r)\n' % (k, v) for k, v in (params or {}).items())
    actions = (
        "begin actions\n"
        "generate_network({overwrite=>1})\n"
        + overrides +
        'simulate({suffix=>"repro",method=>"ode",t_start=>0,t_end=>%d,'
        "n_steps=>%d,print_functions=>1})\n" % (T_END, T_END) +
        "end actions\n"
    )
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "m.bngl")
        with open(path, "w") as fh:
            fh.write(src + "\n" + actions)
        r = subprocess.run(["perl", BNG, path], cwd=d, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError("BNG2.pl failed (rc=%d):\n%s" % (r.returncode, r.stderr[-2000:]))
        return _cols(glob.glob(os.path.join(d, "*_repro.gdat"))[0])


def neg_bin_nll(pred, obs, r=R_DISP):
    """The job's objective, exactly as PyBNF scores it: a negative-binomial -logpmf per point
    with the prediction taken as the MEAN, summed. Two conventions matter and are PyBNF's own
    (pybnf/noise/negative_binomial.py): `prob` is clipped, which is what keeps the pre-t0 rows
    (prediction exactly 0) finite; and a NEGATIVE observation contributes nothing (the
    count-domain guard) -- the NYT series carries a few negative days where a state revised its
    cumulative count downward."""
    keep = obs >= 0
    pred, obs = pred[keep], obs[keep]
    prob = np.clip(r / (r + pred), 1e-10, 1 - 1e-10)
    log_pmf = (gammaln(obs + r) - gammaln(obs + 1) - gammaln(r)
               + r * np.log(prob) + obs * np.log1p(-prob))
    return float(-log_pmf.sum())


def daily(params):
    """Model daily new detected cases on the .exp grid: the cumulative observable differenced
    row-to-row, with row 0 keeping its raw value (PyBNF's `cumulative` convention)."""
    t, cols = simulate(params)
    C = cols[OBS]
    return t, np.concatenate(([C[0]], np.diff(C)))


def smooth7(y):
    """Centered 7-day rolling mean -- the standard smoothing for daily COVID surveillance
    counts, which carry a strong day-of-week reporting cycle plus occasional revision spikes
    and negative-count corrections. The model predicts the smooth MEAN, so the smoothed series
    is the fair target for a shape comparison (the NLL below is scored on the RAW counts, as
    the fit itself is)."""
    k = np.ones(7) / 7.0
    return np.convolve(y, k, mode="same")


def summarize(label, pred, et, ey):
    keep = ey > 0                                  # rel err is undefined at a zero/negative count
    rel = np.abs(pred[keep] - ey[keep]) / ey[keep]
    sy, sp = smooth7(ey), smooth7(pred)
    keep7 = sy > 1
    rel7 = np.abs(sp[keep7] - sy[keep7]) / sy[keep7]
    pk, dpk = int(np.argmax(sp)), int(np.argmax(sy))
    print("  %s" % label)
    print("    neg_bin NLL (the fit objective)  = %.2f" % neg_bin_nll(pred, ey))
    print("    median |rel err| vs raw counts   = %.1f%%   (n = %d days with a positive count)"
          % (100 * np.median(rel), keep.sum()))
    print("    median |rel err| vs 7-day mean   = %.1f%%   (n = %d days)"
          % (100 * np.median(rel7), keep7.sum()))
    print("    peak of the 7-day mean: model %7.0f/day at day %3d (%s)"
          % (sp.max(), et[pk], _date(et[pk])))
    print("                            data  %7.0f/day at day %3d (%s)   -> %+.1f%%"
          % (sy.max(), et[dpk], _date(et[dpk]), 100 * (sp.max() - sy.max()) / sy.max()))
    print("    cumulative reported cases        = %.0f   (data %.0f, %+.1f%%)"
          % (pred.sum(), ey.sum(), 100 * (pred.sum() - ey.sum()) / ey.sum()))


def _date(day):
    import datetime
    return (datetime.date(2020, 1, 21) + datetime.timedelta(days=int(day))).isoformat()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default=None, help='override the fitted "beta fD" pair')
    args = ap.parse_args()
    fit = dict(FIT)
    if args.params:
        fit = dict(zip(FREE, [float(x) for x in args.params.split()]))

    exp = np.loadtxt(os.path.join(HERE, "phoenix.exp"))
    et, ey = exp[:, 0], exp[:, 1]

    t_map, pred_map = daily(MAP)
    t_fit, pred_fit = daily(fit)
    pred_map = np.interp(et, t_map, pred_map)
    pred_fit = np.interp(et, t_fit, pred_fit)

    print("Phoenix MSA -- Mallela et al. 2024, Fig. 5A  (r = %g)" % R_DISP)
    summarize("published MAP   beta=%.6f fD=%.6f" % (MAP["beta"], MAP["fD"]), pred_map, et, ey)
    summarize("L-BFGS-B fit    beta=%.6f fD=%.6f" % (fit["beta"], fit["fD"]), pred_fit, et, ey)
    print("  Delta NLL (fit - MAP) = %+.2f  (negative = the 2-parameter fit improves on the MAP)"
          % (neg_bin_nll(pred_fit, ey) - neg_bin_nll(pred_map, ey)))

    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5))
    for a, logy in zip(ax, (False, True)):
        a.plot(et, ey, "o", ms=2.5, color="#444", alpha=.65,
               label="NYT daily cases (data)", zorder=1)
        a.plot(et, pred_map, "-", color="#888", lw=1.6, label="model @ published MAP", zorder=2)
        a.plot(et, pred_fit, "-", color="#1f77b4", lw=1.8,
               label="model @ L-BFGS-B fit (beta, fD)", zorder=3)
        a.axvline(T0, ls="--", color="#bbb", lw=1)
        a.set(xlabel="model day (day 0 = 2020-01-21)", ylabel="daily new detected cases",
              title=("log scale" if logy else "linear scale"))
        if logy:
            a.set_yscale("log")
            a.set_ylim(bottom=0.5)
        a.legend(frameon=False, fontsize=9)
        a.grid(alpha=0.25)
    fig.suptitle("Phoenix MSA COVID-19 (Mallela et al. 2024, Fig. 5A)  --  "
                 "neg_bin NLL %.0f (MAP) -> %.0f (fit)"
                 % (neg_bin_nll(pred_map, ey), neg_bin_nll(pred_fit, ey)), fontsize=12)
    fig.tight_layout()
    out = os.path.join(HERE, "phoenix_reproduction.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
