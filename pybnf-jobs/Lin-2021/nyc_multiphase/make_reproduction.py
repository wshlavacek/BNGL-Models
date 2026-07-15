#!/usr/bin/env python
"""Reproduction figure for the nyc_multiphase job (Lin et al. 2021, NYC "BigApple").

Simulates nyc_multiphase.bngl over the first COVID-19 wave (model day t = 0..170; day 0 =
2020-01-21) and overlays the model's DAILY new detected cases on the fit data
(nyc_multiphase.exp).

The fit observable fDCs_Cum is the true cumulative detected-case count (integral of the detection
rate); PyBNF's `cumulative` flag differences it row-to-row into daily incidence before scoring.
This script reproduces that as C(t) - C(t-1) on the full daily grid.

By default it simulates at the .bngl nominals, which ARE the authors' published MAP estimate
(BigApple/adaptive_files/MLE_params.txt). UNLIKE the single-phase nyc/ slug, this published MAP
reproduces the data -- a gold-standard reproduction. The two-phase social distancing captures the
broad NYC plateau (days ~68-85). Pass --params "t0 t_delta t_delta2 beta lambda0 p0 lambda1 p1 fD"
to override (e.g. an independent de fit). The neg-bin dispersion r does not affect the mean
trajectory, so this reproduction is deterministic.

Requires BNGPATH set (BNG2.pl) and matplotlib/numpy.
Usage: BNGPATH=... python make_reproduction.py [--params "..."]
"""
import argparse
import glob
import os
import subprocess
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BNG = os.path.join(os.environ["BNGPATH"], "BNG2.pl")
MODEL = os.path.join(HERE, "nyc_multiphase.bngl")
OBS = "fDCs_Cum"                       # the fit observable (a Molecules observable: cumulative counter)
T_END = 170
FREE = ["t0", "t_delta", "t_delta2", "beta", "lambda0", "p0", "lambda1", "p1", "fD"]


def _cols(gdat_path):
    with open(gdat_path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(gdat_path)
    return data[:, 0], {name: data[:, i] for i, name in enumerate(header)}


def simulate(params=None):
    """One BNG2.pl run: generate the finite network and integrate the ODE from t=0 to 170 at every
    integer day. Optional `params` dict overrides the model free-parameter nominals."""
    with open(MODEL) as fh:
        src = fh.read().split("end model")[0] + "end model\n"
    overrides = "".join(f'setParameter("{k}",{v})\n' for k, v in (params or {}).items())
    actions = (
        "begin actions\n"
        "generate_network({overwrite=>1})\n"
        f"{overrides}"
        f'simulate({{suffix=>"repro",method=>"ode",t_start=>0,t_end=>{T_END},'
        f"n_steps=>{T_END},print_functions=>1}})\n"
        "end actions\n"
    )
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "m.bngl")
        with open(path, "w") as fh:
            fh.write(src + "\n" + actions)
        r = subprocess.run(["perl", BNG, path], cwd=d, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"BNG2.pl failed (rc={r.returncode}):\n{r.stderr[-2000:]}")
        return _cols(glob.glob(os.path.join(d, "*_repro.gdat"))[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default=None,
                    help='override "t0 t_delta t_delta2 beta lambda0 p0 lambda1 p1 fD"')
    args = ap.parse_args()
    params, label = None, "model @ published MAP (2-phase)"
    if args.params:
        params = dict(zip(FREE, [float(x) for x in args.params.split()]))
        label = "model @ override params"

    t, cols = simulate(params)
    C = cols[OBS]
    daily, tday = np.diff(C), t[1:]

    exp = np.loadtxt(os.path.join(HERE, "nyc_multiphase.exp"))
    et, ey = exp[:, 0], exp[:, 1]
    keep = ey > 0                                        # skip the pre-epidemic zero rows for rel err
    model_at_data = np.interp(et, tday, daily)
    rel = np.abs(model_at_data[keep] - ey[keep]) / ey[keep]
    med_rel, max_rel = float(np.median(rel)), float(rel.max())
    pk_i, dat_i = int(np.argmax(daily)), int(np.argmax(ey))

    # phase-change days (from the .bngl nominals or overrides)
    p = params or dict(t0=32.831672, t_delta=0.072681, t_delta2=110.222846)
    sigma = p["t0"] + p["t_delta"]
    tau1 = sigma + p["t_delta2"]

    print(f"{label}  (params: {params or 'nominals'})")
    print(f"  fit window (ey>0): median |rel err| = {med_rel*100:.1f}%   max = {max_rel*100:.1f}%")
    print(f"  model peak: {daily.max():.0f}/day at day {int(tday[pk_i])}  ({_date(int(tday[pk_i]))})")
    print(f"  data  peak: {ey.max():.0f}/day at day {int(et[dat_i])}  ({_date(int(et[dat_i]))})")
    print(f"  phase-1 onset sigma={sigma:.1f} ({_date(round(sigma))}); phase-2 tau1={tau1:.1f} ({_date(round(tau1))})")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for a, logy in zip(ax, (False, True)):
        a.plot(et, ey, "o", ms=4, color="#222", label="NYT daily cases (data)", zorder=5)
        a.plot(tday, daily, "-", color="#1f77b4", lw=2, label=label)
        for x, name in [(sigma, "σ (phase 1)"), (tau1, "τ₁ (phase 2)")]:
            a.axvline(x, ls="--", color="#888", lw=1)
        a.set(xlabel="model day (day 0 = 2020-01-21)", ylabel="daily new detected cases",
              title=("log scale" if logy else "linear scale"))
        if logy:
            a.set_yscale("log"); a.set_ylim(bottom=0.5)
        a.legend(frameon=False, fontsize=9); a.grid(alpha=0.25)
    fig.suptitle(f"NYC MSA COVID-19 -- two-phase model @ published MAP  "
                 f"(median rel err {med_rel*100:.0f}%)", fontsize=12)
    fig.tight_layout()
    out = os.path.join(HERE, "nyc_multiphase_reproduction.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


def _date(day):
    import datetime
    return (datetime.date(2020, 1, 21) + datetime.timedelta(days=int(day))).isoformat()


if __name__ == "__main__":
    main()
