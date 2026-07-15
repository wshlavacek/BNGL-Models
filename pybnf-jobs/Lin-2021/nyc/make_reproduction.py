#!/usr/bin/env python
"""Reproduction figure for the nyc job (Mallela et al. 2023 / Lin et al. 2021).

Simulates nyc.bngl over the first COVID-19 wave (model day t = 0..151; day 0 = 2020-01-21)
and overlays the model's DAILY new detected cases on the fit data (nyc.exp, t = 50..151).

The fit observable CumNum_detected_cases_Cum() is a running (cumulative-like) detected-case
count; PyBNF's `cumulative` flag differences it row-to-row into daily incidence before
scoring. This script reproduces that: it differences C(t) - C(t-1) on the full daily grid,
which is the physically meaningful daily incidence. (PyBNF scores on the data's own sample
grid [0, 50, 51, ..., 151], so its FIRST scored point t=50 is differenced against the forced
t=0 baseline -- i.e. C(50) - C(0) = the cumulative through day 50, a single-point edge
artifact -- whereas every point t=51..151 is the true daily increment C(t) - C(t-1). This
figure shows the true daily increment at every point, incl. t=50 = C(50) - C(49).)

By default it simulates at the .bngl nominals, which are an INDEPENDENT PyBNF `de` best fit
(this job's nyc.conf, with S0 corrected to NYC's census population). They are NOT the shipped
DataS1 MLE_params.txt, which does not reproduce the data with the placeholder S0/sigma (it
overshoots ~4x -- see VALIDATION.md). Pass --params "t0 beta lambda0 p0 fD" to override. The
negative-binomial dispersion r does not affect the mean trajectory, so this is deterministic.

Requires BNGPATH set (BNG2.pl) and matplotlib/numpy.
Usage: BNGPATH=... python make_reproduction.py [--params "27.37 0.671 0.105 0.752 0.262"]
"""
import argparse
import glob
import os
import re
import subprocess
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BNG = os.path.join(os.environ["BNGPATH"], "BNG2.pl")
MODEL = os.path.join(HERE, "nyc.bngl")
OBS = "CumNum_detected_cases_Cum"           # the fit observable (a model function)
T_END = 151
FREE = ["t0", "beta", "lambda0", "p0", "fD"]


def _cols(gdat_path):
    with open(gdat_path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(gdat_path)
    return data[:, 0], {name: data[:, i] for i, name in enumerate(header)}


def simulate(params=None):
    """One BNG2.pl run: generate the (finite, ~31-species) network and integrate the ODE
    from t=0 to 151 at every integer day (so C(t) - C(t-1) is the daily increment). Optional
    `params` dict overrides the model free-parameter nominals via setParameter."""
    with open(MODEL) as fh:
        src = fh.read().split("end model")[0] + "end model\n"
    overrides = ""
    if params:
        overrides = "".join(f'setParameter("{k}",{v})\n' for k, v in params.items())
    actions = (
        "begin actions\n"
        "generate_network({overwrite=>1})\n"
        f"{overrides}"
        f'simulate({{suffix=>"repro",method=>"ode",t_start=>0,t_end=>{T_END},'
        f"n_steps=>{T_END},print_functions=>1}})\n"
        "end actions\n"
    )
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "nyc.bngl")
        with open(path, "w") as fh:
            fh.write(src + "\n" + actions)
        r = subprocess.run(["perl", BNG, path], cwd=d, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"BNG2.pl failed (rc={r.returncode}):\n{r.stderr[-2000:]}")
        return _cols(glob.glob(os.path.join(d, "*_repro.gdat"))[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default=None,
                    help='override "t0 beta lambda0 p0 fD" (default: .bngl nominals = published MAP)')
    args = ap.parse_args()
    params = None
    label = "model @ PyBNF de best fit"     # the .bngl nominals ARE the de best fit (see README)
    if args.params:
        vals = [float(x) for x in args.params.split()]
        params = dict(zip(FREE, vals))
        label = "model @ override params"

    t, cols = simulate(params)
    C = cols[OBS]
    daily = np.diff(C)                       # daily new detected = C(t) - C(t-1), t = 1..151
    tday = t[1:]

    exp = np.loadtxt(os.path.join(HERE, "nyc.exp"))
    et, ey = exp[:, 0], exp[:, 1]

    # metrics over the fit window t = 50..151
    model_at_data = np.interp(et, tday, daily)
    rel = np.abs(model_at_data - ey) / ey
    med_rel, max_rel = float(np.median(rel)), float(rel.max())
    pk_i = int(np.argmax(daily)); dat_i = int(np.argmax(ey))
    print(f"{label}  (params: {params or 'nominals'})")
    print(f"  fit window t=50..151: median |rel err| = {med_rel*100:.1f}%   max = {max_rel*100:.1f}%")
    print(f"  model peak: {daily.max():.0f}/day at day {int(tday[pk_i])}  ({_date(int(tday[pk_i]))})")
    print(f"  data  peak: {ey.max():.0f}/day at day {int(et[dat_i])}  ({_date(int(et[dat_i]))})")
    print(f"  cumulative over window:  model {daily[(tday>=50)&(tday<=151)].sum():.0f}  vs  data {ey.sum():.0f}")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for a, logy in zip(ax, (False, True)):
        a.plot(et, ey, "o", ms=4, color="#222", label="NYT daily cases (data)", zorder=5)
        a.plot(tday, daily, "-", color="#d62728", lw=2, label=label)
        a.axvline(63, ls="--", color="#888", lw=1)
        a.text(63, a.get_ylim()[1] if not logy else 1, "  σ=63\n  (~Mar 24: NY PAUSE)",
               fontsize=8, color="#555", va="top")
        a.set(xlabel="model day (day 0 = 2020-01-21)", ylabel="daily new detected cases",
              title=("log scale" if logy else "linear scale"))
        if logy:
            a.set_yscale("log")
        a.legend(frameon=False, fontsize=9); a.grid(alpha=0.25)
    fig.suptitle(f"NYC MSA COVID-19 first wave -- {label}  "
                 f"(median rel err {med_rel*100:.0f}%)", fontsize=12)
    fig.tight_layout()
    out = os.path.join(HERE, "nyc_reproduction.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


def _date(day):
    import datetime
    return (datetime.date(2020, 1, 21) + datetime.timedelta(days=day)).isoformat()


if __name__ == "__main__":
    main()
