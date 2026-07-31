#!/usr/bin/env python
"""Reproduction figure for the Zhang-2023 tumor_growth job.

Simulates tumor_growth.bngl over the 45-day window of Fig. 4D in Zhang et al. (2023) for
each of the four treatment arms and overlays the result on the digitized Bridgeman et al.
(2016) xenograft data (`*.exp`) and on the authors' own plotted model curves.

Three parameter sets are compared:

  * **File S4** -- the authors' deposited nominals (the .bngl values), with the untreated
    Erk = Akt = 1;
  * **Table S1** -- the values Table S1 of the paper reports for the same model, with
    kkill = 0 because neither source gives one;
  * **this fit** -- the differential-evolution result of tumor_growth.conf, which fits all
    six parameters the Methods name as fitted against all four arms.

The headline metric is the job's own objective (PyBNF reports `sos` as half the residual
sum of squares) plus, per arm, the RMSE in fold units and the median relative error. Fold
units are the honest scale for the combination arm, whose data hover near 0.3-fold, where a
relative error divides by a small number.

Simulation goes through BNG2.pl, so the figure reproduces without the bngsim toolchain.

Requires BNGPATH set (BNG2.pl) and matplotlib/numpy.
Usage: BNGPATH=... python make_reproduction.py
"""
import os
import subprocess
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BNG = os.path.join(os.environ["BNGPATH"], "BNG2.pl")
MODEL = os.path.join(HERE, "tumor_growth.bngl")
OBS = "Obs_Cells"
T_END = 45.0

# Erk / Akt per arm -- the normalized peak endothelial pERK1/2 and ppAKT of Fig. 4B and
# Fig. 4C, digitized from the simulation bars (see the model folder's reference/ CSV).
ARMS = {
    "control":    dict(Erk=1.000, Akt=1.000),
    "sunitinib":  dict(Erk=0.736, Akt=0.793),
    "trametinib": dict(Erk=0.580, Akt=1.000),
    "combo":      dict(Erk=0.355, Akt=0.793),
}

SETS = {
    "File S4 (nominals)": {},
    "Table S1 (kkill = 0)": dict(w_OR=0.3, kTD=29.2812, EC50TD=0.5840,
                                 kg=0.2993, klinear=1.9424, kkill=0.0),
    "this fit (all 4 arms)": dict(w_OR=0.3422752725722071, kTD=30.505087975036044,
                                  EC50TD=0.5732582945654124, kg=0.13880486227968328,
                                  klinear=8.943839880074522, kkill=0.026410435081856903),
}
COLORS = {"File S4 (nominals)": "0.55", "Table S1 (kkill = 0)": "tab:orange",
          "this fit (all 4 arms)": "tab:blue"}


def simulate(overrides):
    """Run tumor_growth.bngl with parameter overrides; return (t, Obs_Cells)."""
    text = open(MODEL).read()
    text = text[:text.index("end model") + len("end model")]
    setp = "".join(f'  setParameter("{k}",{v!r})\n' for k, v in overrides.items())
    text += ("\n\nbegin actions\n" + setp +
             "  generate_network({overwrite=>1})\n"
             f'  simulate({{method=>"ode",suffix=>"ode",t_start=>0,t_end=>{T_END},'
             "n_steps=>900,atol=>1e-10,rtol=>1e-8})\nend actions\n")
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "m.bngl")
        open(path, "w").write(text)
        res = subprocess.run(["perl", BNG, "m.bngl"], cwd=tmp,
                             capture_output=True, text=True)
        assert res.returncode == 0, res.stderr or res.stdout[-2000:]
        gdat = os.path.join(tmp, "m_ode.gdat")
        with open(gdat) as fh:
            header = fh.readline().lstrip("#").split()
        d = np.loadtxt(gdat)
    return d[:, 0], d[:, header.index(OBS)]


def load_exp(arm):
    d = np.loadtxt(os.path.join(HERE, f"{arm}.exp"))
    return np.atleast_2d(d)


def main():
    exp = {a: load_exp(a) for a in ARMS}
    curves, metrics = {}, {}
    for sname, over in SETS.items():
        curves[sname] = {}
        sse_tot, n_tot = 0.0, 0
        rows = []
        for arm, sig in ARMS.items():
            t, y = simulate({**sig, **over})
            curves[sname][arm] = (t, y)
            obs = exp[arm]
            pred = np.interp(obs[:, 0], t, y)
            resid = pred - obs[:, 1]
            sse = float((resid ** 2).sum())
            sse_tot += sse
            n_tot += len(obs)
            rows.append((arm, sse, np.sqrt(sse / len(obs)),
                         float(np.median(np.abs(resid) / obs[:, 1]))))
        metrics[sname] = (rows, sse_tot, n_tot)

    print(f"{'parameter set':24s} {'arm':11s} {'SSE':>8s} {'RMSE (fold)':>12s} "
          f"{'median |rel err|':>17s}")
    for sname, (rows, sse_tot, n_tot) in metrics.items():
        for arm, sse, rmse, med in rows:
            print(f"{sname:24s} {arm:11s} {sse:8.3f} {rmse:12.3f} {med:17.3f}")
        print(f"{sname:24s} {'ALL':11s} {sse_tot:8.3f} "
              f"{np.sqrt(sse_tot/n_tot):12.3f}   (pybnf sos = {sse_tot/2:.3f}, "
              f"n = {n_tot})")
        print()

    fig, axes = plt.subplots(1, 4, figsize=(16, 4), sharey=True)
    for ax, arm in zip(axes, ARMS):
        obs = exp[arm]
        for sname in SETS:
            t, y = curves[sname][arm]
            ax.plot(t, y, lw=1.6, color=COLORS[sname], label=sname)
        ax.plot(obs[:, 0], obs[:, 1], "o", mfc="none", mec="k", ms=7,
                label="Bridgeman et al. (2016)")
        ax.set_title(f"{arm}  (Erk = {ARMS[arm]['Erk']:.3f}, "
                     f"Akt = {ARMS[arm]['Akt']:.3f})", fontsize=10)
        ax.set_xlabel("time (day)")
        ax.set_xlim(0, T_END)
        ax.set_ylim(0, 9)
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("tumor growth (fold)")
    axes[0].legend(frameon=False, fontsize=8, loc="upper left")
    fig.suptitle("Angiogenesis-driven tumor growth (Zhang et al. 2023, Fig. 4D) -- "
                 "published parameter values vs. this PyBNF fit", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = os.path.join(HERE, "tumor_growth_reproduction.png")
    fig.savefig(out, dpi=150)
    print(f"Saved {os.path.basename(out)}")


if __name__ == "__main__":
    main()
