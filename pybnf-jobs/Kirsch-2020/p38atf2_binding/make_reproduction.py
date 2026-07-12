#!/usr/bin/env python
"""Reproduction figure for the p38atf2_binding job (PMC7666158).

Runs the model at (a) the authors' published parameters (Supplementary Table 2) and
(b) the PyBNF best-fit parameters, under the equilibrate->stimulate protocol, normalizes
the p38:ATF2 complex observable ``p38ATF2all`` to its basal (t=0) value (= treated/
untreated), and overlays both on the real Fig. 4a WT p38-ATF2 NanoBit binding data
(p38atf2_binding.exp, mean +/- SD). Prints chi-square and relative-error metrics.
Requires BNGPATH set (uses BNG2.pl) and matplotlib.

Usage:  BNGPATH=... python make_reproduction.py
"""
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
MODEL = os.path.join(HERE, "p38atf2_binding.bngl")
OBS_COL = 2   # st.gdat columns: time, pT69pT71, p38ATF2all, ppp38, ppJNK

PUBLISHED = dict(keq7=4.695e-6, kstim7=9.074e-5, keq6=1.74e-5, kstim6=1.16e-4,
                 dp1=6.26e-4, dp2=1.76e-3, dp3=9.54e-3, dp4=4.50e-3)
# PyBNF best fit of the 4 free p38-arm params {keq6, kstim6, dp2, dp4}; from
# output/Results/sorted_params_final.txt (de, obj=16.73, converged in ~19 s). The JNK-arm
# and dp3 stay at their published values (held fixed). All 4 recovered within ~1.3x of
# published (keq6 1.14x, kstim6 1.33x, dp2 1.30x, dp4 1.00x).
BESTFIT = dict(PUBLISHED)
BESTFIT.update(dict(keq6=1.9814e-5, kstim6=1.5449e-4, dp2=2.2902e-3, dp4=4.5215e-3))


def run(params):
    """Run equilibrate(Stim=0)->stimulate(Stim=1); return (t[s], p38ATF2all normalized)."""
    with open(MODEL) as fh:
        src = fh.read()
    src = src.split("end model")[0] + "end model\n"
    for k, v in params.items():
        src = re.sub(rf"(^\s*{k}\s+)[\d.eE+-]+", rf"\g<1>{v}", src, count=1, flags=re.M)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "m.bngl")
        with open(path, "w") as fh:
            fh.write(src)
            fh.write('\ngenerate_network({overwrite=>1})\n')
            fh.write('setParameter("Stim",0)\n')
            fh.write('simulate({prefix=>"eq",method=>"ode",t_start=>-20000,t_end=>0,n_steps=>10})\n')
            fh.write('setParameter("Stim",1)\n')
            fh.write('simulate({prefix=>"st",method=>"ode",t_start=>0,t_end=>2466,n_steps=>2466})\n')
        subprocess.run(["perl", BNG, path], cwd=d, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        g = np.loadtxt(os.path.join(d, "st.gdat"))
    t, p = g[:, 0], g[:, OBS_COL]
    return t, p / p[0]


def main():
    exp = np.loadtxt(os.path.join(HERE, "p38atf2_binding.exp"))
    te, ye, se = exp[:, 0], exp[:, 1], exp[:, 2]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.errorbar(te / 60, ye, yerr=se, fmt="o", ms=5, color="#222", capsize=3, zorder=5,
                label="Fig. 4a WT p38-ATF2 NanoBit (mean +/- SD, n=3)")
    for params, color, lab in [(PUBLISHED, "#1f77b4", "published params (Suppl. Table 2)"),
                               (BESTFIT, "#d62728", "PyBNF best fit (keq6,kstim6,dp2,dp4)")]:
        t, y = run(params)
        ax.plot(t / 60, y, "-", color=color, lw=2, label=lab)
        ym = np.interp(te, t, y)
        # PyBNF's chi_sq objective is 0.5 * sum(((sim-exp)/sigma)^2) (Gaussian -logL up to a constant).
        chisq = 0.5 * float(np.sum(((ym - ye) / se) ** 2))
        relerr = np.abs(ym - ye) / ye
        print(f"{lab:44s}  chi_sq(PyBNF)={chisq:8.2f}  median|relerr|={np.median(relerr):.3f}  "
              f"max|relerr|={relerr.max():.3f}")
    ax.set_xlabel("time after anisomycin (min)")
    ax.set_ylabel("p38-ATF2 binding, treated/untreated\n(model: p38ATF2all normalized to basal)")
    ax.set_title("Kirsch 2020 (PMC7666158) Fig. 4a WT -- p38-ATF2 TAD binding (NanoBit)")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out = os.path.join(HERE, "p38atf2_binding_reproduction.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
