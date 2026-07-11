#!/usr/bin/env python
"""Reproduction figure for the fig7b_timecourse job (PMC7666158).

Runs the model at (a) the authors' published parameters (Supplementary Table 2) and
(b) the PyBNF best-fit parameters, under the equilibrate->stimulate protocol, normalizes
pT69pT71 to its basal (t=0) value (= treated/untreated), and overlays both on the real
Fig. 4a WT pp-ATF2(T69/T71) data (fig7b_timecourse.exp, mean +/- SD). Prints chi-square
and relative-error metrics. Requires BNGPATH set (uses BNG2.pl) and matplotlib.

Usage:  BNGPATH=... python make_reproduction.py
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
MODEL = os.path.join(HERE, "fig7b_timecourse.bngl")

PUBLISHED = dict(keq7=4.695e-6, kstim7=9.074e-5, keq6=1.74e-5, kstim6=1.16e-4,
                 dp1=6.26e-4, dp2=1.76e-3, dp3=9.54e-3, dp4=4.50e-3)
# PyBNF best fit (fit_out/Results/fig7b_timecourse_bestfit.bngl; chi_sq ~= 64.9)
BESTFIT = dict(keq7=3.0272e-5, kstim7=6.0934e-5, keq6=1.0711e-5, kstim6=6.5245e-4,
               dp1=5.0459e-3, dp2=5.5849e-3, dp3=1.0444e-3, dp4=1.2411e-3)


def run(params):
    """Run equilibrate(Stim=0)->stimulate(Stim=1) and return (t[s], pT69pT71 normalized)."""
    with open(MODEL) as fh:
        src = fh.read()
    src = src.split("end model")[0] + "end model\n"
    for k, v in params.items():
        import re
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
    t, p = g[:, 0], g[:, 1]
    return t, p / p[0]


def main():
    exp = np.loadtxt(os.path.join(HERE, "fig7b_timecourse.exp"))
    te, ye, se = exp[:, 0], exp[:, 1], exp[:, 2]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.errorbar(te / 60, ye, yerr=se, fmt="o", ms=5, color="#222", capsize=3, zorder=5,
                label="Fig. 4a WT data (mean +/- SD, n=3)")
    for params, color, lab in [(PUBLISHED, "#1f77b4", "published params (Suppl. Table 2)"),
                               (BESTFIT, "#d62728", "PyBNF best fit")]:
        t, y = run(params)
        ax.plot(t / 60, y, "-", color=color, lw=2, label=lab)
        ym = np.interp(te, t, y)
        # PyBNF's chi_sq objective is 0.5 * sum(((sim-exp)/sigma)^2) (Gaussian -logL up to a constant).
        chisq = 0.5 * float(np.sum(((ym - ye) / se) ** 2))
        relerr = np.abs(ym - ye) / ye
        print(f"{lab:34s}  chi_sq(PyBNF)={chisq:8.2f}  median|relerr|={np.median(relerr):.3f}  "
              f"max|relerr|={relerr.max():.3f}")
    ax.set_xlabel("time after anisomycin (min)")
    ax.set_ylabel("pp-ATF2(T69/T71), treated/untreated\n(model: pT69pT71 normalized to basal)")
    ax.set_title("Kirsch 2020 (PMC7666158) Fig. 4a WT -- JNK/p38/ATF2 phosphoswitch fit")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out = os.path.join(HERE, "fig7b_timecourse_reproduction.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
