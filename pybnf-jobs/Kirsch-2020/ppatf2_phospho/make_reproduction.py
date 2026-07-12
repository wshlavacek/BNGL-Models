#!/usr/bin/env python
"""Reproduction figure for the ppatf2_phospho job (PMC7666158).

Runs the model at (a) the authors' published parameters (Supplementary Table 2) and
(b) the PyBNF best-fit dp3, under the equilibrate->stimulate protocol, in TWO conditions
-- CTR (nominal) and JNK-IN-8 (k1 = k2 = 0) -- and overlays the ABSOLUTE pT69pT71 (no
normalization) on the digitized Fig. 7b pp-ATF2(pT69/pT71) data. Prints sos and
relative-error metrics. Requires BNGPATH set (uses BNG2.pl) and matplotlib.

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
MODEL = os.path.join(HERE, "ppatf2_phospho.bngl")
OBS_COL = 1   # st.gdat columns: time, pT69pT71, p38ATF2all, ppp38, ppJNK

PUBLISHED = dict(keq7=4.695e-6, kstim7=9.074e-5, keq6=1.74e-5, kstim6=1.16e-4,
                 dp1=6.26e-4, dp2=1.76e-3, dp3=9.54e-3, dp4=4.50e-3)
# PyBNF best fit of the single free param dp3 (de, sos=0.0121, ~8 s); from
# output/Results/sorted_params_final.txt. Recovered 1.008e-2 vs published 9.54e-3 = 1.06x.
DP3_BEST = 1.0075e-2


def run(params, jnk_inhibitor=False):
    """equilibrate(Stim=0)->stimulate(Stim=1); return (t[s], pT69pT71 absolute uM).
    jnk_inhibitor=True sets k1=k2=0 (JNK-IN-8) in both phases."""
    p = dict(params)
    if jnk_inhibitor:
        p["k1"] = 0.0
        p["k2"] = 0.0
    with open(MODEL) as fh:
        src = fh.read()
    src = src.split("end model")[0] + "end model\n"
    for k, v in p.items():
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
    return g[:, 0], g[:, OBS_COL]


def main():
    ctr = np.loadtxt(os.path.join(HERE, "ppatf2_phospho_ctr.exp"))
    jnki = np.loadtxt(os.path.join(HERE, "ppatf2_phospho_jnki.exp"))

    best = dict(PUBLISHED); best["dp3"] = DP3_BEST
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for data, jnk, mcol, dcol, dlab in [
            (ctr, False, "#2ca02c", "#1f3fd6", "CTR (digitized Fig. 7b)"),
            (jnki, True, "#ff7f0e", "#d62728", "JNK-IN-8 (digitized Fig. 7b)")]:
        te, ye = data[:, 0], data[:, 1]
        ax.plot(te / 60, ye, "o", ms=7, color=dcol, zorder=5, label=dlab)
        for params, ls, plab in [(PUBLISHED, "-", "published"), (best, "--", "best-fit dp3")]:
            t, y = run(params, jnk_inhibitor=jnk)
            ax.plot(t / 60, y, ls, color=mcol, lw=2,
                    label=f"{'CTR' if not jnk else 'JNK-IN-8'} calc ({plab})")
            ym = np.interp(te, t, y)
            sos = float(np.sum((ym - ye) ** 2))
            rel = np.abs(ym - ye) / ye
            print(f"{'CTR' if not jnk else 'JNKi':4s} {plab:11s}  sos={sos:.4f}  "
                  f"median|relerr|={np.median(rel):.3f}  max|relerr|={rel.max():.3f}")
    ax.set_xlabel("time after anisomycin (min)")
    ax.set_ylabel("pp-ATF2(pT69/pT71) (uM)  [model: pT69pT71, absolute]")
    ax.set_title("Kirsch 2020 (PMC7666158) Fig. 7b -- pp-ATF2(T69/T71) phosphorylation")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out = os.path.join(HERE, "ppatf2_phospho_reproduction.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
