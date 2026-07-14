#!/usr/bin/env python
"""Reproduction figure for cstar_trkab_bmra_exact.

TOP (fold-change reproduction): each receptor model AT THE AUTHORS' PUBLISHED PARAMETERS reproduces
its ligand-stimulation phospho time course (7 readouts at 0/10/45 min) -- pre-equilibrate WITHOUT
ligand, add ligand (Lig_on 0->1), FC(t) = obs(t)/obs(0). Run through BNG2.pl (independent of PyBNF).
The 10-min points are the training set; the 45-min points are the paper's validation set.

BOTTOM (the exact BMRA-CI test): each model's EXACT Eq.14 connection coefficient r_ij at the BASAL
(no-ligand) steady state, vs its 10-min BMRA confidence interval (Table S5). Green = in-band,
red = out; the fit pulls the out-of-band edges toward their CIs. The RTK-incoming (receptor-
dimerization) edges are documented, not constrained (no closed-form Eq.14 C_i; see VALIDATION.md).

Requires BNGPATH set (BNG2.pl). Run from this folder:  python make_reproduction.py
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
MODELS = {"A": ("cstar_trka_bmra_exact.bngl", "cstar_trka.exp", "TrkA / NGF"),
          "B": ("cstar_trkb_bmra_exact.bngl", "cstar_trkb.exp", "TrkB / BDNF")}
READOUTS = ["FC_pTRK", "FC_pAdRTK", "FC_ppERK", "FC_pAKT", "FC_pJNK", "FC_pRSK", "FC_pS6K"]
OBS = {"FC_pTRK": "pTRK", "FC_pAdRTK": "pAdRTK", "FC_ppERK": "ppERK", "FC_pAKT": "pAKT",
       "FC_pJNK": "pJNKt", "FC_pRSK": "pRSK", "FC_pS6K": "pS6Kt"}


def _run(model, actions):
    src = open(os.path.join(HERE, model)).read()
    src = re.sub(r'(?m)^(Lig_on\s*=).*$', r'\1 0', src)
    src += "\n\nbegin actions\ngenerate_network({overwrite=>1})\n" + actions + "\nend actions\n"
    with tempfile.TemporaryDirectory() as d:
        open(os.path.join(d, "m.bngl"), "w").write(src)
        subprocess.run(["perl", BNG, "m.bngl"], cwd=d, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out = {}
        for suff in re.findall(r'suffix=>"(\w+)"', actions):
            g = open(os.path.join(d, f"m_{suff}.gdat")).read().strip().splitlines()
            hdr = g[0].lstrip("#").split()
            out[suff] = [dict(zip(hdr, [float(x) for x in ln.split()])) for ln in g[1:]]
        return out


def timecourse(model):
    """FC(t)=obs(t)/obs(0) at 0/600/2700 s after ligand, from the basal-equilibrated state."""
    acts = ('setParameter("Lig_on",0)\nsimulate({method=>"ode",steady_state=>1,t_end=>1e6,n_steps=>1,suffix=>"pre"})\n'
            'setParameter("Lig_on",1)\nsimulate({method=>"ode",t_start=>0,sample_times=>[0,600,2700],suffix=>"stim",print_functions=>1})\n')
    rows = _run(model, acts)["stim"]
    base = {r: rows[0][OBS[r]] for r in READOUTS}
    return [ {r: rows[k][OBS[r]] / base[r] for r in READOUTS} for k in range(3) ]  # t = 0,600,2700


def basal_r(model):
    acts = 'simulate({method=>"ode",steady_state=>1,t_end=>1e6,n_steps=>1,suffix=>"ss",print_functions=>1})\n'
    v = _run(model, acts)["ss"][-1]
    return {k: val for k, val in v.items() if k.startswith("r_")}


def read_exp(fn):
    lines = [l.split() for l in open(os.path.join(HERE, fn)).read().splitlines() if l and not l.startswith("#")]
    cols = open(os.path.join(HERE, fn)).readline().lstrip("#").split()
    out = []
    for ln in lines:
        d = dict(zip(cols, ln))
        out.append({r: float(d[r + "()"]) for r in READOUTS})
    return out  # rows at t = 0,600,2700


def bands(tag):
    b = {}
    for line in open(os.path.join(HERE, f"bmra_{tag}_rij.prop")):
        m = re.search(r"# (\w+)<-(\w+): BMRA r=([-+0-9.]+)\+/-([0-9.]+)\s+->\s+r in \[([-+0-9.]+),([-+0-9.]+)\]", line)
        if m:
            b[f"r_{m.group(1)}_{m.group(2)}"] = (float(m.group(3)), float(m.group(4)), float(m.group(5)), float(m.group(6)))
    return b


def main():
    T = [0, 10, 45]  # minutes
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.25])
    allerr = []
    for col, tag in enumerate(("A", "B")):
        model, expfn, title = MODELS[tag]
        mtc = timecourse(model)
        dtc = read_exp(expfn)
        ax = fig.add_subplot(gs[0, col])
        cmap = plt.cm.tab10
        errs = []
        for i, r in enumerate(READOUTS):
            m = [mtc[k][r] for k in range(3)]
            d = [dtc[k][r] for k in range(3)]
            errs += [abs(a - b) / abs(b) for a, b in zip(m[1:], d[1:])]  # skip t=0 (==1)
            c = cmap(i)
            ax.plot(T, m, "-", color=c, lw=1.5, label=r[3:])
            ax.plot(T, d, "o", color=c, ms=6, mec="k", mew=0.4)
        allerr += errs
        ax.set_title(f"{title}  (10+45 min, median {np.median(errs)*100:.0f}%)  line=model  o=data", fontsize=10)
        ax.set_xlabel("min after ligand", fontsize=9)
        ax.set_ylabel("fold change vs t=0", fontsize=9)
        ax.set_xticks(T)
        ax.legend(fontsize=7, ncol=2)

        # bottom: exact r_ij vs 10-min BMRA CI
        axr = fig.add_subplot(gs[1, col])
        rv = basal_r(model)
        bd = bands(tag)
        edges = sorted(bd, key=lambda e: bd[e][0])
        y = np.arange(len(edges))
        for yi, e in zip(y, edges):
            mean, std, lo, hi = bd[e]
            axr.plot([lo, hi], [yi, yi], color="#BBBBBB", lw=6, solid_capstyle="butt", zorder=1)
            axr.plot([mean, mean], [yi - 0.3, yi + 0.3], color="#888", lw=1.3, zorder=2)
        inb = [bd[e][2] <= rv[e] <= bd[e][3] for e in edges]
        axr.scatter([rv[e] for e, ok in zip(edges, inb) if ok], [yi for yi, ok in zip(y, inb) if ok],
                    color="#55A868", s=40, zorder=3, label="in-band")
        axr.scatter([rv[e] for e, ok in zip(edges, inb) if not ok], [yi for yi, ok in zip(y, inb) if not ok],
                    color="#C44E52", s=50, marker="D", zorder=3, label="out-of-band")
        axr.axvline(0, color="k", lw=0.6)
        axr.set_yticks(y)
        axr.set_yticklabels([e.replace("r_", "").replace("_", "<-") for e in edges], fontsize=7)
        axr.set_xlabel("connection coefficient r  (grey = BMRA 10-min mean +/- std, Table S5)", fontsize=9)
        axr.set_title(f"Trk{tag}: exact Eq.14 r_ij at basal SS vs 10-min BMRA CI  ({sum(inb)}/{len(edges)} in-band)", fontsize=9)
        axr.legend(fontsize=8, loc="lower right")

    fig.suptitle(f"cstar_trkab_bmra_exact -- model at published params: ligand phospho time courses "
                 f"(overall median {np.median(allerr)*100:.1f}%) + exact BMRA-CI test", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = os.path.join(HERE, "cstar_trkab_bmra_exact_reproduction.png")
    fig.savefig(out, dpi=115)
    print(f"overall median rel err = {np.median(allerr)*100:.1f}%  (n={len(allerr)})")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
