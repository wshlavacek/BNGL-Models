#!/usr/bin/env python
"""
Gate 3a reproduction for reduced_combination: the reduced model at its Table-1 fitted nominals
reproduces the fit target (the original Lipniacki-2004 model's combination-experiment output,
peak-normalized in reduced_combination_<protocol>_{wt,a20ko}.exp). Self-contained -- uses only
committed files (the 12 per-protocol .bngl + the 12 .exp). Reproduces the paper's Fig 2 (WT nuclear
NF-kB across the 6 protocols) and reports the peak-normalized median relative error per protocol.

Run: BNGPATH=$HOME/Simulations/BioNetGen-2.9.3 ~/Code/PyBNF/.venv/bin/python make_reproduction.py
"""
import os, subprocess, tempfile
import numpy as np

BNGPATH = os.environ["BNGPATH"]; BNG2 = os.path.join(BNGPATH, "BNG2.pl")
HERE = os.path.dirname(os.path.abspath(__file__))
PROTOCOLS = ["continuous", "pulse5_60", "pulse5_100", "pulse5_200", "pulse22_5", "pulse45"]
TEND = {"continuous": 14400, "pulse5_60": 43200, "pulse5_100": 43200,
        "pulse5_200": 43200, "pulse22_5": 43200, "pulse45": 43200}
OBS = ["IKK_a", "NFkB_n", "A20", "IkBa_star", "tIkBa"]


def read_exp(path):
    lines = open(path).read().splitlines()
    cols = lines[0].lstrip("#").split()
    a = np.array([[np.nan if v == "NaN" else float(v) for v in ln.split()] for ln in lines[1:]])
    return cols, a


def sim(model_stem, tend):
    body = open(os.path.join(HERE, model_stem + ".bngl")).read().rsplit("end model", 1)[0] + "end model\n"
    ns = max(2000, int(tend / 10))
    act = ("begin actions\ngenerate_network({overwrite=>1})\n"
           f"simulate({{method=>\"ode\",t_start=>0,t_end=>{tend},n_steps=>{ns},print_functions=>1,suffix=>\"o\"}})\nend actions\n")
    with tempfile.TemporaryDirectory() as wd:
        p = os.path.join(wd, "r.bngl"); open(p, "w").write(body + act)
        subprocess.run(["perl", BNG2, p], cwd=wd, capture_output=True, text=True, check=True)
        g = os.path.join(wd, "r_o.gdat")
        hdr = open(g).readline().lstrip("#").split()
        d = np.loadtxt(g, comments="#")
    return {n: d[:, i] for i, n in enumerate(hdr)}


def peaknorm(s, obs, tg):
    y = np.interp(tg, s["time"], s[obs]); return y / np.max(y)


def main():
    print("=== reduced_combination Gate 3a: reduced@Table1 vs target (per-protocol median |rel err|) ===")
    sims = {}
    for name in PROTOCOLS:
        for g in ("wt", "a20ko"):
            sims[(name, g)] = sim(f"reduced_combination_{name}_{g}", TEND[name])
    for name in PROTOCOLS:
        for g, obs_all in (("wt", OBS), ("a20ko", [o for o in OBS if o != "A20"])):
            cols, a = read_exp(os.path.join(HERE, f"reduced_combination_{name}_{g}.exp"))
            s = sims[(name, g)]; errs = []
            for obs in cols[1:]:
                j = cols.index(obs); mask = ~np.isnan(a[:, j]) & (a[:, 0] > 0)
                yn = peaknorm(s, obs, a[mask, 0]); e0 = a[mask, j] - 0.03
                errs.append(np.median(np.abs(yn - e0) / np.clip(np.abs(e0), 0.03, None)))
            print(f"  {name:11s} {g:5s}  median |rel err| = {100*np.median(errs):5.1f}%")

    # Fig 2 reproduction: WT nuclear NF-kB across the 6 protocols
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 3, figsize=(16, 8))
        for ax, name in zip(axes.flat, PROTOCOLS):
            s = sims[(name, "wt")]; tg = np.linspace(0, TEND[name], 600)
            ax.plot(tg / 60, peaknorm(s, "NFkB_n", tg), "r-", lw=1.5, label="reduced@Table1")
            cols, a = read_exp(os.path.join(HERE, f"reduced_combination_{name}_wt.exp"))
            j = cols.index("NFkB_n"); mask = ~np.isnan(a[:, j])
            ax.plot(a[mask, 0] / 60, a[mask, j] - 0.03, "ko", ms=3, label="target (original)")
            ax.set_title(f"WT NFkB_n: {name}"); ax.set_xlabel("min")
        axes.flat[0].legend(fontsize=8)
        fig.suptitle("reduced_combination Gate 3a (Fig 2): WT nuclear NF-kB, reduced@Table1 vs original-model target")
        fig.tight_layout()
        out = os.path.join(HERE, "reduced_combination_reproduction.png"); fig.savefig(out, dpi=90)
        print("wrote", out)
    except Exception as e:
        print("plot skipped:", e)


if __name__ == "__main__":
    main()
