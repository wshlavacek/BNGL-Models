#!/usr/bin/env python
"""
Gate 3a reproduction for reduced_onoff: the reduced model at its Table-1 fitted nominals
(shipped in reduced_onoff.bngl) reproduces the fit target (reduced_onoff_{wt,a20ko}.exp, which
is the original Lipniacki-2004 model's RAW on-off output at the Table-2 times). Self-contained --
uses only committed files (the .bngl + the two .exp). Peak-normalizes the reduced trajectory and
the RAW target per observable, overlays them, and reports the per-observable median relative error.
(Gate 3a is objective-independent -- the switch to the exact Eq-7 objective only changed the .exp
from pre-normalized to raw, so this script now peak-normalizes the target itself for the overlay.)

Run: BNGPATH=$HOME/Simulations/BioNetGen-2.9.3 ~/Code/PyBNF/.venv/bin/python make_reproduction.py
"""
import os, subprocess, tempfile
import numpy as np

BNGPATH = os.environ["BNGPATH"]; BNG2 = os.path.join(BNGPATH, "BNG2.pl")
HERE = os.path.dirname(os.path.abspath(__file__))
RED = os.path.join(HERE, "reduced_onoff.bngl")


def read_exp(path):
    lines = open(path).read().splitlines()
    cols = lines[0].lstrip("#").split()
    a = np.array([[np.nan if v == "NaN" else float(v) for v in ln.split()] for ln in lines[1:]])
    return cols, a


def sim_reduced(t_end=43200, n=4320, c_deg0=False):
    body = open(RED).read().rsplit("end model", 1)[0] + "end model\n"
    ov = 'setParameter("c_deg",0)\n' if c_deg0 else ""
    act = ("begin actions\ngenerate_network({overwrite=>1})\n" + ov +
           f"simulate({{method=>\"ode\",t_start=>0,t_end=>{t_end},n_steps=>{n},print_functions=>1,suffix=>\"o\"}})\nend actions\n")
    with tempfile.TemporaryDirectory() as wd:
        p = os.path.join(wd, "r.bngl"); open(p, "w").write(body + act)
        subprocess.run(["perl", BNG2, p], cwd=wd, capture_output=True, text=True, check=True)
        g = os.path.join(wd, "r_o.gdat")
        hdr = open(g).readline().lstrip("#").split()
        d = np.loadtxt(g, comments="#")
    return {n: d[:, i] for i, n in enumerate(hdr)}


def peaknorm(sim, obs, tgrid):
    y = np.interp(tgrid, sim["time"], sim[obs])
    return y / np.max(y)


def main():
    exps = {"WT": ("reduced_onoff_wt.exp", False), "A20KO": ("reduced_onoff_a20ko.exp", True)}
    sims = {k: sim_reduced(c_deg0=v[1]) for k, v in exps.items()}
    OBS = ["IKK_a", "NFkB_n", "A20", "IkBa_star", "tIkBa"]

    print("=== reduced_onoff Gate 3a: reduced@Table1 vs RAW original target (per-obs median |rel err|, t>0) ===")
    for label, (expf, _) in exps.items():
        cols, a = read_exp(os.path.join(HERE, expf))
        sim = sims[label]
        tsec = a[:, 0]
        for obs in cols[1:]:
            j = cols.index(obs)
            full = ~np.isnan(a[:, j])
            # The .exp now ships RAW original-model output; peak-normalize target AND sim over the
            # observable's measured points, then take the median relative error over t>0.
            tgt = a[full, j] / np.max(a[full, j])
            ys = peaknorm(sim, obs, tsec[full])
            m = tsec[full] > 0
            rel = np.median(np.abs(ys[m] - tgt[m]) / np.clip(np.abs(tgt[m]), 1e-2, None))
            print(f"  {label:6s} {obs:10s} {100*rel:6.1f}%")

    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 5, figsize=(19, 7))
        tgrid = np.linspace(0, 43200, 500)
        for r, (label, (expf, _)) in enumerate(exps.items()):
            cols, a = read_exp(os.path.join(HERE, expf)); sim = sims[label]
            for c, obs in enumerate(OBS):
                ax = axes[r, c]
                if obs not in cols:
                    ax.set_visible(False); continue
                ax.plot(tgrid/60, peaknorm(sim, obs, tgrid), "r-", lw=1.5, label="reduced@Table1")
                j = cols.index(obs); mask = ~np.isnan(a[:, j])
                ax.plot(a[mask, 0]/60, a[mask, j]/np.max(a[mask, j]), "ko", ms=4, label="target (original)")
                ax.set_title(f"{label}: {obs}"); ax.set_xlabel("min")
                if r == 0 and c == 0: ax.legend(fontsize=8)
        fig.suptitle("reduced_onoff Gate 3a: reduced model at Table-1 params reproduces the original-model on-off target")
        fig.tight_layout()
        out = os.path.join(HERE, "reduced_onoff_reproduction.png")
        fig.savefig(out, dpi=90); print("wrote", out)
    except Exception as e:
        print("plot skipped:", e)


if __name__ == "__main__":
    main()
