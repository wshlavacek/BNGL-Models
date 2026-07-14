#!/usr/bin/env python
"""Reproduction figure for cstar_skmel133_bmra_exact.

TOP (fold-change reproduction): the model AT THE AUTHORS' PUBLISHED PARAMETERS reproduces the
24 h single-drug fold-change TRAINING data (11 arms: ERK/AKT/SRC/PKC/mTOR/CDK single-drug panels).
Each fold change is computed the way the fit protocol does it -- a drug-free pre-equilibration to
steady state (the baseline reference), then the inhibitor applied and the network run 24 h --
through BNG2.pl once per arm, so the reproduction is INDEPENDENT of PyBNF. FC = obs(drug, 24 h) /
obs(baseline steady state). ERK/AKT/SRC act via the rate-law factor 1/(1+I_x_conc); PKC/mTOR/CDK
by competitive binding (their inhibitor seed I_x_conc feeds a binding rule).

BOTTOM (the exact BMRA-CI test -- the point of this slug): the model's EXACT Eq.14 connection
coefficient r_ij at the published drug-free steady state, plotted against each BMRA confidence
interval [mean-std, mean+std] (Table S10, with_MYC). 16/20 signaling edges lie in-band; 4 --
AKT<-ERK, AKT<-PKC, AKT<-IRS, MYC<-CDK -- sit just outside (marked red), which is what
bmra_rij.prop makes binding during the constrained fit (job_type=check -> Satisfied 39/43).

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
MODEL = os.path.join(HERE, "cstar_skmel133_bmra_exact.bngl")

READOUTS = ["FC_tIRS", "FC_IRSI", "FC_pERK", "FC_pAKT", "FC_pSRC",
            "FC_pPKC", "FC_pS6K", "FC_pRB", "FC_MYC"]
# fold-change function -> underlying observable emitted to gdat
OBS = {"FC_tIRS": "tIRS", "FC_IRSI": "pIRSI", "FC_pERK": "pERK", "FC_pAKT": "pAKT",
       "FC_pSRC": "pSRC", "FC_pPKC": "pPKCt", "FC_pS6K": "pS6K", "FC_pRB": "pRB",
       "FC_MYC": "MYCt"}

# arm -> (kind, target, dose, exp file). kind 'param' -> setParameter (rate-law inhibitor);
# 'species' -> setConcentration (competitive-binding inhibitor seed). Doses == the conf conditions.
ARMS = [
    ("erk_dose1",  "param",   "I_ERK_conc",  1.3122028733013222, "dose1_ERKinh.exp"),
    ("erk_dose2",  "param",   "I_ERK_conc",  2.6244057466026445, "dose2_ERKinh.exp"),
    ("akt_dose1",  "param",   "I_AKT_conc",  4.32073354133792,   "dose1_AKTinh.exp"),
    ("akt_dose2",  "param",   "I_AKT_conc",  8.64146708267584,   "dose2_AKTinh.exp"),
    ("src_dose1",  "param",   "I_SRC_conc",  0.5173293874300976, "dose1_SRCinh.exp"),
    ("src_dose2",  "param",   "I_SRC_conc",  1.0346587748601952, "dose2_SRCinh.exp"),
    ("pkc_dose1",  "species", "IPKC(PKCBD)", 1.7250107623770197, "dose1_PKCinh.exp"),
    ("pkc_dose2",  "species", "IPKC(PKCBD)", 3.4500215247540393, "dose2_PKCinh.exp"),
    ("mtor_dose1", "species", "ImTOR(mTORBD)", 1.3510883690081217, "dose1_mTORinh.exp"),
    ("mtor_dose2", "species", "ImTOR(mTORBD)", 2.7021767380162434, "dose2_mTORinh.exp"),
    ("cdk_dose1",  "species", "ICDK(CDKBD)", 3.128238440223784,  "dose1_CDKinh.exp"),
]


def _run(actions):
    """Run the model with an appended actions block through BNG2.pl; return list of gdat dicts
    (one per simulate action, in order)."""
    src = open(MODEL).read()
    src += "\n\nbegin actions\ngenerate_network({overwrite=>1})\n" + actions + "\nend actions\n"
    with tempfile.TemporaryDirectory() as d:
        open(os.path.join(d, "m.bngl"), "w").write(src)
        subprocess.run(["perl", BNG, "m.bngl"], cwd=d, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out = []
        for suff in re.findall(r'suffix=>"(\w+)"', actions):
            g = open(os.path.join(d, f"m_{suff}.gdat")).read().strip().splitlines()
            hdr = g[0].lstrip("#").split()
            out.append(dict(zip(hdr, [float(x) for x in g[-1].split()])))
        return out


def arm_fold_change(kind, target, dose):
    """FC per observable via the fit's two-phase protocol: drug-free pre-equilibration to steady
    state, then the drug applied and run 24 h. FC = obs(24 h) / obs(baseline steady state)."""
    if kind == "param":
        pert = f'setParameter("{target}",{dose})\n'
    else:
        pert = f'setConcentration("{target}",{dose})\n'
    acts = (
        'resetConcentrations()\n'
        'setParameter("I_ERK_conc",0)\nsetParameter("I_AKT_conc",0)\nsetParameter("I_SRC_conc",0)\n'
        'setConcentration("IPKC(PKCBD)",0)\nsetConcentration("ImTOR(mTORBD)",0)\nsetConcentration("ICDK(CDKBD)",0)\n'
        'simulate({method=>"ode",steady_state=>1,t_end=>1000000,n_steps=>1,suffix=>"base"})\n'
        + pert +
        'simulate({method=>"ode",t_start=>0,t_end=>86400,n_steps=>2,suffix=>"drug"})\n'
    )
    base, drug = _run(acts)
    return {r: drug[OBS[r]] / base[OBS[r]] for r in READOUTS}


def read_exp(fn):
    lines = [l for l in open(os.path.join(HERE, fn)).read().splitlines() if l and not l.startswith("#")]
    cols = open(os.path.join(HERE, fn)).readline().lstrip("#").split()
    row = dict(zip(cols, lines[-1].split()))  # measured row (t=86400)
    fc = {r: float(row[r + "()"]) for r in READOUTS}
    sd = {r: (float(row[r + "_SD"]) if row[r + "_SD"] != "NaN" else np.nan) for r in READOUTS}
    return fc, sd


def bmra_bands():
    """From bmra_rij.prop: {edge: (mean, std, lo, hi)} for the 20 signaling r_ij."""
    bands = {}
    text = open(os.path.join(HERE, "bmra_rij.prop")).read().splitlines()
    for line in text:
        m = re.search(r"# (\w+)<-(\w+): BMRA r=([-+0-9.]+)\+/-([0-9.]+)\s+->\s+r in \[([-+0-9.]+),([-+0-9.]+)\]", line)
        if m:
            edge = f"r_{m.group(1)}_{m.group(2)}"
            bands[edge] = (float(m.group(3)), float(m.group(4)), float(m.group(5)), float(m.group(6)))
    return bands


def main():
    # --- fold-change panels ---
    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(4, 4, height_ratios=[1, 1, 1, 1.5])
    errs = []
    for i, (arm, kind, target, dose, expfn) in enumerate(ARMS):
        ax = fig.add_subplot(gs[i // 4, i % 4])
        mfc = arm_fold_change(kind, target, dose)
        dfc, dsd = read_exp(expfn)
        m = [mfc[r] for r in READOUTS]
        d = [dfc[r] for r in READOUTS]
        rel = [abs(a - b) / abs(b) for a, b in zip(m, d)]
        errs += rel
        x = np.arange(len(READOUTS))
        ax.bar(x - 0.2, d, 0.4, yerr=[0 if np.isnan(s) else s for s in dsd.values()],
               label="data", color="#4C72B0", capsize=2)
        ax.bar(x + 0.2, m, 0.4, label="model", color="#DD8452")
        ax.axhline(1.0, color="gray", lw=0.6, ls="--")
        ax.set_xticks(x); ax.set_xticklabels([r[3:] for r in READOUTS], rotation=60, fontsize=6)
        note = "  (no SD; not in chi_sq)" if expfn == "dose1_CDKinh.exp" else ""
        ax.set_title(f"{arm}  (median {np.median(rel)*100:.0f}%){note}", fontsize=8)
        ax.set_ylabel("fold change", fontsize=7)
        if i == 0:
            ax.legend(fontsize=7)

    # --- the exact BMRA-CI panel (spans the bottom row) ---
    axr = fig.add_subplot(gs[3, :])
    ss = _run('resetConcentrations()\nsimulate({method=>"ode",steady_state=>1,t_end=>1000000,'
              'n_steps=>1,suffix=>"ss",print_functions=>1})\n')[0]
    bands = bmra_bands()
    edges = list(bands.keys())
    edges.sort(key=lambda e: bands[e][0])
    y = np.arange(len(edges))
    for yi, e in zip(y, edges):
        mean, std, lo, hi = bands[e]
        axr.plot([lo, hi], [yi, yi], color="#BBBBBB", lw=6, solid_capstyle="butt", zorder=1)
        axr.plot([mean, mean], [yi - 0.3, yi + 0.3], color="#888888", lw=1.5, zorder=2)
    rvals = [ss[e] for e in edges]
    inband = [(bands[e][2] <= ss[e] <= bands[e][3]) for e in edges]
    axr.scatter([r for r, ok in zip(rvals, inband) if ok],
                [yi for yi, ok in zip(y, inband) if ok],
                color="#55A868", s=45, zorder=3, label="model r (in-band)")
    axr.scatter([r for r, ok in zip(rvals, inband) if not ok],
                [yi for yi, ok in zip(y, inband) if not ok],
                color="#C44E52", s=55, zorder=3, marker="D", label="model r (out-of-band)")
    axr.axvline(0, color="k", lw=0.6)
    axr.set_yticks(y); axr.set_yticklabels([e.replace("r_", "").replace("_", "<-") for e in edges], fontsize=7)
    axr.set_xlabel("connection coefficient r  (grey bar = BMRA mean +/- std, Table S10 withMyc)", fontsize=9)
    axr.set_title(f"EXACT Eq.14 r_ij at published params vs BMRA CI  "
                  f"({sum(inband)}/{len(edges)} signaling edges in-band; the fit pulls the rest in)", fontsize=10)
    axr.legend(fontsize=8, loc="lower right")

    fig.suptitle(f"cstar_skmel133_bmra_exact -- model at published params: single-drug fold changes "
                 f"(overall median {np.median(errs)*100:.1f}%, n={len(errs)}) + exact BMRA-CI test", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = os.path.join(HERE, "cstar_skmel133_bmra_exact_reproduction.png")
    fig.savefig(out, dpi=110)
    print(f"overall median rel err = {np.median(errs)*100:.1f}%  (n={len(errs)})")
    print(f"signaling edges in-band at published params: {sum(inband)}/{len(edges)}")
    for e, r, ok in zip(edges, rvals, inband):
        if not ok:
            print(f"  OUT: {e.replace('r_','').replace('_','<-')}  r={r:+.3f}  band=[{bands[e][2]:+.3f},{bands[e][3]:+.3f}]")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
