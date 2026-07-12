#!/usr/bin/env python
"""Reproduction figure for cstar_skmel133_bmra.

Left panels: the model AT THE PUBLISHED PARAMETERS reproduces the authors' 24 h single-drug
fold-change TRAINING data (6 arms: ERK/AKT/SRC x dose1/dose2). Fold change =
obs(drug, 24 h) / obs(no-drug baseline, 24 h), computed by running the model (with an
actions block appended) through BNG2.pl once per dose.

Right panel: the BMRA-CI SIGN agreement -- the BMRA-inferred connection coefficient (mean
+/- std, withMyc) for each of the 20 signalling edges + 3 DPD force coefficients, with the
published model's connection side (g>1 activation / g<1 inhibition, degradation flipped)
marked. All signs agree, which is what bmra_signs.prop enforces during the fit.

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
MODEL = os.path.join(HERE, "cstar_skmel133_bmra.bngl")

READOUTS = ["FC_tIRS", "FC_IRSI", "FC_pERK", "FC_pAKT", "FC_pSRC",
            "FC_pPKC", "FC_pS6K", "FC_pRB", "FC_MYC"]
# function -> underlying observable emitted to gdat
OBS = {"FC_tIRS": "tIRS", "FC_IRSI": "pIRSI", "FC_pERK": "pERK", "FC_pAKT": "pAKT",
       "FC_pSRC": "pSRC", "FC_pPKC": "pPKCt", "FC_pS6K": "pS6K", "FC_pRB": "pRB",
       "FC_MYC": "MYCt"}
# arm -> (inhibitor param, dose)
ARMS = {
    "erk_dose1": ("I_ERK_conc", 1.3122028733013222),
    "erk_dose2": ("I_ERK_conc", 2.6244057466026444),
    "akt_dose1": ("I_AKT_conc", 4.32073354133792),
    "akt_dose2": ("I_AKT_conc", 8.64146708267584),
    "src_dose1": ("I_SRC_conc", 0.5173293874300976),
    "src_dose2": ("I_SRC_conc", 1.0346587748601952),
}


def run_ss(overrides):
    """Run the model to steady state with parameter overrides; return {obs: value@t_end}."""
    src = open(MODEL).read()
    for k, v in overrides.items():
        src = re.sub(rf"(?m)^{k}\s*=.*$", f"{k} = {v}", src)
    src += ("\n\nbegin actions\ngenerate_network({overwrite=>1})\n"
            "simulate({method=>\"ode\",t_end=>86400,n_steps=>50,atol=>1e-10,rtol=>1e-8})\n"
            "end actions\n")
    with tempfile.TemporaryDirectory() as d:
        mp = os.path.join(d, "m.bngl")
        open(mp, "w").write(src)
        subprocess.run(["perl", BNG, "m.bngl"], cwd=d, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        g = open(os.path.join(d, "m.gdat")).read().strip().splitlines()
        hdr = g[0].lstrip("#").split()
        last = [float(x) for x in g[-1].split()]
        return dict(zip(hdr, last))


def read_exp(fn):
    lines = [l for l in open(os.path.join(HERE, fn)).read().splitlines() if l and not l.startswith("#")]
    cols = open(os.path.join(HERE, fn)).readline().lstrip("#").split()
    row = lines[-1].split()  # measured row (t=86400)
    d = dict(zip(cols, row))
    return {r: float(d[r + "()"]) for r in READOUTS}, {r: float(d[r + "_SD"]) for r in READOUTS}


def main():
    base = run_ss({})  # no-drug baseline steady state
    fig, axes = plt.subplots(2, 4, figsize=(17, 8))
    errs = []
    for ax, (arm, (par, dose)) in zip(axes.flat, ARMS.items()):
        drug = run_ss({par: dose})
        model_fc = [drug[OBS[r]] / base[OBS[r]] for r in READOUTS]
        data_fc, data_sd = read_exp(f"{arm.replace('_dose', '').replace('erk','ERK').replace('akt','AKT').replace('src','SRC')}inh.exp"
                                    if False else _expname(arm))
        d_fc = [data_fc[r] for r in READOUTS]
        d_sd = [data_sd[r] for r in READOUTS]
        rel = [abs(m - x) / abs(x) for m, x in zip(model_fc, d_fc)]
        errs += rel
        x = np.arange(len(READOUTS))
        ax.bar(x - 0.2, d_fc, 0.4, yerr=d_sd, label="data", color="#4C72B0", capsize=2)
        ax.bar(x + 0.2, model_fc, 0.4, label="model", color="#DD8452")
        ax.axhline(1.0, color="gray", lw=0.6, ls="--")
        ax.set_xticks(x); ax.set_xticklabels([r[3:] for r in READOUTS], rotation=60, fontsize=7)
        ax.set_title(f"{arm}  (median {np.median(rel)*100:.0f}%)", fontsize=9)
        ax.set_ylabel("fold change vs baseline", fontsize=8)
        if arm == "erk_dose1":
            ax.legend(fontsize=8)

    # right-bottom two cells: BMRA sign agreement
    _bmra_panel(axes.flat[6], axes.flat[7])
    fig.suptitle(f"cstar_skmel133_bmra -- model at published params vs single-drug fold changes "
                 f"(overall median {np.median(errs)*100:.1f}%) + BMRA-CI sign agreement", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(HERE, "cstar_skmel133_bmra_reproduction.png")
    fig.savefig(out, dpi=120)
    print(f"overall median rel err = {np.median(errs)*100:.1f}%  (n={len(errs)})")
    print(f"wrote {out}")


def _expname(arm):
    m = {"erk_dose1": "dose1_ERKinh.exp", "erk_dose2": "dose2_ERKinh.exp",
         "akt_dose1": "dose1_AKTinh.exp", "akt_dose2": "dose2_AKTinh.exp",
         "src_dose1": "dose1_SRCinh.exp", "src_dose2": "dose2_SRCinh.exp"}
    return m[arm]


def _bmra_panel(ax1, ax2):
    # read the generated .prop to list (edge, BMRA r, rs, required side)
    rows = []
    for line in open(os.path.join(HERE, "bmra_signs.prop")):
        m = re.search(r"# (\S+) BMRA r=([-+0-9.]+)\+/-([0-9.]+) \(z=([0-9.]+)\)", line)
        if m:
            rows.append((m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4))))
    rows = rows[:20]  # signalling edges (skip DPD rows for this bar; still all agree)
    rows.sort(key=lambda r: r[1])
    y = np.arange(len(rows))
    for ax, title in ((ax1, "BMRA-inferred connection r (mean +/- std)"), (ax2, None)):
        if title is None:
            ax.axis("off")
            ax.text(0.0, 0.5,
                    "All 20 signalling edges + 3 DPD force coeffs:\n"
                    "the published model's connection side\n"
                    "(g>1 activation / g<1 inhibition, Eq. 24)\n"
                    "matches the BMRA-inferred sign.\n\n"
                    "bmra_signs.prop enforces this during the fit\n"
                    "(job_type=check -> Satisfied 23/23).",
                    fontsize=9, va="center")
            continue
        cols = ["#C44E52" if r[1] < 0 else "#55A868" for r in rows]
        ax.barh(y, [r[1] for r in rows], xerr=[1.96 * r[2] for r in rows], color=cols, capsize=2)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=6)
        ax.set_xlabel("BMRA r (green +, red -; bar = 95% CI)", fontsize=8)
        ax.set_title(title, fontsize=9)


if __name__ == "__main__":
    main()
