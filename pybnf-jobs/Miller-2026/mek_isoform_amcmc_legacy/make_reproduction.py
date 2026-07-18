#!/usr/bin/env python
"""Gate-3a reproduction for the Miller-2026 MEK-isoform aMCMC / Bayesian-UQ job.

The aMCMC holds 25 rate constants at Kocieniewski & Lipniacki's original defaults (baked into the
models) and samples d3, u3, sigma + 3 scaling factors. This script evaluates the model at the
POSTERIOR MEDIANS of the two feedback rate constants (Fig 7: d3 ~ 1.0e-3, u3 ~ 0.4e-3) with the
output scaling factors near their posterior center, runs BNG2.pl (c1 = 0.02), and:
  * overlays the scaled WT outputs vs WT.exp (reproduces Fig 6 median trajectories);
  * plots MEK_pRDS(t) for the five cell lines (reproduces Fig 5 median trajectories);
  * reports BPSL constraint satisfaction at the posterior median.

Run (from this folder, with BNGPATH set):  python make_reproduction.py
Writes: mek_isoform_amcmc_reproduction.png
"""
import os, re, glob, subprocess, tempfile
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BNG = os.path.join(os.environ["BNGPATH"], "BNG2.pl")
MODELS = ["WT", "KO", "N78G", "T292A", "T292D"]

# Posterior medians (Fig 7 marginals) + output-scale centers (~ starting_params / posterior).
POST = {"d3":1.0e-3, "u3":0.4e-3, "sigma":1.1,
        "scalepEGFR":5.27e-5, "scalepERK":3.89e-6, "scalepSos1":2.0e-4}

def run(m, wd):
    txt = open(os.path.join(HERE, f"MEK1_{m}.bngl")).read()
    for k, v in POST.items():
        txt = re.sub(rf'\b{re.escape(k)}__FREE\b', repr(v), txt)
    open(os.path.join(wd, f"{m}.bngl"), "w").write(txt)
    subprocess.run(["perl", BNG, f"{m}.bngl"], cwd=wd, capture_output=True, text=True)
    g = glob.glob(os.path.join(wd, f"*_{m}.gdat"))[0]
    L = [l for l in open(g) if l.strip()]; h = L[0].lstrip("#").split(); d = {x: [] for x in h}
    for l in L[1:]:
        for x, val in zip(h, l.split()): d[x].append(float(val))
    return d

tmp = tempfile.mkdtemp()
G = {}
for m in MODELS:
    wd = os.path.join(tmp, m); os.makedirs(wd); G[m] = run(m, wd)
def at(d, c, t):
    i = min(range(len(d["time"])), key=lambda j: abs(d["time"][j]-t)); return d[c][i]

EXP = {}
for ln in open(os.path.join(HERE, "WT.exp")):
    if ln.strip() and not ln.strip().startswith("#"):
        p = ln.split(); EXP[float(p[0])] = {"pSOS1": p[2], "pEGFR": p[3], "pERK": p[4]}
et = sorted(EXP)
FCOL = {"pSOS1":"scaled_pSOS1","pEGFR":"scaled_pEGFR","pERK":"scaled_pERK"}
F_quant = sum((at(G["WT"], FCOL[sp], t) - float(EXP[t][sp]))**2
              for t in et for sp in ("pSOS1","pEGFR","pERK") if EXP[t][sp] != "nan")

def side(s):
    m = re.match(r'([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)\s+at\s+time=(\d+)', s.strip())
    return ("obs", m.group(1), m.group(2), int(m.group(3))) if m else ("num", float(s), None, None)
nsat = ntot = 0
for m in MODELS:
    for raw in open(os.path.join(HERE, f"{m}.prop")):
        raw = raw.strip()
        if not raw or raw.startswith("#"): continue
        line = re.split(r'confidence|weight', raw)[0].strip()
        g = re.match(r'(.+?)([<>]=?)(.+)', line); L, op, R = g.groups()
        lt, rt = side(L), side(R)
        lv = at(G[lt[1]], lt[2], lt[3]) if lt[0]=="obs" else lt[1]
        rv = at(G[rt[1]], rt[2], rt[3]) if rt[0]=="obs" else rt[1]
        ok = (lv > rv) if op.startswith(">") else (lv < rv)
        ntot += 1; nsat += ok
print(f"F_quant (18 WT points) = {F_quant:.2f}   [K&L-baseline model at posterior-median d3,u3]")
print(f"BPSL constraints satisfied at posterior median: {nsat} / {ntot}")

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.6))
colors = {"pEGFR":"#E69F00","pSOS1":"#56B4E9","pERK":"#009E73"}; mk = {"pEGFR":"o","pSOS1":"s","pERK":"^"}
for sp in ("pEGFR","pSOS1","pERK"):
    axA.plot(G["WT"]["time"], G["WT"][FCOL[sp]], "-", color=colors[sp], lw=2, label=sp)
    axA.plot([t for t in et if EXP[t][sp]!="nan"], [float(EXP[t][sp]) for t in et if EXP[t][sp]!="nan"],
             mk[sp], color=colors[sp], ms=7, mec="k", mew=0.5, ls="none")
axA.set_xlabel("Time (s)"); axA.set_ylabel("Phosphorylation (AU)"); axA.set_xlim(0, 3600)
axA.set_title("A  WT posterior-predictive median (Fig 6)"); axA.legend(frameon=False)
vcol = {"WT":"#0072B2","KO":"#E69F00","N78G":"#56B4E9","T292A":"#009E73","T292D":"#CC79A7"}
vls = {"WT":"-","KO":"--","N78G":":","T292A":"-.","T292D":"-"}
for m in MODELS:
    axB.plot(G[m]["time"], [x/1e5 for x in G[m][f"MEK_pRDS_{m}"]], vls[m], color=vcol[m], lw=2, label=m)
for tx in (300, 1800, 3600): axB.axvline(tx, color="k", ls=":", lw=0.7, alpha=0.5)
axB.set_xlabel("Time (s)"); axB.set_ylabel(r"Phosphorylated MEK ($\times10^{5}$)"); axB.set_xlim(0, 3600)
axB.set_title("B  MEK_pRDS median, 5 cell lines (Fig 5)"); axB.legend(frameon=False)
fig.suptitle(f"Miller-2026 aMCMC reproduction at posterior-median d3,u3 (K&L baseline) — "
             f"F_quant={F_quant:.1f}, {nsat}/{ntot} BPSL", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(HERE, "mek_isoform_amcmc_reproduction.png"); fig.savefig(out, dpi=130)
print("wrote", out)
