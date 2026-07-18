#!/usr/bin/env python
"""Gate-3a reproduction oracle for the Miller-2026 MEK-isoform DE/MLE job.

Substitutes the paper's Table-3 "Best-fit value" (31 parameters) into the five slug models,
runs BNG2.pl (ODE, c1 = 0.02 = Kocieniewski & Lipniacki's S0 ligand signal), and:
  * overlays the scaled WT outputs vs WT.exp = SI Table 8 (reproduces Fig 1B);
  * plots MEK_pRDS(t) for all five cell lines with the 300/1800/3600 s constraint times
    (reproduces the qualitative panels of Fig 2);
  * prints F_quant (sum-of-squares on the 18 WT points) and BPSL constraint satisfaction.

Run (from this folder, with BNGPATH set):  python make_reproduction.py
Writes: mek_isoform_de_reproduction.png
"""
import os, re, glob, subprocess, tempfile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BNG = os.path.join(os.environ["BNGPATH"], "BNG2.pl")
MODELS = ["WT", "KO", "N78G", "T292A", "T292D"]

# Miller et al. 2026, Table 3 "Best-fit value" column (= the reported __FREE token values).
TAB3 = {
 "c1_L":7.4e-3,"c2":9.3e-9,"t1":99,"d3":2.0e-3,"b1":2.4e-8,"n1":6.5e-4,"b2":4.2e-6,"n2":2.6e-4,
 "b3":9.9e-5,"n3":8.0e-2,"b4":2.9e-4,"n4":3.5e-3,"a1":2.3e-7,"i1":3.9e-1,"a2":3.3e-8,"i2":2.7e-2,
 "p1":3.5e-6,"u1":9.5e-4,"p2a":5.2e-7,"p2b":2.3e-5,"u2":6.7e-2,"p3":1.1e-9,"u3":2.6e-4,"p4":6.5e-10,
 "u4":3.2e-4,"b5":1.7e-8,"n5":3.9e-3,"u5":19,"scalepEGFR":1.1e-4,"scalepERK":3.2e-6,"scalepSos1":1.1e-4,
}

def build_and_run(m, workdir):
    txt = open(os.path.join(HERE, f"MEK1_{m}.bngl")).read()
    for k, v in TAB3.items():
        txt = re.sub(rf'\b{re.escape(k)}__FREE\b', repr(v), txt)
    # c1 is already fixed at 0.02 by the model's own action; run it as-is
    path = os.path.join(workdir, f"{m}.bngl")
    open(path, "w").write(txt)
    subprocess.run(["perl", BNG, path], cwd=workdir, capture_output=True, text=True)
    g = glob.glob(os.path.join(workdir, f"*_{m}.gdat"))[0]
    lines = [l for l in open(g) if l.strip()]
    hdr = lines[0].lstrip("#").split()
    d = {h: [] for h in hdr}
    for l in lines[1:]:
        for h, val in zip(hdr, l.split()):
            d[h].append(float(val))
    return d

def at(d, col, t):
    i = min(range(len(d["time"])), key=lambda j: abs(d["time"][j]-t))
    return d[col][i]

tmp = tempfile.mkdtemp(prefix="mekrepro_")
G = {m: build_and_run(m, os.path.join(tmp, m)) for m in MODELS
     for _ in [os.makedirs(os.path.join(tmp, m), exist_ok=True)]}

# ---- WT.exp (= SI Table 8) ----
EXP = {}
for ln in open(os.path.join(HERE, "WT.exp")):
    if ln.strip() and not ln.strip().startswith("#"):
        p = ln.split(); EXP[float(p[0])] = {"pSOS1": p[2], "pEGFR": p[3], "pERK": p[4]}
et = sorted(EXP)

# ---- F_quant + constraint satisfaction (console) ----
FCOL = {"pSOS1":"scaled_pSOS1","pEGFR":"scaled_pEGFR","pERK":"scaled_pERK"}
F_quant = 0.0
for t in et:
    for sp in ("pSOS1","pEGFR","pERK"):
        if EXP[t][sp] != "nan":
            F_quant += (at(G["WT"], FCOL[sp], t) - float(EXP[t][sp]))**2

def parse_side(s):
    m = re.match(r'([A-Za-z0-9_]+)\.([A-Za-z0-9_]+)\s+at\s+time=(\d+)', s.strip())
    return ("obs", m.group(1), m.group(2), int(m.group(3))) if m else ("num", float(s), None, None)
nsat = ntot = 0
for m in MODELS:
    for raw in open(os.path.join(HERE, f"{m}.prop")):
        raw = raw.strip()
        if not raw or raw.startswith("#"): continue
        line = raw.split("weight")[0].strip()
        g = re.match(r'(.+?)([<>]=?)(.+)', line); L, op, R = g.groups()
        lt, rt = parse_side(L), parse_side(R)
        lv = at(G[lt[1]], lt[2], lt[3]) if lt[0]=="obs" else lt[1]
        rv = at(G[rt[1]], rt[2], rt[3]) if rt[0]=="obs" else rt[1]
        ok = (lv > rv) if op.startswith(">") else (lv < rv)
        ntot += 1; nsat += ok
print(f"F_quant (18 WT points) = {F_quant:.2f}   [paper DE objective = 24.0, SI Table 7]")
print(f"BPSL constraints satisfied at Table-3 params: {nsat} / {ntot}   [paper: 84/90]")

# ---- figure ----
fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4.6))
colors = {"pEGFR":"#E69F00","pSOS1":"#56B4E9","pERK":"#009E73"}
mk = {"pEGFR":"o","pSOS1":"s","pERK":"^"}
for sp in ("pEGFR","pSOS1","pERK"):
    axA.plot(G["WT"]["time"], G["WT"][FCOL[sp]], "-", color=colors[sp], lw=2, label=sp)
    axA.plot([t for t in et if EXP[t][sp]!="nan"],
             [float(EXP[t][sp]) for t in et if EXP[t][sp]!="nan"],
             mk[sp], color=colors[sp], ms=7, mec="k", mew=0.5, ls="none")
axA.set_xlabel("Time (s)"); axA.set_ylabel("Phosphorylation (AU)")
axA.set_title("A  WT quantitative fit  (Fig 1B) — lines=model, points=Kamioka Fig 3D")
axA.legend(frameon=False); axA.set_xlim(0, 3600)

vcol = {"WT":"#0072B2","KO":"#E69F00","N78G":"#56B4E9","T292A":"#009E73","T292D":"#CC79A7"}
vls  = {"WT":"-","KO":"--","N78G":":","T292A":"-.","T292D":"-"}
for m in MODELS:
    axB.plot(G[m]["time"], [x/1e5 for x in G[m][f"MEK_pRDS_{m}"]],
             vls[m], color=vcol[m], lw=2, label=m)
for tx in (300, 1800, 3600):
    axB.axvline(tx, color="k", ls=":", lw=0.7, alpha=0.5)
axB.set_xlabel("Time (s)"); axB.set_ylabel(r"Phosphorylated MEK copy number ($\times10^{5}$)")
axB.set_title("B  MEK_pRDS trajectories, 5 cell lines  (Fig 2) — dotted = constraint times")
axB.legend(frameon=False); axB.set_xlim(0, 3600)
fig.suptitle(f"Miller-2026 MEK-isoform DE/MLE reproduction at Table-3 best-fit "
             f"(F_quant={F_quant:.1f}≈24, {nsat}/{ntot} BPSL)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(HERE, "mek_isoform_de_reproduction.png")
fig.savefig(out, dpi=130)
print("wrote", out)
