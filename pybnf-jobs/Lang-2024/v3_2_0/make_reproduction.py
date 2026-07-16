#!/usr/bin/env python
"""Reproduction figure for the Lang-2024 v3_2_0 job.

Simulates v3_2_0.bngl over the fit window (t = 0..1.36e5, the Stallaert 2-cycle grid) and
overlays the model's 8 nuclear observation functions on the fit data (v3_2_0.exp).

By default it simulates at the .bngl nominals, which are the PEtab v3.2.0 `nominalValue`
column (paulflang/cell_cycle_petab v3.2.0). These are a *reference point*, NOT the
cluster-converged best fit: the paper fits 177 parameters with cooperative scatter search
(saCeSS) on a cluster, and that optimum is not shipped in the PEtab problem. So this figure
documents the SETUP (model + data + observation model wired correctly), not a reproduction of
the paper's fit quality -- see VALIDATION.md. Pass --params-file <PyBNF best_fit.txt> to
overlay a fitted parameter set instead.

Requires BNGPATH set (BNG2.pl) and matplotlib/numpy.
Usage: BNGPATH=... python make_reproduction.py [--params-file path/to/pybnf_best.txt]
"""
import argparse
import glob
import os
import subprocess
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BNG = os.path.join(os.environ["BNGPATH"], "BNG2.pl")
MODEL = os.path.join(HERE, "v3_2_0.bngl")
EXP = os.path.join(HERE, "v3_2_0.exp")
T_END = 136471.6126
N_STEPS = 1200

# (exp column, model function, pretty label) in paper order
OBS = [
    ("obs_cycA__nuc_median()",  "obs_cycA__nuc_median",  "cyclin A"),
    ("obs_cycB1__nuc_median()", "obs_cycB1__nuc_median", "cyclin B1"),
    ("obs_cycE__nuc_median()",  "obs_cycE__nuc_median",  "cyclin E"),
    ("obs_E2F1__nuc_median()",  "obs_E2F1__nuc_median",  "E2F1"),
    ("obs_pRB__nuc_median()",   "obs_pRB__nuc_median",   "pRB (Ser807/811)"),
    ("obs_Skp2__nuc_median()",  "obs_Skp2__nuc_median",  "Skp2"),
    ("obs_p21__nuc_median()",   "obs_p21__nuc_median",   "p21"),
    ("obs_p27__nuc_median()",   "obs_p27__nuc_median",   "p27"),
]


def _cols(gdat_path):
    with open(gdat_path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(gdat_path)
    return data[:, 0], {name: data[:, i] for i, name in enumerate(header)}


def simulate(overrides=None):
    """One BNG2.pl run: generate the finite 73-species network and integrate the ODE from
    t=0 to T_END, printing the observation functions. `overrides` is a {param: value} dict."""
    with open(MODEL) as fh:
        src = fh.read().split("end model")[0] + "end model\n"
    setp = "".join(f'setParameter("{k}",{v})\n' for k, v in (overrides or {}).items())
    actions = (
        "begin actions\n"
        "generate_network({overwrite=>1})\n"
        f"{setp}"
        f'simulate({{suffix=>"repro",method=>"ode",t_start=>0,t_end=>{T_END},'
        f"n_steps=>{N_STEPS},print_functions=>1,atol=>1e-8,rtol=>1e-8}})\n"
        "end actions\n"
    )
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "v3_2_0.bngl")
        with open(path, "w") as fh:
            fh.write(src + "\n" + actions)
        r = subprocess.run(["perl", BNG, path], cwd=d, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"BNG2.pl failed (rc={r.returncode}):\n{r.stderr[-2000:]}")
        return _cols(glob.glob(os.path.join(d, "*_repro.gdat"))[0])


def load_params_file(path):
    """Read a PyBNF best-fit params file (`name  value` per line) into a dict."""
    ov = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("#")[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    ov[parts[0]] = float(parts[1])
                except ValueError:
                    pass
    return ov


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params-file", default=None,
                    help="PyBNF best_fit params file to overlay (default: .bngl PEtab nominals)")
    args = ap.parse_args()
    overrides = load_params_file(args.params_file) if args.params_file else None
    label = "model @ fitted params" if overrides else "model @ PEtab v3.2.0 nominals"

    t, cols = simulate(overrides)

    exp = np.loadtxt(EXP)
    et = exp[:, 0]

    fig, axes = plt.subplots(2, 4, figsize=(18, 8))
    rels = []
    sos_total = 0.0
    for k, (expcol, fn, pretty) in enumerate(OBS):
        ax = axes.flat[k]
        y = exp[:, k + 1]
        m_full = cols[fn]
        m_at = np.interp(et, t, m_full)
        rel = np.median(np.abs(m_at - y) / np.abs(y))
        rels.append(rel)
        sos_total += float(np.sum((m_at - y) ** 2))
        ax.plot(et, y, "o", ms=2.2, color="#222", label="data", zorder=5)
        ax.plot(t, m_full, "-", color="#d62728", lw=1.6, label=label)
        ax.set_title(f"{pretty}   (med rel err {rel*100:.0f}%)", fontsize=10)
        ax.set_xlabel("time (s)"); ax.grid(alpha=0.25)
        if k == 0:
            ax.legend(frameon=False, fontsize=8)
    med = float(np.median(rels))
    fig.suptitle(f"Lang 2024 RPE-1 cell cycle (v3.2.0) -- 8 nuclear observables, {label}   "
                 f"(overall median rel err {med*100:.0f}%, SOS {sos_total:.0f})", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = os.path.join(HERE, "v3_2_0_reproduction.png")
    fig.savefig(out, dpi=120)
    print(f"{label}")
    print(f"  overall median |rel err| = {med*100:.1f}%   SOS = {sos_total:.1f}   (8 obs x 600 pts)")
    for (expcol, fn, pretty), r in zip(OBS, rels):
        print(f"    {pretty:20s} median rel err {r*100:5.1f}%")
    print("wrote", out)


if __name__ == "__main__":
    main()
