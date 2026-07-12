#!/usr/bin/env python
"""Reproduction figure for the fceri_gamma job (Gupta & Mendes 2018 benchmark).

This is a SYNTHETIC-DATA parameter-recovery problem: fceri_gamma2.exp was generated at the
ground-truth parameters (fceri_gamma2_ground_truth.bngl), and a fit tries to recover them. The
"reproduction" is therefore the recovery check's foundation -- the model AT the ground-truth
parameters must reproduce the synthetic data. This script runs the model at those ground-truth
values (Gillespie SSA) and overlays it on fceri_gamma2.exp for all six observables.

The ~58,000-reaction network is expanded ONCE (generate_network, minutes) and N Gillespie SSA
replicates are run within that single BNG2.pl invocation (each a fast `simulate` reusing the
generated network), then averaged -- exactly what the conf's `smoothing` does. Requires BNGPATH
set (uses BNG2.pl) and matplotlib/numpy.
Usage:  BNGPATH=... python make_reproduction.py [n_reps]   (default 5)
"""
import glob
import os
import re
import subprocess
import sys
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BNG = os.path.join(os.environ["BNGPATH"], "BNG2.pl")
MODEL = os.path.join(HERE, "fceri_gamma2.bngl")
GROUND_TRUTH = os.path.join(HERE, "fceri_gamma2_ground_truth.bngl")

OBS = ["LynFree", "RecMon", "RecPbeta", "RecPgamma", "RecSyk", "RecSykPS"]
T_END, N_STEPS = 100, 10
N_REPS = int(sys.argv[1]) if len(sys.argv) > 1 else 5


def ground_truth_params():
    """Read the 20 recovery-target values from fceri_gamma2_ground_truth.bngl."""
    names = ("kp1 kp2 kpL kmL kpLs kmLs kpS kmS kpSs kmSs "
             "pLb pLbs pLg pLgs pLS pLSs pSS pSSs dm dc").split()
    vals = {}
    with open(GROUND_TRUTH) as fh:
        for line in fh:
            m = re.match(r"\s*(\w+)\s+([\d.eE+-]+)\s*$", line)
            if m and m.group(1) in names:
                vals[m.group(1)] = float(m.group(2))
    return vals


def run_ssa_reps(n_reps):
    """Generate the network once at ground-truth params, run n_reps SSA sims; return {obs: (t, mean)}."""
    truth = ground_truth_params()
    with open(MODEL) as fh:
        src = fh.read().split("end reaction rules")[0] + "end reaction rules\n"
    for k, v in truth.items():          # set nominal -> ground-truth
        src = re.sub(rf"(^\s*{k}\s+)[\d.eE+-]+", rf"\g<1>{v:.12g}", src, count=1, flags=re.M)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "fceri.bngl")
        with open(path, "w") as fh:
            fh.write(src)
            # Generate the network once, snapshot the seed state, then run each SSA replicate
            # from that saved initial state. resetConcentrations() before every simulate is
            # REQUIRED: consecutive simulate() actions otherwise CONTINUE from the previous run's
            # end-state (BNG carry-over), so only the first replicate would start from the seed.
            fh.write("\ngenerate_network({overwrite=>1,max_iter=>100})\n")
            fh.write("saveConcentrations()\n")
            for i in range(n_reps):     # each simulate reuses the generated network from the seed state
                fh.write("resetConcentrations()\n")
                fh.write(f'simulate({{method=>"ssa",t_start=>0,t_end=>{T_END},'
                         f'n_steps=>{N_STEPS},seed=>{1000 + i},suffix=>"rep{i}"}})\n')
        subprocess.run(["perl", BNG, path], cwd=d, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        cols_per_rep = []
        for i in range(n_reps):
            gdat = glob.glob(os.path.join(d, f"*rep{i}.gdat"))[0]
            with open(gdat) as fh:
                header = fh.readline().lstrip("#").split()
            data = np.loadtxt(gdat)
            cols_per_rep.append({name: data[:, j] for j, name in enumerate(header)})
    t = cols_per_rep[0]["time"]
    return t, {o: np.mean([c[o] for c in cols_per_rep], axis=0) for o in OBS}


def main():
    exp = np.genfromtxt(os.path.join(HERE, "fceri_gamma2.exp"), names=True, deletechars="")
    print(f"fceri_gamma SSA reproduction at the GROUND-TRUTH parameters (recovery target), "
          f"reps={N_REPS}:")
    t, model = run_ssa_reps(N_REPS)

    fig, ax = plt.subplots(2, 3, figsize=(13, 7.5))
    for a, o in zip(ax.ravel(), OBS):
        yd = exp[o]
        ym = np.interp(exp["time"], t, model[o])
        nz = yd > 3         # relative error only where counts are meaningfully above SSA noise
        rel = np.abs(ym[nz] - yd[nz]) / yd[nz]
        med = np.median(rel) if nz.any() else float("nan")
        a.plot(exp["time"], yd, "o", color="#222", ms=6, zorder=5, label="synthetic data (ground truth)")
        a.plot(t, model[o], "-", color="#1f77b4", lw=2, label="SSA @ ground truth (mean)")
        a.set(title=f"{o}  (median |relerr|={med:.2f})", xlabel="time (s)", ylabel="molecules/cell")
        a.legend(frameon=False, fontsize=8); a.grid(alpha=0.25)
        print(f"  {o:10s} median|relerr|(y>3) = {med:.3f}")
    fig.suptitle("FceRI gamma-chain (Gupta & Mendes 2018) -- SSA at ground-truth params vs. synthetic data\n"
                 "(recovery target for the 20-parameter fit)", fontsize=12)
    fig.tight_layout()
    out = os.path.join(HERE, "fceri_gamma_reproduction.png")
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
