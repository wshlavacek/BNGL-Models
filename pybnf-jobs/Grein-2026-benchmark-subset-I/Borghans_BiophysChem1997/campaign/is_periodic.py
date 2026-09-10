#!/usr/bin/env python3
"""Is a fitted parameter vector on a periodic attractor? Three checks, far past the record.

PyBNF simulates a fit only to the last measurement (t = 8.98 here), so a fit's score says
nothing about whether the model oscillates. This runs each vector to t = 400 and applies:

1. Peak train on the second half of the horizon, [200, 400]: number of Z peaks, the median
   spacing and its coefficient of variation (a bursting orbit has uneven spacing), and the
   amplitude trend (largest prominence in the last third of the window over the first third;
   a damped transient trends below 1).
2. Boundedness of the STATE, not the observable: the largest state norm over the second half
   of the window over the first half. A periodic orbit gives 1.000; the OG -1.28 vector gives
   1.33 because store calcium Y grows without bound while Ca barely moves.
3. The fixed point: solve rhs = 0 from the late state, take the Jacobian by central
   differences, and read the eigenvalues. A stable fixed point reached with the trajectory
   speed near zero is a settled transient; an unstable one surrounded by a bounded orbit is
   a periodic attractor. For a vector that is still moving, the fixed point is searched from
   the mass-balance estimate Z* = (v0 + beta*v1)/K_par, which is where the -1.28 vector's
   lives (Z* = 42.6, Y* = 48,711, reached at t ~ 20,000).

Run from anywhere with a python that has numpy and scipy (no PyBNF needed):

    python3 campaign/is_periodic.py                       # the four reference vectors
    python3 campaign/is_periodic.py output/Results/sorted_params_final.txt   # any PyBNF result
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import fsolve
from scipy.signal import find_peaks

sys.path.insert(0, str(Path(__file__).resolve().parent))
import borghans_ode as bo  # noqa: E402

T_END = 400.0


def jacobian(p, x):
    J = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3)
        e[j] = 1e-7 * max(1.0, abs(x[j]))
        J[:, j] = (np.array(bo.rhs(0, x + e, p)) - np.array(bo.rhs(0, x - e, p))) / (2 * e[j])
    return J


def fixed_point(p, guess):
    x, info, ier, msg = fsolve(lambda s: bo.rhs(0, s, p), guess, xtol=1e-13, full_output=True)
    if np.linalg.norm(bo.rhs(0, x, p)) > 1e-8:
        return None, None
    return x, np.linalg.eigvals(jacobian(p, x))


def diagnose(p):
    sol = bo.solve(p, T_END)
    if sol is None:
        return {"verdict": "integration failed"}
    tt = np.linspace(T_END / 2, T_END, 400001)
    X = sol.sol(tt)
    Z = X[0]
    rng = Z.max() - Z.min()
    pk, pr = find_peaks(Z, prominence=0.1 * max(rng, 1e-12))
    ipi = np.diff(tt[pk])
    amp = pr["prominences"]
    n3 = max(len(amp) // 3, 1)
    trend = float(amp[-n3:].max() / amp[:n3].max()) if len(amp) >= 2 else float("nan")
    nrm = np.linalg.norm(X, axis=0)
    half = len(tt) // 2
    growth = float(nrm[half:].max() / nrm[:half].max())
    x_end = sol.sol(T_END)
    speed = float(np.linalg.norm(bo.rhs(0, x_end, p)))
    x_fp, lam = fixed_point(p, x_end)
    if x_fp is None:
        z_star = (p["v0"] + p["beta_par"] * p["v1"]) / p["K_par"]
        x_fp, lam = fixed_point(p, [z_star, max(1.0, x_end[1]) * 10, x_end[2]])
    cv = float(np.std(ipi) / np.mean(ipi)) if len(ipi) > 1 else float("nan")
    if growth > 1.1:
        verdict = "NO ATTRACTOR within the horizon: the state envelope keeps growing"
    elif speed < 1e-4 and lam is not None and lam.real.max() < 0:
        verdict = "TRANSIENT, settled on a stable fixed point"
    elif len(pk) >= 3 and 0.9 < trend < 1.1:
        verdict = "PERIODIC " + ("bursting (uneven spike spacing)" if cv > 0.2 else "(one spike per period)")
    else:
        verdict = "ambiguous: extend the horizon"
    return dict(peaks=len(pk), ipi=float(np.median(ipi)) if len(ipi) else float("nan"), cv=cv,
                trend=trend, growth=growth, speed=speed, x_fp=x_fp, eig=lam, verdict=verdict)


def main(argv):
    if argv:
        vecs = {Path(a).name: (bo.params_from_results(a) if a.endswith(".txt")
                               else __import__("json").load(open(a)).get("params", __import__("json").load(open(a))))
                for a in argv}
    else:
        vecs = {"paper Fig. 8 (BioModels)": bo.profile_observation(bo.paper_fig8()),
                "PEtab nominal": bo.nominal(),
                "best_periodic_fit.json": bo.best_periodic(),
                "OG -1.28 vector": bo.best_ever()}
    print(f"{'vector':26s} {'OG':>7s} {'peaks':>5s} {'spacing':>15s} {'amp trend':>9s} {'envelope':>8s} {'|dx/dt|':>8s}  verdict")
    details = []
    for lab, p in vecs.items():
        d = diagnose(p)
        og = bo.optimality_gap(p)
        if "peaks" not in d:
            print(f"{lab:26s} {og:7.2f}  {d['verdict']}")
            continue
        sp = "-" if not np.isfinite(d["cv"]) else f"{d['ipi']:.2f} CV {d['cv']:.2f}"
        print(f"{lab:26s} {og:7.2f} {d['peaks']:5d} {sp:>15s} {d['trend']:9.3f} {d['growth']:8.3f} {d['speed']:8.1e}  {d['verdict']}")
        details.append((lab, d))
    print("\nfixed points (rhs = 0) and Jacobian eigenvalues:")
    for lab, d in details:
        if d["x_fp"] is None:
            print(f"  {lab:26s} none found")
        else:
            x, lam = d["x_fp"], d["eig"]
            print(f"  {lab:26s} Z*={x[0]:.4g} Y*={x[1]:.4g} A*={x[2]:.4g}  eigenvalues {np.round(lam, 4)}  "
                  f"{'stable' if lam.real.max() < 0 else 'unstable'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
