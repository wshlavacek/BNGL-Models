#!/usr/bin/env python3
"""The Borghans 1997 (Model 2) ODEs, observable and objective, independently of PyBNF.

Everything this slug's README says about the *dynamics* of a fit -- periodic or not, where a
trajectory goes after the record ends -- needs simulations far past t = 9, which PyBNF never
runs (a fit simulates to the last measurement and no further). This is a plain scipy
implementation of the same model, so those statements can be checked by anyone with numpy and
scipy and no PyBNF install.

It is validated against PyBNF's own objective. PyBNF's reduced objective for this problem is
``n * ln(sigma_hat) + const`` with sigma_hat the RMS log10 residual (sigma profiled), so the
*difference* between two points is what is comparable; on that scale this file reproduces
PyBNF exactly:

    point                      this file      PyBNF (score_point.py, noise_profiling = 1)
    nominal - flat line        -32.22         -32.22
    OG -1.28 vector - flat     -82.09         -82.09

(`reduced_objective()` adds the constant so its absolute values match PyBNF too: flat line
-165.982, nominal -198.207, the -1.28 vector -248.069.)

Model: Z cytosolic Ca, Y store Ca, A IP3; observable Ca = scale*Z + offset; lognormal noise,
i.e. Gaussian in log10(Ca). Parameter names are the PEtab problem's.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent.parent          # the job directory
NAMES = ["K2", "K_par", "Ka", "Kd", "Kf", "Kp", "Ky", "Kz", "Vd", "Vm2", "Vm3", "Vp", "beta_par",
         "epsilon_par", "init_A_state", "init_Y_state", "init_Z_state", "n_par", "offset", "scale",
         "sigma", "v0", "v1"]
#: Every parameter that carries one factor of time: scaling all nine by alpha gives Z(alpha t).
RATES = ["v0", "v1", "Vm2", "Vm3", "Kf", "K_par", "Vp", "Vd", "epsilon_par"]
FLAT_PYBNF = -165.9821          # PyBNF's reduced objective of the best horizontal line, sigma profiled


def load_data(path=HERE / "experiment1.exp"):
    d = np.loadtxt(path, comments="#")
    return d[:, 0], d[:, 1]


T_DATA, Y_DATA = load_data()
LOG_Y = np.log10(Y_DATA)
N = len(Y_DATA)


def rhs(t, s, p):
    Z, Y, A = s
    v2 = p["Vm2"] * Z**2 / (p["K2"]**2 + Z**2)
    v3 = (p["Vm3"] * A**4 / (p["Ka"]**4 + A**4) * Y**2 / (p["Ky"]**2 + Y**2)
          * Z**4 / (p["Kz"]**4 + Z**4))
    n = p["n_par"]
    return [p["v0"] + p["beta_par"] * p["v1"] - v2 + v3 + p["Kf"] * Y - p["K_par"] * Z,
            v2 - v3 - p["Kf"] * Y,
            p["Vp"] * p["beta_par"] - p["Vd"] * A**2 / (p["Kp"]**2 + A**2) * Z**n / (p["Kd"]**n + Z**n)
            - p["epsilon_par"] * A]


def solve(p, t_end, rtol=1e-9, atol=1e-12):
    """Dense-output solution on [0, t_end], or None if the integration fails."""
    sol = solve_ivp(rhs, (0.0, float(t_end)),
                    [p["init_Z_state"], p["init_Y_state"], p["init_A_state"]],
                    args=(p,), method="LSODA", rtol=rtol, atol=atol, dense_output=True)
    return sol if sol.success else None


def observe(p, sol, t):
    return p["scale"] * sol.sol(t)[0] + p["offset"]


def sigma_hat(p):
    """RMS log10 residual against the data with the trajectory from the fit's own initial values."""
    sol = solve(p, float(T_DATA[-1]))
    if sol is None:
        return math.nan
    ca = observe(p, sol, T_DATA)
    if np.any(ca <= 0) or not np.all(np.isfinite(ca)):
        return math.nan
    return float(np.sqrt(np.mean((np.log10(ca) - LOG_Y) ** 2)))


FLAT_SIGMA = float(np.std(LOG_Y))       # the best horizontal line, in log space


def reduced_objective(p):
    """PyBNF's reduced objective (sigma profiled) for this problem, matched through the flat line."""
    s = sigma_hat(p)
    return math.nan if not np.isfinite(s) else FLAT_PYBNF + N * math.log(s / FLAT_SIGMA)


def optimality_gap(p):
    """Grein's OG on this slug's scale: solved iff OG < 1.92 (see score.py / README)."""
    return reduced_objective(p) + 246.78635


def rescale_time(p, alpha):
    q = dict(p)
    for k in RATES:
        q[k] *= alpha
    return q


# --- reference vectors ---------------------------------------------------------------------

def params_from_results(path):
    """A parameter dict from a PyBNF ``sorted_params_final.txt`` / ``best_fit_params.txt``."""
    lines = Path(path).read_text().splitlines()
    hdr = [l for l in lines if l.startswith("#")][0].lstrip("#").split()
    row = [l for l in lines if l.strip() and not l.startswith("#")][0].split()
    p = {k: float(v) for k, v in zip(hdr, row) if k in NAMES}
    p.setdefault("sigma", 0.1)
    return p


def nominal():
    """The PEtab nominal point, read from campaign/Borghans_gntr_nominal.conf (its mean: fields)."""
    p = {}
    for line in (HERE / "campaign" / "Borghans_gntr_nominal.conf").read_text().splitlines():
        m = re.match(r"parameter:\s*(\w+),.*?mean:\s*([-+0-9.eE]+)", line)
        if m:
            k, v = m.group(1), float(m.group(2))
            p[k] = v if k.startswith("init_") else 10 ** v
    return p


def paper_fig8():
    """Borghans et al. 1997 Fig. 8 values, as curated in BioModels BIOMD0000000044 (Model 2).
    Initial values are not part of the paper's parameter set; these are the PEtab nominal
    initial values, which only set the phase. scale/offset are placeholders to be profiled."""
    return dict(v0=2.0, v1=1.0, beta_par=0.5, K_par=10.0, K2=0.1, Vm2=6.5, Ka=0.2, Ky=0.2, Vm3=19.5,
                Kz=0.3, Kf=1.0, Vp=2.5, Vd=80.0, Kp=1.0, Kd=0.4, n_par=4.0, epsilon_par=0.1,
                init_Z_state=0.0879205244255038, init_Y_state=0.999348084438687,
                init_A_state=0.99999999999996, scale=1.0, offset=0.3, sigma=0.1)


def best_ever():
    return params_from_results(HERE / "multiple_shooting_prototype" / "verify_best_fit_params.txt")


def best_periodic():
    return json.load(open(HERE / "best_periodic_fit.json"))["params"]


def profile_observation(p):
    """Re-fit only scale and offset (the observation map) to the data; the dynamics untouched."""
    from scipy.optimize import minimize
    sol = solve(p, float(T_DATA[-1]))
    Z = sol.sol(T_DATA)[0]
    A = np.vstack([Z, np.ones_like(Z)]).T
    (sc, off), *_ = np.linalg.lstsq(A, Y_DATA, rcond=None)
    sc, off = max(sc, 1e-3), max(off, 1e-3)

    def f(x):
        s, o = 10 ** x
        ca = s * Z + o
        return 1e9 if np.any(ca <= 0) else np.sqrt(np.mean((np.log10(ca) - LOG_Y) ** 2))

    r = minimize(f, [np.log10(sc), np.log10(off)], method="Nelder-Mead",
                 options={"xatol": 1e-6, "fatol": 1e-9, "maxiter": 2000})
    q = dict(p)
    q["scale"], q["offset"] = 10 ** r.x[0], 10 ** r.x[1]
    return q


if __name__ == "__main__":
    print("validation against PyBNF (differences from the flat line, sigma profiled):")
    for lab, p, ref in [("nominal", nominal(), -32.22), ("OG -1.28 vector", best_ever(), -82.09)]:
        d = reduced_objective(p) - FLAT_PYBNF
        print(f"  {lab:16s} {d:+8.2f}   PyBNF {ref:+.2f}")
