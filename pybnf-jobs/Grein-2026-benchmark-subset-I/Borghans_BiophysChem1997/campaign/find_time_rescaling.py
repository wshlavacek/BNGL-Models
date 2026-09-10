#!/usr/bin/env python3
"""Which parameters carry the model's time-rescaling symmetry? Determined, not assumed.

ADR-0109 records that Borghans has an exact time-rescaling direction -- "multiply its 9 rate
constants by alpha and Z(t) -> Z(alpha t)". If scaling a set S by alpha is a time rescaling,
then simulating with S scaled and reading at time t must reproduce the unscaled model read at
time alpha*t:

    y(t ; theta with S*alpha)  ==  y(alpha*t ; theta)

The set follows from the equations: every term of every right-hand side must pick up one
factor of alpha, so S is the six velocities (Vd, Vm2, Vm3, Vp, v0, v1) plus the three
first-order rate constants (Kf, K_par, epsilon_par). Concentration constants (K2, Ka, Kd, Kp,
Ky, Kz), the Hill coefficient, beta (it only multiplies v1 and Vp) and the initial values are
not in S. This script VERIFIES that set numerically and shows the 7-parameter set the first
version reported (no Kf, no epsilon_par) is not a symmetry.

Why the first version got it wrong, recorded so the mistake is not repeated: it tested at
three random draws from the prior box. Draws from this box are flat -- settled to a steady
state -- by t ~ 0.03 in all but ~1 in 8,000 cases, and a constant trajectory satisfies
y(t) == y(alpha*t) for ANY set, so the test was vacuous and stopped at the first set that
happened to pass its tolerance. A symmetry test needs test points with dynamics on the grid;
this version takes the model file's own parameter values (the SBML defaults, which oscillate)
and perturbations of them, and rejects any test point whose trajectory is nearly constant.

Nothing here touches the data or any fitted point -- the symmetry is a property of the ODE
system, derivable by anyone holding the model file and before ever looking at an observation.

Running
-------
Needs PyBNF's own interpreter, not this repo's. The `pybnf` imports here reach
pybnf.pset / pybnf.parse, which pull in roadrunner and distributed; those live in
PyBNF's venv and are not installed in BNGL-Models'. `.envrc.local` exports PYBNF_PY for this:

    "$PYBNF_PY" find_time_rescaling.py
"""

from __future__ import annotations

import itertools
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

CAMPAIGN = Path(__file__).resolve().parent
# The job directory. This script lives in campaign/, but the model, the .exp data and
# every run happen one level up, and the confs' paths are relative to it.
HERE = CAMPAIGN.parent
os.chdir(HERE)
# Run with the PyBNF environment's interpreter, per this collection's convention.
# Set PYBNF_SRC to prepend a source checkout instead.
_pybnf_src = os.environ.get("PYBNF_SRC")
if _pybnf_src:
    sys.path.insert(0, _pybnf_src)

CONF = "Borghans_BiophysChem1997.conf"
MODEL = "model_Borghans_BiophysChem1997.xml"
HORIZON = 9.0
ALPHAS = (0.5, 1.37)   # not near 1, so a false positive cannot hide in the tolerance
TOL = 1e-4             # the floor is bngsim's own integration tolerance (~4e-6 here); the wrong set is at ~1e0
MIN_DYNAMICS = 0.2     # a test point must move by >= 20 % of its max over the grid

VELOCITIES = ["Vd", "Vm2", "Vm3", "Vp", "v0", "v1"]
FIRST_ORDER = ["Kf", "K_par", "epsilon_par"]
EXPECTED = sorted(VELOCITIES + FIRST_ORDER)
REPORTED_2026_08 = sorted(VELOCITIES + ["K_par"])       # what the first version printed
MAYBE_RATE = ["Kf", "K_par", "epsilon_par", "Kd", "K2", "Ka", "Ky", "Kz", "Kp",
              "beta_par", "n_par"]


def simulate(backend, variables, params, times):
    """Z_state read AT the requested times (the segment backend accepts explicit sample
    times; interpolating an oscillator from the model's own grid costs ~1e-3, the order of
    the mismatch being measured)."""
    from pybnf.pset import PSet
    from pybnf.shooting.backend import SegmentSimulationFailed
    pset = PSet([v.set_value(params[v.name]) for v in variables])
    try:
        data = backend.simulate(pset, np.asarray(times, dtype=float), None)
    except (SegmentSimulationFailed, Exception):
        return None
    arr = np.asarray(data.data, dtype=float)
    idx = data.cols["Z_state"]
    grid = arr[:, data.cols[data.indvar]]
    if len(grid) != len(times) or not np.allclose(grid, times, rtol=0, atol=1e-9):
        return None
    return arr[:, idx]


def sbml_defaults():
    """The model file's own parameter values -- model knowledge, not a fitted point."""
    root = ET.parse(MODEL).getroot()
    ns = {"s": root.tag.split("}")[0].strip("{")}
    vals = {p.get("id"): float(p.get("value")) for p in root.iterfind(".//s:parameter", ns)}
    # scale/offset/sigma do not enter Z_state; any in-box value will do.
    vals.setdefault("scale", 1.0); vals.setdefault("offset", 0.3); vals.setdefault("sigma", 0.1)
    return vals


def main() -> int:
    from pybnf.parse import load_config
    from pybnf.pset import MutationSet
    from pybnf.shooting import BngsimSegmentBackend

    config = load_config(CONF)
    variables = config.variables
    model = list(config.models.values())[0]
    action = model.actions[0]
    backend = BngsimSegmentBackend(model, action, MutationSet(), action.suffix, timeout=10)
    rng = np.random.default_rng(7)

    base_grid = np.linspace(0.2, HORIZON / max(ALPHAS), 60)

    def dynamic(params):
        y = simulate(backend, variables, params, base_grid * max(ALPHAS))
        if y is None:
            return False
        return (np.max(y) - np.min(y)) / max(np.max(np.abs(y)), 1e-12) >= MIN_DYNAMICS

    # Test points WITH dynamics: the SBML defaults plus two log-normal perturbations of them.
    defaults = sbml_defaults()
    if not dynamic(defaults):
        print("the SBML default point does not oscillate on the grid; cannot test"); return 2
    points = [defaults]
    tries = 0
    while len(points) < 3 and tries < 200:
        tries += 1
        p = {k: (v * 10 ** rng.normal(0.0, 0.3) if not k.startswith("init_") else v)
             for k, v in defaults.items()}
        if dynamic(p):
            points.append(p)
    print(f"{len(points)} test points with dynamics (>= {MIN_DYNAMICS:.0%} of max over the grid)\n")

    def mismatch(subset, params, alpha):
        grid = np.linspace(0.2, HORIZON / alpha, 60)
        base = simulate(backend, variables, params, grid * alpha)
        scaled = dict(params)
        for name in subset:
            scaled[name] = params[name] * alpha
        got = simulate(backend, variables, scaled, grid)
        if base is None or got is None:
            return np.inf
        denom = max(float(np.max(base) - np.min(base)), 1e-12)
        return float(np.max(np.abs(got - base)) / denom)

    def worst(subset):
        return max(mismatch(subset, p, a) for p in points for a in ALPHAS)

    print("y(t; S*alpha) vs y(alpha*t; theta), max |diff| / range, over all points and alphas")
    w_exp = worst(EXPECTED); w_old = worst(REPORTED_2026_08); w_vel = worst(VELOCITIES)
    print(f"  expected 9-rate set {EXPECTED}: {w_exp:.2e}")
    print(f"  set reported 2026-08 (no Kf, no epsilon_par): {w_old:.2e}")
    print(f"  velocities only: {w_vel:.2e}\n")

    # The blind search the first version ran, now at points where it can fail.
    found = None
    for k in range(0, 4):
        for extra in itertools.combinations(MAYBE_RATE, k):
            subset = VELOCITIES + list(extra)
            if worst(subset) < TOL:
                found = sorted(subset)
                break
        if found:
            break
    print(f"  blind search over velocities + up to 3 of {MAYBE_RATE}: "
          f"{'EXACT ' + str(found) if found else 'nothing exact'}")

    ok = w_exp < TOL and w_old > 1e-2 and found == EXPECTED
    print("\nRESULT:", "the 9-rate set is the exact time rescaling; the 7-set is not"
          if ok else "UNEXPECTED -- inspect the numbers above")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
