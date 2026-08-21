#!/usr/bin/env python3
"""Which parameters carry the model's time-rescaling symmetry? Determined, not assumed.

ADR-0109 records that Borghans has an exact time-rescaling direction -- "multiply its 9 rate
constants by alpha and Z(t) -> Z(alpha t)" -- and that only a -4.5 %/+2.3 % window in period
beats a flat line. That makes period a one-dimensional, analytically-known direction through
the parameter space, which is a far better handle on this problem than any search.

But the ADR does not say WHICH nine, and guessing from parameter names would make every
result downstream meaningless. So this determines the set empirically, from the model alone:
if scaling a set S by alpha is a time rescaling, then simulating with S scaled and reading at
time t must reproduce the unscaled model read at time alpha*t.

    y(t ; theta with S*alpha)  ==  y(alpha*t ; theta)

Nothing here touches the data or any fitted point -- the symmetry is a property of the ODE
system, derivable by anyone holding the model file and before ever looking at an observation.
That is what makes exploiting it legitimate rather than privileged.

Running
-------
Needs PyBNF's own interpreter, not this repo's. The `pybnf` imports here reach
pybnf.pset / pybnf.parse, which pull in roadrunner and distributed; those live in
PyBNF's venv and are not installed in BNGL-Models'. Under plain python3 the import
succeeds at the top level and fails later, inside the run. `.envrc.local` exports
PYBNF_PY for this:

    "$PYBNF_PY" find_time_rescaling.py
"""

from __future__ import annotations

import itertools
import os
import sys
from pathlib import Path

import numpy as np

CAMPAIGN = Path(__file__).resolve().parent
# The job directory. This script lives in campaign/, but the model, the .exp data and
# every run happen one level up, and the confs' paths are relative to it.
HERE = CAMPAIGN.parent
os.chdir(HERE)
# Run with the PyBNF environment's interpreter, per this collection's convention.
# Set PYBNF_SRC to prepend a source checkout instead (that is how these runs were made).
_pybnf_src = os.environ.get("PYBNF_SRC")
if _pybnf_src:
    sys.path.insert(0, _pybnf_src)

CONF = "Borghans_bench_a1_cmaes.conf"
HORIZON = 9.0
ALPHA = 1.37          # not near 1, so a false positive cannot hide in the tolerance

# Candidates, coarsest first. The names are the model's own; which of them are *rates* is
# exactly the question, so every set is tested rather than argued for.
VELOCITIES = ["Vd", "Vm2", "Vm3", "Vp", "v0", "v1"]
MAYBE_RATE = ["Kf", "K_par", "epsilon_par", "Kd", "K2", "Ka", "Ky", "Kz", "Kp",
              "beta_par", "n_par"]


def simulate(backend, variables, params, times):
    """Z_state read AT the requested times.

    Via the multiple-shooting segment backend rather than model.execute, because it accepts
    explicit sample times. The first version of this script interpolated from the model's own
    output grid onto the target times, and linear interpolation of an oscillator costs ~1e-3 --
    the same order as the mismatch being measured, so the test could not resolve its own
    question. Simulating at the times removes that error entirely.
    """
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


def main() -> int:
    from pybnf.parse import load_config
    from pybnf.pset import PSet

    config = load_config(CONF)
    variables = config.variables
    model = list(config.models.values())[0]
    from pybnf.shooting import BngsimSegmentBackend
    from pybnf.pset import MutationSet
    action = model.actions[0]
    backend = BngsimSegmentBackend(model, action, MutationSet(),
                                   action.suffix, timeout=10)
    rng = np.random.default_rng(7)

    # Three unrelated points, so a set that happens to work at one is not mistaken for a
    # symmetry of the system.
    points = []
    while len(points) < 3:
        p = {v.name: v.value_from_quantile(rng.random()).value for v in variables}
        if simulate(backend, variables, p, np.linspace(0.1, 1.0, 5)) is not None:
            points.append(p)

    grid = np.linspace(0.2, HORIZON / ALPHA, 40)     # inside the horizon after rescaling

    def mismatch(subset, params):
        base = simulate(backend, variables, params, grid * ALPHA)
        scaled = dict(params)
        for name in subset:
            scaled[name] = params[name] * ALPHA
        got = simulate(backend, variables, scaled, grid)
        if base is None or got is None:
            return np.inf
        denom = max(float(np.max(np.abs(base))), 1e-12)
        return float(np.max(np.abs(got - base)) / denom)

    print(f"testing y(t; S*alpha) == y(alpha*t; theta)  at alpha = {ALPHA}\n")
    tested = []
    for k in range(0, len(MAYBE_RATE) + 1):
        for extra in itertools.combinations(MAYBE_RATE, k):
            subset = VELOCITIES + list(extra)
            worst = max(mismatch(subset, p) for p in points)
            tested.append((worst, tuple(subset)))
            if worst < 1e-6:
                print(f"  EXACT (max rel mismatch {worst:.2e}, {len(subset)} params): "
                      f"{sorted(subset)}")
                return 0
        best = min(tested)
        print(f"  no exact set with {k} extra; best so far {best[0]:.3e} "
              f"({len(best[1])} params: {sorted(best[1])})")
        if k >= 3:
            break
    print("\nno exact time-rescaling set found among the candidates tested")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
