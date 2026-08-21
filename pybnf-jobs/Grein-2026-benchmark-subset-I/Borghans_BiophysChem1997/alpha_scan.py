#!/usr/bin/env python3
"""Exploit the model's exact time-rescaling symmetry: a 1-D scan in period.

Borghans is hard for one reason (ADR-0109): only a -4.5 %/+2.3 % window in *period* beats a
flat line, and everywhere else a correctly-shaped oscillator scores WORSE than fitting no
dynamics at all. So the objective actively steers a search away from the region a fit has to
pass through, and finding that window by search in a 20-dimensional box spanning 8 decades
per axis is what nobody has reliably done.

But period is not really 20-dimensional here. `find_time_rescaling.py` establishes -- from the
model's equations alone, numerically, at three unrelated points -- that

    y(t ; theta with S*alpha) == y(alpha*t ; theta)   to 1.3e-08,  S = RESCALE below

so alpha is an exact one-dimensional direction that moves period at fixed trajectory shape.
That turns "find a 3 % window in 20-D" into "scan a line", which costs ~100 simulations
instead of ~33,000.

WHAT THIS DOES AND DOES NOT USE
-------------------------------
The symmetry is a property of the ODE system, derivable by anyone holding the model file
before ever looking at an observation. It is model knowledge, not solution knowledge, and is
categorically different from the PEtab nominal point (which ADR-0110 correctly calls
privileged).

So, deliberately:

* starts are prior-box draws -- the nominal point is never used;
* the alpha grid is wide and answer-agnostic (`--lo`/`--hi` decades either way), not a window
  centred where the answer is known to be;
* the winner is chosen by the fit's OWN objective on the fit's OWN data, which is fitting;
* no fitted solution enters at any point.

Usage:  python alpha_scan.py --starts oscillating_starts.json --points 121
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
os.chdir(HERE)
# Run with the PyBNF environment's interpreter, per this collection's convention.
# Set PYBNF_SRC to prepend a source checkout instead (that is how these runs were made).
_pybnf_src = os.environ.get("PYBNF_SRC")
if _pybnf_src:
    sys.path.insert(0, _pybnf_src)

CONF = "Borghans_alpha_scoring.conf"

#: The exact time-rescaling set, determined numerically by find_time_rescaling.py rather than
#: guessed from parameter names -- the six velocities alone give a 0.37 mismatch, so K_par is
#: essential and easy to miss.
RESCALE = ("K_par", "Vd", "Vm2", "Vm3", "Vp", "v0", "v1")

BOX = {"log": (1e-3, 1e5), "lin": (0.0, 1.0)}
LINEAR = {"init_A_state", "init_Y_state", "init_Z_state"}


def build():
    from pybnf.parse import load_config
    config = load_config(CONF)
    model = list(config.models.values())[0]
    return config, model


def scored(config, model, params):
    """The fit's own reduced objective at one point, or None."""
    from pybnf.pset import PSet
    pset = PSet([v.set_value(params[v.name]) for v in config.variables])
    model.param_set = pset
    try:
        sim = model.execute(str(HERE / "_alpha"), "alpha", config.config["wall_time_sim"])
    except Exception:
        return None
    name = model.name
    try:
        return config.obj.evaluate_multiple({name: sim}, {name: config.exp_data[name]},
                                            pset, show_warnings=False)
    except Exception:
        return None


def rescaled(params, alpha):
    """``params`` with the time-rescaling set multiplied by ``alpha``, clipped to the box.

    A point pushed outside the fit's own box by the scan is not a point the fit could report,
    so it is clipped rather than silently searched out of bounds.
    """
    out = dict(params)
    lo, hi = BOX["log"]
    for name in RESCALE:
        out[name] = float(np.clip(params[name] * alpha, lo, hi))
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--starts", default="oscillating_starts.json")
    parser.add_argument("--points", type=int, default=121)
    parser.add_argument("--lo", type=float, default=-3.0, help="log10 alpha lower bound")
    parser.add_argument("--hi", type=float, default=3.0, help="log10 alpha upper bound")
    parser.add_argument("--out", default="alpha_scan_results.json")
    args = parser.parse_args()

    config, model = build()
    starts = json.loads(Path(args.starts).read_text())
    grid = 10.0 ** np.linspace(args.lo, args.hi, args.points)

    print(f"alpha grid: {args.points} points over 10^{args.lo:g} .. 10^{args.hi:g}")
    print(f"rescaling set ({len(RESCALE)}): {', '.join(RESCALE)}")
    print("reference: flat line = -165.98 ; PEtab nominal = -198.10 ; target ~ -244.87\n")

    results = []
    for i, item in enumerate(starts, 1):
        base = item["params"]
        j0 = scored(config, model, base)
        curve = []
        for alpha in grid:
            j = scored(config, model, rescaled(base, alpha))
            curve.append(None if j is None or not np.isfinite(j) else float(j))
        finite = [(j, a) for j, a in zip(curve, grid) if j is not None]
        if not finite:
            print(f"start {i:2d}: nothing integrated across the grid")
            continue
        best_j, best_a = min(finite)
        results.append({"start": i, "start_objective": j0, "best_alpha": float(best_a),
                        "best_objective": best_j,
                        "params": rescaled(base, best_a)})
        flag = ""
        if best_j < -198.10:
            flag = "   <-- better than the PEtab nominal point"
        elif best_j < -165.98:
            flag = "   <-- beats a flat line"
        print(f"start {i:2d}: start J = {j0 if j0 is None else f'{j0:9.3f}'}   "
              f"best alpha = {best_a:9.4g}   best J = {best_j:9.3f}{flag}")

    Path(args.out).write_text(json.dumps(results, indent=2) + "\n")
    if results:
        best = min(results, key=lambda r: r["best_objective"])
        print(f"\nbest over all starts: J = {best['best_objective']:.4f} "
              f"at alpha = {best['best_alpha']:.4g} (start {best['start']})")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
