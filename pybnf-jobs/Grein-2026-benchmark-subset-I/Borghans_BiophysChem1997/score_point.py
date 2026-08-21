#!/usr/bin/env python3
"""Score one parameter point directly, with and without noise profiling.

Validating a start is not the same as running a fit from it: a fit only tells you where it
*ended*. This evaluates the objective at the point itself, which is what says whether a
proposed anchor is the point you think it is.

Usage:  python score_point.py --starts radius_starts.json [--index 0]
        python score_point.py --nominal            # the means of Borghans_gntr_nominal.conf
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


def score(params: dict, conf_name: str) -> float | None:
    from pybnf.parse import load_config
    from pybnf.pset import PSet

    config = load_config(conf_name)
    pset = PSet([v.set_value(params[v.name]) for v in config.variables])
    model = list(config.models.values())[0]
    model.param_set = pset
    try:
        sim = model.execute(str(HERE / "_score"), "score", config.config["wall_time_sim"])
    except Exception as exc:
        print("   simulation failed:", type(exc).__name__, exc)
        return None
    name = model.name
    return config.obj.evaluate_multiple({name: sim}, {name: config.exp_data[name]}, pset,
                                        show_warnings=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--starts")
    parser.add_argument("--index", type=int, default=0)
    parser.add_argument("--nominal", action="store_true")
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()

    if args.nominal:
        from make_radius_starts import nominal_point
        points = [nominal_point()]
    else:
        data = json.loads(Path(args.starts).read_text())
        points = [d["params"] for d in data] if args.all else [data[args.index]["params"]]

    # Two confs differing only in noise_profiling, so the comparison is like-for-like.
    for i, params in enumerate(points):
        plain = score(params, "Borghans_bench_a1_cmaes.conf")
        print(f"point {i}: reduced objective (sigma searched, at its start value) = "
              f"{'failed' if plain is None else f'{plain:.6f}'}")
    print("\nreference: PEtab nominal reduced_objective = -198.1017 (nominal_check.json)")
    print("           flat line = -165.98 ; target ~ -244.87 ; prototype solve = -248.07")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
