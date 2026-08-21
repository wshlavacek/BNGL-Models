#!/usr/bin/env python3
"""Perturbations of the PEtab nominal point at a fixed radius -- the regime the #563
prototype's one solve actually came from.

Everything measured on Borghans so far compares methods where neither works: uninformed box
draws (0/24 either way) and oscillating box draws (a wash at 32x the cost). The single result
that motivates this whole feature is different in kind -- the prototype reached
``OG = -1.282656`` from a **radius 0.4 perturbation of the PEtab nominal point**, and it is
the only run that has ever solved this problem by any method in this repository.

That number is doing all the work and it has never been replicated. It is one run. The
prototype's own paired sweeps at convergence-region radii were 24-24 with medians tied, which
is not obviously compatible with "multiple shooting solves this and single shooting does not".
So the decisive experiment is the narrow one: at exactly that radius, from the same kind of
start, does ``ms`` solve where ``gntr`` does not -- reproducibly, across seeds?

This uses the nominal point, which is privileged information. That is the POINT: the claim
under test is a claim about the convergence region, and you cannot test a convergence-region
claim from starts that are not in it. It is not an acceptance-benchmark result and must never
be reported as one.

The perturbation is isotropic in each parameter's own sampling space (log10 for the
loguniform parameters, linear for the three bounded initials), which is the space the
prototype perturbed in and the space every optimizer here searches.

Usage:  python make_radius_starts.py --radius 0.4 --n 10
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NOMINAL = HERE / "Borghans_gntr_nominal.conf"
LINEAR = {"init_A_state", "init_Y_state", "init_Z_state"}
LOG_BOUNDS = (1e-3, 1e5)
LIN_BOUNDS = (0.0, 1.0)

LINE = re.compile(
    r"^parameter:\s*(\w+),\s*prior:\s*normal,\s*(?:parameter_scale:\s*(\w+),\s*)?"
    r"mean:\s*([-\d.eE+]+)")


def nominal_point() -> dict:
    """``{name: value}`` at the PEtab nominal point, read off the shipped conf."""
    out = {}
    for line in NOMINAL.read_text().splitlines():
        m = LINE.match(line.strip())
        if not m:
            continue
        name, scale, mean = m.group(1), m.group(2), float(m.group(3))
        out[name] = 10.0 ** mean if scale == "log10" else mean
    return out


def perturb(point: dict, radius: float, rng) -> dict:
    """One isotropic perturbation of ``point`` at ``radius``, in sampling space.

    The direction is uniform on the sphere and the magnitude is exactly ``radius``, so every
    start sits at the same distance -- a shell, not a ball. A ball would confound "how far
    out does this still work" with "how often", which is the question this is asking.
    """
    names = sorted(point)
    direction = rng.normal(size=len(names))
    direction /= np.linalg.norm(direction)
    out = {}
    for name, step in zip(names, direction * radius):
        if name in LINEAR:
            value = float(np.clip(point[name] + step, *LIN_BOUNDS))
        else:
            value = float(np.clip(10.0 ** (math.log10(point[name]) + step), *LOG_BOUNDS))
        out[name] = value
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", type=float, default=0.4)
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--out", default="radius_starts.json")
    args = parser.parse_args()

    point = nominal_point()
    if len(point) != 23:
        raise SystemExit(f"expected 23 nominal parameters, parsed {len(point)}")
    rng = np.random.default_rng(args.seed)
    starts = [{"params": perturb(point, args.radius, rng),
               "peaks": -1, "amplitude": float("nan"),
               "relative_residual": float("nan"), "radius": args.radius}
              for _ in range(args.n)]
    (HERE / args.out).write_text(json.dumps(starts, indent=2) + "\n")
    print(f"{args.n} starts at radius {args.radius} around the PEtab nominal point "
          f"({len(point)} parameters) -> {args.out}")
    print("NOTE: these use the nominal point, which is privileged information. This is a "
          "convergence-region measurement, NOT an acceptance-benchmark result.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
