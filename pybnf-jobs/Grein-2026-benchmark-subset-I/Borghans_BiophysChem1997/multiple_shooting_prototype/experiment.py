"""Convergence-region sweep: single shooting vs multiple shooting from the same starts.

The claim under test is issue #563's motivation -- that multiple shooting "can enlarge
the useful convergence region without changing the final scientific objective".  The
measurement is therefore basin size, not raw multistart luck: draw starts by perturbing
the PEtab nominal point at increasing radius, run both transcriptions from the SAME
start, and score both through the SAME certified single-shoot reconstruction.

Perturbing a known reference point is a basin measurement, not a benchmark result; the
acceptance benchmark in #563 (no fitted solution supplied) is a separate run.

    python experiment.py <out.jsonl> <shard> <nshards> [radii] [seeds]
"""
import json
import os
import sys
import time

import numpy as np

import msproto as M

SOLVED = M.JSTAR + M.SOLVED_OG - M.PAPER_OFFSET
#: "Kept the dynamics": better than the best no-dynamics (flat-line) solution the
#: global searches all collapse to, which #563 measures at -165.98.
FLAT = -165.98


def start_point(seed, radius):
    """A start drawn by perturbing the PEtab nominal point.  Log10-scaled parameters
    get a N(0, radius) kick in their own sampling space; the linear init_* states get
    a proportionally smaller one (their box is one unit wide, not eight)."""
    theta = M.nominal_theta()
    rng = np.random.default_rng(seed)
    x = theta.copy()
    x[:19] += rng.normal(scale=radius, size=19)
    x[19:22] = np.clip(x[19:22] + rng.normal(scale=radius * 0.3, size=3), 0.01, 0.99)
    lo, hi = M.Layout(1).bounds()
    return np.clip(x, lo, hi)


#: Wall-clock seconds each arm may spend on one start.  Both arms get the same cap;
#: single shooting converges and stops well inside it, which is reported as its
#: actual cost rather than hidden.
BUDGET = 240.0


def run_one(method, seed, radius, **kw):
    x0 = start_point(seed, radius)
    t0 = time.time()
    if method == 'ss':
        _, info = M.solve_single_shooting(
            x0, deadline=time.monotonic() + BUDGET)
        stages = []
    else:
        schedule = kw.pop('schedule', (8, 4, 2, 1))
        _, info = M.solve_homotopy(x0, schedule=schedule, budget=BUDGET, **kw)
        stages = info.get('stages', [])
    J = info['certified']
    return {
        'method': method, 'seed': seed, 'radius': radius,
        'start_J': M.certify(x0),
        'certified': None if not np.isfinite(J) else J,
        'OG': None if not np.isfinite(J) else M.optimality_gap(J),
        'solved': bool(np.isfinite(J) and J < SOLVED),
        'kept_dynamics': bool(np.isfinite(J) and J < FLAT),
        'n_sim': info['n_sim'], 'seconds': time.time() - t0,
        'stages': stages,
    }


def main():
    out_path = sys.argv[1]
    shard = int(sys.argv[2])
    nshards = int(sys.argv[3])
    radii = [float(r) for r in (sys.argv[4] if len(sys.argv) > 4
                                else '0.1,0.2,0.4,0.8').split(',')]
    seeds = list(range(int(sys.argv[5]) if len(sys.argv) > 5 else 8))

    jobs = [(meth, s, r) for r in radii for s in seeds for meth in ('ss', 'ms')]
    jobs = [j for i, j in enumerate(jobs) if i % nshards == shard]
    with open(out_path, 'a', buffering=1) as fh:
        for meth, s, r in jobs:
            rec = run_one(meth, s, r)
            fh.write(json.dumps(rec) + '\n')
            print('[%d] %-3s seed=%-3d r=%.2f  start %8.2f -> %8.2f  '
                  '%s  %5.1fs' % (shard, meth, s, r, rec['start_J'],
                                  rec['certified'] if rec['certified'] is not None
                                  else float('nan'),
                                  'SOLVED' if rec['solved'] else
                                  ('dyn' if rec['kept_dynamics'] else '   '),
                                  rec['seconds']), flush=True)


if __name__ == '__main__':
    main()
