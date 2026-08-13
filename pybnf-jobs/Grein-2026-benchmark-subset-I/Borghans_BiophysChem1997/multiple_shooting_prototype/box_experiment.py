"""The acceptance-style row: multistart from the PEtab fit box, no fitted solution
supplied, no burst-specific knots.

This is the honest benchmark framing issue #563 asks for -- starts drawn from the
job's own box (`loguniform_var  0.001 100000` on 19 parameters, `uniform_var 0 1` on
the three initial states), not perturbations of a known point.  It produces the
"multistart multiple-shooting GNTR" row of the acceptance table, with plain multistart
single shooting under an identical budget as its control.

    python box_experiment.py <out.jsonl> <shard> <nshards> [nstarts]
"""
import json
import sys
import time

import numpy as np

import msproto as M
from experiment import BUDGET, SOLVED, FLAT


def box_start(seed):
    """A Latin-hypercube-free, plain uniform draw from the fit box in sampling space
    -- the same box the shipped `Borghans_BiophysChem1997.conf` searches."""
    rng = np.random.default_rng(10_000 + seed)
    lo, hi = M.Layout(1).bounds()
    return lo + rng.random(len(lo)) * (hi - lo)


def main():
    out_path, shard, nshards = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    nstarts = int(sys.argv[4]) if len(sys.argv) > 4 else 24
    jobs = [(meth, s) for s in range(nstarts) for meth in ('ss', 'ms')]
    jobs = [j for i, j in enumerate(jobs) if i % nshards == shard]
    with open(out_path, 'a', buffering=1) as fh:
        for meth, seed in jobs:
            x0 = box_start(seed)
            t0 = time.time()
            if meth == 'ss':
                _, info = M.solve_single_shooting(
                    x0, deadline=time.monotonic() + BUDGET)
                stages = []
            else:
                # 4-2-1, the best schedule in the segment-count comparison: the
                # many-short-segments end of the homotopy is under-determined on this
                # problem (one observed state of three, ~14 points per segment).
                _, info = M.solve_homotopy(x0, schedule=(4, 2, 1), budget=BUDGET)
                stages = info['stages']
            J = info['certified']
            rec = {'method': meth, 'seed': seed, 'radius': -1.0,
                   'start_J': M.certify(x0),
                   'certified': None if not np.isfinite(J) else J,
                   'OG': None if not np.isfinite(J) else M.optimality_gap(J),
                   'solved': bool(np.isfinite(J) and J < SOLVED),
                   'kept_dynamics': bool(np.isfinite(J) and J < FLAT),
                   'n_sim': info['n_sim'], 'seconds': time.time() - t0,
                   'stages': stages}
            fh.write(json.dumps(rec) + '\n')
            print('[%d] %-3s seed=%-3d  start %9.2f -> %9.2f  %s  %5.0fs'
                  % (shard, meth, seed, rec['start_J'],
                     rec['certified'] if rec['certified'] is not None else float('nan'),
                     'SOLVED' if rec['solved'] else
                     ('dyn' if rec['kept_dynamics'] else '   '),
                     rec['seconds']), flush=True)


if __name__ == '__main__':
    main()
