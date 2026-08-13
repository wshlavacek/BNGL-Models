"""Is the segment count the confound?

The first sweep ran one homotopy schedule (8-4-2-1).  Its stage traces show the m = 8
stage certifying WORSE than its own start -- consistent with the segmented problem
being under-determined here: 111 points, one observed state of three, ~14 points per
segment, and 21 free auxiliary states, so the data term can be satisfied without
correct dynamics and continuity carries nearly all the information.  If that is the
story, a coarser transcription should do better, and this measures it against the
single-shooting results already in hand for the same starts.

    python schedules.py <out.jsonl> <shard> <nshards> [radius] [nseeds]
"""
import json
import sys
import time

import numpy as np

import msproto as M
from experiment import BUDGET, SOLVED, FLAT, start_point

SCHEDULES = {
    '2-1': (2, 1),
    '3-1': (3, 1),
    '4-2-1': (4, 2, 1),
    '8-4-2-1': (8, 4, 2, 1),
}


def main():
    out_path, shard, nshards = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    radius = float(sys.argv[4]) if len(sys.argv) > 4 else 0.2
    nseeds = int(sys.argv[5]) if len(sys.argv) > 5 else 8
    jobs = [(name, s) for s in range(nseeds) for name in SCHEDULES]
    jobs = [j for i, j in enumerate(jobs) if i % nshards == shard]
    with open(out_path, 'a', buffering=1) as fh:
        for name, seed in jobs:
            x0 = start_point(seed, radius)
            t0 = time.time()
            _, info = M.solve_homotopy(x0, schedule=SCHEDULES[name], budget=BUDGET)
            J = info['certified']
            rec = {'schedule': name, 'seed': seed, 'radius': radius,
                   'start_J': M.certify(x0),
                   'certified': None if not np.isfinite(J) else J,
                   'solved': bool(np.isfinite(J) and J < SOLVED),
                   'kept_dynamics': bool(np.isfinite(J) and J < FLAT),
                   'n_sim': info['n_sim'], 'seconds': time.time() - t0,
                   'stages': info['stages']}
            fh.write(json.dumps(rec) + '\n')
            print('[%d] %-9s seed=%-3d -> %9.2f  %5.0fs  %s'
                  % (shard, name, seed,
                     rec['certified'] if rec['certified'] is not None else float('nan'),
                     rec['seconds'],
                     ' '.join('%d:%.1f' % (s['m'], s['certified'])
                              for s in info['stages'])), flush=True)


if __name__ == '__main__':
    main()
