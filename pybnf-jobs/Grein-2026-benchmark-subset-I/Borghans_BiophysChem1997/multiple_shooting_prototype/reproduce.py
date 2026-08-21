"""Re-run the solved start and save its parameter vector for independent verification.

Writes rerun_seed<seed>_r<radius>.json, NOT the committed solved_*.json beside it, and
diffs the two when the recorded one exists. The committed file is the recorded result;
writing there would overwrite it with whatever this run produced, and on a rebuilt ODE
library that is routinely a much worse point -- the solve does not survive a rebuild
(see the slug README). An earlier version of this script did overwrite it.
"""
import json
import os
import sys
import time

import numpy as np

import msproto as M
from experiment import start_point

seed = int(sys.argv[1]) if len(sys.argv) > 1 else 3
radius = float(sys.argv[2]) if len(sys.argv) > 2 else 0.4
budget = float(sys.argv[3]) if len(sys.argv) > 3 else 600.0
schedule = tuple(int(s) for s in (sys.argv[4] if len(sys.argv) > 4 else '8,4,2,1').split(','))

x0 = start_point(seed, radius)
print('seed %d radius %.2f schedule %s  start J = %.4f'
      % (seed, radius, schedule, M.certify(x0)), flush=True)

t0 = time.time()
theta, info = M.solve_homotopy(x0, schedule=schedule, budget=budget, verbose=True)
print('certified J = %.6f   OG = %.6f   %.0fs  nsim=%d'
      % (info['certified'], M.optimality_gap(info['certified']),
         time.time() - t0, info['n_sim']))

L = M.Layout(1)
names = L.theta_names
lin = {}
for i, n in enumerate(names):
    lin[n] = float(10.0 ** theta[i]) if i < len(M.DYN) + len(M.OBS) else float(theta[i])

# independent recomputation from a fresh simulation
tr = M.Transcription(1)
out = tr.evaluate(theta, want_jac=False)
S = float(np.sum(out['r'] ** 2))
sigma = float(np.sqrt(S / M.N_OBS))
J = tr.reduced_objective(out['r'])
print()
print('independent recomputation:')
print('  profiled sigma      = %.10f' % sigma)
print('  reduced objective   = %.6f' % J)
print('  J_paper (-log L)    = %.6f' % M.paper_nll(J))
print('  J* (Grein)          = %.6f' % M.JSTAR)
print('  optimality gap OG   = %.6f   -> %s'
      % (M.optimality_gap(J), 'SOLVED' if M.optimality_gap(J) < 1.92 else 'not solved'))

lin['sigma'] = sigma
recorded = 'solved_seed%d_r%s.json' % (seed, radius)
out_path = 'rerun_seed%d_r%s.json' % (seed, radius)
with open(out_path, 'w') as fh:
    json.dump({'seed': seed, 'radius': radius, 'schedule': list(schedule),
               'reduced_objective': J, 'J_paper': M.paper_nll(J),
               'OG': M.optimality_gap(J), 'sigma_profiled': sigma,
               'params_linear': lin,
               'theta_sampling_space': [float(v) for v in theta],
               'theta_names': names}, fh, indent=2)
print()
print('wrote %s' % out_path)
if os.path.exists(recorded):
    with open(recorded) as fh:
        ref = json.load(fh)
    dOG = M.optimality_gap(J) - ref['OG']
    print('recorded %s: OG = %.6f   this run: OG = %.6f   delta = %+.6f'
          % (recorded, ref['OG'], M.optimality_gap(J), dOG))
    if abs(dOG) > 1e-6:
        print('  the recorded result did NOT reproduce here. That is a known property of\n'
              '  this problem, not necessarily a bug: the path to the solve does not\n'
              '  survive a rebuild of the ODE library. %s is left untouched.' % recorded)
print()
print('parameters (linear scale):')
for n in sorted(lin):
    print('  %-14s %.10g' % (n, lin[n]))
