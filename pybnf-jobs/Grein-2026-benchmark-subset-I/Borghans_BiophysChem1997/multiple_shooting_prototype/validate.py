"""Validate the prototype against PyBNF's own numbers before trusting any result.

  V1  the objective reproduces nominal_check.json's reduced objective / J_paper / OG
  V2  a multi-segment transcription seeded by single-shooting theta has c == 0 and
      the SAME data residuals as single shooting (the transcription is exact)
  V3  the analytic Jacobian of (r, c) matches central finite differences
"""
import time

import numpy as np

import msproto as M

np.set_printoptions(precision=6, suppress=False, linewidth=130)

print('=' * 78)
print('V1  objective calibration at the PEtab nominal point')
print('=' * 78)
theta = M.nominal_theta()
tr1 = M.Transcription(1)
t0 = time.time()
out = tr1.evaluate(theta, want_jac=False)
print('    one single-shoot evaluation (no jac): %.3f s' % (time.time() - t0))
S = float(np.sum(out['r'] ** 2))
sigma_nom = 10 ** -0.9791296359262935
J_fixed = M.N_OBS * np.log(sigma_nom) + S / (2 * sigma_nom ** 2)
print('    reduced objective at nominal sigma : %.6f   (expected -198.101699)' % J_fixed)
print('    delta                              : %.3e' % abs(J_fixed - (-198.1016990896733)))
J_prof = tr1.reduced_objective(out['r'])
print('    reduced objective, sigma profiled  : %.6f' % J_prof)
print('    profiled sigma                     : %.6f  (nominal %.6f)'
      % (tr1.sigma_profiled(out['r']), sigma_nom))
print('    J_paper at nominal sigma           : %.6f   (expected  -83.323678)'
      % M.paper_nll(J_fixed))
print('    OG at nominal sigma                : %.6f   (expected   48.684799)'
      % M.optimality_gap(J_fixed))
print()
print('    solved threshold, reduced scale    : %.4f' % (M.JSTAR + M.SOLVED_OG - M.PAPER_OFFSET))

print()
print('=' * 78)
print('V2  the m-segment transcription is exact at a continuous point')
print('=' * 78)
for m in (2, 4, 8, 16):
    x = M.seed_aux(theta, m)
    tr = M.Transcription(m)
    o = tr.evaluate(x, want_jac=False)
    dr = np.abs(o['r'] - out['r']).max()
    print('    m=%-3d  ||c||_inf = %.3e   max |r_m - r_1| = %.3e   J = %.6f'
          % (m, np.abs(o['c']).max() if len(o['c']) else 0.0, dr,
             tr.reduced_objective(o['r'])))

print()
print('=' * 78)
print('V3  analytic Jacobian vs central finite differences  (m = 4)')
print('=' * 78)
m = 4
x = M.seed_aux(theta, m)
# perturb off the continuous point so the continuity block is not identically zero
rng = np.random.default_rng(0)
x = x + rng.normal(scale=0.02, size=x.shape)
tr = M.Transcription(m)
t0 = time.time()
o = tr.evaluate(x, want_jac=True)
print('    one %d-segment evaluation WITH jac: %.3f s' % (m, time.time() - t0))
Jr, Jc = o['Jr'], o['Jc']

fd_r = np.zeros_like(Jr)
fd_c = np.zeros_like(Jc)
for k in range(tr.L.n):
    h = 1e-6 * max(abs(x[k]), 1.0)
    xp, xm = x.copy(), x.copy()
    xp[k] += h
    xm[k] -= h
    op = tr.evaluate(xp, want_jac=False)
    om = tr.evaluate(xm, want_jac=False)
    fd_r[:, k] = (op['r'] - om['r']) / (2 * h)
    fd_c[:, k] = (op['c'] - om['c']) / (2 * h)

for label, A, B in (('data Jacobian  dr/dx', Jr, fd_r),
                    ('continuity Jac dc/dx', Jc, fd_c)):
    denom = np.maximum(np.abs(A), 1.0)
    err = np.abs(A - B) / denom
    k = int(np.unravel_index(np.argmax(err), err.shape)[1])
    print('    %s: worst relative disagreement %.3e  (column %d = %s)'
          % (label, err.max(), k,
             (tr.L.theta_names + ['z%d_%s' % (1 + i // 3, M.STATES[i % 3])
                                  for i in range(tr.L.n_aux)])[k]))
