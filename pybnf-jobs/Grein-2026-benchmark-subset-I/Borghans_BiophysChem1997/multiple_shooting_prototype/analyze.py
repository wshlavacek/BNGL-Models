"""Summarize a convergence-region sweep: python analyze.py results.jsonl"""
import collections
import json
import sys

import numpy as np

import msproto as M

SOLVED = M.JSTAR + M.SOLVED_OG - M.PAPER_OFFSET
FLAT = -165.98

recs = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
by = collections.defaultdict(list)
for r in recs:
    by[(r['radius'], r['method'])].append(r)

print('solved threshold (reduced scale) %.3f ; flat-line ceiling %.2f ; nominal -198.21'
      % (SOLVED, FLAT))
print()
print('%-7s %-4s %4s %10s %10s %10s %7s %7s %9s' % (
    'radius', 'meth', 'n', 'median J', 'best J', 'median st', 'kept', 'solved', 'med sim'))
print('-' * 84)
for radius in sorted({r['radius'] for r in recs}):
    for meth in ('ss', 'ms'):
        rs = by.get((radius, meth), [])
        if not rs:
            continue
        J = np.array([r['certified'] if r['certified'] is not None else np.inf
                      for r in rs], float)
        st = np.array([r['start_J'] for r in rs], float)
        print('%-7.2f %-4s %4d %10.2f %10.2f %10.2f %6.0f%% %6.0f%% %9.0f' % (
            radius, meth, len(rs), np.median(J), J.min(), np.median(st),
            100 * np.mean([r['kept_dynamics'] for r in rs]),
            100 * np.mean([r['solved'] for r in rs]),
            np.median([r['n_sim'] for r in rs])))
    print()

print('Paired per start (same seed, same radius): does MS beat SS?')
print('%-7s %5s %10s %10s   %s' % ('radius', 'seed', 'SS', 'MS', 'winner'))
print('-' * 60)
wins = collections.Counter()
for radius in sorted({r['radius'] for r in recs}):
    for seed in sorted({r['seed'] for r in recs}):
        s = [r for r in by[(radius, 'ss')] if r['seed'] == seed]
        m = [r for r in by[(radius, 'ms')] if r['seed'] == seed]
        if not (s and m):
            continue
        a = s[0]['certified'] if s[0]['certified'] is not None else np.inf
        b = m[0]['certified'] if m[0]['certified'] is not None else np.inf
        w = 'MS' if b < a - 1e-6 else ('SS' if a < b - 1e-6 else 'tie')
        wins[w] += 1
        print('%-7.2f %5d %10.2f %10.2f   %s' % (radius, seed, a, b, w))
print()
print('paired wins:', dict(wins))
