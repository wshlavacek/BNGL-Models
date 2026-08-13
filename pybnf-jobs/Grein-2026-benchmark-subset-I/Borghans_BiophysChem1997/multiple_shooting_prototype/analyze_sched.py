"""Segment-count comparison, paired against the single-shooting results already in
hand for the same starts: python analyze_sched.py sched.jsonl sweep.jsonl"""
import collections
import json
import sys

import numpy as np

import msproto as M

SOLVED = M.JSTAR + M.SOLVED_OG - M.PAPER_OFFSET

sched = [json.loads(l) for l in open(sys.argv[1]) if l.strip()]
sweep = [json.loads(l) for l in open(sys.argv[2]) if l.strip()]

radius = sched[0]['radius']
ss = {r['seed']: r for r in sweep if r['method'] == 'ss' and r['radius'] == radius}

by = collections.defaultdict(dict)
for r in sched:
    by[r['schedule']][r['seed']] = r

names = ['2-1', '3-1', '4-2-1', '8-4-2-1']
seeds = sorted(set.intersection(*[set(by[n]) for n in names if by[n]]) & set(ss))

print('radius %.2f, %d paired starts.  Lower is better; solved needs %.1f'
      % (radius, len(seeds), SOLVED))
print()
hdr = '%5s %10s' % ('seed', 'SS')
for n in names:
    hdr += ' %10s' % ('MS ' + n)
print(hdr)
print('-' * len(hdr))
for s in seeds:
    row = '%5d %10.2f' % (s, ss[s]['certified'])
    for n in names:
        row += ' %10.2f' % by[n][s]['certified']
    print(row)
print('-' * len(hdr))

row = '%5s %10.2f' % ('med', np.median([ss[s]['certified'] for s in seeds]))
for n in names:
    row += ' %10.2f' % np.median([by[n][s]['certified'] for s in seeds])
print(row)
row = '%5s %10.2f' % ('best', min(ss[s]['certified'] for s in seeds))
for n in names:
    row += ' %10.2f' % min(by[n][s]['certified'] for s in seeds)
print(row)
row = '%5s %10.0f' % ('wins', sum(1 for s in seeds
                                  if all(ss[s]['certified'] <= by[n][s]['certified']
                                         for n in names)))
for n in names:
    row += ' %10.0f' % sum(1 for s in seeds if by[n][s]['certified'] < ss[s]['certified'])
print(row + '   (SS column: starts where SS beat every MS schedule)')
row = '%5s %10.0f' % ('sims', np.median([ss[s]['n_sim'] for s in seeds]))
for n in names:
    row += ' %10.0f' % np.median([by[n][s]['n_sim'] for s in seeds])
print(row)
