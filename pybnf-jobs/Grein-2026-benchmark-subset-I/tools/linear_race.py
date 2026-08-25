#!/usr/bin/env python
"""Head-to-head at a matched simulation budget: does removing the linear coefficients help a FIT?

`linear_profile.py` measures the landscape. That is a proxy, and lanl/PyBNF#572's item 3 asks for
the thing it is a proxy for: whether a global optimizer, given the same number of simulations, ends
up somewhere better when the scale and offset are solved for at every point instead of searched.

Both sides run the same optimizer -- a small rand/1/bin differential evolution written out here so
the budget is exact and identical -- from the same seeds, on the same conf. The only difference is
the one under test:

    searched   the optimizer carries every free parameter, the linear ones included
    profiled   the optimizer carries the rest, and each candidate is scored with the objective
               already minimized over the linear ones

Both sides report a value of the SAME objective, so the two numbers compare directly: the profiled
side's answer is that objective with its linear coefficients at their conditional optimum.

Usage:
    linear_race.py <slug-dir> [--params a,b] [--budgets 60,120,250] [--seeds 6]
                              [--truth FILE.json] [--over a,b] [--swap a,b]
                              [--noise-profiling] [--closed-form] [--popsize 12]

Gotchas:

* **Match simulations, not evaluations and not wall clock.** #572's whole proposal is to replace
  the inner minimization with one small linear solve, so charging its cost to the profiled side
  prices a tool artefact rather than the feature. `--closed-form` makes the artefact go away
  entirely by using `linear_profile.py::_varpro` as the inner solve, which is what the feature
  would compute. Wall clock is reported either way so the size of the artefact stays visible.

* **`--closed-form` with `--noise-profiling` is only valid where one sigma covers every scored
  point.** `_varpro` weights by `1/sigma^2`, and a profiled sigma moves with every evaluation; with
  a single sigma that is a scalar multiple of the whole weight matrix and the argmin is unchanged,
  but with two estimated sigmas on different observables it is not, and the closed form would be
  solving a different weighted problem than the search converges to.

* **A budget spent is not a budget used.** A point that fails to integrate still costs a simulation,
  and the two sides do not fail at the same rate, because the profiled side draws from a smaller
  box. Both counts are reported.

* **An easy problem cannot discriminate.** Sweep `--budgets` downward until the two sides separate;
  that is the regime a real model with a seventy-second simulation lives in. At a large enough
  budget both sides land on the optimum and the race says nothing.

* **This is a proxy for a real fit, not a real fit.** PyBNF's own optimizers have restarts,
  refinement and convergence tests this does not. The comparison is fair between the two sides and
  is not a prediction of what `job_type = de` would do with either.
"""
import argparse
import json
import logging
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from linear_profile import (Scorer, _distances, _load, _profile,  # noqa: E402
                            _resolve, _simulate, _varpro)


class Budget:
    """A simulation counter. The cap becomes an exhausted flag rather than an exception."""

    def __init__(self, cap):
        self.cap, self.used, self.failed = cap, 0, 0

    def spend(self):
        self.used += 1
        return self.used <= self.cap

    @property
    def exhausted(self):
        return self.used >= self.cap


def evaluator(config, searched, profiled, budget, tag, closed_form, maxiter, grid):
    """`f(u) -> objective` over the searched parameters' sampling-space vector `u`."""
    counter = [0]

    def evaluate(u):
        if not budget.spend():
            return np.inf
        counter[0] += 1
        values = {v.name: v.from_sampling_space(x) for v, x in zip(searched, u)}
        for v in profiled:
            # A starting value for the inner solve, and the value the simulation runs at. The
            # linear coefficients do not enter the model, so which one is used does not change
            # the trajectory; the box midpoint is as good as any.
            values.setdefault(v.name, v.from_sampling_space(
                0.5 * (v.to_sampling_space(v.lower_bound) + v.to_sampling_space(v.upper_bound))))
        for v in config.variables:
            values.setdefault(v.name, v.value)
        simdata = _simulate(config, {v.name: max(min(values[v.name], v.upper_bound),
                                                 v.lower_bound) for v in config.variables},
                            '%s%d' % (tag, counter[0]))
        if simdata is None:
            budget.failed += 1
            return np.inf
        scorer = Scorer(config, simdata, values)
        if not profiled:
            return scorer.score({})
        if closed_form:
            value, _at, note = _varpro(config, scorer, profiled)
            if value is not None and np.isfinite(value):
                return float(value)
            raise SystemExit('--closed-form is not applicable to this conf: %s' % note)
        value, _at, _n = _profile(scorer, profiled, values, False, grid=grid, maxiter=maxiter)
        return value if np.isfinite(value) else np.inf

    return evaluate


def differential_evolution(evaluate, lo, hi, budget, seed, popsize=12, f=0.7, cr=0.9):
    """rand/1/bin over a box in sampling space, stopping when `budget` runs out."""
    rng = np.random.default_rng(seed)
    n = max(4, popsize)
    pop = rng.uniform(lo, hi, size=(n, len(lo)))
    fit = np.array([evaluate(u) for u in pop])
    while not budget.exhausted:
        for i in range(n):
            if budget.exhausted:
                break
            a, b, c = rng.choice([j for j in range(n) if j != i], 3, replace=False)
            trial = np.clip(pop[a] + f * (pop[b] - pop[c]), lo, hi)
            mask = rng.random(len(lo)) < cr
            mask[rng.integers(len(lo))] = True
            candidate = np.where(mask, trial, pop[i])
            value = evaluate(candidate)
            if value <= fit[i]:
                pop[i], fit[i] = candidate, value
    best = int(np.argmin(fit))
    return float(fit[best]), pop[best]


def run_side(config, profiled_names, size, seed, popsize, tag, closed_form, maxiter, grid):
    searched = [v for v in config.variables if v.name not in profiled_names]
    profiled = [v for v in config.variables if v.name in profiled_names]
    lo = np.array([v.to_sampling_space(v.lower_bound) for v in searched])
    hi = np.array([v.to_sampling_space(v.upper_bound) for v in searched])
    budget = Budget(size)
    evaluate = evaluator(config, searched, profiled, budget, tag, closed_form, maxiter, grid)
    started = time.time()
    value, u = differential_evolution(evaluate, lo, hi, budget, seed, popsize=popsize)
    theta = {v.name: float(v.from_sampling_space(x)) for v, x in zip(searched, u)}
    return dict(objective=value, theta=theta, simulations=budget.used, failed=budget.failed,
                seconds=time.time() - started)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('slug_dir')
    ap.add_argument('--params', default=None, help='comma-separated names to profile')
    ap.add_argument('--budgets', default='60,120,250,500,1000',
                    help='comma-separated simulation budgets, applied to both sides')
    ap.add_argument('--seeds', type=int, default=6)
    ap.add_argument('--popsize', type=int, default=12)
    ap.add_argument('--grid', type=int, default=15, help='inner prescan grid, without --closed-form')
    ap.add_argument('--maxiter', type=int, default=200, help='inner Nelder-Mead cap')
    ap.add_argument('--closed-form', action='store_true',
                    help="use #572's variable projection as the inner solve (exact and free)")
    ap.add_argument('--noise-profiling', action='store_true')
    ap.add_argument('--truth', default=None, metavar='FILE.json')
    ap.add_argument('--over', default=None)
    ap.add_argument('--swap', action='append', default=[], metavar='A,B')
    ap.add_argument('--conf', default=None,
                    help='which .conf in the slug directory to use (default: the one named '
                         'after the directory, else the first alphabetically)')
    ap.add_argument('--json', default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.ERROR)

    slug_dir = os.path.abspath(args.slug_dir)
    base = os.path.basename(slug_dir)
    conf = os.path.join(slug_dir, args.conf or (base + '.conf'))
    if not os.path.exists(conf):
        conf = os.path.join(slug_dir, sorted(f for f in os.listdir(slug_dir)
                                             if f.endswith('.conf'))[0])
    truth = json.load(open(_resolve(args.truth, slug_dir))) if args.truth else None
    os.chdir(slug_dir)
    config = _load(conf, ['noise_profiling = 1'] if args.noise_profiling else [])

    by_name = {v.name: v for v in config.variables}
    if args.params:
        want = [p.strip() for p in args.params.split(',') if p.strip()]
    else:
        from linear_scope import scan
        want = [p['name'] for p in scan(conf)['params']]
    want = [w for w in want if w in by_name]
    if not want:
        raise SystemExit('nothing to profile in this conf')
    over = ([n.strip() for n in args.over.split(',')] if args.over
            else [v.name for v in config.variables if v.name not in want])
    swaps = [[n.strip() for n in s.split(',')] for s in args.swap]
    budgets = [int(b) for b in args.budgets.split(',')]

    print('slug        : %s' % base)
    print('profiling   : %s   (%d of %d free parameters)'
          % (', '.join(want), len(want), len(config.variables)))
    print('budgets     : %s simulations, %d seed(s), popsize %d'
          % (', '.join(str(b) for b in budgets), args.seeds, args.popsize))
    print('inner solve : %s' % ('#572 variable projection (exact)' if args.closed_form
                                else 'numerical, grid %d then Nelder-Mead' % args.grid))
    print('noise_profiling: %s' % bool(args.noise_profiling))
    if truth:
        print('distance    : over %s%s' % (', '.join(over),
              '' if not swaps else '  (swaps: %s)' % '; '.join(','.join(g) for g in swaps)))
    print()

    rows = []
    for size in budgets:
        for seed in range(args.seeds):
            for side, profiled_names in (('searched', set()), ('profiled', set(want))):
                out = run_side(config, profiled_names, size, 1000 + seed, args.popsize,
                               '%s%d_%d_' % (side[:4], size, seed), args.closed_form,
                               args.maxiter, args.grid)
                out.update(side=side, seed=seed, budget=size)
                rows.append(out)
                print('  %5d sims  seed %d  %-8s best %13.6f   (%d failed)  %.0fs'
                      % (size, seed, side, out['objective'], out['failed'], out['seconds']))
                sys.stdout.flush()

    distance = (_distances(config, rows, truth, over, swaps) if truth
                else np.full(len(rows), np.nan))
    print()
    print('%-7s %-9s %13s %13s %10s' % ('sims', 'side', 'best', 'median', 'distance'))
    for size in budgets:
        for side in ('searched', 'profiled'):
            index = [i for i, r in enumerate(rows) if r['budget'] == size and r['side'] == side]
            values = np.array([rows[i]['objective'] for i in index])
            print('%-7d %-9s %13.6f %13.6f %10.3f'
                  % (size, side, values.min(), np.median(values),
                     np.nanmedian(distance[index])))
    if args.json:
        for row, d in zip(rows, distance):
            row['distance'] = None if not np.isfinite(d) else float(d)
        with open(args.json, 'w') as fh:
            json.dump({'slug': base, 'params': want, 'budgets': budgets,
                       'noise_profiling': bool(args.noise_profiling),
                       'closed_form': bool(args.closed_form), 'rows': rows},
                      fh, indent=2, sort_keys=True, default=float)
        print('wrote', args.json)


if __name__ == '__main__':
    main()
