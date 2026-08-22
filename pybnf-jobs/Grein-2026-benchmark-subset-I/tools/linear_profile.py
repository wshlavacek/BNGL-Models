#!/usr/bin/env python
"""What does the objective landscape look like once the observable-layer linear parameters
are profiled out instead of searched?

The evaluation lanl/PyBNF#572 asks for, before anything is built. #572 proposes removing an
observable's `scale`/`offset` from the search and replacing them, at each theta, with their
conditional optimum. It states a real counter-hypothesis: an offset is a cheap descent
direction, `Borghans`'s known no-dynamics attractor is reachable *purely* through the linear
parameters (`a = 0, b = mean(d)`), and profiling hands that solution to every theta for free
-- which could **compress** the ranking a global sampler sees rather than sharpen it.

This measures it. For each of N points drawn from the conf's own sampling distribution:

    searched   PyBNF's objective at the drawn linear-parameter values -- what a sampler ranks today
    profiled   the same objective minimized over the linear parameters, theta held fixed

The profile is taken **numerically, through PyBNF's own `evaluate_multiple`**, not through a
re-derived closed form. That is deliberate and it is the point:

* #572's closed form (variable projection over `Phi = [s, 1]`) is a **linear least-squares**
  identity, and four of the six slugs it names score their residual on a **log** scale, where
  it does not hold. A numeric profile is the exact conditional optimum whatever the family is,
  so the landscape question can be answered separately from the "is there a closed form"
  question. Run `linear_scope.py` for the second one.
* Re-deriving the algebra and measuring *it* would confirm the algebra, not the landscape.
  ADR-0108 pinned its sigma closed form against a numeric minimization of PyBNF's reported
  objective for the same reason.

Usage:
    linear_profile.py <slug-dir> [--params scale,offset] [-n N] [--seed S]
                      [--noise-profiling] [--linear-space] [--closed-form] [--maxiter M]
                      [--point LABEL=FILE.json ...] [--json out.json]

`--params` defaults to every parameter `linear_scope.py` classifies as affine in the
prediction. `--noise-profiling` adds `noise_profiling = 1` to the conf, so sigma is profiled
by PyBNF's own ADR-0108 closed form and the numeric search covers only the linear parameters
-- that is #572's evaluation item 2, "scale, offset and sigma all profiled". `--point` scores
a named parameter vector before the box draws: the reference optimum, the PEtab nominal point,
a known oscillating start, or a ladder of interpolated points to walk a section.

`--closed-form` additionally computes #572's own variable-projection solution and prints its gap
to the numeric optimum. It applies only to a linear-scale Gaussian objective, and it must be run
with sigma **searched** rather than profiled, so that `W` is the same matrix on both sides. Where
both can reach the same point the two agree to ~1e-14; where they differ, the projection has left
the parameter's declared support -- see `_varpro` and
`pybnf-jobs/Synthetic-2026-linear-observable/README.md`.

The **flat** column is the objective of a *constant* prediction at its best constant -- the
no-dynamics score, and the number the whole counter-hypothesis is about. Read it carefully,
because the obvious reading is wrong on the slugs that matter:

* On a **linear-scale** family the constant column is in the span at every theta, so the
  profiled objective can never be worse than the flat line -- it is an upper bound.
* On a **log** family it is neither a bound nor a floor. `log(a*s + b)` is not a projection,
  and a wrong-shaped `s` makes the fit *worse* than no `s` at all, so the conditional optimum
  is frequently `a = 0` exactly. On `Borghans` that happens at 39 of 76 box draws: the flat
  line is not what the reduced surface is bounded by, it is what the reduced surface **is**.
  The `FLAT` marker in each row is that test.

Gotchas, each of which produced a plausible wrong answer first:

* **A Result can only be scored once.** The measurement layer materializes its observable
  columns onto the Data *in place* (ADR-0036), and `_add_column` raises rather than
  overwriting, so the second scoring of the same Result dies. The inner minimization scores
  one Result hundreds of times, so every scoring gets a `deepcopy` of the simulated data.
  Skip that and the tool cannot profile at all -- it fails loudly, which is the good case.

* **`FreeParameter.set_value` clamps to the box.** Building a `PSet` to carry a trial value
  silently returns the bound instead whenever the profile wants to leave the box -- and
  leaving the box is the entire premise of profiling. The inner loop therefore hands
  `evaluate_multiple` a duck-typed `(name, value)` sequence, which is all it reads
  (`{p.name: p.value for p in pset}`), and never a PSet.

* **Profile in the parameter's own sampling space.** A `loguniform_var` offset is declared
  positive; its unconstrained least-squares optimum need not be. Optimizing `log10(value)`
  keeps the profile on the support the fit declared, so `searched` and `profiled` are the
  same quantity with and without a box rather than two different models. `--linear-space`
  drops that and profiles on the real line, which is what VarPro would do -- the difference
  is worth knowing and is not assumed away here.

* **Normalize before the deepcopy, not after.** `res.normalize` transforms the raw simulation
  columns; the measurement layer runs later, at scoring time. Normalizing per inner
  evaluation would re-normalize an already-normalized column.
"""
import argparse
import copy
import json
import logging
import os
import sys
import tempfile
import time

import numpy as np

logging.basicConfig(level=logging.ERROR)


class _Var:
    """A duck-typed pset entry: `evaluate_multiple` reads only `.name` / `.value`.

    Deliberately not a `FreeParameter`, whose `set_value` clamps to the declared box -- see
    the module docstring. A profiled coefficient has no box, and a clamped trial value would
    make the measured "profiled" landscape a bounded-search landscape wearing its name.
    """
    __slots__ = ('name', 'value')

    def __init__(self, name, value):
        self.name = name
        self.value = value


def _load(conf_path, extra_lines=()):
    """Load a conf, optionally with extra keys appended.

    A key is read when `load_config` BUILDS the objective, so an override has to be in the
    conf text -- assigning `config.config[...]` afterwards is a silent no-op (the same trap
    `box_probe.py --atol` documents).
    """
    from pybnf.parse import load_config
    if not extra_lines:
        return load_config(conf_path)
    with open(conf_path) as fh:
        text = fh.read()
    drop = {ln.split('=')[0].strip() for ln in extra_lines}
    keep = [ln for ln in text.splitlines() if ln.split('=')[0].strip() not in drop]
    tmp = conf_path + '.linear_profile.conf'
    with open(tmp, 'w') as fh:
        fh.write('\n'.join(keep) + '\n' + '\n'.join(extra_lines) + '\n')
    try:
        return load_config(tmp)
    finally:
        os.unlink(tmp)


def _simulate(config, values, tag):
    """One fresh, normalized `simdata` at the given `{name: value}` point, or `None`."""
    from pybnf.algorithms import core
    from pybnf.pset import PSet
    pset = PSet([v.set_value(values[v.name]) for v in config.variables])
    with tempfile.TemporaryDirectory() as sim_dir:
        job = core.Job(list(config.models.values()), pset, tag, sim_dir,
                       config.config['wall_time_sim'], None,
                       config.config['normalization'], config.postprocessing, True,
                       stochastic_seed_policy=config.config['stochastic_seed'])
        res = core.run_job(job)
    if getattr(res, 'failed', False) or res.simdata is None:
        return None
    res.normalize(config.config['normalization'])
    res.postprocess_data(config.postprocessing)
    return res.simdata


class Scorer:
    """Scores one fixed simulation repeatedly, at varying linear-parameter values."""

    def __init__(self, config, simdata, base_values):
        self.config = config
        self.simdata = simdata
        self.base = dict(base_values)
        self.calls = 0

    def score(self, overrides):
        self.calls += 1
        vals = dict(self.base)
        vals.update(overrides)
        pset = [_Var(n, v) for n, v in vals.items()]
        try:
            out = self.config.obj.evaluate_multiple(
                copy.deepcopy(self.simdata), self.config.exp_data, pset,
                self.config.constraints, show_warnings=False)
        except Exception:
            return np.inf
        if out is None or not np.isfinite(out):
            return np.inf
        return float(out)


def _space(var, linear_space):
    """`(to_u, from_u)` for one free parameter -- its own sampling space, or the real line."""
    if linear_space:
        return (lambda x: x), (lambda u: u)
    return var.to_sampling_space, var.from_sampling_space


def _profile(scorer, prof_vars, start, linear_space, grid=21, maxiter=2000):
    """Minimize the objective over `prof_vars`, everything else fixed. Returns
    `(best_score, {name: value}, n_calls)`.

    For one or two coefficients a coarse grid over the declared box precedes the local
    polish: the profiled surface of a *log*-family offset is not convex in `log10(offset)`,
    and a bare local solver started at the drawn value converges to whichever side of the
    data it started on. On the corpus that difference is tens of objective units, i.e. larger
    than the whole effect being measured.
    """
    from scipy.optimize import minimize

    names = [v.name for v in prof_vars]
    spaces = [_space(v, linear_space) for v in prof_vars]

    def at(u):
        return scorer.score({n: fr(x) for n, (_to, fr), x in zip(names, spaces, u)})

    u0 = [to(start[v.name]) for v, (to, _fr) in zip(prof_vars, spaces)]
    starts = [np.array(u0, dtype=float)]

    if len(prof_vars) <= 2:
        axes = []
        for v, (to, _fr) in zip(prof_vars, spaces):
            lo, hi = to(v.lower_bound), to(v.upper_bound)
            if not np.isfinite(lo) or not np.isfinite(hi):
                lo, hi = u0[len(axes)] - 5.0, u0[len(axes)] + 5.0
            axes.append(np.linspace(lo, hi, grid))
        mesh = np.stack([m.ravel() for m in np.meshgrid(*axes, indexing='ij')], axis=1)
        scores = np.array([at(pt) for pt in mesh])
        order = np.argsort(scores)[:3]
        starts.extend(mesh[i] for i in order if np.isfinite(scores[i]))

    best, best_u = np.inf, np.array(u0, dtype=float)
    for s in starts:
        r = minimize(at, s, method='Nelder-Mead',
                     options={'xatol': 1e-8, 'fatol': 1e-10, 'maxiter': maxiter})
        if np.isfinite(r.fun) and r.fun < best:
            best, best_u = float(r.fun), r.x
    return best, {n: fr(x) for n, (_to, fr), x in zip(names, spaces, best_u)}, scorer.calls


def _flat_reference(config, scorer, prof_vars, linear_space, maxiter=2000):
    """The objective of a prediction held *constant*, at its best constant -- the
    no-dynamics score.

    Realized by driving every non-intercept linear coefficient to (effectively) zero and
    profiling only the intercept, so the prediction column is a constant. Returns `None`
    when the slug has no intercept-like coefficient, in which case a constant prediction is
    not in the model's reach and the whole counter-hypothesis does not apply to it.

    **This number is only a property of the data when sigma is profiled too.** With a searched
    free sigma it is evaluated at whatever sigma the current draw happened to carry, so it
    moves from draw to draw by hundreds of objective units and is not a reference line at all
    -- it is a different number each row wearing the same name. `main` therefore recomputes it
    at every point and reports the spread; a spread that is not ~0 means the run was made
    without `--noise-profiling` and the comparison against it is meaningless.
    """
    intercepts = [v for v in prof_vars if 'offset' in v.name or 'background' in v.name]
    scales = [v for v in prof_vars if v not in intercepts]
    if not intercepts:
        return None, {}
    # Exactly zero, not the declared lower bound. A `loguniform_var` scale bottoms out at
    # 1e-3, and on a trajectory spanning ten decades `1e-3 * s` is not remotely constant --
    # the "flat" line then moves with the draw by hundreds of objective units and reads as
    # sigma contamination. `_Var` does not clamp, so 0.0 is reachable and the prediction
    # column really is the intercept.
    tiny = {v.name: 0.0 for v in scales}
    sub = Scorer(config, scorer.simdata, {**scorer.base, **tiny})
    start = {v.name: scorer.base[v.name] for v in intercepts}
    val, at, _n = _profile(sub, intercepts, start, linear_space, maxiter=maxiter)
    return val, {**tiny, **at}


def _varpro(config, scorer, prof_vars):
    """#572's own closed form, computed and scored: `c* = (Phi^T W Phi)^-1 Phi^T W d`.

    Returns `(score_at_c_star, {name: value}, note)`, or `(None, {}, reason)` when the
    construction does not apply to this slug.

    `Phi` is built by *evaluating the model's own prediction* at basis coefficient vectors
    rather than by parsing the formula: with every profiled coefficient at 0 the aligned
    prediction is `Phi.0`, and setting coefficient `j` to 1 gives `Phi.e_j`, so
    `Phi_j = pred(e_j) - pred(0)`. That is formula-agnostic, it needs no new sensitivity, and
    it doubles as a check -- if `pred(0)` is not ~0 the observable is not linear in these
    coefficients and the whole construction is inapplicable, which is reported rather than
    silently absorbed into the intercept.

    `LikelihoodObjective.aligned_prediction_data` supplies `(prediction, observation,
    variance)` for exactly the scored points, and returns `None` unless **every** scored
    observable is a linear-scale Gaussian -- which is precisely #572's precondition, so the
    inapplicable case reports itself.

    **Run this with sigma SEARCHED, not profiled.** `W = diag(1/sigma^2)` has to be the same
    matrix on both sides of the comparison; under `noise_profiling` sigma moves with every
    evaluation, so the closed form would be solving a different weighted problem than the one
    the numeric profile converges to and the two would disagree for a reason that is not an
    error in either.

    Per-point fit weights are not on this seam. A fixture that declares none (weight 1
    everywhere) is unaffected; a job with weights would need them folded into `W`.
    """
    names = [v.name for v in prof_vars]

    def aligned(overrides):
        vals = dict(scorer.base)
        vals.update(overrides)
        pset = [_Var(n, x) for n, x in vals.items()]
        return config.obj.aligned_prediction_data(
            copy.deepcopy(scorer.simdata), config.exp_data, pset)

    base = aligned({n: 0.0 for n in names})
    if base is None:
        return None, {}, 'not a linear-scale Gaussian at every scored point'
    p0, d, var = base
    if np.max(np.abs(p0)) > 1e-8 * max(1.0, float(np.max(np.abs(d)))):
        return None, {}, 'prediction at c = 0 is not zero; the observable is not linear here'

    cols = []
    for j, n in enumerate(names):
        step = aligned({m: (1.0 if m == n else 0.0) for m in names})
        if step is None:
            return None, {}, 'aligned prediction unavailable at basis vector %d' % j
        cols.append(step[0] - p0)
    phi = np.column_stack(cols)

    w = 1.0 / np.asarray(var, dtype=float)
    if not np.all(np.isfinite(w)):
        return None, {}, 'a scored point has non-finite variance'
    a = phi.T @ (w[:, None] * phi)
    b = phi.T @ (w * d)
    try:                                  # pinv, not solve: Phi is singular wherever s is
        c = np.linalg.solve(a, b)         # constant over the series, which a sampler visits
    except np.linalg.LinAlgError:
        c = np.linalg.pinv(a) @ b
    at = {n: float(c[j]) for j, n in enumerate(names)}
    return scorer.score(at), at, 'ok'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('slug_dir')
    ap.add_argument('--params', default=None,
                    help='comma-separated names to profile (default: linear_scope.py\'s set)')
    ap.add_argument('-n', type=int, default=25)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--noise-profiling', action='store_true',
                    help='add noise_profiling = 1, so sigma is profiled by ADR-0108 too')
    ap.add_argument('--linear-space', action='store_true',
                    help='profile on the real line rather than the parameter\'s sampling space')
    ap.add_argument('--point', action='append', default=[], metavar='LABEL=FILE.json',
                    help='a named {name: value} point to score before the box draws '
                         '(repeatable) -- the reference optimum, the PEtab nominal point, '
                         'a known oscillating start')
    ap.add_argument('--closed-form', action='store_true',
                    help="also compute #572's variable-projection solution and report its gap "
                         "to the numeric optimum; requires a linear-scale Gaussian objective "
                         "and a SEARCHED sigma (see _varpro)")
    ap.add_argument('--maxiter', type=int, default=2000,
                    help='Nelder-Mead iteration cap for the inner solve; lower it only for a\n'
                         'slug whose simulation is slow, and say so with the numbers')
    ap.add_argument('--json', default=None)
    args = ap.parse_args()

    slug_dir = os.path.abspath(args.slug_dir)
    base = os.path.basename(slug_dir)
    conf = os.path.join(slug_dir, base + '.conf')
    if not os.path.exists(conf):
        conf = os.path.join(slug_dir, sorted(f for f in os.listdir(slug_dir)
                                             if f.endswith('.conf'))[0])
    os.chdir(slug_dir)

    extra = ['noise_profiling = 1'] if args.noise_profiling else []
    config = _load(conf, extra)

    if args.params:
        want = [p.strip() for p in args.params.split(',') if p.strip()]
    else:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from linear_scope import scan
        want = [p['name'] for p in scan(conf)['params']]
    by_name = {v.name: v for v in config.variables}
    missing = [w for w in want if w not in by_name]
    if missing:
        raise SystemExit('not searched free parameters in this conf: %s' % ', '.join(missing))
    prof_vars = [by_name[w] for w in want]

    print('slug        : %s' % base)
    print('conf        : %s' % os.path.basename(conf))
    print('objective   : %s' % type(config.obj).__name__)
    print('searched k  : %d   profiled here: %d  (%s)'
          % (len(config.variables), len(prof_vars), ', '.join(want)))
    print('noise_profiling: %s   space: %s'
          % (bool(extra), 'linear' if args.linear_space else 'sampling'))
    if getattr(config, 'profiled_noise_params', None):
        print('profiled sigma : %s' % ', '.join(sorted(config.profiled_noise_params)))
    print()

    rng = np.random.default_rng(args.seed)
    points = []
    for spec in args.point:
        label, _, path = spec.partition('=')
        with open(path) as fh:
            given = json.load(fh)
        missing_v = [v.name for v in config.variables if v.name not in given]
        if missing_v:
            raise SystemExit('point %r has no value for: %s' % (label, ', '.join(missing_v)))
        points.append((label, {v.name: float(given[v.name]) for v in config.variables}))
    for i in range(args.n):
        vals = {}
        for v in config.variables:
            lo, hi = v.to_sampling_space(v.lower_bound), v.to_sampling_space(v.upper_bound)
            vals[v.name] = v.from_sampling_space(rng.uniform(lo, hi))
        points.append(('box%03d' % i, vals))

    rows = []
    flats = []
    print('%-9s %13s %13s %10s %13s %-4s  %s'
          % ('point', 'searched', 'profiled', 'gain', 'flat', '', 'at'))
    t0 = time.time()
    for label, vals in points:
        simdata = _simulate(config, vals, label)
        if simdata is None:
            print('%-9s %13s' % (label, 'SIM FAIL'))
            rows.append({'label': label, 'failed': True})
            continue
        scorer = Scorer(config, simdata, vals)
        searched = scorer.score({})
        flat, _flat_at = _flat_reference(config, scorer, prof_vars, args.linear_space,
                                         maxiter=args.maxiter)
        if flat is not None and np.isfinite(flat):
            flats.append(flat)
        profiled, at, ncalls = _profile(scorer, prof_vars, vals, args.linear_space,
                                        maxiter=args.maxiter)
        gain = searched - profiled
        # Did the inner solve just choose the no-dynamics answer? That is the counter-
        # hypothesis made measurable: not "the flat line is available" but "the flat line
        # IS the conditional optimum here", which makes the reduced objective constant.
        collapsed = (flat is not None and np.isfinite(flat)
                     and abs(profiled - flat) < 1e-6)
        varpro = varpro_at = varpro_note = None
        if args.closed_form:
            varpro, varpro_at, varpro_note = _varpro(config, scorer, prof_vars)
        print('%-9s %13.6f %13.6f %10.4g %13s %-4s  %s'
              % (label, searched, profiled, gain,
                 '-' if flat is None else '%.6f' % flat,
                 'FLAT' if collapsed else '',
                 ' '.join('%s=%.4g' % kv for kv in sorted(at.items()))))
        if args.closed_form:
            print('%-9s %13s %13s %10s %13s %-4s  %s'
                  % ('  varpro', '',
                     '-' if varpro is None else '%.6f' % varpro,
                     '' if varpro is None else '%+.3g' % (varpro - profiled),
                     '', '',
                     varpro_note if varpro is None
                     else ' '.join('%s=%.4g' % kv for kv in sorted(varpro_at.items()))))
        sys.stdout.flush()
        rows.append({'label': label, 'failed': False, 'searched': searched,
                     'profiled': profiled, 'gain': gain, 'at': at, 'flat': flat,
                     'collapsed': bool(collapsed), 'inner_calls': ncalls,
                     'varpro': varpro, 'varpro_at': varpro_at, 'varpro_note': varpro_note,
                     'theta': {k: v for k, v in vals.items() if k not in want}})

    ok = [r for r in rows if not r['failed'] and np.isfinite(r['searched'])
          and np.isfinite(r['profiled'])]
    print('\n%d/%d points scored in %.1f s' % (len(ok), len(rows), time.time() - t0))
    flat_ref = float(np.median(flats)) if flats else None
    if flats:
        spread = max(flats) - min(flats)
        print('flat-line reference recomputed at every point: median %.6f, spread %.3g'
              % (flat_ref, spread))
        if spread > 1e-3:
            print('  ^ NOT a constant. The flat line is being evaluated at each draw\'s own '
                  'sigma, so it is not a reference line. Re-run with --noise-profiling.')
    summary = {}
    if ok and flat_ref is not None:
        s = np.array([r['searched'] for r in ok])
        p = np.array([r['profiled'] for r in ok])
        summary = {
            'flat': flat_ref,
            'flat_spread': float(max(flats) - min(flats)),
            'n': len(ok),
            'searched_better_than_flat': int((s < flat_ref).sum()),
            'profiled_better_than_flat': int((p < flat_ref).sum()),
            'searched_spread': float(s.max() - s.min()),
            'profiled_spread': float(p.max() - p.min()),
            'searched_median_gap_to_flat': float(np.median(s) - flat_ref),
            'profiled_median_gap_to_flat': float(np.median(p) - flat_ref),
            'rank_corr': float(np.corrcoef(np.argsort(np.argsort(s)),
                                           np.argsort(np.argsort(p)))[0, 1]),
            'collapsed_to_flat': int(sum(bool(r.get('collapsed')) for r in ok)),
        }
        print('draws beating flat, searched  : %d / %d' % (summary['searched_better_than_flat'], len(ok)))
        print('draws beating flat, profiled  : %d / %d' % (summary['profiled_better_than_flat'], len(ok)))
        print('spread (max-min)  searched    : %.4f' % summary['searched_spread'])
        print('spread (max-min)  profiled    : %.4f' % summary['profiled_spread'])
        print('Spearman rank corr searched~profiled: %.4f' % summary['rank_corr'])
        print('draws whose profile IS the flat line : %d / %d'
              % (summary['collapsed_to_flat'], len(ok)))

    if args.json:
        with open(args.json, 'w') as fh:
            json.dump({'slug': base, 'conf': os.path.basename(conf), 'params': want,
                       'noise_profiling': bool(extra),
                       'linear_space': bool(args.linear_space),
                       'seed': args.seed, 'rows': rows, 'summary': summary},
                      fh, indent=2, sort_keys=True, default=float)
        print('wrote', args.json)


if __name__ == '__main__':
    main()
