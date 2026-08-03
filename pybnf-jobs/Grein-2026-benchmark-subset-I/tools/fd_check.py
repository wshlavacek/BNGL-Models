#!/usr/bin/env python
"""Finite-difference check of a real PyBNF job's assembled gradient.

Mirrors GradientOptimizer.gradient_at: build one ExperimentRouting per scored
(model, suffix), apply the union sensitivity request, simulate once, assemble the
gradient in sampling space u -- then compare against central differences of the same
objective in the same space. Every free parameter is reported, so a column that is
structurally zero stands out from one that is merely small.

This found lanl/PyBNF#534 the first time it was pointed at a real model: a free
parameter bound by id that reached the trajectory only by seeding another entity
assembled to exactly 0 while central differences said -10.4.

Five things to know before trusting a red row:

* **Step size.** On a stiff model 1e-7 is pure roundoff and 1e-6 is marginal; 3e-4
  (the default here) gave the best agreement on Laske. Sweep `-h` before concluding
  anything -- a REAL defect does not move with the step size, FD noise does. Sweep
  UPWARD too: the small-h end is the roundoff end, so a column drifting away from the
  assembled value as h shrinks is converging at the other end.
* **The default point.** With no param-values.json this evaluates at each parameter's
  median-quantile value, which may be a silly point for the model. Pass the PEtab
  nominal values (see tools/README.md) for a meaningful comparison.
* **A gradient check is only a test where the gradient is large.** Most slugs here have
  `OG_nominal` ~ 0 -- their PEtab nominal point IS the optimum, where the true gradient
  vanishes and both sides are noise. Every such slug flags red at the nominal point and
  every one of those flags is an artefact. Use `--disp F` to re-evaluate F of each box
  width away, which is what separated Fiedler's real defect from eleven artefacts
  (lanl/PyBNF#535).
* **Bounds.** `FreeParameter.set_value` clamps an out-of-box value, so a parameter on a
  bound has f(u+h) == f(u) and its central difference is a half-step -- or exactly zero
  when both sides clamp, which looks like a dead column. The evaluation point is kept
  `8*h` clear of both bounds and any parameter that had to be moved is reported.
* **Both sides must model the same objective.** The two ways this script used to compare
  unlike quantities are fixed below and described in tools/README.md: each experiment now
  travels with its `data_key` (without it `objective._scale_factors` keeps the last-scored
  experiment's series factors, silently), and the assembled side now adds the constraint
  penalty gradient that `loss_at` has always included. Both found while auditing this
  script for lanl/PyBNF#537.

Usage:  fd_check.py <slug-dir | job.conf> [param-values.json] [-h STEP] [--disp F]
"""
import json
import os
import sys
import tempfile

import numpy as np
from pybnf.algorithms import core
from pybnf.gradient import (
    apply_routings, assemble_constraint_gradient, assemble_gaussian_gradient, route_for_model,
)
from pybnf.parse import load_config
from pybnf.pset import PSet


def _condition_for_suffix(model, suffix):
    best = None
    for mut in getattr(model, 'mutants', []) or []:
        ms = getattr(mut, 'suffix', '')
        if ms and suffix.endswith(ms) and (best is None or len(ms) > len(best.suffix)):
            best = mut
    return best


def simulate(config, pset, tag, routings=None):
    models = list(config.models.values())
    if routings is not None:
        for m in models:
            m.set_scored_suffixes(config.exp_data.get(m.name, {}))
            apply_routings(m, routings[m.name])
    with tempfile.TemporaryDirectory() as sim_dir:
        job = core.Job(models, pset, tag, sim_dir, config.config['wall_time_sim'], None,
                       config.config['normalization'], config.postprocessing, True,
                       stochastic_seed_policy=config.config['stochastic_seed'])
        res = core.run_job(job)
    if getattr(res, 'failed', False) or res.simdata is None:
        raise SystemExit('simulation failed at this point')
    res.normalize(config.config['normalization'])
    res.postprocess_data(config.postprocessing)
    return res


def main():
    target = os.path.abspath(sys.argv[1])
    has_values = len(sys.argv) > 2 and not sys.argv[2].startswith('-')
    values = json.load(open(sys.argv[2])) if has_values else None
    step = float(sys.argv[sys.argv.index('-h') + 1]) if '-h' in sys.argv else 3e-4
    disp = float(sys.argv[sys.argv.index('--disp') + 1]) if '--disp' in sys.argv else 0.0
    if os.path.isdir(target):
        slug_dir = target
        confs = [os.path.join(slug_dir, f) for f in sorted(os.listdir(slug_dir))
                 if f.endswith('.conf')]
        if not confs:
            raise SystemExit(f'no .conf in {slug_dir}')
        conf = confs[0]
    else:
        conf, slug_dir = target, os.path.dirname(target)

    os.chdir(slug_dir)
    config = load_config(conf)
    free = list(config.variables)
    if values:
        free = [v.set_value(values[v.name]) for v in free]
    else:
        free = [v.set_value(v.initial_value_from_quantile(0.5)
                            if hasattr(v, 'initial_value_from_quantile') else v.value)
                for v in free]
    # Move the evaluation point off the optimum if asked, and keep every parameter clear of its
    # box either way: set_value CLAMPS out of bounds, so a parameter sitting ON a bound has
    # f(u+h) == f(u) and reads as a dead column that is nothing of the sort.
    lo = np.array([v.to_sampling_space(v.lower_bound) for v in free])
    hi = np.array([v.to_sampling_space(v.upper_bound) for v in free])
    u = np.array([v.to_sampling_space(v.value) for v in free])
    if disp:
        u = u + disp * (hi - lo) * np.array([1.0 if j % 2 == 0 else -1.0 for j in range(len(free))])
    u = np.clip(u, lo + 8.0 * step, hi - 8.0 * step)
    moved = [v.name for v, a, b in zip(free, u, [v.to_sampling_space(v.value) for v in free])
             if a != b]
    if moved:
        print('moved off the optimum / away from a bound:', moved)
    free = [v.set_value(v.from_sampling_space(x)) for v, x in zip(free, u)]
    names = [v.name for v in free]
    pset = PSet(free)

    # Routings, exactly as _setup_gradient_path builds them.
    routings, per_model = {}, {}
    for m in config.models.values():
        wt = route_for_model(m, names, None)
        mine = {}
        for suffix in config.exp_data.get(m.name, {}):
            mine[suffix] = route_for_model(m, names, _condition_for_suffix(m, suffix))
        per_model[m.name] = [wt, *mine.values()]
        for suffix, r in mine.items():
            routings[(m.name, suffix)] = r
    print('point-dependent routings:',
          sum(1 for r in routings.values() if r.is_point_dependent), '/', len(routings))

    def loss_at(u_vec):
        pv = [v.from_sampling_space(u) for v, u in zip(free, u_vec, strict=False)]
        ps = PSet([v.set_value(x) for v, x in zip(free, pv, strict=False)])
        res = simulate(config, ps, 'fd')
        return config.obj.evaluate_multiple(res.simdata, config.exp_data, ps, config.constraints)

    u0 = np.array([v.to_sampling_space(v.value) for v in free])
    grad_fd = np.zeros(len(free))
    for j in range(len(free)):
        up, dn = u0.copy(), u0.copy()
        up[j] += step
        dn[j] -= step
        grad_fd[j] = (loss_at(up) - loss_at(dn)) / (2.0 * step)

    res = simulate(config, pset, 'grad', routings=per_model)
    # Score this Result first, exactly as the fit does before gradient_at: that binds the
    # objective's _pset_values and materializes any measurement-model observable column onto
    # the Data. Assembling without it silently yields zero columns for every free parameter
    # that lives in an observable formula or a noise scale.
    config.obj.evaluate_multiple(res.simdata, config.exp_data, pset, config.constraints)
    values_at = {p.name: p.value for p in pset}
    resolved = {key: r.at_point(values_at) for key, r in routings.items()}
    experiments = []
    for mname, by_suffix in res.simdata.items():
        me = config.exp_data.get(mname, {})
        for suffix, sim in by_suffix.items():
            if suffix in me:
                # The suffix travels as the experiment's data_key, exactly as gradient_at passes
                # it. Omitting it does not merely refuse an analytic 'scale' column -- it also
                # leaves objective._scale_factors pointing at whichever experiment
                # evaluate_multiple scored LAST, so every experiment's residual is scored with
                # the wrong series factor. Silent, and a ratio of two profiled c* reads as a
                # clean factor (lanl/PyBNF#537).
                experiments.append((sim, me[suffix], resolved[(mname, suffix)], suffix))
    g = assemble_gaussian_gradient(config.obj, experiments, free)
    if config.constraints:
        # loss_at scores the constraint penalties (evaluate_multiple takes config.constraints),
        # so the assembled side must carry their gradient too or every constraint-touched column
        # reads red as an artefact. This is the second half of gradient_at (lanl/PyBNF#537).
        g.gradient = g.gradient + assemble_constraint_gradient(
            config.constraints, res.simdata, resolved, free)

    print(f'\n{"param":26s} {"assembled":>15s} {"central diff":>15s} {"rel err":>10s}')
    scale = max(np.max(np.abs(grad_fd)), 1e-30)
    worst = 0.0
    for n, a, f in zip(g.param_names, g.gradient, grad_fd, strict=False):
        rel = abs(a - f) / max(abs(f), scale * 1e-3)
        worst = max(worst, rel)
        flag = '  <-- ' if rel > 1e-2 else ''
        print(f'{n:26s} {a:15.6g} {f:15.6g} {rel:10.2e}{flag}')
    print(f'\nworst relative error: {worst:.2e}')


if __name__ == '__main__':
    main()
