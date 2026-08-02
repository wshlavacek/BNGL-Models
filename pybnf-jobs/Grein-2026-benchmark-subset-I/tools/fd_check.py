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

Two things to know before trusting a red row:

* **Step size.** On a stiff model 1e-7 is pure roundoff and 1e-6 is marginal; 3e-4
  (the default here) gave the best agreement on Laske. Sweep `-h` before concluding
  anything -- a REAL defect does not move with the step size, FD noise does.
* **The default point.** With no param-values.json this evaluates at each parameter's
  median-quantile value, which may be a silly point for the model. Pass the PEtab
  nominal values (see tools/README.md) for a meaningful comparison.

Usage:  fd_check.py <slug-dir | job.conf> [param-values.json] [-h STEP]
"""
import json
import os
import sys
import tempfile

import numpy as np
from pybnf.algorithms import core
from pybnf.gradient import apply_routings, assemble_gaussian_gradient, route_for_model
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
    experiments = []
    for mname, by_suffix in res.simdata.items():
        me = config.exp_data.get(mname, {})
        for suffix, sim in by_suffix.items():
            if suffix in me:
                experiments.append((sim, me[suffix], routings[(mname, suffix)].at_point(values_at)))
    g = assemble_gaussian_gradient(config.obj, experiments, free)

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
