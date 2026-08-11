#!/usr/bin/env python
"""How many points sampled from a slug's fit box actually integrate — and which solve fails.

A multistart fit is only as large as the number of start points whose forward model runs.
This draws N points from the conf's own sampling distribution (the box, in sampling space)
and simulates each, reporting how many survive. The point of the tool is the `--sens` flag:
run it twice, with and without, and the difference isolates the *forward-sensitivity* solve
from the *state* solve.

That distinction is the whole answer on two slugs. Both had been recorded as
lanl/bngsim#196 cases — "a scalar atol cannot serve a model spanning ten decades":

    Weber_BMC2015      plain 7/11 in 0.6 s   with sensitivities 2/11 in 80 s
    Brannmark_JBC2010  plain 30/30 in 1.0 s  with sensitivities 19/30 in 210 s

In both, the state solve is fine and the sensitivity solve is what dies — because CVODES
derives its sensitivity tolerances from the state ones. Measuring the aggregate, or
measuring through `fd_check.py` (which always applies the sensitivity request), cannot see
this and reads as an unconditional failure.

Usage:  box_probe.py <slug-dir> [-n N] [--sens] [--atol A] [--rtol R] [--wall S]
                                [--seed S] [--nominal UPSTREAM_PETAB_DIR]

Gotchas, all of which produced a plausible wrong answer first:

* **`--atol` writes a temporary conf; it does not poke `config.config`.** The tolerance keys
  are read when `load_config` BUILDS the model, so setting `config.config['sbml_atol']`
  afterwards is a silent no-op that leaves `_config_atol` at None. The symptom is a sweep
  whose rows are all identical, which reads as "the tolerance makes no difference".
  The header line prints `_config_atol` for exactly this reason — check it moved.

* **An explicit `--atol` turns the ADR-0105 per-species vector OFF.** `sbml_atol` is the
  documented off-switch for that whole mechanism, so `--atol <the same scalar the derivation
  produces>` is NOT a no-op comparison: it is vector-vs-scalar. On `Brannmark` both give
  19/30, but the vector is 1.4x faster (210 s vs 303 s). Compare against the no-`--atol` run
  to see what the vector is worth.

* **`--sens` is not free and not optional for a gradient job.** `gntr`/`trf` and friends
  apply the union sensitivity request on every evaluation, so the `--sens` number is the one
  that predicts start mortality. Without it you are measuring a job nobody is running.
  Predictions have held: 22/30 forecast 24 dead starts out of 100 on Weber's real run.

* **A dead point here is not necessarily a tolerance problem.** Points that fail both with
  and without `--sens`, and fail fast, are bad parameter points. Weber has two of those in
  eleven, and they stay dead at every tolerance.
"""
import argparse
import csv
import logging
import os
import sys
import tempfile
import time

import numpy as np
from pybnf.algorithms import core
from pybnf.gradient import apply_routings, route_for_model
from pybnf.parse import load_config
from pybnf.pset import PSet


def nominal_values(petab_dir):
    path = [os.path.join(petab_dir, f) for f in sorted(os.listdir(petab_dir))
            if f.startswith('parameters') and f.endswith('.tsv')][0]
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if str(row.get('estimate', '1')).strip() in ('1', 'true', 'True'):
                out[row['parameterId']] = float(row['nominalValue'])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('slug_dir')
    ap.add_argument('-n', type=int, default=20)
    ap.add_argument('--atol', type=float, default=None)
    ap.add_argument('--rtol', type=float, default=None)
    ap.add_argument('--wall', type=float, default=None)
    ap.add_argument('--seed', type=int, default=7)
    ap.add_argument('--nominal', default=None,
                    help='upstream PEtab dir; adds the nominal point as point 0')
    ap.add_argument('--sens', action='store_true',
                    help='apply the gradient sensitivity request (what gntr asks for)')
    args = ap.parse_args()

    logging.basicConfig(level=logging.ERROR)
    slug_dir = os.path.abspath(args.slug_dir)
    nominals = nominal_values(os.path.abspath(args.nominal)) if args.nominal else None
    conf = [os.path.join(slug_dir, f) for f in sorted(os.listdir(slug_dir))
            if f.endswith('.conf')][0]
    os.chdir(slug_dir)
    # The tolerance keys are read off the config when load_config BUILDS the model, so an
    # override has to be in the conf text -- setting config.config afterwards is a no-op that
    # leaves _config_atol None and reads as "the override made no difference".
    probe_conf = conf
    if args.atol is not None or args.rtol is not None or args.wall is not None:
        extra = []
        if args.atol is not None:
            extra.append('sbml_atol = %r' % args.atol)
        if args.rtol is not None:
            extra.append('sbml_rtol = %r' % args.rtol)
        if args.wall is not None:
            extra.append('wall_time_sim = %r' % args.wall)
        probe_conf = os.path.join(slug_dir, '_box_probe.conf')
        with open(conf) as fh:
            text = fh.read()
        # A repeated key is a parse error, so drop any the override replaces.
        keep = [ln for ln in text.splitlines()
                if ln.split('=')[0].strip() not in
                {e.split('=')[0].strip() for e in extra}]
        with open(probe_conf, 'w') as fh:
            fh.write('\n'.join(keep) + '\n' + '\n'.join(extra) + '\n')
    try:
        config = load_config(probe_conf)
    finally:
        if probe_conf != conf:
            os.unlink(probe_conf)
    models = list(config.models.values())
    print('conf:', os.path.basename(conf))
    for m in models:
        print('  model %s: _config_atol=%r _config_rtol=%r nominal_state_scale=%r wall=%r'
              % (m.name, getattr(m, '_config_atol', 'n/a'),
                 getattr(m, '_config_rtol', 'n/a'),
                 getattr(m, '_nominal_state_scale', 'n/a'),
                 config.config['wall_time_sim']))

    rng = np.random.default_rng(args.seed)
    free = list(config.variables)
    if args.sens:
        # The union sensitivity request a gntr fit applies, exactly as fd_check.py builds it.
        names = [v.name for v in free]
        for m in models:
            m.set_scored_suffixes(config.exp_data.get(m.name, {}))
            wt = route_for_model(m, names, None)
            mine = []
            for suffix in config.exp_data.get(m.name, {}):
                best = None
                for mut in getattr(m, 'mutants', []) or []:
                    ms = getattr(mut, 'suffix', '')
                    if ms and suffix.endswith(ms) and (best is None or len(ms) > len(best.suffix)):
                        best = mut
                mine.append(route_for_model(m, names, best))
            apply_routings(m, [wt, *mine])
    psets = []
    labels = []
    if nominals:
        psets.append(PSet([v.set_value(nominals[v.name]) for v in free]))
        labels.append('nominal')
    for i in range(args.n):
        vals = []
        for v in free:
            lo = v.to_sampling_space(v.lower_bound)
            hi = v.to_sampling_space(v.upper_bound)
            u = rng.uniform(lo, hi)
            vals.append(v.set_value(v.from_sampling_space(u)))
        psets.append(PSet(vals))
        labels.append('box%02d' % i)

    ok = 0
    t0 = time.time()
    for label, pset in zip(labels, psets):
        with tempfile.TemporaryDirectory() as sim_dir:
            job = core.Job(models, pset, label, sim_dir,
                           config.config['wall_time_sim'], None,
                           config.config['normalization'], config.postprocessing, True,
                           stochastic_seed_policy=config.config['stochastic_seed'])
            t1 = time.time()
            res = core.run_job(job)
        failed = bool(getattr(res, 'failed', False)) or res.simdata is None
        if not failed:
            ok += 1
        print('  %-9s %-6s %6.2fs' % (label, 'FAIL' if failed else 'ok', time.time() - t1))
        sys.stdout.flush()
    print('SUMMARY  %d/%d integrated in %.1fs' % (ok, len(psets), time.time() - t0))


if __name__ == '__main__':
    main()
