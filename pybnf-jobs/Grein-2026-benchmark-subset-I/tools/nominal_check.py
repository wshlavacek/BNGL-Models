#!/usr/bin/env python
"""Recompute a subset-I slug's nominal_check.json.

Evaluates PyBNF's objective at the PEtab `nominalValue` point and puts it on the paper's
Eq. 6 scale the same way a real run does -- through `likelihood_information_criteria`,
which is exactly what writes `information_criteria.txt` at the end of a fit.

Usage:  nominal_check.py <slug-dir> <upstream-petab-dir> [--write]
"""
import csv
import json
import os
import sys
import tempfile

from pybnf.algorithms import core
from pybnf.objective import likelihood_information_criteria
from pybnf.parse import load_config
from pybnf.pset import PSet


def nominal_values(petab_dir):
    """{parameterId: nominalValue} for the estimated parameters (linear scale, per spec)."""
    path = os.path.join(petab_dir, 'parameters.tsv')
    if not os.path.exists(path):
        # PEtab v1 layout: parameters_<slug>.tsv
        path = [os.path.join(petab_dir, f) for f in sorted(os.listdir(petab_dir))
                if f.startswith('parameters') and f.endswith('.tsv')][0]
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if str(row.get('estimate', '1')).strip() in ('1', 'true', 'True'):
                out[row['parameterId']] = float(row['nominalValue'])
    return out


def evaluate(conf_path, nominals):
    config = load_config(conf_path)
    free = [v.name for v in config.variables]
    missing = [n for n in free if n not in nominals]
    if missing:
        raise SystemExit('no nominalValue for: {}'.format(', '.join(missing)))

    pset = PSet([v.set_value(nominals[v.name]) for v in config.variables])
    exp_data = config.exp_data

    def simulate(tag):
        """One fresh Result at the nominal point.

        A measurement-model layer (ADR-0036) materializes its observable columns onto the
        Data in place, so a second scoring pass over the *same* Result collides with the
        columns the first one added. A real fit re-simulates for the information criteria
        for exactly this reason (`_compute_information_criteria`), so do the same here.
        """
        with tempfile.TemporaryDirectory() as sim_dir:
            job = core.Job(list(config.models.values()), pset, tag, sim_dir,
                           config.config['wall_time_sim'], None,
                           config.config['normalization'], config.postprocessing, True,
                           stochastic_seed_policy=config.config['stochastic_seed'])
            res = core.run_job(job)
        if getattr(res, 'failed', False) or res.simdata is None:
            raise SystemExit('the nominal point failed to simulate')
        res.normalize(config.config['normalization'])
        res.postprocess_data(config.postprocessing)
        return res

    objective = config.obj
    reduced = objective.evaluate_multiple(
        simulate('nominal_obj').simdata, exp_data, pset, config.constraints)
    ic = likelihood_information_criteria(
        objective, simulate('nominal_ic').simdata, exp_data, pset, len(config.variables))
    return reduced, ic


def main():
    slug_dir = os.path.abspath(sys.argv[1])
    petab_dir = os.path.abspath(sys.argv[2])
    write = '--write' in sys.argv
    slug = os.path.basename(slug_dir)
    conf = os.path.join(slug_dir, f'{slug}.conf')
    if not os.path.exists(conf):
        conf = [os.path.join(slug_dir, f) for f in os.listdir(slug_dir) if f.endswith('.conf')][0]

    jstar = float(open(os.path.join(slug_dir, 'jstar.txt')).read().split()[0])
    old_path = os.path.join(slug_dir, 'nominal_check.json')
    old = json.load(open(old_path)) if os.path.exists(old_path) else {}

    cwd = os.getcwd()
    os.chdir(slug_dir)
    try:
        reduced, ic = evaluate(conf, nominal_values(petab_dir))
    finally:
        os.chdir(cwd)

    j_paper = -ic.log_likelihood
    new = dict(old)
    new.update({
        'jstar': jstar,
        'J_paper': j_paper,
        'reduced_objective': reduced,
        'n_scored': ic.n,
        'k': ic.k,
        'OG_nominal': j_paper - jstar,
    })
    print(json.dumps({k: new[k] for k in
                      ('slug', 'jstar', 'J_paper', 'reduced_objective', 'n_scored', 'k',
                       'OG_nominal') if k in new}, indent=2))
    if old:
        print('was: J_paper={!r}  OG_nominal={!r}  n={!r}'.format(
            old.get('J_paper'), old.get('OG_nominal'), old.get('n_scored')))
    if write:
        with open(old_path, 'w') as fh:
            json.dump(new, fh, indent=2)
        print('wrote', old_path)


if __name__ == '__main__':
    main()
