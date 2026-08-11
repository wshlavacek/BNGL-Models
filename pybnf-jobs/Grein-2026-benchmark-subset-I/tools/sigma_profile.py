#!/usr/bin/env python
"""How far is a slug's PEtab nominal point *really*, once its estimated sigmas are free?

`OG_nominal` (nominal_check.json) evaluates the objective at the PEtab `nominalValue` vector --
every parameter, including the noise parameters. For the many slugs here whose nominal sigmas are
placeholders (typically 1), that number is dominated by `sum n_j log sigma_j` being ~0 instead of
its MLE value, and it is NOT a measure of how far the *dynamics* are from the optimum. Issue #38
orders the remaining tuning candidates "roughly by nominal-point distance, i.e. plausibly by
difficulty", so the distortion is load-bearing.

This computes the honest version: hold every non-noise parameter at nominal, set every ESTIMATED
sigma to its MLE, and report the resulting OG. That is the best OG reachable without moving the
dynamics at all -- a lower bound on what the nominal point is worth, and the right difficulty proxy.

No PyBNF and no simulation: the reference trajectory is upstream's own `simulatedData` table, which
is the same oracle 2c uses. So this inherits 2c's coverage -- it works exactly for the slugs whose
`simulatedData` joins to `measurementData`.

SELF-CHECK. With the nominal sigmas substituted, the same code path must reproduce
nominal_check.json's `J_paper` -- the paper-scale Eq. 6 NLL. If it does not, this script is wrong,
not the slug, and the profiled number it prints is not to be believed.

Check against `J_paper`, NOT against `reduced_objective`. PyBNF's reduced objective drops only the
*parameter-independent* per-point constants; `sum n_j log sigma_j` depends on a fitted sigma and so
stays inside it. Comparing to `reduced_objective` reports a spurious failure on every estimated-sigma
slug -- which is exactly what the first version of this script did, on Brannmark, Laske and Zhao.

Usage:  sigma_profile.py <upstream-root> <slug-dir>...
"""
import glob
import json
import math
import os
import sys

import numpy as np
import pandas as pd

LN2PI = math.log(2 * math.pi)


def _find(d, *patterns):
    for p in patterns:
        hits = sorted(glob.glob(os.path.join(d, p)))
        if hits:
            return hits[0]
    return None


def _upstream_dir(root, slug, slug_dir=None):
    """Locate the upstream problem directory for a local slug.

    Usually the two names are the same. They are not for `Schwen_PONE2015`, which is renamed
    locally (the paper is 2015; upstream's directory, the model filename and Grein et al.'s own
    tables all say 2014). The local slug is therefore NOT a reliable join key, and guessing is
    what would silently skip such a slug -- so consult the `upstream_slug` the slug records in
    its own nominal_check.json before falling back to the directory name.
    """
    names = [slug]
    if slug_dir:
        try:
            with open(os.path.join(slug_dir, 'nominal_check.json')) as fh:
                declared = json.load(fh).get('upstream_slug')
            if declared and declared != slug:
                names.insert(0, declared)
        except (OSError, ValueError):
            pass
    for name in names:
        for cand in (os.path.join(root, 'Benchmark-Models', name), os.path.join(root, name)):
            if os.path.isdir(cand):
                return cand
    return None


def _read_tables(d):
    """measurement / simulated / observable / parameter tables, tolerating both layouts."""
    meas_p = _find(d, 'measurementData_training_*.tsv', 'measurementData_*.tsv')
    # measurementData_test_* is a held-out split; never the training objective.
    if meas_p and '_test_' in os.path.basename(meas_p):
        meas_p = None
    sim_p = _find(d, 'simulatedData_*.tsv')
    obs_p = _find(d, 'observables_*.tsv', 'observable_*.tsv')
    par_p = _find(d, 'parameters_*.tsv', 'parameters.tsv')
    if not all((meas_p, sim_p, obs_p, par_p)):
        missing = [n for n, p in (('measurement', meas_p), ('simulated', sim_p),
                                  ('observable', obs_p), ('parameters', par_p)) if not p]
        raise FileNotFoundError('no ' + '/'.join(missing) + ' table')
    return (pd.read_csv(meas_p, sep='\t'), pd.read_csv(sim_p, sep='\t'),
            pd.read_csv(obs_p, sep='\t'), pd.read_csv(par_p, sep='\t'))


def _transform(obs_row):
    for col in ('observableTransformation', 'observableTransform'):
        if col in obs_row and isinstance(obs_row[col], str):
            return obs_row[col].strip().lower()
    return 'lin'


def _apply(scale, y):
    if scale in ('lin', ''):
        return y
    if scale == 'log10':
        return np.log10(y)
    if scale in ('log', 'ln'):
        return np.log(y)
    raise ValueError('unsupported observableTransformation: ' + scale)


def _jacobian(scale, y):
    """sum log|d(transformed)/dy| -- the change-of-variables term the NLL carries."""
    if scale in ('lin', ''):
        return 0.0
    if scale == 'log10':
        return float(-np.sum(np.log(np.abs(y) * math.log(10))))
    if scale in ('log', 'ln'):
        return float(-np.sum(np.log(np.abs(y))))
    raise ValueError(scale)


def analyse(slug_dir, upstream_root, profile=True):
    slug = os.path.basename(os.path.abspath(slug_dir))
    up = _upstream_dir(upstream_root, slug, slug_dir)
    if up is None:
        raise FileNotFoundError('no upstream directory for ' + slug)
    meas, sim, obs, par = _read_tables(up)

    # `datasetId` is deliberately NOT a join key. PEtab defines it as a visualization
    # grouping label, not part of a measurement's identity, and the two tables are free to
    # disagree on it. `Smith_BMCSystBiol2013` is where that bit: 13 of its 62 rows -- all
    # ten of figure 2C and all three of figure 2D -- are tagged `fig2C`/`fig2c`/`fig2D` in
    # measurementData and `fig2A` in simulatedData, so including the column joined 49 of 62
    # and the slug was recorded as having no oracle at all. On the identity keys it joins
    # 62 of 62, one-to-one. Dropping the column cannot silently over-match: the
    # `len(j) != len(meas)` guard below rejects any key that stops being unique.
    keys = [c for c in ('observableId', 'simulationConditionId', 'preequilibrationConditionId',
                        'time', 'observableParameters', 'noiseParameters')
            if c in meas.columns and c in sim.columns]
    # The tables disagree on dtype wherever a column is sometimes blank (NaN -> float64 on one
    # side, str on the other), and on formatting where one side writes `1` and the other `1.0`.
    # Normalize both the same way: a numeric cell to its float repr, anything else to a bare
    # string, so the two sides are comparable without smuggling in a false match.
    def _key(v):
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return ''
        s = str(v).strip()
        if s in ('', 'nan'):
            return ''
        try:
            return repr(float(s))
        except ValueError:
            return s

    meas, sim = meas.copy(), sim.copy()
    for c in keys:
        if c != 'time':
            meas[c] = meas[c].map(_key)
            sim[c] = sim[c].map(_key)
    j = meas.merge(sim, on=keys, suffixes=('_meas', '_sim'))
    if len(j) != len(meas):
        # A duplicated join means the keys do not identify a row; a short one means rows are
        # missing. Either way the residuals would be wrong, so refuse rather than guess.
        raise ValueError(f'joined {len(j)} of {len(meas)} rows on {keys}')
    simcol = 'simulation' if 'simulation' in j.columns else next(
        c for c in j.columns if c.lower().startswith('simulat') and c not in keys)

    par = par.set_index('parameterId')
    estimated = {p for p in par.index if int(par.loc[p].get('estimate', 1)) == 1}
    nominal = par['nominalValue'].to_dict() if 'nominalValue' in par.columns else {}

    obs_idx = obs.set_index('observableId')

    # Each row's sigma is named by its noiseParameters cell (a parameter id, or a literal).
    def sigma_key(row):
        raw = row.get('noiseParameters')
        if raw is None or (isinstance(raw, float) and math.isnan(raw)):
            return ('formula', row['observableId'])
        return ('name', str(raw).strip())

    j['_sig'] = j.apply(sigma_key, axis=1)
    j['_scale'] = j['observableId'].map(lambda o: _transform(obs_idx.loc[o]))

    n = len(j)
    nll = n / 2 * LN2PI
    reduced = 0.0
    rows = []
    for key, g in j.groupby('_sig'):
        scales = set(g['_scale'])
        if len(scales) != 1:
            raise ValueError('mixed transformations under one sigma: ' + str(key))
        scale = scales.pop()
        r = (_apply(scale, g['measurement'].to_numpy(float))
             - _apply(scale, g[simcol].to_numpy(float)))
        S, nj = float(np.sum(r ** 2)), len(g)

        kind, name = key
        is_est = kind == 'name' and name in estimated
        if is_est and profile:
            sigma = math.sqrt(S / nj)
            src = 'profiled'
        elif kind == 'name' and name in nominal:
            sigma = float(nominal[name])
            src = 'nominal' if is_est else 'fixed'
        else:
            candidate = name
            if kind == 'formula':
                # No `noiseParameters` cell, so sigma is the observable's own `noiseFormula`.
                # `sigma_key` already routes these here; this finishes the lookup it names.
                # A literal number there is fully determined -- `Smith_BMCSystBiol2013`
                # carries `noiseFormula = 1.0` on all nine observables, i.e. sigma == 1 --
                # while a formula naming parameters needs the per-measurement binding whose
                # absence put the row on this branch, so it still refuses below.
                try:
                    candidate = str(obs_idx.loc[name]['noiseFormula']).strip()
                except KeyError:
                    candidate = name
            try:
                sigma = float(candidate)
                src = 'literal'
            except (TypeError, ValueError):
                raise ValueError('cannot resolve sigma for ' + str(key))
        nll += nj * math.log(sigma) + S / (2 * sigma ** 2)
        reduced += S / (2 * sigma ** 2)
        nll += _jacobian(scale, g['measurement'].to_numpy(float))
        rows.append((str(name), nj, S, sigma, src, scale))

    jstar = float(open(os.path.join(slug_dir, 'jstar.txt')).read().split()[0])
    return {'slug': slug, 'n': n, 'nll': nll, 'jstar': jstar, 'og': nll - jstar,
            'reduced': reduced, 'groups': rows,
            'n_estimated': sum(1 for r in rows if r[4] in ('profiled', 'nominal'))}


def main():
    root, slug_dirs = sys.argv[1], sys.argv[2:]
    print(f'{"slug":36s} {"n":>4s} {"estσ":>4s} {"OG_nominal":>12s} {"OG_σ-profiled":>14s} '
          f'{"inflation":>10s} {"self-check":>11s}')
    print('-' * 100)
    for sd in slug_dirs:
        slug = os.path.basename(os.path.abspath(sd))
        try:
            nc = json.load(open(os.path.join(sd, 'nominal_check.json')))
        except FileNotFoundError:
            print(f'{slug:36s} {"":>4s} -- no nominal_check.json')
            continue
        try:
            prof = analyse(sd, root, profile=True)
            chk = analyse(sd, root, profile=False)
        except Exception as exc:
            print(f'{slug:36s} {nc["n_scored"]:4d} {"":>4s} {nc["OG_nominal"]:12.4g} '
                  f'{"skipped":>14s}   {type(exc).__name__}: {exc}')
            continue
        # Self-check: with nominal sigmas this must reproduce nominal_check.json's J_paper.
        # Reported, not thresholded -- read it against the inflation column. A self-check of the
        # same order as the inflation means the trajectory is not the nominal one (or is written
        # to too few digits), and the profiled number says nothing.
        d = abs(chk['nll'] - nc['J_paper'])
        note = '  (no estimated σ -- profiling is a no-op)' if prof['n_estimated'] == 0 else ''
        print(f'{slug:36s} {prof["n"]:4d} {prof["n_estimated"]:4d} {nc["OG_nominal"]:12.4g} '
              f'{prof["og"]:14.4f} {nc["OG_nominal"] - prof["og"]:10.4g} {d:11.2e}{note}')


if __name__ == '__main__':
    main()
