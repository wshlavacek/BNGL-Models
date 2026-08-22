#!/usr/bin/env python
"""Which declared parameters could a linear-observable profile actually remove -- and on
what residual scale is each one linear?

Issue lanl/PyBNF#572 proposes profiling an observable's additive **offset**, and the coupled
`(scale, offset)` pair, out of the search by variable projection: stack `Phi = [s, 1]`, solve
`c* = (Phi^T W Phi)^-1 Phi^T W d`, score the projection residual. That construction is a
**linear least-squares** identity. It holds when the family's residual is `d - y` and `y` is
affine in the profiled coefficients.

A PyBNF noise family declares the space its residual lives in (`family.additive_on.ln_base`:
0 for a linear-scale Gaussian/Laplace, `ln 10` for `lognormal`, `1` for `lnnormal`). On a
**log** family the residual is `log d - log(a*s + b)`. That is affine in `log a` when `b = 0`
-- which is exactly why ADR-0066's geometric-mean ratio is a closed form for a pure scale --
and it is **not** affine in `b` for any `b`. So the linear-profiling question is not one
question but two, and this script answers them separately per parameter:

    linear-in-prediction   does the observable formula depend affinely on this parameter?
    linear-in-residual     ... AND does its observable score on a linear-scale family?

Only a parameter that is both is a candidate for a closed-form profile.

Usage:  linear_scope.py <slug-dir>...            (or a corpus root: every slug under it)
        linear_scope.py --json out.json <dirs>

Reported per candidate parameter:

* `observables` -- every scored observable whose formula names it. A parameter naming more
  than one observable, or one observable measured in more than one experiment, is **tied**:
  its profile is a solve over the stacked series, not the per-series solve #572 writes down.
* `series` -- how many (experiment, observable) series it is tied across.
* `ln_base` -- the residual scale(s) of those observables. Mixed scales on one tied parameter
  mean there is no single space its profile could be taken in.
* `also_noise` -- the parameter is *simultaneously* bound as a noise parameter (Fiedler's
  `s_pErk_*`). It is not a free linear coefficient: moving it moves sigma too, so the
  projection would be solving the wrong problem. #572 asks for this to be refused by name.

Gotchas, each of which produced a plausible wrong answer first:

* **A row-varying token is a parameter too.** `Brannmark`'s `k_IRSiP_1Step` never appears in
  the conf's `observable:` line -- it is bound per data row through the measurement-params
  table, to `observableParameter1_IRS1_P`. Reading only `MeasurementModel.formula` finds two
  of Brannmark's scales and misses the two the issue is about. Both binding routes are walked
  here, and a per-row token is credited to the observable it is bound in.

* **Affine in the declared name is not the same as affine in the span.** `Schwen` declares
  `scale*(IR1 + IR1in + offset)`, which is *quadratic* in `(scale, offset)` jointly -- the
  cross term `scale*offset` -- while spanning the same two columns as `[s, 1]`. So the
  per-parameter affinity test below is reported alongside a **joint** test (is the formula
  affine in all of the observable's candidate parameters *at once*?), and a formula that is
  separately-but-not-jointly affine is flagged `reparam`. A profile solves for span
  coefficients; mapping back to the declared names is a separate, and at `scale -> 0`
  ill-posed, step.

* **`estimate = 0` is not searched.** Only parameters the conf declares as free variables are
  counted; a fixed observable coefficient is not a search dimension and removing it would
  change the fit.
"""
import argparse
import json
import os
import sys

import logging

logging.basicConfig(level=logging.ERROR)


def _sympy_expr(formula):
    """The parsed PEtab math for one observableFormula, as a sympy expression."""
    from pybnf.petab.formula import _parse, _require_petab_math
    return _parse(_require_petab_math(), formula, source='observableFormula')


def _symbol(expr, name):
    """The symbol object named ``name`` inside ``expr``, or ``None``.

    **Gotcha, and it silently inverts every answer this script gives.** PEtab's math parser
    builds symbols carrying assumptions (``Symbol('scale', real=True)``), which do **not**
    compare equal to a bare ``sympy.Symbol('scale')``. So `sp.Symbol(name) in
    expr.free_symbols` is `False` for a symbol that is plainly there, `diff` against the bare
    symbol is identically 0, and every affinity test passes vacuously -- reporting the whole
    corpus as linear, including formulas that are nothing of the kind. Resolve by *name*
    against the expression's own symbols instead of constructing one.
    """
    for s in expr.free_symbols:
        if str(s) == name:
            return s
    return None


def _role_of(expr, name):
    """How one affine parameter enters the formula: ``'scale'``, ``'offset'`` or ``'affine'``.

    An affine dependence is ``expr = A(rest)*p + B(rest)``. The distinction matters because a
    **log** family has a closed form for one of the three and not the others:

    * ``scale``  -- ``B == 0`` identically, so ``expr = A*p`` and ``log expr = log A + log p``
      is affine in ``log p``. This is ADR-0066's geometric-mean ratio, and it survives on a
      log family. Note ``B == 0`` is a property of the whole formula: Borghans's ``scale`` in
      ``Z_state*scale + offset`` is **not** a pure scale, because ``B = offset``.
    * ``offset`` -- ``A == 1``.
    * ``affine`` -- everything else (``A`` non-unit and ``B`` nonzero).
    """
    import sympy as sp
    p = _symbol(expr, name)
    if p is None:
        return None
    a = sp.simplify(sp.diff(expr, p))
    b = sp.simplify(expr.subs(p, 0))
    if b == 0:
        return 'scale'
    if a == 1:
        return 'offset'
    return 'affine'


def _is_affine_in(expr, names):
    """Is ``expr`` jointly affine in every symbol of ``names``?

    Affine means every second partial (including the cross terms) vanishes identically --
    the exact condition under which ``expr`` can be written ``Phi(rest) . c + const`` with
    ``c`` the named symbols. Checking only ``d2/dp2 == 0`` per parameter admits Schwen's
    ``scale*(x + offset)``, which is affine in each of ``scale``/``offset`` separately and
    not in the pair.
    """
    import sympy as sp
    present = [s for s in (_symbol(expr, n) for n in names) if s is not None]
    if not present:
        return True
    for i, a in enumerate(present):
        for b in present[i:]:
            if sp.simplify(sp.diff(expr, a, b)) != 0:
                return False
    return True


def _observable_formulas(obj):
    """``{observable_id: (formula, kind)}`` for every measurement-model observable --
    both the pre-materialized constant-per-observable models and the row-varying ones."""
    out = {}
    for mm in (obj.measurement.models if obj.measurement else []):
        out[mm.observable_id] = (mm.formula, 'column')
    for oid, pm in getattr(obj, '_per_measurement_models', {}).items():
        out[oid] = (pm.formula, 'per-row')
    return out


def _row_bound_params(config, observable_id, placeholder_names):
    """The declared parameter ids a row-varying observable's placeholders resolve to, keyed
    by placeholder: ``{placeholder: {param_name: set(data_key)}}``.

    A per-measurement binding table maps a placeholder (``observableParameter1_x``) to a
    token per data row; a numeric token inlines and is not a parameter.

    **Keyed by placeholder, not flattened.** Schwen's ``observable_Insulin`` is
    ``op1 + op2 * g(.; op3, op4)`` -- two placeholders affine, two inside a Michaelis-Menten
    denominator and not. Flattening the bindings loses which token was which, so the only
    available affinity test is "are they all affine", which rejects the whole observable and
    silently drops the two parameters that are exactly what #572 calls partial separability.
    """
    found = {ph: {} for ph in placeholder_names}
    for _model, by_suffix in config.exp_data.items():
        for suffix, data in by_suffix.items():
            table = getattr(data, 'measurement_params', None)
            if not table or observable_id not in table:
                continue
            for ph in placeholder_names:
                for token in table[observable_id].get(ph, []):
                    try:
                        float(token)
                    except (TypeError, ValueError):
                        found[ph].setdefault(str(token), set()).add(suffix)
    return found


def _noise_param_names(obj):
    """Every parameter name any noise source reads -- a free sigma, or a name appearing in
    a sigma *formula* (Fiedler binds the same token to both layers, so the formula path is
    the one that catches it)."""
    import re
    names = set()
    specs = [obj._default_sources()] + [s for _f, s in obj.overrides.values()]
    for sources in specs:
        for src in sources.values():
            for attr in ('name', 'formula'):
                val = getattr(src, attr, None)
                if isinstance(val, str):
                    if attr == 'name':
                        names.add(val)
                    else:
                        names.update(re.findall(r'[A-Za-z_]\w*', val))
    # A row-varying sigma resolves its placeholder from the binding table, like an
    # observable's; walk it the same way.
    return names


def _row_bound_noise_params(config, obj):
    """Parameter ids a *noise* placeholder resolves to per row, keyed by observable."""
    import re
    out = {}
    for col, (_family, sources) in obj.overrides.items():
        phs = set()
        for src in sources.values():
            formula = getattr(src, 'formula', None)
            if isinstance(formula, str):
                phs.update(re.findall(r'noiseParameter\d+_?\w*', formula))
        if not phs:
            continue
        for _model, by_suffix in config.exp_data.items():
            for _suffix, data in by_suffix.items():
                table = getattr(data, 'measurement_params', None)
                if not table or col not in table:
                    continue
                for ph in phs:
                    for token in table[col].get(ph, []):
                        try:
                            float(token)
                        except (TypeError, ValueError):
                            out.setdefault(str(token), set()).add(col)
    return out


def scan(conf_path):
    """Classify one slug's observable-layer parameters. Returns a plain-data dict."""
    from pybnf.parse import load_config

    slug_dir = os.path.dirname(os.path.abspath(conf_path))
    cwd = os.getcwd()
    os.chdir(slug_dir)
    try:
        config = load_config(os.path.abspath(conf_path))
    finally:
        os.chdir(cwd)

    obj = config.obj
    free = {v.name for v in config.variables}
    formulas = _observable_formulas(obj)
    noise_names = _noise_param_names(obj)
    row_noise = _row_bound_noise_params(config, obj)

    # Which (experiment, observable) series exist -- an observable measured in N experiments
    # is N series, and a parameter naming it is tied across all of them.
    series_of = {}
    for _model, by_suffix in config.exp_data.items():
        for suffix, data in by_suffix.items():
            for col in data.cols:
                series_of.setdefault(col, set()).add(suffix)

    params = {}

    def note(name, observable, kind, tied_series, role):
        rec = params.setdefault(name, {
            'name': name, 'observables': {}, 'series': set(), 'roles': set(),
            'searched': name in free, 'also_noise': False, 'binding': kind})
        rec['observables'].setdefault(observable, kind)
        rec['series'] |= set(tied_series)
        rec['roles'].add(role)

    for oid, (formula, kind) in sorted(formulas.items()):
        expr = _sympy_expr(formula)
        symbols = {str(s) for s in expr.free_symbols}
        direct = sorted(symbols & free)
        placeholders = sorted(s for s in symbols if s.startswith('observableParameter'))
        row_bound = _row_bound_params(config, oid, placeholders)

        # Which of this observable's parameters enter the *prediction* affinely.
        candidates = []
        for p in direct:
            if _is_affine_in(expr, [p]):
                candidates.append((p, p, series_of.get(oid, set())))
        # A row-varying token's affinity is its *placeholder's* -- tested per placeholder, so
        # a partially separable observable contributes its affine tokens and not the others.
        for ph, bound in sorted(row_bound.items()):
            if not _is_affine_in(expr, [ph]):
                continue
            for p, suffixes in sorted(bound.items()):
                candidates.append((p, ph, suffixes))

        joint_names = sorted({sym for _p, sym, _s in candidates})
        joint_affine = _is_affine_in(expr, joint_names)

        for p, sym, suffixes in candidates:
            if p not in free:
                continue
            note(p, oid, kind, suffixes, _role_of(expr, sym) or 'affine')
            params[p]['reparam'] = not joint_affine

    # Family / residual scale per observable, and the double-binding flag.
    fam_of = {}
    for oid in formulas:
        family, _sources = obj._spec_for(oid)
        fam_of[oid] = {'family': type(family).__name__,
                       'ln_base': float(family.additive_on.ln_base)}
    for name, rec in params.items():
        rec['also_noise'] = bool(name in noise_names or name in row_noise)
        rec['ln_bases'] = sorted({fam_of[o]['ln_base'] for o in rec['observables']})
        rec['families'] = sorted({fam_of[o]['family'] for o in rec['observables']})
        rec['n_series'] = len(rec['series'])
        rec['series'] = sorted(rec['series'])
        rec['observables'] = sorted(rec['observables'])
        rec['roles'] = sorted(rec['roles'])
        rec['linear_in_prediction'] = True
        # Which closed form, if any, actually exists for this parameter:
        #   linear-lsq   -- a linear-scale family; variable projection over `[s, 1]` applies
        #                   whatever the role (this is #572's construction)
        #   log-geomean  -- a log family, but the parameter is a PURE multiplicative factor
        #                   of the whole formula, so it is affine in log space (ADR-0066)
        #   none         -- a log family and the formula is not homogeneous in it
        # A parameter also bound as a noise parameter is `none` regardless: moving it moves
        # sigma, so it is not a free linear coefficient (Fiedler).
        if rec['also_noise']:
            rec['closed_form'] = 'none'
        elif rec['ln_bases'] == [0.0]:
            rec['closed_form'] = 'linear-lsq'
        elif rec['roles'] == ['scale']:
            rec['closed_form'] = 'log-geomean'
        else:
            rec['closed_form'] = 'none'
        rec['linear_in_residual'] = rec['closed_form'] == 'linear-lsq'

    return {
        'slug': os.path.basename(slug_dir),
        'conf': os.path.basename(conf_path),
        # `sos` / `norm_sos` are not likelihoods: they carry no sigma, so `_spec_for`'s
        # family is a structural default and only its `ln_base` (the residual's space) is
        # meaningful. `Smith_BMCSystBiol2013` is the corpus's `sos` slug, and its residual
        # IS `sim - exp` on the linear scale, so its classification below stands -- but the
        # sigma-weighting question #572 raises simply does not arise there.
        'objective': type(obj).__name__,
        'n_searched': len(free),
        'observables': fam_of,
        'params': [params[k] for k in sorted(params)],
    }


def _default_conf(slug_dir):
    """A slug's canonical conf: `<slug>.conf`, else the only/first `.conf` present."""
    base = os.path.basename(os.path.abspath(slug_dir))
    exact = os.path.join(slug_dir, base + '.conf')
    if os.path.exists(exact):
        return exact
    confs = sorted(f for f in os.listdir(slug_dir) if f.endswith('.conf'))
    return os.path.join(slug_dir, confs[0]) if confs else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dirs', nargs='+')
    ap.add_argument('--json', default=None)
    args = ap.parse_args()

    slug_dirs = []
    for d in args.dirs:
        d = os.path.abspath(d)
        if _default_conf(d):
            slug_dirs.append(d)
        else:
            slug_dirs.extend(sorted(os.path.join(d, x) for x in os.listdir(d)
                                    if os.path.isdir(os.path.join(d, x))
                                    and _default_conf(os.path.join(d, x))))

    out = []
    for d in slug_dirs:
        conf = _default_conf(d)
        try:
            rec = scan(conf)
        except Exception as exc:                      # a slug that will not load is data
            rec = {'slug': os.path.basename(d), 'conf': os.path.basename(conf),
                   'error': '%s: %s' % (type(exc).__name__, exc), 'params': []}
        out.append(rec)
        if 'error' in rec:
            print('%-34s  LOAD FAILED  %s' % (rec['slug'], rec['error']))
            continue
        cands = rec['params']
        if not cands:
            print('%-34s  k=%-3d  no observable-layer linear parameters' % (rec['slug'], rec['n_searched']))
            continue
        n_lsq = sum(p['closed_form'] == 'linear-lsq' for p in cands)
        n_geo = sum(p['closed_form'] == 'log-geomean' for p in cands)
        n_none = sum(p['closed_form'] == 'none' for p in cands)
        print('%-34s  k=%-3d  %d affine obs params:  %d linear-lsq, %d log-geomean, %d NO closed form'
              % (rec['slug'], rec['n_searched'], len(cands), n_lsq, n_geo, n_none))
        for p in cands:
            flags = []
            if p['also_noise']:
                flags.append('ALSO-NOISE')
            if p.get('reparam'):
                flags.append('reparam')
            if p['n_series'] > 1:
                flags.append('tied x%d' % p['n_series'])
            print('    %-32s %-7s ln_base=%-8s %-12s %s'
                  % (p['name'], '/'.join(p['roles']),
                     ','.join('%g' % b for b in p['ln_bases']),
                     p['closed_form'],
                     ' '.join(flags)))
        sys.stdout.flush()

    if args.json:
        with open(args.json, 'w') as fh:
            json.dump(out, fh, indent=2, sort_keys=True)
        print('wrote', args.json)


if __name__ == '__main__':
    main()
