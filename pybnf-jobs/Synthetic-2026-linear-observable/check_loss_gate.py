#!/usr/bin/env python
"""Which noise families may profile a linear coefficient: the residual's space, or the loss?

ADR-0123's finding-1 table sorts families by the space their residual lives in
(`family.additive_on.ln_base`: 0 linear, nonzero log) and puts `laplace` in the least-squares
column alongside `gaussian` and `sos`. PyBNF's own `LikelihoodObjective.is_linear_gaussian()`
disagrees: it is `True` for this fixture's `gaussian` conf and `False` for the `laplace` one.

They cannot both be the right gate for lanl/PyBNF#572, and this decides it by measurement.

Variable projection minimizes a SUM OF SQUARES. A Laplace likelihood is a sum of ABSOLUTE
residuals, so its best `(scale, offset)` is a least-absolute-deviations fit. The residual being on
a linear scale is necessary and not sufficient. This builds the design matrix once, computes both
candidate solutions, and scores each of them under each conf. Under `gaussian` the least squares
point wins; under `laplace` the L1 point wins, by more than the whole effect #572 is about.

The design matrix is built from the `gaussian` conf, by the same basis-evaluation trick
`linear_profile.py::_varpro` uses -- with every profiled coefficient at 0 the aligned prediction is
`Phi.0`, and setting coefficient j to 1 gives `Phi.e_j`. That is legitimate for both confs because
the model, the observable formula and the data are identical between them; only the loss differs.
`aligned_prediction_data` refuses the `laplace` conf outright, which is the same refusal stated a
second way.

Run from this directory:  python check_loss_gate.py
"""
import json
import logging
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'Grein-2026-benchmark-subset-I', 'tools'))
import linear_profile as lp  # noqa: E402

NAMES = ('scale', 'offset')
CONFS = ('linear_observable.conf', 'linear_observable_laplace.conf')


def design(config, scorer):
    """`(Phi, d, w)` over the scored points, from the objective's own aligned prediction."""
    def aligned(overrides):
        values = dict(scorer.base)
        values.update(overrides)
        return config.obj.aligned_prediction_data(
            __import__('copy').deepcopy(scorer.simdata), config.exp_data,
            [lp._Var(n, x) for n, x in values.items()])

    base = aligned({n: 0.0 for n in NAMES})
    if base is None:
        raise SystemExit('aligned_prediction_data refused the conf the design is built from')
    p0, d, var = base
    columns = [aligned({m: (1.0 if m == n else 0.0) for m in NAMES})[0] - p0 for n in NAMES]
    return np.column_stack(columns), np.asarray(d, float), 1.0 / np.asarray(var, float)


def main():
    logging.basicConfig(level=logging.ERROR)
    os.chdir(HERE)
    truth = json.load(open('truth.json'))

    configs = {name: lp._load(name, []) for name in CONFS}
    for name, config in configs.items():
        print('%-34s is_linear_gaussian() = %s   additive_on.ln_base = %g'
              % (name, config.obj.is_linear_gaussian(),
                 config.obj.noise.additive_on.ln_base))
    print()

    scorers, simdata = {}, None
    for name, config in configs.items():
        simdata = lp._simulate(config, truth, 'gate')
        if simdata is None:
            raise SystemExit('the truth point did not integrate')
        scorers[name] = lp.Scorer(config, simdata, truth)

    gaussian = configs[CONFS[0]]
    phi, d, w = design(gaussian, scorers[CONFS[0]])
    squares = np.linalg.solve(phi.T @ (w[:, None] * phi), phi.T @ (w * d))

    from scipy.optimize import minimize
    absolute = minimize(lambda c: float(np.sum(np.sqrt(w) * np.abs(d - phi @ c))),
                        np.array(squares), method='Nelder-Mead',
                        options={'xatol': 1e-10, 'fatol': 1e-12,
                                 'maxiter': 20000, 'maxfev': 20000}).x

    print('%d scored points; the truth is scale = %g, offset = %g'
          % (phi.shape[0], truth['scale'], truth['offset']))
    print('least squares       (scale, offset) = (%.6f, %.6f)' % tuple(squares))
    print('least abs deviation (scale, offset) = (%.6f, %.6f)' % tuple(absolute))
    print()
    print('%-34s %16s %16s %s' % ('conf', 'least squares', 'by L1', 'which wins'))
    for name in CONFS:
        score = scorers[name].score
        a, b = score(dict(zip(NAMES, squares))), score(dict(zip(NAMES, absolute)))
        print('%-34s %16.6f %16.6f %s by %.4f'
              % (name, a, b, 'least squares' if a < b else 'L1', abs(a - b)))
    print()
    print('The residual space is the same in both rows. The loss is not, and the loss is what')
    print('decides which solve is the conditional optimum.')


if __name__ == '__main__':
    main()
