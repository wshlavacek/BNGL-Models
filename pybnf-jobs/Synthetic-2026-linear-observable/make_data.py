#!/usr/bin/env python
"""Generate `experiment1.exp` for the linear-observable fixture, from a closed form.

The model is `A -> B -> 0` with `A(0) = 1`, `B(0) = 0`, so

    B(t) = k1 / (k2 - k1) * (exp(-k1 t) - exp(-k2 t))

and the observable is `y = scale * B + offset`. Committed so the data is reproducible and
so the fixture's truth is auditable rather than folklore.

**Generated from the closed form, not from a PyBNF simulation, on purpose.** If the SBML
in this directory did not encode the ODE above, the two would disagree and the profile at
the truth would not recover `scale = 3` / `offset = 1.5`. That recovery is the fixture's
own correctness test, and it is only a test because the data came from somewhere else.

Usage:  python make_data.py [--write]
"""
import argparse
import json
import math
import os

TRUTH = {'k1': 0.8, 'k2': 0.25, 'scale': 3.0, 'offset': 1.5, 'sigma': 0.05}
TIMES = [0.25 * i for i in range(1, 81, 3)]          # 0.25 .. 19.75, 27 points
SEED = 20260821


def b_of_t(t, k1, k2):
    return k1 / (k2 - k1) * (math.exp(-k1 * t) - math.exp(-k2 * t))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()

    import numpy as np
    rng = np.random.default_rng(SEED)
    here = os.path.dirname(os.path.abspath(__file__))

    rows = []
    for t in TIMES:
        clean = TRUTH['scale'] * b_of_t(t, TRUTH['k1'], TRUTH['k2']) + TRUTH['offset']
        rows.append((t, clean + rng.normal(0.0, TRUTH['sigma'])))

    text = '# time\ty\n' + ''.join('%.10g\t%.10g\n' % r for r in rows)
    if args.write:
        with open(os.path.join(here, 'experiment1.exp'), 'w') as fh:
            fh.write(text)
        with open(os.path.join(here, 'truth.json'), 'w') as fh:
            json.dump(TRUTH, fh, indent=1, sort_keys=True)
        print('wrote experiment1.exp (%d points) and truth.json' % len(rows))
    else:
        print(text, end='')

    # The k1 <-> k2 mirror, stated numerically rather than asserted. Swapping the two rate
    # constants multiplies B by k2/k1, which a free `scale` absorbs exactly -- so this
    # fixture has TWO global optima, and a fit that reports k1 = 0.25 has not failed.
    mirror = TRUTH['scale'] * TRUTH['k1'] / TRUTH['k2']
    worst = max(abs(TRUTH['scale'] * b_of_t(t, TRUTH['k1'], TRUTH['k2'])
                    - mirror * b_of_t(t, TRUTH['k2'], TRUTH['k1'])) for t in TIMES)
    print('k1<->k2 mirror: scale %.6g reproduces the same curve; max deviation %.3g'
          % (mirror, worst))


if __name__ == '__main__':
    main()
