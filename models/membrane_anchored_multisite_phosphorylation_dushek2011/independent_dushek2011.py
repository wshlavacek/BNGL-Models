"""Independent implementation of the Dushek et al. (2011) two-step model.

Level-1 verification for ``models/membrane_anchored_multisite_phosphorylation_dushek2011``.

The check available here is unusually strong, and it is worth saying why. Every reaction in
the model is a first-order transition of a single substrate molecule -- the kinase and
phosphatase are held in excess and enter only through the pseudo-first-order encounter rates
``ET*k_plus`` and ``FT*k_plus``, so no reaction has two reactants. The generated network has
**zero bimolecular reactions**, which makes the whole model a linear continuous-time Markov
chain on 7N+3 states. Its steady state is therefore not something to integrate towards: it is
the null space of the generator, available exactly from linear algebra with no solver
tolerance at all. That is the "exact witness" `stochastic-verification.md` §2 asks for, and it
beats integrating the same ODEs twice.

The chain is assembled here from Eq. 1 and Table 1 of the paper, state by state. It shares no
code with BioNetGen and never reads the generated ``.net`` -- which is the point: a rule that
matches the wrong context, or a site multiplicity attached to the wrong transition, changes
the chain and shows up in the comparison.

The state space, derived by hand:

===============  ==========================  =========================
state            meaning                     count
===============  ==========================  =========================
``free_j``       no enzyme in the encounter  N+1
                 complex, j sites phosphorylated
``E0_j``         kinase present and active,  N+1
                 not bound to a site
``Eb_j``         kinase bound to one of the  N     (j = 0..N-1)
                 N-j unphosphorylated sites
``E1_j``         kinase refractory           N     (j = 1..N; j=0 is
                                                   unreachable, since
                                                   only catalysis makes
                                                   an enzyme refractory
                                                   and it raises j)
``F0_j``         phosphatase present, active N+1
``Fb_j``         phosphatase bound           N     (j = 1..N)
``F1_j``         phosphatase refractory      N     (j = 0..N-1)
===============  ==========================  =========================

which totals 7N+3 = 143 at N = 20 -- the "143 coupled ODEs" the paper reports, and the count
BioNetGen independently produces. Dropping the two refractory blocks gives 5N+3 = 103, the
non-refractory variant.

Run this file directly for a self-test.
"""

from __future__ import annotations

import numpy as np

# Published parameter values. Fig. 2 caption and Table 1, except as noted.
KOFF = 1.0      # /s     unbinding rate
KR = 0.1        # /s     modification rate
MU = 1.0        # /s     reactivation rate, Fig. 2 C,D caption
A_ENC = 1e-4    # um^2   encounter area; caption prints A = 0.01^2, main text pi*s^2 = 1e-4
EF_TOT = 1000.0  # um^-2  total enzyme; stated in the SI Fig. S1 caption

# The two regimes of Fig. 2, as (k_plus, k_on) in um^2/s.
#
# The reaction-limited k_on is NOT the main text's k+/100. That value puts the N = 1 curve of
# Fig. 2 A at 0.313 -> 0.687 across the plotted range, where the published panel measures
# 0.476 at log10([E]/[F]) = -2; k_on = 1 gives 0.471. See verify_dushek2011.ipynb, which
# makes that measurement.
REGIMES = {
    "diffusion": (0.1, 10.0),
    "reaction": (10.0, 1.0),
}


def states(n_sites: int, refractory: bool = True):
    """Every state of the chain, as (block, j) pairs, in a fixed order."""
    s = [("free", j) for j in range(n_sites + 1)]
    s += [("E0", j) for j in range(n_sites + 1)] + [("Eb", j) for j in range(n_sites)]
    s += [("F0", j) for j in range(n_sites + 1)] + [("Fb", j) for j in range(1, n_sites + 1)]
    if refractory:
        s += [("E1", j) for j in range(1, n_sites + 1)]
        s += [("F1", j) for j in range(n_sites)]
    return s


def generator(n_sites: int, k_plus: float, k_on: float, ET: float, FT: float, *,
              koff: float = KOFF, kr: float = KR, mu: float = MU, area: float = A_ENC,
              refractory: bool = True):
    """The CTMC generator Q, built from Eq. 1. Rows are sources, columns destinations."""
    k_minus = k_plus / area          # local diffusion rate, Table 1: k- = k+/A
    k_on_star = k_on / area          # local single-site on-rate, Table 1: k*on = k_on/A
    S = states(n_sites, refractory)
    idx = {s: i for i, s in enumerate(S)}
    Q = np.zeros((len(S), len(S)))

    def add(a, b, rate):
        Q[idx[a], idx[b]] += rate
        Q[idx[a], idx[a]] -= rate

    for j in range(n_sites + 1):
        # forming the encounter complex, at X*k+ ; and leaving it, at k-
        add(("free", j), ("E0", j), ET * k_plus)
        add(("free", j), ("F0", j), FT * k_plus)
        add(("E0", j), ("free", j), k_minus)
        add(("F0", j), ("free", j), k_minus)
        # binding inside it. lambda = N-j free sites for the kinase, j for the phosphatase.
        if j < n_sites:
            add(("E0", j), ("Eb", j), (n_sites - j) * k_on_star)
        if j > 0:
            add(("F0", j), ("Fb", j), j * k_on_star)
    for j in range(n_sites):
        add(("Eb", j), ("E0", j), koff)
        add(("Eb", j), ("E1", j + 1) if refractory else ("E0", j + 1), kr)
    for j in range(1, n_sites + 1):
        add(("Fb", j), ("F0", j), koff)
        add(("Fb", j), ("F1", j - 1) if refractory else ("F0", j - 1), kr)
    if refractory:
        for j in range(1, n_sites + 1):
            add(("E1", j), ("E0", j), mu)         # reactivates in place
            add(("E1", j), ("free", j), k_minus)  # or diffuses away and is replaced
        for j in range(n_sites):
            add(("F1", j), ("F0", j), mu)
            add(("F1", j), ("free", j), k_minus)
    return Q, S


def steady_state(n_sites: int, k_plus: float, k_on: float, ET: float, FT: float, **kw):
    """Exact stationary distribution: the null space of Q^T with the probabilities summing to 1."""
    Q, S = generator(n_sites, k_plus, k_on, ET, FT, **kw)
    n = len(S)
    M = np.vstack([Q.T, np.ones(n)])
    b = np.zeros(n + 1)
    b[-1] = 1.0
    p, *_ = np.linalg.lstsq(M, b, rcond=None)
    return p, S


def mean_phosphorylation(n_sites: int, regime: str = "diffusion", r: float = 1.0, *,
                         ef_tot: float = EF_TOT, refractory: bool = True, **kw) -> float:
    """The paper's Eq. 2: total phosphorylation normalised by N*S_T, at [E]/[F] = `r`."""
    k_plus, k_on = REGIMES[regime]
    ET = ef_tot * r / (1.0 + r)
    FT = ef_tot / (1.0 + r)
    p, S = steady_state(n_sites, k_plus, k_on, ET, FT, refractory=refractory, **kw)
    return float(sum(p[i] * s[1] for i, s in enumerate(S)) / n_sites)


def dose_response(n_sites: int, regime: str = "diffusion", *, refractory: bool = True,
                  lo: float = -2.0, hi: float = 2.0, n: int = 81, **kw):
    """(log10 [E]/[F], <S>) over the range Fig. 2 plots."""
    x = np.linspace(lo, hi, n)
    return x, np.array([mean_phosphorylation(n_sites, regime, 10.0**t,
                                             refractory=refractory, **kw) for t in x])


def hill_number(n_sites: int, regime: str = "diffusion", *, refractory: bool = True,
                **kw) -> float:
    """Hill number of the dose-response, by a free-asymptote Hill fit.

    The asymptotes are fitted rather than taken from the window because none of these curves
    saturates within the four decades plotted. The fit reproduces the published insets for
    N >= 2; at N = 1 it is not meaningful, because the N = 1 dose-response is not a Hill
    function at all -- a Hill fit leaves a residual of 0.021 on a curve whose total span is
    0.057 -- so this returns NaN there rather than a number that would look like a result.
    """
    from scipy.optimize import curve_fit
    if n_sites < 2:
        return float("nan")
    x, y = dose_response(n_sites, regime, refractory=refractory, n=161, **kw)
    r = 10.0**x

    def f(rr, lo_, hi_, K, nh):
        return lo_ + (hi_ - lo_) * rr**nh / (K**nh + rr**nh)

    p, _ = curve_fit(f, r, y, p0=[y[0], y[-1], 1.0, 1.0],
                     bounds=([-1, -1, 1e-6, 1e-3], [2, 2, 1e6, 50]), maxfev=400000)
    return float(p[3])


if __name__ == "__main__":
    print(f"{'N':>3} {'states':>7}  " + "  ".join(f"{k:>12}" for k in
                                                  ("2A", "2B", "2C", "2D")))
    cfg = [("reaction", False), ("diffusion", False), ("reaction", True), ("diffusion", True)]
    for N in (1, 2, 5, 10, 20):
        vals = [mean_phosphorylation(N, reg, 100.0, refractory=refr) for reg, refr in cfg]
        print(f"{N:>3} {len(states(N)):>7}  " + "  ".join(f"{v:12.6f}" for v in vals))
    print("\n<S> at [E]/[F] = 100 per panel; Hill numbers (N >= 2):")
    for N in (2, 5, 10, 20):
        h = [hill_number(N, reg, refractory=refr) for reg, refr in cfg]
        print(f"  N={N:>2}  " + "  ".join(f"{v:6.3f}" for v in h)
              + f"   published Fig. 2D fit {0.58 + 0.24 * N:6.3f}")
