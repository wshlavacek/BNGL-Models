#!/usr/bin/env python3
"""Independent isotopomer balance equations for the four-flux network of Mu et al. (2010).

This module deliberately does NOT go through BioNetGen or through the BNGL carbon
fate maps. It implements the *other* method the chapter reviews, the isotopomer
mapping matrix (IMM) formalism of Sec. 15.2.1 and Sec. 15.3.2, so that the curated
BNGL model can be checked against a construction that shares nothing with it beyond
the published network.

The state of a two-carbon metabolite is its isotopomer distribution vector (IDV) of
Eq. 15.14, ordered (00, 01, 10, 11) where the first digit is the labeling state of C1
and the second that of C2, so the index of an isotopomer is 2*c1 + c2. For a reaction
A -> B carrying flux r, Eq. 15.1 gives

    V_B dI_B/dt = r * (IMM_{A>B} . I_A)

and a metabolite's balance is the sum of such terms less its own consumption. The
four IMMs of the network are Eqs. 15.15 to 15.18; only IMM_{S>M2} is not the identity,
because the S -> M2 step transposes the two carbons. Note that Eqs. 15.27 to 15.30 as
printed omit that transposition, while Eq. 15.16, Fig. 15.3, and the BNGL listing of
Sec. 15.3.3 all carry it; without it a feed labeled at C1 could never label C2 of the
product and the whole example would be empty, so Eq. 15.16 is the one implemented here.

Writing the balances on fractions rather than on amounts is exact for this network:
each pool's inflow equals its outflow for every flux set used in Fig. 15.4, so the
pool sizes are constant and amount = pool size * fraction throughout.

The initial condition follows the chapter's description of the experiment as a step
change of the input feed at t = 0. S is that input, and Eq. 15.4 makes v_S = v1 + v2,
so Eqs. 15.19 to 15.22 reduce to dS/dt = 0: S carries the feed composition from t = 0
onward while the internal pools M1, M2, and P start unlabeled. `simulate_lagged_S`
integrates the other reading, in which S starts unlabeled and relaxes to the feed, so
the two can be compared against the digitized figure.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

# Isotopomer index of a labeling pattern, Eq. 15.14 ordering.
I00, I01, I10, I11 = 0, 1, 2, 3

# Eq. 15.15, S -> M1: carbon order preserved.
IMM_S_M1 = np.eye(4)

# Eq. 15.16, S -> M2: the two carbons are transposed, so S01 <-> M2_10.
IMM_S_M2 = np.zeros((4, 4))
IMM_S_M2[I00, I00] = 1.0
IMM_S_M2[I10, I01] = 1.0
IMM_S_M2[I01, I10] = 1.0
IMM_S_M2[I11, I11] = 1.0

# Eqs. 15.17 and 15.18, M1 -> P and M2 -> P: carbon order preserved.
IMM_M1_P = np.eye(4)
IMM_M2_P = np.eye(4)

# Flux and pool settings of the four curves of Fig. 15.4. The caption fixes
# V_S = V_M1 = V_P = 1 throughout and varies V_M2; the dotted curve instead
# redistributes flux between the two branches at uniform pool sizes.
BASE = dict(
    V_S=1.0, V_M1=1.0, V_M2=1.0, V_P=1.0,
    v_vs=1.0, v_v1=0.5, v_v2=0.5, v_v3=0.4, v_v4=0.4,
    v_vm1=0.1, v_vm2=0.1, v_vp=0.8, f_label=0.1,
)
CONDITIONS = {
    "ode": BASE,
    "ode_m2low": {**BASE, "V_M2": 0.5},
    "ode_m2high": {**BASE, "V_M2": 2.0},
    "ode_altflux": {**BASE, "v_v1": 0.6, "v_v2": 0.4, "v_v3": 0.5, "v_v4": 0.3},
}


def check_balances(p: dict) -> None:
    """Assert the stoichiometric balances of Eqs. 15.4 to 15.7."""
    assert np.isclose(p["v_vs"], p["v_v1"] + p["v_v2"]), "Eq. 15.4"
    assert np.isclose(p["v_v1"], p["v_vm1"] + p["v_v3"]), "Eq. 15.5"
    assert np.isclose(p["v_v2"], p["v_vm2"] + p["v_v4"]), "Eq. 15.6"
    assert np.isclose(p["v_vp"], p["v_v3"] + p["v_v4"]), "Eq. 15.7"


def feed_idv(f_label: float) -> np.ndarray:
    """Feed composition after the t = 0 switch: 1 - f_label unlabeled, f_label at C1."""
    idv = np.zeros(4)
    idv[I00] = 1.0 - f_label
    idv[I10] = f_label
    return idv


def rhs(_t: float, y: np.ndarray, p: dict) -> np.ndarray:
    """Isotopomer balances of Eqs. 15.19 to 15.34, written on fractions."""
    S, M1, M2, P = y[0:4], y[4:8], y[8:12], y[12:16]
    F = feed_idv(p["f_label"])
    dS = (p["v_vs"] * F - (p["v_v1"] + p["v_v2"]) * S) / p["V_S"]
    dM1 = (p["v_v1"] * (IMM_S_M1 @ S) - (p["v_v3"] + p["v_vm1"]) * M1) / p["V_M1"]
    dM2 = (p["v_v2"] * (IMM_S_M2 @ S) - (p["v_v4"] + p["v_vm2"]) * M2) / p["V_M2"]
    dP = (p["v_v3"] * (IMM_M1_P @ M1) + p["v_v4"] * (IMM_M2_P @ M2)
          - p["v_vp"] * P) / p["V_P"]
    return np.concatenate([dS, dM1, dM2, dP])


def _integrate(condition: str, t_eval: np.ndarray, y0: np.ndarray) -> dict[str, np.ndarray]:
    p = CONDITIONS[condition]
    check_balances(p)
    # Always integrate from t = 0, the instant of the feed switch, whatever the
    # first requested output time is.
    assert t_eval[0] >= 0.0
    sol = solve_ivp(
        rhs, (0.0, float(t_eval[-1])), y0, t_eval=t_eval,
        args=(p,), method="LSODA", rtol=1e-11, atol=1e-13,
    )
    assert sol.success, sol.message
    y = sol.y
    return {
        "t": sol.t,
        "S": y[0:4], "M1": y[4:8], "M2": y[8:12], "P": y[12:16],
        "P00": y[12 + I00], "P01": y[12 + I01],
        "P10": y[12 + I10], "P11": y[12 + I11],
    }


def simulate(condition: str, t_eval: np.ndarray) -> dict[str, np.ndarray]:
    """Integrate one Fig. 15.4 condition and return the four IDVs plus P's fractions.

    S carries the feed composition from t = 0; M1, M2, and P start unlabeled.
    """
    unlabeled = np.array([1.0, 0.0, 0.0, 0.0])
    f = feed_idv(CONDITIONS[condition]["f_label"])
    return _integrate(condition, t_eval, np.concatenate([f, unlabeled, unlabeled, unlabeled]))


def simulate_lagged_S(condition: str, t_eval: np.ndarray) -> dict[str, np.ndarray]:
    """The alternative reading, in which S also starts unlabeled.

    Everything else is identical, so the difference between this and `simulate` is
    exactly the first-order lag of V_S/(v1 + v2) that S then imposes on the branches.
    """
    return _integrate(condition, t_eval, np.tile(np.array([1.0, 0.0, 0.0, 0.0]), 4))


def steady_state_P(condition: str) -> np.ndarray:
    """Closed-form isotopic steady state of the product IDV.

    With the balances of Eqs. 15.4 to 15.7 satisfied, S relaxes to the feed
    composition, each intermediate inherits it through its own IMM, and
    v_p * P = v_3 * M1 + v_4 * M2. Pool sizes cancel, which is the chapter's
    statement that pool size does not affect the steady-state distribution.
    """
    p = CONDITIONS[condition]
    check_balances(p)
    S = feed_idv(p["f_label"])
    M1 = IMM_S_M1 @ S
    M2 = IMM_S_M2 @ S
    return (p["v_v3"] * (IMM_M1_P @ M1) + p["v_v4"] * (IMM_M2_P @ M2)) / p["v_vp"]


if __name__ == "__main__":
    t = np.linspace(0.0, 20.0, 41)
    for name in CONDITIONS:
        out = simulate(name, t)
        ss = steady_state_P(name)
        s_drift = np.max(np.abs(out["S"] - out["S"][:, :1]))
        print(f"{name:12s} P01(20) = {out['P01'][-1]:.6f}   "
              f"P01(inf) = {ss[I01]:.6f}   sum(P) = {out['P'][:, -1].sum():.12f}   "
              f"max|S(t)-S(0)| = {s_drift:.1e}")
