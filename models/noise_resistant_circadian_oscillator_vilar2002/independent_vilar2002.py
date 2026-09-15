#!/usr/bin/env python3
"""Independent implementation of the Vilar et al. (2002) circadian oscillator.

Transcribed from the paper: the rate equations are Eq. 1 as printed, and the
reaction list is read off the Fig. 1 network diagram and the terms of Eq. 1.
Nothing here parses or reuses the BioNetGen-generated network, and it shares no
code with BioNetGen. That is the whole point of verification level 1: a check
that starts from the `.net` file verifies the integrator, not the model, and
would agree with a mis-transcribed rule just as happily as a correct one.

Two checks live here, per `references/stochastic-verification.md` §2:

* `rhs` -- the deterministic limit, integrated with SciPy. Cheap, and it catches
  transcription errors in the rules, the rate laws and the initial conditions
  independently of any noise.
* `gillespie` -- a direct-method simulator built from the sixteen reactions. This
  is the strongest check the model admits, because it tests the stochastic
  mechanism itself rather than its mean.

Propensity conventions, stated because they are not all obvious by inspection:

* Every bimolecular step here is heterodimeric -- A with a promoter, A with R --
  so each propensity is k*x*y with no symmetry factor. There is no homodimeric
  step, and so none of the factor-of-two ambiguity that Samoilov's driver has.
* `delta_A` appears on two reactions. Free A degrades (reaction 13) and the
  complex C decays by degradation of the A inside it, releasing R (reaction 16).
  Eq. 1 writes the second as the `+delta_A*C` term in dR/dt and `-delta_A*C` in
  dC/dt.
* Promoter-bound activator does not degrade: Eq. 1's dD_A'/dt has no delta term.
  Repressor inside the complex does not degrade at delta_R either.
"""

from __future__ import annotations

import numpy as np

# Fig. 1 caption. Rates are per hour; the two gammas and gamma_C are per molecule
# per hour. The cell volume is unity, so counts and concentrations coincide.
NOMINAL = dict(
    alpha_A=50.0,
    alpha_A_act=500.0,
    alpha_R=0.01,
    alpha_R_act=50.0,
    beta_A=50.0,
    beta_R=5.0,
    delta_MA=10.0,
    delta_MR=0.5,
    delta_A=1.0,
    delta_R=0.2,
    gamma_A=1.0,
    gamma_R=1.0,
    gamma_C=2.0,
    theta_A=50.0,
    theta_R=100.0,
)

# Fig. 5: the only change is a fourfold slower repressor degradation, which turns
# the trace of Eq. 3 negative and makes the fixed point stable.
STABLE_FIXED_POINT = dict(NOMINAL, delta_R=0.05)

# Fig. 7: both translation rates and both transcript degradation rates times 100.
LOW_MRNA = dict(NOMINAL, beta_A=5000.0, beta_R=500.0, delta_MA=1000.0, delta_MR=50.0)

PARAMETERIZATIONS = {
    "nominal": NOMINAL,
    "stable_fixed_point": STABLE_FIXED_POINT,
    "low_mrna": LOW_MRNA,
}

# State order, following Eq. 1
SPECIES = ("D_A", "D_R", "D_Ap", "D_Rp", "M_A", "A", "M_R", "R", "C")
IDX = {s: i for i, s in enumerate(SPECIES)}


def initial_state() -> np.ndarray:
    """D_A = D_R = 1 molecule, everything else zero (Fig. 1 caption)."""
    y = np.zeros(len(SPECIES))
    y[IDX["D_A"]] = 1.0
    y[IDX["D_R"]] = 1.0
    return y


def rhs(t: float, y: np.ndarray, p: dict) -> np.ndarray:
    """Eq. 1 of Vilar et al. (2002), transcribed term by term."""
    D_A, D_R, D_Ap, D_Rp, M_A, A, M_R, R, C = y
    dD_A = p["theta_A"] * D_Ap - p["gamma_A"] * D_A * A
    dD_R = p["theta_R"] * D_Rp - p["gamma_R"] * D_R * A
    dD_Ap = p["gamma_A"] * D_A * A - p["theta_A"] * D_Ap
    dD_Rp = p["gamma_R"] * D_R * A - p["theta_R"] * D_Rp
    dM_A = p["alpha_A_act"] * D_Ap + p["alpha_A"] * D_A - p["delta_MA"] * M_A
    dA = (
        p["beta_A"] * M_A
        + p["theta_A"] * D_Ap
        + p["theta_R"] * D_Rp
        - A * (p["gamma_A"] * D_A + p["gamma_R"] * D_R + p["gamma_C"] * R + p["delta_A"])
    )
    dM_R = p["alpha_R_act"] * D_Rp + p["alpha_R"] * D_R - p["delta_MR"] * M_R
    dR = p["beta_R"] * M_R - p["gamma_C"] * A * R + p["delta_A"] * C - p["delta_R"] * R
    dC = p["gamma_C"] * A * R - p["delta_A"] * C
    return np.array([dD_A, dD_R, dD_Ap, dD_Rp, dM_A, dA, dM_R, dR, dC])


def integrate(p: dict, t_end: float, n_steps: int):
    """Deterministic solution on the same grid the BNGL protocol samples."""
    from scipy.integrate import solve_ivp

    t_eval = np.linspace(0.0, t_end, n_steps + 1)
    sol = solve_ivp(
        rhs,
        (0.0, t_end),
        initial_state(),
        t_eval=t_eval,
        args=(p,),
        method="LSODA",
        rtol=1e-10,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    return sol.t, sol.y.T


# Sixteen reactions as (stoichiometry, propensity). Stoichiometry is a tuple of
# (species index, change); the propensity closure reads the parameter dict.
_A, _R, _C = IDX["A"], IDX["R"], IDX["C"]
_DA, _DR, _DAp, _DRp = IDX["D_A"], IDX["D_R"], IDX["D_Ap"], IDX["D_Rp"]
_MA, _MR = IDX["M_A"], IDX["M_R"]

REACTIONS = [
    # 1-2: activator binds and unbinds its own promoter
    (((_A, -1), (_DA, -1), (_DAp, +1)), lambda y, p: p["gamma_A"] * y[_DA] * y[_A]),
    (((_A, +1), (_DA, +1), (_DAp, -1)), lambda y, p: p["theta_A"] * y[_DAp]),
    # 3-4: activator binds and unbinds the repressor promoter
    (((_A, -1), (_DR, -1), (_DRp, +1)), lambda y, p: p["gamma_R"] * y[_DR] * y[_A]),
    (((_A, +1), (_DR, +1), (_DRp, -1)), lambda y, p: p["theta_R"] * y[_DRp]),
    # 5-8: basal and activated transcription of both genes
    (((_MA, +1),), lambda y, p: p["alpha_A"] * y[_DA]),
    (((_MA, +1),), lambda y, p: p["alpha_A_act"] * y[_DAp]),
    (((_MR, +1),), lambda y, p: p["alpha_R"] * y[_DR]),
    (((_MR, +1),), lambda y, p: p["alpha_R_act"] * y[_DRp]),
    # 9-10: translation
    (((_A, +1),), lambda y, p: p["beta_A"] * y[_MA]),
    (((_R, +1),), lambda y, p: p["beta_R"] * y[_MR]),
    # 11-12: transcript degradation
    (((_MA, -1),), lambda y, p: p["delta_MA"] * y[_MA]),
    (((_MR, -1),), lambda y, p: p["delta_MR"] * y[_MR]),
    # 13-14: degradation of free protein
    (((_A, -1),), lambda y, p: p["delta_A"] * y[_A]),
    (((_R, -1),), lambda y, p: p["delta_R"] * y[_R]),
    # 15: sequestration of the activator, irreversible in Eq. 1
    (((_A, -1), (_R, -1), (_C, +1)), lambda y, p: p["gamma_C"] * y[_A] * y[_R]),
    # 16: the complex decays because the activator in it degrades, freeing R
    (((_C, -1), (_R, +1)), lambda y, p: p["delta_A"] * y[_C]),
]


def gillespie(p: dict, t_end: float, n_steps: int, seed: int):
    """Direct-method SSA on the sixteen reactions, sampled on a uniform grid."""
    rng = np.random.default_rng(seed)
    grid = np.linspace(0.0, t_end, n_steps + 1)
    out = np.zeros((grid.size, len(SPECIES)))

    y = initial_state()
    t = 0.0
    nxt = 0
    props = np.zeros(len(REACTIONS))
    while True:
        for i, (_, f) in enumerate(REACTIONS):
            props[i] = f(y, p)
        total = props.sum()
        if total <= 0.0:
            break
        t_new = t + rng.exponential(1.0 / total)
        while nxt < grid.size and grid[nxt] < t_new:
            out[nxt] = y
            nxt += 1
        if nxt >= grid.size:
            return grid, out
        t = t_new
        k = int(np.searchsorted(np.cumsum(props), rng.random() * total))
        for idx, change in REACTIONS[k][0]:
            y[idx] += change
    while nxt < grid.size:
        out[nxt] = y
        nxt += 1
    return grid, out


def peaks(t: np.ndarray, y: np.ndarray, t_from: float = 100.0):
    """Peak times and heights of an oscillation, after a burn-in.

    A peak is a local maximum above the midpoint of the post-burn-in range, taken
    one per excursion: the series is thresholded at the midpoint and the maximum
    within each contiguous run above it is returned. Thresholding rather than
    derivative-testing is what makes this usable on a stochastic trace, where
    every cycle carries dozens of one-molecule reversals that a naive local
    maximum test would count as peaks.
    """
    m = t >= t_from
    tt, yy = t[m], y[m]
    if yy.size == 0 or yy.max() <= yy.min():
        return np.array([]), np.array([])
    thr = 0.5 * (yy.max() + yy.min())
    above = yy > thr
    edges = np.diff(above.astype(int))
    starts = np.where(edges == 1)[0] + 1
    ends = np.where(edges == -1)[0] + 1
    if above[0]:
        starts = np.r_[0, starts]
    if above[-1]:
        ends = np.r_[ends, yy.size]
    n = min(starts.size, ends.size)
    pt, ph = [], []
    for s, e in zip(starts[:n], ends[:n]):
        j = s + int(np.argmax(yy[s:e]))
        pt.append(tt[j])
        ph.append(yy[j])
    return np.array(pt), np.array(ph)


def period(t: np.ndarray, y: np.ndarray, t_from: float = 100.0) -> float:
    pt, _ = peaks(t, y, t_from)
    return float(np.mean(np.diff(pt))) if pt.size > 1 else float("nan")
