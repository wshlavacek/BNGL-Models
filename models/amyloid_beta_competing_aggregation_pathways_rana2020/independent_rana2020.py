"""Independent NumPy/SciPy implementation of the Rana et al. (2020) competing-pathway model.

Written from the paper's equations, not from BioNetGen's generated network, so that
comparing it against BioNetGen output independently checks the network generation, the
molecule-deletion and molecule-synthesis rules, the statistical factor BioNetGen applies to
the four identical A_1 patterns of the condensation rule, and the observable definitions.

Reference
---------
Rana P, Bose P, Vaidya A, Rangachari V, Ghosh P (2020). Global fitting and parameter
identifiability for amyloid-beta aggregation with competing pathways. 2020 IEEE 20th
International Conference on BioInformatics and BioEngineering (BIBE), pp. 73-78.
doi:10.1109/BIBE50027.2020.00020

Equations transcribed
---------------------
Eq. 4, on pathway
    A_i + A_1 <-> A_i+1                 i = 1..11
    A_i + F   <-> F                     i = 1..11,  F == A_12
Eq. 4-II, off pathway
    4 A_1 + L    <-> A'_4
    A'_i + A_1   <-> A'_i+1             i = 4..11
    A'_12 + A'_i <-> F'_1               i = 4..11
Eq. 5, switching
    A_4 <-> A'_4
Fluxes (the paper's own rate laws)
    H_i   = k_nuon [A_i][A_1] - k_nuon_ [A_i+1]                  i = 1..11
    I_i   = k_fbon [A_i][F]   - k_fbon_ [F]                      i = 1..11
    G'_1  = k_con  [A_1]^4[L] - k_con_  [A'_4]
    H'_i  = k_nuoff [A'_3+i][A_1]   - k_nuoff_ [A'_4+i]          i = 1..8
    I'_i  = k_fboff [A'_3+i][A'_12] - k_fboff_ [F'_1]            i = 1..8
    J     = k_swi [A_4] - k_swi_ [A'_4]

Stoichiometry follows directly: A_1 is consumed twice by H_1 (both partners are monomer)
and four times by G'_1; F is not consumed by I_i, because the published reaction has the
fibril on both sides; and total peptide is therefore not conserved, by construction.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

# state layout: A_1..A_12 (0..11), A'_4..A'_12 (12..20), F'_1 (21), L (22)
N_STATE = 23
I_FP1 = 21
I_L = 22
_W = np.arange(1, 12)

PARAM_NAMES = ("k_nuon", "k_nuon_", "k_fbon", "k_fbon_", "k_con", "k_con_",
               "k_nuoff", "k_nuoff_", "k_fboff", "k_fboff_", "k_swi", "k_swi_")


def rhs(t, y, p, off_pathway=True):
    A = y[0:12]
    Ap = y[12:21]
    A1 = A[0]
    F = A[11]

    H = p["k_nuon"] * A[0:11] * A1 - p["k_nuon_"] * A[1:12]
    I = p["k_fbon"] * A[0:11] * F - p["k_fbon_"] * F

    d = np.zeros(N_STATE)
    d[0] = -H.sum() - H[0] - I[0]
    d[1:11] = H[0:10] - H[1:11] - I[1:11]
    d[11] = H[10]

    if off_pathway:
        Fp1, L = y[I_FP1], y[I_L]
        G1 = p["k_con"] * A1**4 * L - p["k_con_"] * Ap[0]
        Hp = p["k_nuoff"] * Ap[0:8] * A1 - p["k_nuoff_"] * Ap[1:9]
        Ip = p["k_fboff"] * Ap[0:8] * Ap[8] - p["k_fboff_"] * Fp1
        J = p["k_swi"] * A[3] - p["k_swi_"] * Ap[0]

        d[0] += -4.0 * G1 - Hp.sum()
        d[3] -= J
        d[12] = G1 - Hp[0] - Ip[0] + J
        d[13:20] = Hp[0:7] - Hp[1:8] - Ip[1:8]
        d[20] = Hp[7] - Ip.sum()
        d[I_FP1] = Ip.sum()
        d[I_L] = -G1
    return d


def y0(abeta=25.0, micelle=0.0):
    y = np.zeros(N_STATE)
    y[0] = abeta
    y[I_L] = micelle
    return y


def integrate(p, y_init, t_eval, *, off_pathway=True, rtol=1e-8, atol=1e-14):
    """Integrate on the given grid. t_eval must start at 0."""
    t_eval = np.asarray(t_eval, float)
    t_end = float(max(t_eval[-1], 1e-9))
    sol = solve_ivp(rhs, (0.0, t_end), np.asarray(y_init, float), t_eval=t_eval,
                    args=(p, off_pathway), method="LSODA", rtol=rtol, atol=atol)
    if not sol.success or sol.y.shape[1] != len(t_eval):
        raise RuntimeError(sol.message)
    return sol.y.T


def protocol(p, t_meas, *, micelle_before, micelle_after, t_event=0.0,
             abeta=25.0, off_pathway=True, **kw):
    """The two-phase design the PyBNF confs express as preequilibrate + condition.

    Phase 1 runs for exactly ``t_event`` hours with ``micelle_before`` pseudo-micelle, is
    not measured, and hands its end state to phase 2, whose clock restarts at zero and
    whose free micelle pool is reset to ``micelle_after``. ``t_event = 0`` degenerates to a
    single measured phase, which is the no-event protocol.
    """
    y = y0(abeta, micelle_before)
    if t_event > 0:
        y = integrate(p, y, np.array([0.0, t_event]), off_pathway=off_pathway, **kw)[-1]
    if off_pathway:
        y = y.copy()
        y[I_L] = micelle_after
    grid = np.asarray(t_meas, float)
    if grid[0] > 0:
        grid = np.concatenate([[0.0], grid])
        return integrate(p, y, grid, off_pathway=off_pathway, **kw)[1:]
    return integrate(p, y, grid, off_pathway=off_pathway, **kw)


def tht(states, map_on, map_off=0.0):
    """map_on*[F] + map_off*[F'_1] — the paper's signal_on + signal_off."""
    return map_on * states[:, 11] + map_off * states[:, I_FP1]
