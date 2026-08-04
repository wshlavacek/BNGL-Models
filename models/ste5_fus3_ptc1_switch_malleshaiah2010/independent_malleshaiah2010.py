#!/usr/bin/env python3
"""Independent ODE implementation of the Malleshaiah et al. model equations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Layout:
    S: slice
    K0: slice
    K1: slice
    K2: slice
    P0: slice
    P1: slice
    P2: slice
    A: slice
    Kfree: int
    Pfree: int
    Afree: int
    size: int


def layout(n_sites: int) -> Layout:
    cursor = 0

    def take(length: int) -> slice:
        nonlocal cursor
        result = slice(cursor, cursor + length)
        cursor += length
        return result

    S = take(n_sites + 1)
    K0 = take(n_sites + 1)
    K1 = take(n_sites)
    K2 = take(n_sites)
    P0 = take(n_sites + 1)
    P1 = take(n_sites)
    P2 = take(n_sites)
    A = take(n_sites)
    Kfree, Pfree, Afree = cursor, cursor + 1, cursor + 2
    return Layout(S, K0, K1, K2, P0, P1, P2, A, Kfree, Pfree, Afree, cursor + 3)


def inputs(alpha_nM: float) -> tuple[float, float, float]:
    ptc1 = 1.2 + 39.0 * alpha_nM**2.3 / (alpha_nM**2.3 + 240.0**2.3)
    active = 5.8 * alpha_nM**1.3 / (alpha_nM**1.3 + 1680.0**1.3)
    inactive = 197.0 - active
    return inactive, ptc1, active


def rhs_factory(n_sites: int, molar_association_units: bool):
    ix = layout(n_sites)
    scale = 1e-9 if molar_association_units else 1.0
    f1_P, f1_K, f4_K = 186000.0 * scale, 12000.0 * scale, 109000.0 * scale
    f2_P, f3_P, b1_P, b2_P, k_P = 327.0, 0.3, 22.0, 0.12, 0.5
    f2_K, f3_K, b6_K, k_K = 850.0, 0.1, 24.0, 1.13
    b_K = np.array([99.0, 42.0, 21.0, 13.0, 10.0])

    def rhs(_time: float, y: np.ndarray) -> np.ndarray:
        dy = np.zeros_like(y)
        S, K0, K1, K2 = y[ix.S], y[ix.K0], y[ix.K1], y[ix.K2]
        P0, P1, P2, A = y[ix.P0], y[ix.P1], y[ix.P2], y[ix.A]

        for n in range(n_sites + 1):
            flux = f1_K * y[ix.Kfree] * S[n]
            dy[ix.Kfree] -= flux
            dy[ix.S.start + n] -= flux
            dy[ix.K0.start + n] += flux

            flux = b_K[n] * K0[n]
            dy[ix.Kfree] += flux
            dy[ix.S.start + n] += flux
            dy[ix.K0.start + n] -= flux

            flux = f1_P * y[ix.Pfree] * S[n]
            dy[ix.Pfree] -= flux
            dy[ix.S.start + n] -= flux
            dy[ix.P0.start + n] += flux

            flux = b1_P * P0[n]
            dy[ix.Pfree] += flux
            dy[ix.S.start + n] += flux
            dy[ix.P0.start + n] -= flux

        for n in range(n_sites):
            flux = (n_sites - n) * f2_K * K0[n]
            dy[ix.K0.start + n] -= flux
            dy[ix.K1.start + n] += flux

            flux = b6_K * K1[n]
            dy[ix.K0.start + n] += flux
            dy[ix.K1.start + n] -= flux

            flux = b_K[n] * K1[n]
            dy[ix.K1.start + n] -= flux
            dy[ix.K2.start + n] += flux

            flux = f3_K * K2[n]
            dy[ix.K1.start + n] += flux
            dy[ix.K2.start + n] -= flux

            flux = k_K * K1[n]
            dy[ix.K1.start + n] -= flux
            dy[ix.K0.start + n + 1] += flux

            flux = k_K * K2[n]
            dy[ix.K2.start + n] -= flux
            dy[ix.Kfree] += flux
            dy[ix.S.start + n + 1] += flux

            flux = (n_sites - n) * f4_K * y[ix.Afree] * S[n]
            dy[ix.Afree] -= flux
            dy[ix.S.start + n] -= flux
            dy[ix.A.start + n] += flux

            flux = b6_K * A[n]
            dy[ix.Afree] += flux
            dy[ix.S.start + n] += flux
            dy[ix.A.start + n] -= flux

            flux = k_K * A[n]
            dy[ix.Afree] += flux
            dy[ix.A.start + n] -= flux
            dy[ix.S.start + n + 1] += flux

        for n in range(1, n_sites + 1):
            j = n - 1
            flux = n * f2_P * P0[n]
            dy[ix.P0.start + n] -= flux
            dy[ix.P1.start + j] += flux

            flux = b2_P * P1[j]
            dy[ix.P0.start + n] += flux
            dy[ix.P1.start + j] -= flux

            flux = b1_P * P1[j]
            dy[ix.P1.start + j] -= flux
            dy[ix.P2.start + j] += flux

            flux = f3_P * P2[j]
            dy[ix.P1.start + j] += flux
            dy[ix.P2.start + j] -= flux

            flux = k_P * P1[j]
            dy[ix.P1.start + j] -= flux
            dy[ix.P0.start + n - 1] += flux

            flux = k_P * P2[j]
            dy[ix.P2.start + j] -= flux
            dy[ix.Pfree] += flux
            dy[ix.S.start + n - 1] += flux

        return dy

    return ix, rhs


def conserved_totals(y: np.ndarray, n_sites: int) -> np.ndarray:
    ix = layout(n_sites)
    ste5 = sum(y[part].sum() for part in (ix.S, ix.K0, ix.K1, ix.K2, ix.P0,
                                         ix.P1, ix.P2, ix.A))
    kinase = y[ix.Kfree] + y[ix.K0].sum() + y[ix.K1].sum() + y[ix.K2].sum()
    ptc1 = y[ix.Pfree] + y[ix.P0].sum() + y[ix.P1].sum() + y[ix.P2].sum()
    active = y[ix.Afree] + y[ix.A].sum()
    return np.array([ste5, kinase, ptc1, active])


def mean_phosphorylation(y: np.ndarray, n_sites: int) -> float:
    ix = layout(n_sites)
    by_state = np.zeros(n_sites + 1)
    by_state += y[ix.S] + y[ix.K0] + y[ix.P0]
    by_state[:-1] += y[ix.K1] + y[ix.K2] + y[ix.A]
    by_state[1:] += y[ix.P1] + y[ix.P2]
    return float(np.arange(n_sites + 1) @ by_state / 52.0)


def simulate(alpha_nM: float, n_sites: int, molar_association_units: bool):
    ix, rhs = rhs_factory(n_sites, molar_association_units)
    y0 = np.zeros(ix.size)
    inactive, ptc1, active = inputs(alpha_nM)
    y0[ix.S.start] = 52.0
    y0[ix.Kfree], y0[ix.Pfree], y0[ix.Afree] = inactive, ptc1, active
    solution = solve_ivp(
        rhs,
        (0.0, 100000.0),
        y0,
        method="LSODA",
        rtol=2e-10,
        atol=1e-11,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    final = solution.y[:, -1]
    conservation_error = np.max(
        np.abs(conserved_totals(final, n_sites) - conserved_totals(y0, n_sites))
    )
    return mean_phosphorylation(final, n_sites), conservation_error
