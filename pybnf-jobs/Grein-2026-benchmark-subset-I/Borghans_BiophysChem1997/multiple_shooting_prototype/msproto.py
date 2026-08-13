"""Multiple-shooting prototype for PyBNF issue #563, on Borghans_BiophysChem1997.

Deliberately outside PyBNF's fit machinery: this answers "does multiple shooting
enlarge the useful convergence region on this problem" before any of it is built.

Transcription
-------------
Split [0, T] at m-1 interior knots.  Segment j is integrated from its own start
state; segment 0's start state is the model's own (init_Z/Y/A, which are fitted
parameters of the PEtab problem), and each interior knot carries an auxiliary
state z_j.  Continuity is imposed as

    c_j = Phi_{j-1}(z_{j-1}, theta) - z_j = 0,      j = 1 .. m-1

and enforced by an augmented Lagrangian whose inner problem is a bound-constrained
Gauss-Newton least-squares solve (scipy `trf` with an analytic Jacobian -- the same
step math PyBNF's `gntr` runs, so the prototype's inner optimizer is not a
confound).

Objective
---------
The PEtab problem's own: lognormal (log10-additive Gaussian) noise on
Ca = Z_state * scale + offset, with the noise scale sigma ESTIMATED.  Following
issue #563's point 4 (and PyBNF #562 / ADR-0108) sigma is profiled out
analytically FROM THE DATA TERM ONLY, so continuity violation can never be
absorbed into the reported noise scale:

    S(x)      = sum_i (log10 y_i - log10 mu_i)^2
    sigma^2   = S / n
    J_reduced = n/2 * log(S/n) + n/2
    J_paper   = J_reduced + n/2 log(2 pi) + sum_i log(y_i ln 10)   (= -log L)
    OG        = J_paper - J*,     solved iff OG < 1.92

Certified reconstruction
------------------------
Every reported number comes from discarding z and re-simulating theta with
ordinary single shooting.  A multiple-shooting run that leaves continuity
unconverged therefore scores as what it actually is.
"""

from __future__ import annotations

import os
import time
import warnings

import numpy as np
from scipy.optimize import least_squares

import bngsim

# --------------------------------------------------------------------------- #
# Problem constants
# --------------------------------------------------------------------------- #
JOB = ('/Users/l119605/Code/BNGL-Models/pybnf-jobs/Grein-2026-benchmark-subset-I/'
       'Borghans_BiophysChem1997')
MODEL_XML = os.path.join(JOB, 'model_Borghans_BiophysChem1997.xml')
DATA_EXP = os.path.join(JOB, 'experiment1.exp')
JSTAR = -132.00847649739424
SOLVED_OG = 1.92

STATES = ['Z_state', 'Y_state', 'A_state']

#: The 17 parameters the ODE right-hand side reads (log10-scaled in the fit box).
DYN = ['K2', 'K_par', 'Ka', 'Kd', 'Kf', 'Kp', 'Ky', 'Kz', 'Vd', 'Vm2', 'Vm3',
       'Vp', 'beta_par', 'epsilon_par', 'n_par', 'v0', 'v1']
#: Observation parameters (log10-scaled); sigma is profiled, not searched.
OBS = ['scale', 'offset']

#: The PEtab/PyBNF fit box, in each parameter's own sampling space.
BOX_LOG = (-3.0, 5.0)          # loguniform_var  0.001 .. 100000
BOX_INIT = (0.0, 1.0)          # uniform_var     init_*  0 .. 1
BOX_Z_LOG = (-6.0, 3.0)        # auxiliary segment-start states (log10)

#: PEtab nominal point, in sampling space (log10 for DYN/OBS, linear for init_*).
NOMINAL = {
    'K2': -1.0002234671823191, 'K_par': 1.0573648249340193,
    'Ka': -0.7042263782210675, 'Kd': -0.4061630917376330,
    'Kf': 0.05612305941176864, 'Kp': -0.0014932957226511624,
    'Ky': -0.6981800184435972, 'Kz': -0.5180509225644002,
    'Vd': 1.9678938235635597, 'Vm2': 0.8725713912648427,
    'Vm3': 1.3560741603461504, 'Vp': 0.4404143724192836,
    'beta_par': 0.05074788136550951, 'epsilon_par': -0.7871941785182649,
    'n_par': 0.6130501154695961, 'v0': 0.3650735522312953,
    'v1': 0.0021174686013911382,
    'scale': -0.1817258330925958, 'offset': -0.5217746358806918,
    'init_Z_state': 0.087920524425504, 'init_Y_state': 0.999348084438687,
    'init_A_state': 0.99999999999996,
}

LN10 = np.log(10.0)
#: ``timeout`` bounds a pathological parameter point the way PyBNF's own
#: ``wall_time_sim = 10`` does on this job -- without it a stiff point can spend
#: minutes inside CVODE's analytical-Jacobian failure + finite-difference retry, and
#: a multistart sweep is dominated by points that were never going to score.
_TOL = dict(rtol=1e-9, atol=1e-11, max_steps=200000, timeout=0.5)


def load_data():
    d = np.loadtxt(DATA_EXP, skiprows=1)
    return d[:, 0].copy(), d[:, 1].copy()


TIMES, YOBS = load_data()
N_OBS = len(TIMES)
T_END = float(TIMES[-1])
#: The parameter-independent constants PyBNF's reduced objective drops.
PAPER_OFFSET = N_OBS * 0.5 * np.log(2 * np.pi) + float(np.sum(np.log(YOBS * LN10)))


class SimulationFailed(Exception):
    pass


class OutOfTime(Exception):
    """Raised from inside the residual callback to unwind an inner solve that has
    run past its wall-clock share.  A truncated inner solve is a normal outcome
    here -- the outer loop keeps the iterate it reached."""


# --------------------------------------------------------------------------- #
# Variable layout
# --------------------------------------------------------------------------- #
class Layout:
    """The (theta, z) vector layout for an m-segment transcription.

    ``theta`` is the biological block -- 17 log10 dynamic parameters, log10 scale,
    log10 offset, and the three linear ``init_*`` that ARE segment 0's start state.
    ``z`` is the auxiliary block: 3 log10 states per interior knot.  A layout with
    ``m == 1`` has an empty auxiliary block and is exactly single shooting.
    """

    def __init__(self, m):
        self.m = int(m)
        self.theta_names = list(DYN) + list(OBS) + ['init_%s' % s for s in STATES]
        self.n_theta = len(self.theta_names)                       # 22
        self.n_aux = 3 * (self.m - 1)
        self.n = self.n_theta + self.n_aux
        self.i_dyn = slice(0, len(DYN))
        self.i_scale = len(DYN)
        self.i_offset = len(DYN) + 1
        self.i_init = slice(len(DYN) + 2, self.n_theta)
        self.knots = np.linspace(0.0, T_END, self.m + 1)

    def bounds(self):
        lo = np.empty(self.n)
        hi = np.empty(self.n)
        lo[:len(DYN) + 2] = BOX_LOG[0]
        hi[:len(DYN) + 2] = BOX_LOG[1]
        lo[self.i_init] = BOX_INIT[0]
        hi[self.i_init] = BOX_INIT[1]
        if self.n_aux:
            lo[self.n_theta:] = BOX_Z_LOG[0]
            hi[self.n_theta:] = BOX_Z_LOG[1]
        return lo, hi

    def z_of(self, x, j):
        """Segment ``j``'s start state, in model units (Z, Y, A order)."""
        if j == 0:
            return np.asarray(x[self.i_init], dtype=float)
        k = self.n_theta + 3 * (j - 1)
        return 10.0 ** np.asarray(x[k:k + 3], dtype=float)

    def z_slice(self, j):
        """Where segment ``j``'s auxiliary block lives in the vector (j >= 1)."""
        k = self.n_theta + 3 * (j - 1)
        return slice(k, k + 3)

    def sample_times(self, j):
        """Output times for segment ``j``: its start knot, the data inside it, its
        end knot.  The end knot is always present so the continuity defect can be
        read off the same run that produced the data rows."""
        t0, t1 = self.knots[j], self.knots[j + 1]
        mask = (TIMES > t0) & (TIMES <= t1)
        pts = np.concatenate(([t0], TIMES[mask], [t1]))
        pts = np.unique(pts)
        return pts, np.flatnonzero(mask)


# --------------------------------------------------------------------------- #
# Simulation
# --------------------------------------------------------------------------- #
_TEMPLATE = None
_SIMS = {}


def _template():
    global _TEMPLATE
    if _TEMPLATE is None:
        _TEMPLATE = bngsim.Model.from_sbml(MODEL_XML)
    return _TEMPLATE


def _sim_for(want_jac):
    """One long-lived Simulator per (sensitivity on/off), reused across every
    segment and every evaluation.  Constructing a sensitivity-bearing Simulator
    costs ~17 ms even on a warm codegen cache, which at m segments per residual
    evaluation dominates the 50 ms the integration itself takes; mutating the
    model behind an existing Simulator and ``reset()``-ing gives bit-identical
    states AND sensitivities (verified against a fresh Simulator)."""
    sim = _SIMS.get(want_jac)
    if sim is None:
        kwargs = dict(sensitivity_params=DYN, sensitivity_ic=STATES) if want_jac else {}
        sim = bngsim.Simulator(_template().clone(), method='ode', **kwargs)
        _SIMS[want_jac] = sim
    return sim


def _engine(theta_lin, z):
    m = _template().clone()
    for name in DYN:
        m.set_param(name, float(theta_lin[name]))
    for name, value in zip(STATES, z):
        m.set_concentration(name, float(value))
    m.save_concentrations()
    m.reset()
    return m


def simulate_segment(theta_lin, z, sample_times, *, want_jac):
    """One segment.  Returns ``(x, dx_dp, dx_dz)`` at the requested times, with
    ``x`` in (Z, Y, A) order, ``dx_dp`` of shape ``(n_times, 3, len(DYN))`` and
    ``dx_dz`` of shape ``(n_times, 3, 3)``."""
    sim = _sim_for(want_jac)
    model = sim.model
    for name in DYN:
        model.set_param(name, float(theta_lin[name]))
    for name, value in zip(STATES, z):
        model.set_concentration(name, float(value))
    model.save_concentrations()
    model.reset()
    pts = [float(p) for p in sample_times]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = sim.run(t_span=(pts[0], pts[-1]), n_points=len(pts),
                          sample_times=pts, **_TOL)
    except Exception as exc:                       # a non-integrable point
        raise SimulationFailed(str(exc)) from exc

    names = list(res.species_names)
    idx = [names.index(s) for s in STATES]
    x = np.asarray(res.species, dtype=float)[:, idx]
    if not np.all(np.isfinite(x)):
        raise SimulationFailed('non-finite state')
    if not want_jac:
        return x, None, None

    p_axis = list(res.sensitivity_params)
    ic_axis = list(res.sensitivity_ic_species)
    sp = np.asarray(res.sensitivities, dtype=float)
    si = np.asarray(res.sensitivities_ic, dtype=float)
    dx_dp = sp[:, idx, :][:, :, [p_axis.index(p) for p in DYN]]
    dx_dz = si[:, idx, :][:, :, [ic_axis.index(s) for s in STATES]]
    if not (np.all(np.isfinite(dx_dp)) and np.all(np.isfinite(dx_dz))):
        raise SimulationFailed('non-finite sensitivities')
    return x, dx_dp, dx_dz


# --------------------------------------------------------------------------- #
# Residual model
# --------------------------------------------------------------------------- #
class Transcription:
    """The m-segment residual model: data residuals + continuity defects, with
    the analytic Jacobian of both with respect to the full (theta, z) vector."""

    def __init__(self, m, *, c_scale=1.0):
        self.L = Layout(m)
        self.c_scale = float(c_scale)

    # -- pieces ------------------------------------------------------------- #
    def _theta_lin(self, x):
        return {name: 10.0 ** float(x[i]) for i, name in enumerate(DYN)}

    def evaluate(self, x, *, want_jac=True):
        """Simulate every segment.  Returns a dict with the data residuals ``r``
        (unweighted log10 residuals), the continuity defects ``c``, and their
        Jacobians in sampling space."""
        L = self.L
        x = np.asarray(x, dtype=float)
        theta_lin = self._theta_lin(x)
        scale = 10.0 ** float(x[L.i_scale])
        offset = 10.0 ** float(x[L.i_offset])

        r = np.zeros(N_OBS)
        Jr = np.zeros((N_OBS, L.n)) if want_jac else None
        c = np.zeros(L.n_aux)
        Jc = np.zeros((L.n_aux, L.n)) if want_jac else None

        for j in range(L.m):
            pts, data_idx = L.sample_times(j)
            z = L.z_of(x, j)
            xs, dx_dp, dx_dz = simulate_segment(theta_lin, z, pts, want_jac=want_jac)

            # d z_j / d (its own block in sampling space)
            if j == 0:
                dz_du = np.eye(3)                      # init_* are linear
                z_cols = np.arange(L.i_init.start, L.i_init.stop)
            else:
                dz_du = np.diag(z * LN10)              # z_j stored as log10
                sl = L.z_slice(j)
                z_cols = np.arange(sl.start, sl.stop)

            # --- data rows ---------------------------------------------------
            if len(data_idx):
                rows = [int(np.searchsorted(pts, TIMES[i])) for i in data_idx]
                Zi = xs[rows, 0]
                mu = scale * Zi + offset
                if np.any(mu <= 0):
                    raise SimulationFailed('non-positive prediction')
                r[data_idx] = np.log10(YOBS[data_idx]) - np.log10(mu)
                if want_jac:
                    # d r / d q = -(1/(ln10 * mu)) * d mu / d q
                    pref = -1.0 / (LN10 * mu)
                    dmu_dp = scale * dx_dp[rows, 0, :]          # (n_i, n_dyn)
                    Jr[np.ix_(data_idx, np.arange(len(DYN)))] = (
                        pref[:, None] * dmu_dp * np.array(
                            [theta_lin[p] * LN10 for p in DYN])[None, :])
                    Jr[data_idx, L.i_scale] = pref * Zi * scale * LN10
                    Jr[data_idx, L.i_offset] = pref * offset * LN10
                    dmu_dz = scale * dx_dz[rows, 0, :] @ dz_du   # (n_i, 3)
                    Jr[np.ix_(data_idx, z_cols)] = pref[:, None] * dmu_dz

            # --- continuity row block (segment j's end feeds knot j+1) -------
            if j < L.m - 1:
                blk = slice(3 * j, 3 * j + 3)
                z_next = L.z_of(x, j + 1)
                c[blk] = (xs[-1, :] - z_next) / self.c_scale
                if want_jac:
                    Jc[blk, :len(DYN)] = (
                        dx_dp[-1, :, :] * np.array([theta_lin[p] * LN10 for p in DYN])[None, :]
                    ) / self.c_scale
                    Jc[np.ix_(np.arange(blk.start, blk.stop), z_cols)] = (
                        dx_dz[-1, :, :] @ dz_du) / self.c_scale
                    nxt = L.z_slice(j + 1)
                    Jc[np.ix_(np.arange(blk.start, blk.stop),
                              np.arange(nxt.start, nxt.stop))] -= (
                        np.diag(z_next * LN10) / self.c_scale)

        out = {'r': r, 'c': c}
        if want_jac:
            out['Jr'] = Jr
            out['Jc'] = Jc
        return out

    # -- objective ---------------------------------------------------------- #
    @staticmethod
    def sigma_profiled(r):
        return float(np.sqrt(np.sum(r ** 2) / N_OBS))

    @staticmethod
    def reduced_objective(r):
        """PyBNF's reduced objective at the profiled sigma."""
        S = float(np.sum(r ** 2))
        if S <= 0:
            return -np.inf
        return 0.5 * N_OBS * np.log(S / N_OBS) + 0.5 * N_OBS


def paper_nll(reduced):
    return reduced + PAPER_OFFSET


def optimality_gap(reduced):
    return paper_nll(reduced) - JSTAR


# --------------------------------------------------------------------------- #
# Single-shoot certification
# --------------------------------------------------------------------------- #
def certify(x_theta):
    """Discard z, re-simulate theta with ordinary single shooting, and return the
    reduced objective.  Every reported score comes through here."""
    t = Transcription(1)
    try:
        out = t.evaluate(np.asarray(x_theta, dtype=float)[:t.L.n_theta], want_jac=False)
    except SimulationFailed:
        return np.inf
    return t.reduced_objective(out['r'])


# --------------------------------------------------------------------------- #
# Augmented-Lagrangian outer loop
# --------------------------------------------------------------------------- #
BIG = 1e3


def _packed(tr, x, sigma, lam, rho, *, want_jac):
    """The augmented-Lagrangian residual (and Jacobian) scipy's `trf` minimizes:

        f = [ r / sigma ; sqrt(rho) * (c + lam / rho) ]

    so that 0.5 ||f||^2 == data term + lam^T c + (rho/2) ||c||^2 + const."""
    out = tr.evaluate(x, want_jac=want_jac)
    r, c = out['r'], out['c']
    f = [r / sigma]
    if len(c):
        f.append(np.sqrt(rho) * (c + lam / rho))
    f = np.concatenate(f)
    if not want_jac:
        return f, out
    J = [out['Jr'] / sigma]
    if len(c):
        J.append(np.sqrt(rho) * out['Jc'])
    return f, np.vstack(J), out


def solve(m, x0, *, rho0=10.0, outer=12, inner_max=80, c_tol=1e-7,
          eta=0.5, gamma=5.0, verbose=False, c_scale=1.0, best=None,
          deadline=None):
    """Augmented-Lagrangian multiple shooting from ``x0`` (a full (theta, z) vector).

    The penalty starts *small* and escalates slowly (``gamma = 3``), on purpose.
    Multiple shooting's convergence-region win comes from being allowed to hold a
    discontinuous trajectory while the parameters move (Balsa-Canto et al. 2008,
    section "Multiple-shooting"); a penalty that forces continuity in the first few
    outer iterations throws that away and leaves nothing but a more expensive single
    shooting.  The multipliers -- not the penalty -- are what eventually close the
    defects.

    ``best`` is an optional shared ``{'J': ..., 'theta': ...}`` record: every outer
    iterate is certified through ordinary single shooting and the best one kept, so a
    good parameter vector passed through mid-run is not lost if a later outer
    iteration wanders.  Returns ``(x, info)``.
    """
    tr = Transcription(m, c_scale=c_scale)
    lo, hi = tr.L.bounds()
    x = np.clip(np.asarray(x0, dtype=float), lo, hi)
    lam = np.zeros(tr.L.n_aux)
    rho = float(rho0)
    c_prev = np.inf
    n_sim = 0
    history = []
    c_norm = np.inf

    def _record(xx):
        """Certify this iterate's theta through single shooting and keep the best."""
        J = certify(xx)
        if best is not None and J < best['J']:
            best['J'] = J
            best['theta'] = np.asarray(xx[:tr.L.n_theta], dtype=float).copy()
        return J

    for k in range(outer):
        if deadline is not None and time.monotonic() > deadline:
            break
        try:
            base = tr.evaluate(x, want_jac=False)
        except SimulationFailed:
            return x, {'certified': np.inf, 'c_norm': np.inf, 'outer': k,
                       'status': 'start failed', 'n_sim': n_sim, 'history': history}
        sigma = max(tr.sigma_profiled(base['r']), 1e-8)

        def fun(xx):
            nonlocal n_sim
            if deadline is not None and time.monotonic() > deadline:
                raise OutOfTime
            n_sim += 1
            try:
                f, _ = _packed(tr, xx, sigma, lam, rho, want_jac=False)
            except SimulationFailed:
                f = np.full(N_OBS + tr.L.n_aux, BIG)
            return f

        def jac(xx):
            try:
                _, J, _ = _packed(tr, xx, sigma, lam, rho, want_jac=True)
            except SimulationFailed:
                J = np.zeros((N_OBS + tr.L.n_aux, tr.L.n))
            return J

        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                sol = least_squares(fun, x, jac=jac, bounds=(lo, hi), method='trf',
                                    x_scale='jac', max_nfev=inner_max,
                                    xtol=1e-12, ftol=1e-12, gtol=1e-12)
            x = np.clip(sol.x, lo, hi)
        except OutOfTime:
            _record(x)
            break

        try:
            out = tr.evaluate(x, want_jac=False)
        except SimulationFailed:
            return x, {'certified': np.inf, 'c_norm': np.inf, 'outer': k + 1,
                       'status': 'diverged', 'n_sim': n_sim, 'history': history}
        c_norm = float(np.linalg.norm(out['c'])) if tr.L.n_aux else 0.0
        cert = _record(x)
        history.append({'outer': k, 'rho': rho, 'c_norm': c_norm,
                        'reduced': tr.reduced_objective(out['r']), 'certified': cert})
        if verbose:
            print('    outer %2d: rho=%-8.3g ||c||=%.3e  segmented J=%10.4f  '
                  'certified J=%10.4f' % (k, rho, c_norm,
                                          tr.reduced_objective(out['r']), cert))
        if c_norm <= c_tol:
            break
        if c_norm > eta * c_prev:
            rho *= gamma
        else:
            lam = lam + rho * out['c']
        c_prev = c_norm

    return x, {'certified': certify(x), 'c_norm': c_norm, 'outer': k + 1,
               'status': 'ok', 'n_sim': n_sim, 'history': history}


def coarsen(x, m_from, m_to):
    """Re-seed an ``m_to``-segment vector from an ``m_from``-segment solution.

    Equal-time knots nest whenever ``m_from`` is a multiple of ``m_to``, so every
    knot the coarser transcription keeps already carries a solved auxiliary state
    and is copied straight across; nothing is re-simulated and nothing is guessed.
    """
    Lf, Lt = Layout(m_from), Layout(m_to)
    out = np.zeros(Lt.n)
    out[:Lt.n_theta] = x[:Lf.n_theta]
    if Lt.n_aux == 0:
        return out
    if m_from % m_to:
        raise ValueError('knots only nest when m_from is a multiple of m_to')
    step = m_from // m_to
    for j in range(1, m_to):
        out[Lt.z_slice(j)] = x[Lf.z_slice(j * step)]
    return out


def solve_homotopy(x_theta, schedule=(8, 4, 2, 1), *, verbose=False,
                   budget=None, **kw):
    """The MVP recipe: many short segments first (the easiest landscape), coarsening
    toward one (issue #563, formulation point 3).  The final ``m = 1`` stage IS the
    certified reconstruction -- ordinary single shooting, re-solved from the
    multiple-shooting parameter vector -- so the reported score is never a segmented
    one.

    ``budget`` is the total wall-clock seconds for the whole homotopy, divided evenly
    across stages; a stage that overruns is truncated at whatever iterate it reached,
    which is a normal outcome rather than a failure.  Returns ``(theta, info)``."""
    best = {'J': np.inf, 'theta': np.asarray(x_theta, dtype=float)[:Layout(1).n_theta].copy()}
    schedule = [int(s) for s in schedule]
    x = seed_aux(x_theta, schedule[0])
    n_sim = 0
    stages = []
    per_stage = None if budget is None else float(budget) / len(schedule)
    for i, m in enumerate(schedule):
        if i:
            x = coarsen(x, schedule[i - 1], m)
        deadline = None if per_stage is None else time.monotonic() + per_stage
        if verbose:
            print('  -- stage m = %d --' % m)
        if m == 1:
            # The certified reconstruction polishes the best theta the run produced,
            # not merely the last one.  Every outer iterate has already been certified
            # through single shooting, and on a landscape this multimodal the coarsest
            # segmented stage is not reliably the one that held the best parameters.
            seed_theta = (best['theta'] if np.isfinite(best['J'])
                          else x[:Layout(1).n_theta])
            x_t, info = solve_single_shooting(seed_theta, deadline=deadline)
            x = x_t
            if info['certified'] < best['J']:
                best['J'] = info['certified']
                best['theta'] = np.asarray(x_t, dtype=float).copy()
        else:
            x, info = solve(m, x, verbose=verbose, best=best, deadline=deadline, **kw)
        n_sim += info['n_sim']
        stages.append({'m': m, 'certified': info['certified'],
                       'c_norm': info.get('c_norm', 0.0), 'n_sim': info['n_sim']})
        if verbose:
            print('     stage m=%d -> certified J = %.4f (||c|| = %.2e)'
                  % (m, info['certified'], info.get('c_norm', 0.0)))
    return best['theta'], {'certified': best['J'], 'n_sim': n_sim, 'stages': stages}


def solve_single_shooting(x0, *, inner_max=400, rounds=6, deadline=None):
    """The control: the identical inner solver, identical settings, identical
    objective -- with no segments, no auxiliary states, and no outer loop.

    Every knob that could otherwise explain a difference is held equal to the
    multiple-shooting arm: the same ``trf`` driver with ``x_scale='jac'`` and the
    same tolerances, the same profiled sigma, the same wall-clock deadline, and
    enough sigma-reprofiling rounds that the control is never the truncated arm.
    """
    tr = Transcription(1)
    lo, hi = tr.L.bounds()
    x = np.clip(np.asarray(x0, dtype=float)[:tr.L.n_theta], lo, hi)
    n_sim = 0
    try:
        base = tr.evaluate(x, want_jac=False)
    except SimulationFailed:
        return x, {'certified': np.inf, 'c_norm': 0.0, 'status': 'start failed',
                   'n_sim': 0}
    sigma = max(tr.sigma_profiled(base['r']), 1e-8)

    for _ in range(rounds):                  # re-profile sigma between rounds
        if deadline is not None and time.monotonic() > deadline:
            break

        def fun(xx):
            nonlocal n_sim
            if deadline is not None and time.monotonic() > deadline:
                raise OutOfTime
            n_sim += 1
            try:
                return tr.evaluate(xx, want_jac=False)['r'] / sigma
            except SimulationFailed:
                return np.full(N_OBS, BIG)

        def jac(xx):
            try:
                return tr.evaluate(xx, want_jac=True)['Jr'] / sigma
            except SimulationFailed:
                return np.zeros((N_OBS, tr.L.n))

        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                sol = least_squares(fun, x, jac=jac, bounds=(lo, hi), method='trf',
                                    x_scale='jac', max_nfev=inner_max, xtol=1e-12,
                                    ftol=1e-12, gtol=1e-12)
            x_new = np.clip(sol.x, lo, hi)
        except OutOfTime:
            break
        moved = float(np.max(np.abs(x_new - x)))
        x = x_new
        try:
            sigma = max(tr.sigma_profiled(tr.evaluate(x, want_jac=False)['r']), 1e-8)
        except SimulationFailed:
            break
        if moved < 1e-9:                     # reprofiling no longer moves the point
            break
    return x, {'certified': certify(x), 'c_norm': 0.0, 'status': 'ok', 'n_sim': n_sim}


# --------------------------------------------------------------------------- #
# Start points
# --------------------------------------------------------------------------- #
def nominal_theta():
    L = Layout(1)
    return np.array([NOMINAL[n] if n in NOMINAL else NOMINAL[n]
                     for n in L.theta_names], dtype=float)


def seed_aux(x_theta, m):
    """Initialize the auxiliary states by single-shooting theta once and reading
    the state at each knot -- the generic cold start (no fitted solution, no
    burst-specific knots).  Continuity holds exactly at the start point, so the
    outer loop begins feasible and every subsequent discontinuity is the
    optimizer's own choice."""
    L = Layout(m)
    x = np.zeros(L.n)
    x[:L.n_theta] = np.asarray(x_theta, dtype=float)[:L.n_theta]
    if L.n_aux == 0:
        return x
    tr = Transcription(1)
    theta_lin = tr._theta_lin(x)
    z0 = np.asarray(x[L.i_init], dtype=float)
    pts = np.unique(np.concatenate((L.knots, TIMES)))
    try:
        xs, _, _ = simulate_segment(theta_lin, z0, pts, want_jac=False)
    except SimulationFailed:
        x[L.n_theta:] = np.log10(np.tile(np.maximum(z0, 1e-3), L.m - 1))
        return x
    for j in range(1, L.m):
        idx = int(np.searchsorted(pts, L.knots[j]))
        state = np.maximum(xs[idx, :], 1e-8)
        x[L.z_slice(j)] = np.log10(state)
    return x
