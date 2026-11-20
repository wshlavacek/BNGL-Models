"""Independent SciPy integration of the generated MARM2.0 network (Fröhlich et al., 2023).

The curated arm of `run_frohlich2023.py`. It reads the `.net` that BioNetGen generates from
the curated `.bngl`, rebuilds `dy/dt = S r(y)` from scratch -- constant expressions, sparse
stoichiometry, observables, mass-action fluxes -- and integrates it with SciPy's implicit BDF
method, on a sparse finite-difference Jacobian or on the analytic one (`make_jacobian`).
Nothing here calls BioNetGen's `run_network`.

Fixed (`$`) species are held constant. Four rate laws are state-dependent functions: the three
ERK-induced transcription rates (DUSP, EGFR, SPRY) depend on pERK, and EGFR internalization
depends on plasma-membrane pEGFR. They are written out explicitly below from the `.bngl`
`begin functions` block rather than evaluated from the `.net` text, so a change to either side
shows up as a disagreement. Every other function in the model is an output mapping evaluated
after integration (`perk_if`).

`parse_net`, `dynamic_values`, `make_rhs`, `jacobian_sparsity` and `load_gdat` are the
level-1 implementation in `verify_frohlich2023.ipynb`, factored out unchanged in behaviour so
the driver integrates the network the same way; the notebook still carries its inline copy.
"""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import csc_matrix

DYNAMIC_NAMES = [
    "msynthesize_pERK_DUSP_ksyn",
    "catalyze_EGFR_kendo",
    "msynthesize_pERK_EGFR_ksyn",
    "msynthesize_pERK_SPRY_ksyn",
]


def parse_net(path):
    """Parse a BioNetGen `.net` into the arrays the right-hand side needs."""
    env = {"ln": math.log, "log": math.log, "exp": math.exp, "sqrt": math.sqrt}
    params, y0, patterns, reactions, groups = {}, [], [], [], []
    block = None
    for raw in open(path):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("begin "):
            block = line.split()[1]
            continue
        if line.startswith("end "):
            block = None
            continue
        body = line.split("#", 1)[0].strip()
        fields = body.split()
        if block == "parameters":
            name = fields[1]
            expr = "".join(fields[2:]).replace("^", "**")
            params[name] = float(eval(expr, {"__builtins__": {}}, env | params))
        elif block == "species":
            patterns.append(fields[1])
            y0.append(float(eval("".join(fields[2:]), {"__builtins__": {}}, env | params)))
        elif block == "reactions":
            reactants = [int(x) - 1 for x in fields[1].split(",") if x != "0"]
            products = [int(x) - 1 for x in fields[2].split(",") if x != "0"]
            reactions.append((reactants, products, fields[3]))
        elif block == "groups":
            group = body.split(None, 2)
            members = [] if len(group) < 3 else [
                (int(token.split("*")[-1]) - 1,
                 float(token.split("*")[0]) if "*" in token else 1.0)
                for token in group[2].split(",") if token
            ]
            groups.append((group[1], members))

    n, m = len(y0), len(reactions)
    group_map = {
        name: (np.array([i for i, _ in members], dtype=int),
               np.array([w for _, w in members], dtype=float))
        for name, members in groups
    }
    rows, cols, values = [], [], []
    ra = np.full(m, -1, dtype=int)
    rb = np.full(m, -1, dtype=int)
    constants = np.empty(m)
    dynamic_kind = np.full(m, -1, dtype=int)
    for j, (reactants, products, expr) in enumerate(reactions):
        if reactants:
            ra[j] = reactants[0]
        if len(reactants) > 1:
            rb[j] = reactants[1]
        if len(reactants) > 2:
            raise ValueError(f"reaction {j + 1} has more than two reactants")
        for species in reactants:
            rows.append(species)
            cols.append(j)
            values.append(-1.0)
        for species in products:
            rows.append(species)
            cols.append(j)
            values.append(1.0)
        for kind, name in enumerate(DYNAMIC_NAMES):
            if name in expr:
                dynamic_kind[j] = kind
                expr = expr.replace(name, "1")
        constants[j] = float(eval(expr, {"__builtins__": {}}, env | params))
    S = csc_matrix((values, (rows, cols)), shape=(n, m))
    fixed = np.array(["$" in pattern for pattern in patterns])
    return {"params": params, "y0": np.array(y0), "patterns": patterns,
            "reactions": reactions, "groups": group_map, "S": S,
            "ra": ra, "rb": rb, "constants": constants,
            "dynamic_kind": dynamic_kind, "fixed": fixed,
            "first": np.flatnonzero((ra >= 0) & (rb < 0)),
            "second": np.flatnonzero(rb >= 0),
            "dynamic": np.flatnonzero(dynamic_kind >= 0)}


def group_value(net, name, y):
    """Value of a `.net` group (observable) on a state vector or a (species x time) matrix."""
    indices, weights = net["groups"][name]
    return np.dot(weights, y[indices])


def dynamic_values(net, y):
    """The four state-dependent rate laws, from the `.bngl` `begin functions` block."""
    p = net["params"]
    perk = group_value(net, "modulation_synthesize_pERK_DUSP_ERK", y)
    return np.array([
        p["msynthesize_pERK_DUSP_kbase"] *
        (1 + perk * p["synthesize_pERK_DUSP_ERK_gexpslope"] /
         (perk + p["synthesize_pERK_DUSP_ERK_kM"])),
        p["catalyze_EGFR_endo_kcat"] /
        (group_value(net, "pEGFR_PM", y) + p["catalyze_EGFR_deg_kM"]),
        p["msynthesize_pERK_EGFR_kbase"] *
        (1 + perk * p["synthesize_pERK_EGFR_ERK_gexpslope"] /
         (perk + p["synthesize_pERK_EGFR_ERK_kM"])) * p["EGFR_crispr"],
        p["msynthesize_pERK_SPRY_kbase"] *
        (1 + perk * p["synthesize_pERK_SPRY_ERK_gexpslope"] /
         (perk + p["synthesize_pERK_SPRY_ERK_kM"]))
    ])


def fluxes(net, y):
    """Reaction fluxes r(y), one per `.net` reaction."""
    ra, rb = net["ra"], net["rb"]
    first, second, dynamic = net["first"], net["second"], net["dynamic"]
    rates = net["constants"].copy()
    rates[dynamic] *= dynamic_values(net, y)[net["dynamic_kind"][dynamic]]
    rates[first] *= y[ra[first]]
    rates[second] *= y[ra[second]] * y[rb[second]]
    return rates


def make_rhs(net):
    """dy/dt = S r(y), with fixed species held constant."""
    S, fixed = net["S"], net["fixed"]

    def rhs(_time, y):
        derivative = S @ fluxes(net, y)
        derivative[fixed] = 0.0
        return derivative
    return rhs


def dynamic_derivatives(net, y):
    """d/dO of the four state-dependent rate laws, O being the observable each one reads."""
    p = net["params"]
    perk = group_value(net, "modulation_synthesize_pERK_DUSP_ERK", y)
    pegfr = group_value(net, "pEGFR_PM", y)

    def saturating(kbase, slope, km):
        return kbase * slope * km / (perk + km) ** 2

    return np.array([
        saturating(p["msynthesize_pERK_DUSP_kbase"], p["synthesize_pERK_DUSP_ERK_gexpslope"],
                   p["synthesize_pERK_DUSP_ERK_kM"]),
        -p["catalyze_EGFR_endo_kcat"] / (pegfr + p["catalyze_EGFR_deg_kM"]) ** 2,
        saturating(p["msynthesize_pERK_EGFR_kbase"], p["synthesize_pERK_EGFR_ERK_gexpslope"],
                   p["synthesize_pERK_EGFR_ERK_kM"]) * p["EGFR_crispr"],
        saturating(p["msynthesize_pERK_SPRY_kbase"], p["synthesize_pERK_SPRY_ERK_gexpslope"],
                   p["synthesize_pERK_SPRY_ERK_kM"]),
    ])


DYNAMIC_OBSERVABLES = ["modulation_synthesize_pERK_DUSP_ERK", "pEGFR_PM",
                       "modulation_synthesize_pERK_EGFR_ERK", "modulation_synthesize_pERK_SPRY_ERK"]


def make_flux_jacobian(net):
    """dr/dy as a sparse (reactions x species) matrix, analytically."""
    n, m = len(net["y0"]), len(net["reactions"])
    ra, rb, c = net["ra"], net["rb"], net["constants"]
    kind = net["dynamic_kind"]
    first, second, dynamic = net["first"], net["second"], net["dynamic"]
    unary = np.flatnonzero(ra >= 0)
    obs = [net["groups"][name] for name in DYNAMIC_OBSERVABLES]
    dyn_rows, dyn_cols, dyn_w = [np.zeros(0, int)], [np.zeros(0, int)], [np.zeros(0)]
    for j in dynamic:
        idx, w = obs[kind[j]]
        dyn_rows.append(np.full(len(idx), j))
        dyn_cols.append(idx)
        dyn_w.append(w)
    dyn_rows, dyn_cols, dyn_w = map(np.concatenate, (dyn_rows, dyn_cols, dyn_w))
    rows = np.concatenate([unary, second, dyn_rows])
    cols = np.concatenate([ra[unary], rb[second], dyn_cols])

    def flux_jacobian(y):
        scale = c.copy()
        scale[dynamic] *= dynamic_values(net, y)[kind[dynamic]]
        # mass action: d(scale * y_a [* y_b]) / dy_a and / dy_b (A + A sums to 2 scale y_a)
        va = scale[unary] * np.where(rb[unary] >= 0, y[np.maximum(rb[unary], 0)], 1.0)
        vb = scale[second] * y[ra[second]]
        # state-dependent rate laws: c * (mass-action part) * dD/dO * w_i over O's species
        mass = np.ones(m)
        mass[first] = y[ra[first]]
        mass[second] = y[ra[second]] * y[rb[second]]
        vd = (c[dyn_rows] * mass[dyn_rows]
              * dynamic_derivatives(net, y)[kind[dyn_rows]] * dyn_w)
        return csc_matrix((np.concatenate([va, vb, vd]), (rows, cols)), shape=(m, n))
    return flux_jacobian


def make_jacobian(net):
    """Analytic sparse Jacobian of make_rhs, S @ dr/dy with fixed species' rows zeroed.

    Only the Newton iterations of the implicit integrator use it, so it changes how fast the
    integrator converges, never what it converges to. Long horizons need it: with the
    finite-difference Jacobian, steps of tens of hours fail Newton convergence on the stiff,
    near-zero species and BDF crawls.
    """
    S = net["S"].tocsr()
    flux_jacobian = make_flux_jacobian(net)
    mask = np.where(net["fixed"], 0.0, 1.0)[:, None]

    def jac(_time, y):
        return csc_matrix((S @ flux_jacobian(y)).multiply(mask))
    return jac


def jacobian_sparsity(net):
    """Jacobian pattern: mass-action coupling plus the observables behind the four functions."""
    rows, cols = [], []
    extra = [net["groups"]["modulation_synthesize_pERK_DUSP_ERK"][0],
             net["groups"]["pEGFR_PM"][0],
             net["groups"]["modulation_synthesize_pERK_EGFR_ERK"][0],
             net["groups"]["modulation_synthesize_pERK_SPRY_ERK"][0]]
    for j, (reactants, products, _expr) in enumerate(net["reactions"]):
        dependencies = list(set(reactants))
        kind = net["dynamic_kind"][j]
        if kind >= 0:
            dependencies += list(extra[kind])
        for affected in set(reactants + products):
            if not net["fixed"][affected]:
                for dependency in dependencies:
                    rows.append(affected)
                    cols.append(dependency)
    return csc_matrix((np.ones(len(rows), dtype=bool), (rows, cols)),
                      shape=(len(net["y0"]), len(net["y0"])))


def integrate(net, times, y0=None, rtol=1e-9, atol=1e-11, analytic_jacobian=False):
    """Integrate from `y0` (default: the `.net` seed values) and return the solve_ivp result.

    By default BDF forms its Jacobian by finite differences on the sparsity pattern;
    `analytic_jacobian=True` hands it make_jacobian instead.
    """
    times = np.asarray(times, dtype=float)
    y0 = net["y0"] if y0 is None else np.asarray(y0, dtype=float)
    jac = {"jac": make_jacobian(net)} if analytic_jacobian else {
        "jac_sparsity": jacobian_sparsity(net)}
    solution = solve_ivp(make_rhs(net), (times[0], times[-1]), y0, method="BDF",
                         t_eval=times, rtol=rtol, atol=atol, **jac)
    if not solution.success:
        raise RuntimeError(solution.message)
    return solution


def perk_if(net, y):
    """pERK_IF_obs = pERK/tERK * pERK_IF_scale + pERK_IF_offset, as in the `.bngl`."""
    p = net["params"]
    return (group_value(net, "pERK", y) / group_value(net, "tERK", y)
            * p["pERK_IF_scale"] + p["pERK_IF_offset"])


def load_gdat(path):
    """Header and data of a BioNetGen `.gdat`."""
    header = open(path).readline().lstrip("#").split()
    return header, np.atleast_2d(np.loadtxt(path))
