"""Independent SciPy integration of a BioNetGen reaction network.

Reads a generated `.net` file and rebuilds the ODE system from scratch — species,
stoichiometry, mass-action fluxes, observables and functional rate laws — without
using BioNetGen's own `run_network` integrator. This is the level-1 independent
implementation for `models/ksr1_raf_mek_inhibitor_response_imoto2026`: BioNetGen
supplies the network, SciPy supplies the solution.

The Imoto et al. (2026) network has 2709 species and 27457 unidirectional
reactions and is stiff (rate constants span eight decades), so the right-hand
side is assembled as sparse matrices and integrated with an implicit BDF method
using a supplied Jacobian sparsity pattern.
"""

import math
import re

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix, lil_matrix

_SAFE = {
    "ln": math.log, "log": math.log10, "log10": math.log10, "exp": math.exp,
    "sqrt": math.sqrt, "abs": abs, "pow": pow,
    "__builtins__": {},
}


class Network:
    """A BioNetGen `.net` file rebuilt as a sparse ODE system."""

    def __init__(self, path):
        self.params, self.species, self.reactions = {}, [], []
        self.groups, self.functions = {}, []
        self._read(path)
        self._compile()

    # ------------------------------------------------------------------ parsing

    def _read(self, path):
        section = None
        for raw in open(path):
            line = raw.split("#")[0].rstrip()
            if not line.strip():
                continue
            if line.startswith("begin "):
                section = line[6:].strip()
                continue
            if line.startswith("end "):
                section = None
                continue
            f = line.split()
            indexed = f[0].isdigit()
            if section == "parameters":
                self.params[f[1] if indexed else f[0]] = f[2] if indexed else f[1]
            elif section == "species":
                name, val = (f[1], f[2]) if indexed else (f[0], f[1])
                self.species.append((name, val))
            elif section == "reactions":
                self.reactions.append(
                    ([int(x) for x in f[1].split(",")],
                     [int(x) for x in f[2].split(",")],
                     f[3])
                )
            elif section == "groups":
                off = 1 if indexed else 0
                terms = f[off + 1] if len(f) > off + 1 else ""
                d = {}
                for t in terms.split(","):
                    if not t:
                        continue
                    c, i = t.split("*") if "*" in t else ("1", t)
                    d[int(i) - 1] = d.get(int(i) - 1, 0.0) + float(c)
                self.groups[f[off]] = d
            elif section == "functions":
                # `.net` writes `<i> name() <expr>`; an evaluated `.net` writes
                # `name() = <expr>`.
                body = line.strip()
                if indexed:
                    body = body.split(None, 1)[1]
                name, eq, expr = body.partition("=")
                if not eq:
                    name, _, expr = body.partition(" ")
                self.functions.append((name.strip().rstrip("()"), expr.strip()))

    # ---------------------------------------------------------------- compiling

    def _compile(self):
        env = dict(_SAFE)
        for k, v in self.params.items():
            env[k] = float(v) if _is_number(v) else _ev(v, env)
        self.env = env

        n = len(self.species)
        self.names = [s for s, _ in self.species]
        self.constant = np.array([s.lstrip("@a-z:").startswith("$") or "::$" in s
                                  for s in self.names])
        self.x0 = np.array([float(v) if _is_number(v) else _ev(v, env)
                            for _, v in self.species])

        # observables as a sparse matrix
        gm = lil_matrix((len(self.groups), n))
        for r, (_, d) in enumerate(self.groups.items()):
            for i, c in d.items():
                gm[r, i] = c
        self.group_names = list(self.groups)
        self.G = csr_matrix(gm)

        # functions, in dependency order as written
        self.func_names = [f for f, _ in self.functions]
        # a function reference inside another function is written `name()`
        self.func_code = [
            compile(_ev_src(re.sub(r"\b(\w+)\(\)", r"\1", e)), "<fn>", "eval")
            for _, e in self.functions
        ]

        # reactions: constant part and function-valued part
        m = len(self.reactions)
        self.k = np.zeros(m)
        self.kfun = {}  # reaction index -> (multiplier, function name)
        for j, (_, _, rate) in enumerate(self.reactions):
            mult, sym = _split_rate(rate)
            if sym is None or sym in env:
                self.k[j] = _ev(rate, env)
            else:
                self.kfun[j] = (mult, sym)

        # reactant index arrays (at most 2 reactants in this network)
        self.r1 = np.array([r[0] - 1 for r, _, _ in self.reactions])
        self.r2 = np.array([(r[1] - 1 if len(r) > 1 else -1)
                            for r, _, _ in self.reactions])
        self.has_r2 = self.r2 >= 0

        # stoichiometry
        sm = lil_matrix((n, m))
        for j, (reac, prod, _) in enumerate(self.reactions):
            for i in reac:
                sm[i - 1, j] -= 1
            for i in prod:
                if i:  # index 0 means "nothing"
                    sm[i - 1, j] += 1
        sm[self.constant, :] = 0
        self.S = csr_matrix(sm)

        # Jacobian sparsity. Mass-action coupling first: species i influences
        # species l when some reaction consumes i and changes l.
        reactant = lil_matrix((n, m))
        for j, (reac, _, _) in enumerate(self.reactions):
            for i in reac:
                reactant[i - 1, j] = 1
        pat = ((abs(self.S) != 0).astype(np.int8) @ csr_matrix(reactant).T.astype(np.int8))
        pat = (pat != 0).tolil()
        # Functional rate laws depend on every species entering an observable, so
        # the rows they move are dense over that whole set.
        if self.kfun:
            cols = sorted({i for name in self.group_names for i in self.groups[name]})
            rows = sorted({i for j in self.kfun
                           for i in self.S[:, j].nonzero()[0]})
            for i in rows:
                pat[i, cols] = 1
        self.jac_sparsity = csr_matrix(pat.astype(np.int8))

    # ------------------------------------------------------------------ dynamics

    def rhs(self, _t, x):
        env = self.env
        obs = self.G @ x
        for name, val in zip(self.group_names, obs, strict=False):
            env[name] = val
        for name, code in zip(self.func_names, self.func_code, strict=False):
            env[name] = eval(code, env)
        k = self.k.copy()
        for j, (mult, sym) in self.kfun.items():
            k[j] = mult * env[sym]
        v = k * x[self.r1]
        v[self.has_r2] *= x[self.r2[self.has_r2]]
        return self.S @ v

    def run(self, t_eval, x0=None, atol=1e-8, rtol=1e-8, t0=0.0):
        x0 = self.x0 if x0 is None else x0
        sol = solve_ivp(
            self.rhs, (t0, float(np.max(t_eval))), x0, t_eval=np.asarray(t_eval),
            method="BDF", atol=atol, rtol=rtol, jac_sparsity=self.jac_sparsity,
        )
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol

    def observe(self, name, x):
        """Evaluate an observable (or the state itself) on a state vector or matrix."""
        i = self.group_names.index(name)
        return np.asarray(self.G[i] @ x).ravel()

    def evaluate(self, name, x):
        """Value of an observable or a function on one state vector or a state matrix."""
        if name in self.group_names:
            return self.observe(name, x)
        x = np.atleast_2d(np.asarray(x).T).T
        out = []
        for col in range(x.shape[1]):
            for gname, val in zip(self.group_names, self.G @ x[:, col], strict=False):
                self.env[gname] = val
            for fname, code in zip(self.func_names, self.func_code, strict=False):
                self.env[fname] = eval(code, self.env)
            out.append(self.env[name])
        return np.asarray(out)

    def steady_state(self, x0=None, t_end=1.0e5, atol=1e-8, rtol=1e-8):
        """Integrate far enough that the state stops moving; return the final state."""
        sol = self.run([t_end], x0=x0, atol=atol, rtol=rtol)
        return sol.y[:, -1]

    def set_param(self, **kw):
        """Change parameters and re-derive rate constants and initial amounts."""
        self.env.update({k: float(v) for k, v in kw.items()})
        for k, v in self.params.items():
            if k not in kw:
                self.env[k] = float(v) if _is_number(v) else _ev(v, self.env)
        self.x0 = np.array([float(v) if _is_number(v) else _ev(v, self.env)
                            for _, v in self.species])
        for j, (_, _, rate) in enumerate(self.reactions):
            if j not in self.kfun:
                self.k[j] = _ev(rate, self.env)


# ---------------------------------------------------------------------- helpers


def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def _ev_src(expr):
    return expr.replace("^", "**")


def _ev(expr, env):
    return eval(_ev_src(expr), {"__builtins__": {}}, env)


def _split_rate(rate):
    """'2*MM_uMEK' -> (2.0, 'MM_uMEK'); '0.5*k_local1' -> (0.5, 'k_local1')."""
    m = re.fullmatch(r"\s*(?:([\d.eE+-]+)\s*\*\s*)?([A-Za-z_]\w*)\s*", rate)
    if not m:
        return 1.0, None
    return (float(m.group(1)) if m.group(1) else 1.0), m.group(2)
