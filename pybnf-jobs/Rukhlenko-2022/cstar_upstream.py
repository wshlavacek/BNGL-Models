"""Check a Rukhlenko-2022 slug's `.bngl` against the authors' own file, with BioNetGen alone.

Shared by the per-slug `compare_upstream.py` drivers (cstar_trka, cstar_trkb,
cstar_skmel133_bmra_exact). Nothing here imports bngsim: the networks come from BNG2.pl, the
trajectories from BioNetGen's run_network, and the network comparison from the `.net` text,
evaluated by the small reader below.

The authors' repository, github.com/OleksiiR/cSTAR_Nature, carries no licence ("All rights
reserved. No licence is granted, express or implied, to use, copy, modify, or distribute the
software code and/or data in this repository"), so none of its files is vendored. Each slug's
`upstream.json` pins the path and sha256 of the one file it is checked against, at the commit
below; `authors_file` reads it from a local checkout of that commit (`--upstream DIR`) or from
a git-ignored cache under `dev/` that `--fetch` fills from the pinned raw URL, and refuses any
file whose sha256 differs.

The comparison has four parts, each written to the slug's `upstream_comparison.json`; the
first two are asserted, the last two reported:

  network       both files expanded by BNG2.pl; species matched by canonical label; every
                reaction matched by its canonical reactants and products; rate coefficients,
                functions, observables, parameters and dy/dt compared at sampled states.
  trajectories  the authors' file run by BNG2.pl with its own actions block (the published
                protocol, verbatim), against the slug's file run with the same actions
                translated mechanically onto its parameterization (`translate_actions`).
  solver        the slug's file under the same protocol at atol = rtol = 1e-12: the error of
                BioNetGen's default CVODE tolerances (1e-8).
  protocol      the slug's file with the pre-equilibration run to steady state instead of
                the authors' finite phases: how much the published protocol matters.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
import re
import subprocess
import sys
import urllib.request
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent

REPOSITORY = "https://github.com/OleksiiR/cSTAR_Nature"
COMMIT = "bfb092f2e4bffcfc43b3de01c301303fc6d4fd36"
RAW_URL = "https://raw.githubusercontent.com/OleksiiR/cSTAR_Nature/{commit}/{path}"
# git-ignored (`dev/` in .gitignore): the authors' files are fetched, never committed.
CACHE = REPO_ROOT / "dev" / "upstream" / f"cSTAR_Nature-{COMMIT[:7]}"

# BioNetGen's CVODE defaults (simulate() without atol/rtol), and the tight replicate.
BNG_DEFAULT_TOL = 1e-8
TIGHT_TOL = 1e-12
# Sampled states per network-comparison condition.
N_STATES = 32


# --------------------------------------------------------------------------- provenance


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_pin(slug_dir):
    pin = json.loads((Path(slug_dir) / "upstream.json").read_text(encoding="utf-8"))
    if pin["pinned_commit"] != COMMIT:
        sys.exit(f"upstream.json pins {pin['pinned_commit']}, this module {COMMIT}")
    return pin


def authors_file(pin, upstream=None, fetch=False):
    """The pinned authors' file, sha256-checked; from a checkout, the cache, or the web."""
    info = pin["file"]
    if upstream is not None:
        path = Path(upstream) / info["upstream_path"]
    else:
        path = CACHE / info["upstream_path"]
        if not path.exists():
            if not fetch:
                sys.exit(
                    f"{path} is missing: run with --fetch (needs internet), or pass "
                    f"--upstream DIR, a checkout of {REPOSITORY} at {COMMIT}"
                )
            path.parent.mkdir(parents=True, exist_ok=True)
            url = RAW_URL.format(commit=COMMIT, path=info["upstream_path"])
            with urllib.request.urlopen(url, timeout=60) as response:
                data = response.read()
            if hashlib.sha256(data).hexdigest() != info["sha256"]:
                sys.exit(
                    f"{url}: sha256 {hashlib.sha256(data).hexdigest()} != pinned "
                    f"{info['sha256']}; not written"
                )
            path.write_bytes(data)
    if not path.exists():
        sys.exit(f"{path} does not exist")
    digest = sha256(path)
    if digest != info["sha256"]:
        sys.exit(f"{path}: sha256 {digest} != pinned {info['sha256']} ({REPOSITORY} @ {COMMIT})")
    return path


# ------------------------------------------------------------------------------ BNG2.pl


def find_bng2(bngpath=None):
    root = bngpath or os.environ.get("BNGPATH")
    if not root:
        sys.exit("set BNGPATH to the BioNetGen folder that contains BNG2.pl")
    bng2 = Path(root) / "BNG2.pl"
    if not bng2.exists():
        sys.exit(f"{bng2} not found")
    os.environ["BNGPATH"] = str(root)  # BNG2.pl finds run_network through it
    return bng2


def split_actions(text):
    """(model text through `end model`, the actions that follow it)."""
    match = re.search(r"^\s*end\s+model\s*$", text, flags=re.MULTILINE)
    if match is None:
        sys.exit("no `end model` line")
    return text[: match.end()], text[match.end() :]


def run_bng2(bng2, text, stem, run_dir):
    """Write `text` as `<stem>.bngl` in a fresh `run_dir` and run BNG2.pl on it."""
    run_dir = Path(run_dir)
    if run_dir.exists():
        for old in run_dir.iterdir():
            old.unlink()
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / f"{stem}.bngl").write_text(text, encoding="utf-8")
    proc = subprocess.run(
        ["perl", str(bng2), f"{stem}.bngl"], cwd=run_dir, capture_output=True, text=True
    )
    (run_dir / "BNG2.log").write_text(proc.stdout + proc.stderr, encoding="utf-8")
    if proc.returncode:
        sys.exit(f"BNG2.pl failed in {run_dir}:\n{proc.stdout[-3000:]}\n{proc.stderr[-2000:]}")
    version = re.search(r"BioNetGen version ([\d.]+)", proc.stdout)
    return run_dir, (version.group(1) if version else "unknown")


# -------------------------------------------------------------------------- actions


ACTION_RE = re.compile(r"^\s*(\w+)\s*\((.*)\)\s*$")


def parse_actions(text):
    """The action calls of a BNGL actions section, comments dropped: [(name, argument text)]."""
    calls = []
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line or line.startswith("begin ") or line.startswith("end "):
            continue
        match = ACTION_RE.match(line)
        if match is None:
            sys.exit(f"cannot parse action line {raw!r}")
        calls.append((match.group(1), match.group(2).strip()))
    return calls


def set_concentration_args(arg):
    """`"pattern", value` -> (pattern, value text, whether the value was quoted)."""
    match = re.match(r'^"([^"]+)"\s*,\s*(.+)$', arg)
    if match is None:
        sys.exit(f"cannot parse setConcentration({arg})")
    value = match.group(2).strip()
    quoted = value.startswith('"') and value.endswith('"')
    return match.group(1), value.strip('"'), quoted


def tighten(simulate_arg, tol):
    """Add atol/rtol to a simulate({...}) argument."""
    if "atol" in simulate_arg or "rtol" in simulate_arg:
        sys.exit(f"simulate already sets a tolerance: {simulate_arg}")
    return simulate_arg.replace("{", f"{{atol=>{tol:g},rtol=>{tol:g},", 1)


def render_actions(calls):
    lines = ["", "begin actions"]
    lines += [f"{name}({arg})" for name, arg in calls]
    lines += ["end actions", ""]
    return "\n".join(lines)


def simulate_suffix(arg):
    match = re.search(r'suffix\s*=>\s*"([^"]+)"', arg)
    return match.group(1) if match else None


# ------------------------------------------------------------------------ expressions


FUNCTIONS = {
    "exp": math.exp,
    "ln": math.log,
    "log10": math.log10,
    "log2": math.log2,
    "sqrt": math.sqrt,
    "abs": abs,
    "min": min,
    "max": max,
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "_pi": math.pi,
    "_e": math.e,
}
NAME_RE = re.compile(r"[A-Za-z_]\w*")


def _if(condition, a, b):
    return a if condition else b


def to_python(expr):
    """Translate a BioNetGen math expression into Python source."""
    out = expr.replace("^", "**").replace("&&", " and ").replace("||", " or ")
    out = re.sub(r"(?<![<>=!])!(?!=)", " not ", out)
    return NAME_RE.sub(lambda m: "_if" if m.group(0) == "if" else m.group(0), out)


def evaluate(expr, names):
    return eval(to_python(expr), {"__builtins__": {}, "_if": _if, **FUNCTIONS}, names)


# ------------------------------------------------------------------- species labels


MOLECULE_RE = re.compile(r"[A-Za-z_]\w*\([^)]*\)")


def canonical_species(text):
    """Canonical label of a BioNetGen species string.

    Invariant to molecule order, component order and bond numbering: molecules are grouped by
    their bond-free signature, every ordering within each group is tried with bonds renumbered
    in order of appearance, and the lexicographically smallest rendering wins. The complexes
    in these models hold at most four molecules, so the search is small. A leading `$` (fixed
    species) is kept; compartments do not occur and are rejected.
    """
    fixed = text.startswith("$")
    body = text.lstrip("$")
    if "@" in body:
        raise ValueError(f"compartments are not supported: {text!r}")
    tokens = MOLECULE_RE.findall(body)
    if ".".join(tokens) != body:
        raise ValueError(f"cannot parse species {text!r}")
    molecules = []
    for token in tokens:
        name, inner = re.match(r"^(\w+)\(([^)]*)\)$", token).groups()
        components = []
        for comp in inner.split(",") if inner else []:
            comp, _, bond = comp.partition("!")
            comp, _, state = comp.partition("~")
            components.append((comp, state, bond))
        components.sort(key=lambda c: (c[0], c[1]))
        molecules.append((name, components))

    def signature(mol):
        return mol[0], tuple((c, s, bool(b)) for c, s, b in mol[1])

    order = sorted(range(len(molecules)), key=lambda i: signature(molecules[i]))
    groups = [list(g) for _, g in itertools.groupby(order, key=lambda i: signature(molecules[i]))]
    best = None
    for arrangement in itertools.product(*(itertools.permutations(g) for g in groups)):
        relabel, parts = {}, []
        for i in (i for group in arrangement for i in group):
            name, components = molecules[i]
            rendered = []
            for comp, state, bond in components:
                if bond:
                    bond = relabel.setdefault(bond, str(len(relabel) + 1))
                rendered.append(
                    comp + (f"~{state}" if state else "") + (f"!{bond}" if bond else "")
                )
            parts.append(f"{name}({','.join(rendered)})")
        label = ("$" if fixed else "") + ".".join(parts)
        if best is None or label < best:
            best = label
    return best


# ------------------------------------------------------------------------------ .net


@dataclass
class Net:
    path: Path
    param_expr: dict  # name -> expression text, file order
    species: list  # species strings, file order
    species_expr: list  # initial-amount expression text
    fixed: np.ndarray
    labels: list  # canonical labels
    reactions: list  # (reactants, products, factor, token, rule)
    groups: dict  # name -> [(species index, weight)]
    functions: dict  # name -> expression text, file order
    params: dict = field(default_factory=dict)

    def index(self):
        return {label: i for i, label in enumerate(self.labels)}


def parse_net(path, text=None):
    """Read a BioNetGen `.net` (from `path`, or from `text` with `path` as its name)."""
    params, species, species_expr, reactions, groups, functions = {}, [], [], [], {}, {}
    block = None
    if text is None:
        text = Path(path).read_text(encoding="utf-8")
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("begin "):
            block = line.split()[1]
            continue
        if line.startswith("end "):
            block = None
            continue
        body, _, comment = line.partition("#")
        fields = body.split()
        if not fields:
            continue
        if block == "parameters":
            params[fields[1]] = "".join(fields[2:])
        elif block == "species":
            species.append(fields[1])
            species_expr.append("".join(fields[2:]))
        elif block == "reactions":
            reactants = [int(x) - 1 for x in fields[1].split(",") if x != "0"]
            products = [int(x) - 1 for x in fields[2].split(",") if x != "0"]
            rate = "".join(fields[3:])
            match = re.match(r"^([0-9.eE+-]+)\*(\w+)$", rate)
            factor, token = (float(match.group(1)), match.group(2)) if match else (1.0, rate)
            if not NAME_RE.fullmatch(token):
                sys.exit(f"{path}: unsupported rate law {rate!r}")
            reactions.append((reactants, products, factor, token, comment.strip()))
        elif block == "groups":
            members = []
            for item in fields[2].split(",") if len(fields) > 2 else []:
                weight, _, idx = item.rpartition("*")
                members.append((int(idx) - 1, float(weight) if weight else 1.0))
            groups[fields[1]] = members
        elif block == "functions":
            name = fields[1]
            if not name.endswith("()"):
                sys.exit(f"{path}: function {name!r} takes arguments")
            functions[name[:-2]] = "".join(fields[2:])
    net = Net(
        Path(path),
        params,
        species,
        species_expr,
        np.array([s.startswith("$") for s in species]),
        [canonical_species(s) for s in species],
        reactions,
        groups,
        functions,
    )
    net.params = parameter_values(net)
    return net


def parameter_values(net, overrides=None):
    """Evaluate the parameters in file order; an override replaces the value of its name."""
    overrides = overrides or {}
    unknown = set(overrides) - set(net.param_expr)
    if unknown:
        sys.exit(f"{net.path.name}: no parameters {sorted(unknown)}")
    values = {}
    for name, expr in net.param_expr.items():
        values[name] = (
            float(overrides[name]) if name in overrides else float(evaluate(expr, values))
        )
    return values


def initial_state(net, params=None):
    params = params or net.params
    return np.array([float(evaluate(e, params)) for e in net.species_expr])


class State:
    """Observables, functions and rate coefficients of a network at one state."""

    def __init__(self, net, y, params=None):
        self.net, self.y = net, np.asarray(y, dtype=float)
        self.names = dict(params or net.params)
        for name, members in net.groups.items():
            self.names[name] = sum(w * self.y[i] for i, w in members)
        self.cache = {}
        for name in net.functions:
            self.names[name] = self._caller(name)

    def _caller(self, name):
        def call():
            if name not in self.cache:
                self.cache[name] = float(evaluate(self.net.functions[name], self.names))
            return self.cache[name]

        return call

    def function(self, name):
        return self.names[name]()

    def coefficient(self, reaction):
        _, _, factor, token, _ = reaction
        if token in self.net.functions:
            return factor * self.function(token)
        return factor * float(self.names[token])

    def coefficients(self):
        return np.array([self.coefficient(r) for r in self.net.reactions])

    def fluxes(self):
        k = self.coefficients()
        return np.array(
            [k[j] * np.prod([self.y[i] for i in r[0]]) for j, r in enumerate(self.net.reactions)]
        )

    def rhs(self):
        """dy/dt and each species' gross turnover sum_j |S_ij r_j|."""
        dydt, gross = np.zeros(len(self.y)), np.zeros(len(self.y))
        for flux, (reactants, products, *_rest) in zip(
            self.fluxes(), self.net.reactions, strict=True
        ):
            for i in reactants:
                dydt[i] -= flux
                gross[i] += abs(flux)
            for i in products:
                dydt[i] += flux
                gross[i] += abs(flux)
        dydt[self.net.fixed] = 0.0
        gross[self.net.fixed] = 0.0
        return dydt, gross


def function_dependencies(net):
    """Functions reachable from the reactions' rate laws (the ones that drive dy/dt)."""
    used, stack = set(), [r[3] for r in net.reactions if r[3] in net.functions]
    while stack:
        name = stack.pop()
        if name in used:
            continue
        used.add(name)
        stack += [
            n for n in NAME_RE.findall(net.functions[name]) if n in net.functions and n not in used
        ]
    return used


def rel(a, b, floor=0.0):
    if a == b:
        return 0.0
    return abs(a - b) / max(abs(a), abs(b), floor)


def reaction_keys(net):
    """Reactions grouped by (canonical reactants, canonical products)."""
    table = defaultdict(list)
    for k, (reactants, products, *_rest) in enumerate(net.reactions):
        key = (
            tuple(sorted(net.labels[i] for i in reactants)),
            tuple(sorted(net.labels[i] for i in products)),
        )
        table[key].append(k)
    return table


def compare_networks(authors, job):
    """Species and reaction bijections, and observables, functions, parameters by name."""
    out = {"species": {"authors": len(authors.species), "job": len(job.species)}}
    a_index, j_index = authors.index(), job.index()
    only_a = sorted(set(a_index) - set(j_index))
    only_j = sorted(set(j_index) - set(a_index))
    out["species"].update(
        only_authors=only_a,
        only_job=only_j,
        bijection=not only_a and not only_j,
        same_order=authors.labels == job.labels,
    )
    out["reactions"] = {"authors": len(authors.reactions), "job": len(job.reactions)}
    if only_a or only_j:
        return out
    a_keys, j_keys = reaction_keys(authors), reaction_keys(job)
    unmatched = sorted(set(a_keys) ^ set(j_keys))
    counts = [k for k in a_keys if k in j_keys and len(a_keys[k]) != len(j_keys[k])]
    out["reactions"].update(
        unmatched=[list(map(list, k)) for k in unmatched],
        multiplicity_mismatch=len(counts),
        bijection=not unmatched and not counts,
        same_order_and_rate_law_text=[r[:4] for r in authors.reactions]
        == [r[:4] for r in job.reactions],
    )

    def group_set(net, name):
        return sorted((net.labels[i], w) for i, w in net.groups[name])

    common_g = [g for g in authors.groups if g in job.groups]
    out["observables"] = {
        "only_authors": sorted(set(authors.groups) - set(job.groups)),
        "only_job": sorted(set(job.groups) - set(authors.groups)),
        "differing_definition": sorted(
            g for g in common_g if group_set(authors, g) != group_set(job, g)
        ),
    }
    common_f = [f for f in authors.functions if f in job.functions]
    only_jf = sorted(set(job.functions) - set(authors.functions))
    out["functions"] = {
        "only_authors": sorted(set(authors.functions) - set(job.functions)),
        "only_job": only_jf,
        "only_job_used_by_rate_laws": sorted(set(only_jf) & function_dependencies(job)),
        "differing_text": {
            f: {"authors": authors.functions[f], "job": job.functions[f]}
            for f in common_f
            if authors.functions[f] != job.functions[f]
        },
    }
    common_p = [p for p in authors.params if p in job.params]
    out["parameters"] = {
        "only_authors": sorted(set(authors.params) - set(job.params)),
        "only_job": {p: job.params[p] for p in job.params if p not in authors.params},
        "differing_value": {
            p: {"authors": authors.params[p], "job": job.params[p]}
            for p in common_p
            if authors.params[p] != job.params[p]
        },
        "differing_text": sorted(p for p in common_p if authors.param_expr[p] != job.param_expr[p]),
    }
    return out


# ------------------------------------------------------------------------ trajectories


def load_dat(path):
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().lstrip("#").split()
    return header, np.loadtxt(path, ndmin=2)


def compare_columns(a, b, floor):
    """max_t |a - b| / max(peak, floor) per column, peak over both runs."""
    peak = np.maximum(np.abs(a).max(axis=0), np.abs(b).max(axis=0))
    diff = np.abs(a - b).max(axis=0)
    return diff / np.maximum(peak, floor), diff


def compare_runs(a_dir, a_stem, j_dir, j_stem, suffix, perm, floor, skip_species=()):
    """Species (.cdat, through the bijection) and .gdat columns of one simulate() suffix."""
    out = {}
    ha, A = load_dat(Path(a_dir) / f"{a_stem}_{suffix}.cdat")
    hj, J = load_dat(Path(j_dir) / f"{j_stem}_{suffix}.cdat")
    if A.shape[0] != J.shape[0] or not np.array_equal(A[:, 0], J[:, 0]):
        sys.exit(f"{suffix}: the two runs report different time points")
    J = J[:, 1:][:, perm]  # job species in authors' order
    A = A[:, 1:]
    keep = [i for i in range(A.shape[1]) if i not in skip_species]
    r, d = compare_columns(A[:, keep], J[:, keep], floor)
    worst = int(np.argmax(r))
    out["species"] = {
        "compared": len(keep),
        "max_rel_to_peak": float(r.max()),
        "max_abs": float(d.max()),
        "worst": f"S{keep[worst] + 1}",
        "identical": bool(np.array_equal(A[:, keep], J[:, keep])),
    }
    ga, GA = load_dat(Path(a_dir) / f"{a_stem}_{suffix}.gdat")
    gj, GJ = load_dat(Path(j_dir) / f"{j_stem}_{suffix}.gdat")
    missing = [c for c in ga if c not in gj]
    if missing:
        sys.exit(f"{suffix}: job .gdat lacks {missing}")
    cols = ga[1:]
    GA = GA[:, 1:]
    GJ = GJ[:, [gj.index(c) for c in cols]]
    r, d = compare_columns(GA, GJ, floor)
    worst = int(np.argmax(r))
    out["observables_and_functions"] = {
        "compared": len(cols),
        "max_rel_to_peak": float(r.max()),
        "max_abs": float(d.max()),
        "worst": cols[worst],
        "identical": bool(np.array_equal(GA, GJ)),
    }
    out["time_points"] = int(A.shape[0])
    return out


def verdict(value, tol):
    return {"value": value, "tolerance": tol, "pass": bool(value <= tol)}


# ----------------------------------------------------------------------- translation


@dataclass
class Gate:
    """How the slug's file expresses an authors' setConcentration on one species.

    `params(value, text)` gives the slug parameters (name -> expression text) equivalent to
    the species at `value`, whose exact BioNetGen expression is `text`; `keep_species` says
    whether the species itself is still set as the authors set it.
    """

    params: object
    keep_species: bool


def species_index(net, pattern):
    index = net.index()
    for label in (canonical_species(pattern), "$" + canonical_species(pattern)):
        if label in index:
            return index[label]
    sys.exit(f"{net.path.name}: no species {pattern}")


def expand(expr, net):
    """`expr` with every parameter of `net` replaced, recursively, by its own expression.

    The result is a closed expression in literals that BioNetGen's run_network evaluates to
    the same double as `expr` in `net`.
    """

    def sub(match):
        name = match.group(0)
        if name in net.param_expr:
            return f"({expand(net.param_expr[name], net)})"
        return name

    return NAME_RE.sub(sub, expr)


def translate_actions(calls, authors, job, gates):
    """The authors' action list rewritten for the slug's file.

    BNG2.pl keeps a quoted setConcentration value as an expression that run_network evaluates
    in full precision, but writes a bare number, or any setParameter number, to 15 significant
    digits. So each value keeps the authors' form: a bare number is copied as it is; a quoted
    expression is copied as it is when every parameter it names has the same value in both
    files, and is otherwise expanded into the authors' literals (`expand`). A gated species is
    tracked through setConcentration / saveConcentrations / resetConcentrations as BioNetGen
    tracks it, and before each simulate the gate parameters are set, as quoted expressions, to
    match its current value. Returns the new call list and, per simulate suffix, the authors'
    value of each gated species at its start.
    """
    y0 = initial_state(authors)
    current = {
        p: (
            float(y0[species_index(authors, p)]),
            expand(authors.species_expr[species_index(authors, p)], authors),
        )
        for p in gates
    }
    saved = dict(current)
    out, gated_at = [], {}
    for name, arg in calls:
        if name == "writeSBML":
            continue
        if name == "setConcentration":
            pattern, text, quoted = set_concentration_args(arg)
            value = float(evaluate(text, authors.params))
            names = [n for n in NAME_RE.findall(text) if n in authors.params]
            same = all(n in job.params and job.params[n] == authors.params[n] for n in names)
            if quoted and not same:
                text = expand(text, authors)
            if pattern in gates:
                current[pattern] = (value, expand(text, authors) if quoted else text)
                if not gates[pattern].keep_species:
                    continue
            out.append((name, f'"{pattern}","{text}"' if quoted else f'"{pattern}",{text}'))
        elif name == "saveConcentrations":
            saved = dict(current)
            out.append((name, arg))
        elif name == "resetConcentrations":
            current = dict(saved)
            out.append((name, arg))
        elif name == "simulate":
            for pattern, (value, text) in current.items():
                for param, ptext in gates[pattern].params(value, text).items():
                    out.append(("setParameter", f'"{param}","{ptext}"'))
            gated_at[simulate_suffix(arg)] = {p: v for p, (v, _t) in current.items()}
            out.append((name, arg))
        elif name == "generate_network":
            out.append((name, arg))
        else:
            sys.exit(f"no translation for action {name}({arg})")
    return out, gated_at


def with_tolerance(calls, tol):
    return [(n, tighten(a, tol) if n == "simulate" else a) for n, a in calls]


def set_t_end(calls, suffix, t_end):
    """A copy of `calls` with simulate(suffix) run to `t_end` instead."""
    out, hit = [], False
    for name, arg in calls:
        if name == "simulate" and simulate_suffix(arg) == suffix:
            arg, n = re.subn(r"t_end\s*=>\s*[0-9.eE+-]+", f"t_end=>{t_end:g}", arg)
            hit = n == 1
        out.append((name, arg))
    if not hit:
        sys.exit(f"no simulate suffix {suffix}")
    return out


def drop_calls(calls, predicate):
    return [(n, a) for n, a in calls if not predicate(n, a)]


# --------------------------------------------------------------------------- driver


@dataclass
class Spec:
    """What a slug's driver supplies.

    gates       authors' species pattern -> Gate (see `translate_actions`)
    conditions  [(name, sampler, compare_coefficients)]; `sampler(rng, authors, job)` returns
                the authors' state, the slug's state in the authors' species order, and the
                slug's parameter overrides that put it in the same condition
    expected    the documented differences between the two files (see
                `expected_differences`), asserted to be the only ones
    probes      name -> (f(calls) -> calls, [(authors suffix, probe suffix)] whose end states
                are compared); reported, not asserted
    watch       .gdat columns whose end values the probes report
    """

    slug: str
    job_bngl: Path
    gates: dict
    conditions: list
    expected: dict
    probes: dict
    notes: dict = field(default_factory=dict)
    watch: list = field(default_factory=list)


TOLERANCES = {
    "rate_coefficient_max_rel": (
        1e-12,
        "both sides are the same IEEE doubles: the slug's "
        "file only multiplies some constants by a gate equal to 1 or "
        "replaces a fixed species by a parameter of the same value",
    ),
    "rhs_max_rel_to_gross_flux": (
        1e-12,
        "dy/dt of both networks at the same states, "
        "relative to each species' gross turnover sum_j |S_ij r_j|",
    ),
    "common_functions_max_rel": (
        1e-12,
        "every function both files define, at the same "
        "states (the observable and function columns of the .gdat)",
    ),
    "common_observables_max_rel": (1e-12, "every observable both files define, at the same states"),
    "trajectory_max_rel_to_peak": (
        1e-9,
        "max_t |authors - slug| / max(peak, 1e-8) for every "
        "species, observable and function of a phase. The two runs "
        "give the same run_network binary the same equations, "
        "parameter doubles and initial state (the translation keeps "
        "the authors' expressions, which BNG2.pl passes on in full "
        "precision), and CVODE is deterministic, so they must agree "
        "to the 13 significant digits BioNetGen prints; 1e-9 allows "
        "only last-digit rounding. 1e-8 is the solver's absolute "
        "tolerance, below which values are solver noise. The bound "
        "is far below the solver's own error (solver_error), so even "
        "an equivalent reformulation that changed CVODE's step "
        "sequence would be flagged: the claim is identity",
    ),
}


def sample_lognormal(rng, net, spread=(1e-3, 3.0), signed=()):
    """A positive log-uniform state; species whose label is in `signed` get U(-10, 10)."""
    lo, hi = np.log(spread[0]), np.log(spread[1])
    y = np.exp(rng.uniform(lo, hi, len(net.species)))
    for i, label in enumerate(net.labels):
        if label in signed:
            y[i] = rng.uniform(-10.0, 10.0)
    return y


def network_check(spec, authors, job):
    """Structure, then equations under each of the spec's sampled conditions."""
    struct = compare_networks(authors, job)
    if not (struct["species"]["bijection"] and struct["reactions"]["bijection"]):
        return struct, False
    perm = np.array([job.index()[label] for label in authors.labels])
    rng = np.random.default_rng(20221005)
    common_f = [f for f in authors.functions if f in job.functions]
    common_g = [g for g in authors.groups if g in job.groups]
    struct["conditions"] = {}
    ok = True
    for name, sampler, coefficients in spec.conditions:
        worst = defaultdict(float)
        for _ in range(N_STATES):
            y_a, y_j_authors_order, overrides = sampler(rng, authors, job)
            y_j = np.empty_like(y_j_authors_order)
            y_j[perm] = y_j_authors_order
            s_a = State(authors, y_a)
            s_j = State(job, y_j, parameter_values(job, overrides))
            if coefficients:
                worst["rate_coefficient_max_rel"] = max(
                    worst["rate_coefficient_max_rel"],
                    coefficient_gap(authors, job, s_a.coefficients(), s_j.coefficients()),
                )
            d_a, g_a = s_a.rhs()
            d_j, g_j = s_j.rhs()
            d_j, g_j = d_j[perm], g_j[perm]
            scale = np.maximum(np.maximum(g_a, g_j), 1e-300)
            worst["rhs_max_rel_to_gross_flux"] = max(
                worst["rhs_max_rel_to_gross_flux"], float(np.max(np.abs(d_a - d_j) / scale))
            )
            for f in common_f:
                worst["common_functions_max_rel"] = max(
                    worst["common_functions_max_rel"], rel(s_a.function(f), s_j.function(f))
                )
            for g in common_g:
                worst["common_observables_max_rel"] = max(
                    worst["common_observables_max_rel"], rel(s_a.names[g], s_j.names[g])
                )
        entry = {"states": N_STATES}
        for key, value in worst.items():
            tol = TOLERANCES[key][0]
            entry[key] = verdict(value, tol)
            ok &= entry[key]["pass"]
        struct["conditions"][name] = entry
    return struct, ok


def coefficient_gap(authors, job, k_a, k_j):
    """Largest relative gap between matched reactions' coefficient multisets."""
    ta, tj = reaction_keys(authors), reaction_keys(job)
    return max(
        rel(x, z)
        for key, idx in ta.items()
        for x, z in zip(sorted(k_a[idx]), sorted(k_j[tj[key]]), strict=True)
    )


def expected_differences(spec, struct):
    """The documented differences between the two files, asserted to be the only ones."""
    found = {
        "parameters_only_job": sorted(struct["parameters"]["only_job"]),
        "parameters_only_authors": struct["parameters"]["only_authors"],
        "parameters_differing_value": sorted(struct["parameters"]["differing_value"]),
        "observables_only_job": struct["observables"]["only_job"],
        "observables_only_authors": struct["observables"]["only_authors"],
        "observables_differing_definition": struct["observables"]["differing_definition"],
        "functions_only_job_used_by_rate_laws": struct["functions"]["only_job_used_by_rate_laws"],
        "functions_only_authors": struct["functions"]["only_authors"],
        "functions_differing_text": sorted(struct["functions"]["differing_text"]),
    }
    checks = {}
    for key, value in found.items():
        want = sorted(spec.expected.get(key, []))
        checks[key] = {"found": value, "documented": want, "pass": value == want}
    return checks


def main(spec, argv=None):
    import argparse

    parser = argparse.ArgumentParser(description=f"Check {spec.slug} against the authors' file")
    parser.add_argument("--upstream", type=Path, help=f"a checkout of {REPOSITORY} at {COMMIT}")
    parser.add_argument(
        "--fetch",
        action="store_true",
        help="download the pinned authors' file into the git-ignored cache",
    )
    parser.add_argument("--bngpath", help="BioNetGen folder with BNG2.pl (default $BNGPATH)")
    parser.add_argument(
        "--workdir", type=Path, help="where BNG2.pl runs (default dev/runs/<slug>, git-ignored)"
    )
    args = parser.parse_args(argv)

    slug_dir = spec.job_bngl.parent
    pin = load_pin(slug_dir)
    authors_path = authors_file(pin, args.upstream, args.fetch)
    bng2 = find_bng2(args.bngpath)
    work = args.workdir or (REPO_ROOT / "dev" / "runs" / spec.slug)
    work.mkdir(parents=True, exist_ok=True)

    authors_text = authors_path.read_text(encoding="utf-8")
    job_text = spec.job_bngl.read_text(encoding="utf-8")
    a_model, a_actions = split_actions(authors_text)
    j_model, j_actions = split_actions(job_text)
    if parse_actions(j_actions):
        sys.exit(f"{spec.job_bngl.name} has actions after `end model`")
    a_stem, j_stem = authors_path.stem, spec.job_bngl.stem

    # 1. The authors' file, verbatim: its own actions are the published protocol.
    a_dir, version = run_bng2(bng2, authors_text, a_stem, work / "authors")
    authors = parse_net(a_dir / f"{a_stem}.net")
    calls = parse_actions(a_actions)

    # 2. The slug's file under the same protocol, translated; and the tight replicate.
    # The slug's own network first: the translation needs its parameter values.
    n_dir, _ = run_bng2(
        bng2,
        j_model + "\n" + render_actions([("generate_network", "{overwrite=>1}")]),
        j_stem,
        work / "job_network",
    )
    job = parse_net(n_dir / f"{j_stem}.net")
    translated, gated_at = translate_actions(calls, authors, job, spec.gates)
    j_dir, _ = run_bng2(bng2, j_model + "\n" + render_actions(translated), j_stem, work / "job")
    t_dir, _ = run_bng2(
        bng2,
        j_model + "\n" + render_actions(with_tolerance(translated, TIGHT_TOL)),
        j_stem,
        work / "job_tight",
    )

    report = {
        "slug": spec.slug,
        "generated_by": f"pybnf-jobs/Rukhlenko-2022/{spec.slug}/compare_upstream.py",
        "bionetgen_version": version,
        "job_file": {"path": spec.job_bngl.name, "sha256": sha256(spec.job_bngl)},
        "authors_file": {
            "repository": REPOSITORY,
            "commit": COMMIT,
            "path": pin["file"]["upstream_path"],
            "sha256": sha256(authors_path),
        },
        "protocol": {
            "authors_actions": [f"{n}({a})" for n, a in calls],
            "job_actions": [f"{n}({a})" for n, a in translated],
        },
        "tolerances": {k: {"value": v, "reason": r} for k, (v, r) in TOLERANCES.items()},
        "notes": spec.notes,
    }

    # 3. Network identity.
    struct, ok = network_check(spec, authors, job)
    if "parameters" in struct:
        struct["documented_differences"] = expected_differences(spec, struct)
        ok &= all(c["pass"] for c in struct["documented_differences"].values())
    else:
        ok = False
    report["network"] = struct

    # 4. Trajectories, phase by phase.
    if struct["species"]["bijection"]:
        perm = [job.index()[label] for label in authors.labels]
        tol = TOLERANCES["trajectory_max_rel_to_peak"][0]
        report["trajectories"], report["solver_error"] = {}, {}
        for suffix, gated in gated_at.items():
            skip, designed = [], {}
            for pattern, value in gated.items():
                gate = spec.gates[pattern]
                if gate.keep_species:
                    continue
                j_seed = float(initial_state(job)[species_index(job, pattern)])
                if value != j_seed:
                    i = species_index(authors, pattern)
                    skip.append(i)
                    designed[pattern] = designed_species(
                        a_dir, a_stem, j_dir, j_stem, suffix, i, perm[i], value, j_seed
                    )
            res = compare_runs(a_dir, a_stem, j_dir, j_stem, suffix, perm, BNG_DEFAULT_TOL, skip)
            res["species_max_rel_to_peak"] = verdict(res["species"]["max_rel_to_peak"], tol)
            res["observables_max_rel_to_peak"] = verdict(
                res["observables_and_functions"]["max_rel_to_peak"], tol
            )
            if designed:
                res["gated_species"] = designed
                ok &= all(d["pass"] for d in designed.values())
            ok &= (
                res["species_max_rel_to_peak"]["pass"]
                and res["observables_max_rel_to_peak"]["pass"]
            )
            report["trajectories"][suffix] = res
            ident = list(range(len(job.species)))
            err = compare_runs(t_dir, j_stem, j_dir, j_stem, suffix, ident, BNG_DEFAULT_TOL)
            report["solver_error"][suffix] = {
                "species_max_rel_to_peak": err["species"]["max_rel_to_peak"],
                "species_max_abs": err["species"]["max_abs"],
                "species_worst": err["species"]["worst"],
                "observables_max_rel_to_peak": err["observables_and_functions"]["max_rel_to_peak"],
                "observables_max_abs": err["observables_and_functions"]["max_abs"],
                "observables_worst": err["observables_and_functions"]["worst"],
            }

    # 5. Protocol sensitivity: the same protocol with modified pre-equilibration (reported).
    report["protocol_sensitivity"] = {}
    if struct["species"]["bijection"]:
        for name, (modify, pairs) in spec.probes.items():
            p_dir, _ = run_bng2(
                bng2,
                j_model + "\n" + render_actions(modify(translated)),
                j_stem,
                work / f"probe_{name}",
            )
            report["protocol_sensitivity"][name] = probe_report(
                spec, authors, a_dir, a_stem, p_dir, j_stem, perm, gated_at, pairs
            )

    report["verdict"] = "PASS" if ok else "FAIL"
    out = slug_dir / "upstream_comparison.json"
    out.write_text(json.dumps(report, indent=1, sort_keys=False) + "\n", encoding="utf-8")
    summarize(report)
    return 0 if ok else 1


def designed_species(a_dir, a_stem, j_dir, j_stem, suffix, ia, ij, a_value, j_value):
    """A gated species is constant at its designed value on each side."""
    _, A = load_dat(Path(a_dir) / f"{a_stem}_{suffix}.cdat")
    _, J = load_dat(Path(j_dir) / f"{j_stem}_{suffix}.cdat")
    a, j = A[:, 1 + ia], J[:, 1 + ij]
    return {
        "authors_value": a_value,
        "job_value": j_value,
        "authors_max_dev": float(np.abs(a - a_value).max()),
        "job_max_dev": float(np.abs(j - j_value).max()),
        "pass": bool(
            np.abs(a - a_value).max() <= BNG_DEFAULT_TOL
            and np.abs(j - j_value).max() <= BNG_DEFAULT_TOL
        ),
    }


def probe_report(spec, authors, a_dir, a_stem, p_dir, j_stem, perm, gated_at, pairs):
    out = {"states": {}, "trajectories": {}}
    for a_suffix, p_suffix in pairs:
        _, A = load_dat(Path(a_dir) / f"{a_stem}_{a_suffix}.cdat")
        _, P = load_dat(Path(p_dir) / f"{j_stem}_{p_suffix}.cdat")
        a_end, p_end = A[-1, 1:], P[-1, 1:][perm]
        drift = np.abs(P[-1, 1:] - P[-2, 1:]) / np.maximum(np.abs(P[-1, 1:]), BNG_DEFAULT_TOL)
        r = np.abs(a_end - p_end) / np.maximum(
            np.maximum(np.abs(a_end), np.abs(p_end)), BNG_DEFAULT_TOL
        )
        gated = {
            species_index(authors, pat) for pat in spec.gates if not spec.gates[pat].keep_species
        }
        r[list(gated)] = 0.0
        worst = int(np.argmax(r))
        ga, GA = load_dat(Path(a_dir) / f"{a_stem}_{a_suffix}.gdat")
        gp, GP = load_dat(Path(p_dir) / f"{j_stem}_{p_suffix}.gdat")
        out["states"][f"{a_suffix} vs {p_suffix}"] = {
            "end_state_max_rel": float(r.max()),
            "worst": f"S{worst + 1}",
            "worst_species": authors.species[worst],
            "probe_last_interval_drift_max_rel": float(drift.max()),
            "probe_end_time": float(P[-1, 0]),
            "authors_end_time": float(A[-1, 0]),
            "watched_at_end": {
                c: {"authors": float(GA[-1, ga.index(c)]), "probe": float(GP[-1, gp.index(c)])}
                for c in spec.watch
            },
        }
    for suffix in gated_at:
        path = Path(p_dir) / f"{j_stem}_{suffix}.gdat"
        if not path.exists():
            continue
        ga, GA = load_dat(Path(a_dir) / f"{a_stem}_{suffix}.gdat")
        gp, GP = load_dat(path)
        if GA.shape[0] != GP.shape[0] or not np.array_equal(GA[:, 0], GP[:, 0]):
            continue  # a phase whose duration the probe changed
        cols = ga[1:]
        r, _ = compare_columns(GA[:, 1:], GP[:, [gp.index(c) for c in cols]], BNG_DEFAULT_TOL)
        order = np.argsort(-r)[:3]
        out["trajectories"][suffix] = {
            "max_rel_to_peak": float(r.max()),
            "worst": {cols[i]: float(r[i]) for i in order},
            "watched_at_end": {
                c: {"authors": float(GA[-1, ga.index(c)]), "probe": float(GP[-1, gp.index(c)])}
                for c in spec.watch
            },
        }
    return out


def summarize(report):
    net = report["network"]
    print(
        f"{report['slug']}: BioNetGen {report['bionetgen_version']}; authors' "
        f"{report['authors_file']['path']} @ {report['authors_file']['commit'][:7]}"
    )
    print(
        f"  species {net['species']['authors']}/{net['species']['job']} bijection "
        f"{net['species']['bijection']}; reactions {net['reactions']['authors']}/"
        f"{net['reactions']['job']} bijection {net['reactions']['bijection']}"
    )
    for name, cond in net.get("conditions", {}).items():
        parts = [
            f"{k} {v['value']:.2e} (tol {v['tolerance']:g})"
            for k, v in cond.items()
            if isinstance(v, dict)
        ]
        print(f"  [{name}] " + "; ".join(parts))
    for key, check in net.get("documented_differences", {}).items():
        if check["found"] or not check["pass"]:
            print(f"  {key}: {check['found']} {'ok' if check['pass'] else 'UNDOCUMENTED'}")
    for suffix, res in report.get("trajectories", {}).items():
        err = report["solver_error"][suffix]
        same = res["species"]["identical"] and res["observables_and_functions"]["identical"]
        solver = max(err["species_max_rel_to_peak"], err["observables_max_rel_to_peak"])
        print(
            f"  {suffix}: species {res['species_max_rel_to_peak']['value']:.2e}, "
            f"observables {res['observables_max_rel_to_peak']['value']:.2e} "
            f"(identical {same}); solver error {solver:.2e}"
        )
    for name, probe in report.get("protocol_sensitivity", {}).items():
        for key, st in probe["states"].items():
            print(
                f"  probe {name}: {key} end-state max rel {st['end_state_max_rel']:.3g} "
                f"({st['worst_species']}); drift {st['probe_last_interval_drift_max_rel']:.1e}; "
                f"watched {st['watched_at_end']}"
            )
        for suffix, t in probe["trajectories"].items():
            print(f"    {suffix}: max rel {t['max_rel_to_peak']:.3g} {t['worst']}")
    print(f"  verdict: {report['verdict']}")
