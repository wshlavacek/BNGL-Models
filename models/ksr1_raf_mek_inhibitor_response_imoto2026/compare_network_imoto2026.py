#!/usr/bin/env python
"""Check the BioNetGen-generated network against the deposited Imoto2026 SBML.

The authors publish their model only as a flattened SBML export of the PySB
network (`dev/papers/Imoto2026/{mcf7,psn1}_sbml/`), not as a rule set, so the
`.bngl` in this folder is a reconstruction. This module is the test of that
reconstruction: it puts both networks into the same canonical form and compares
them exhaustively —

    species        every species, matched up to molecule order and bond numbering
    reactions      every unidirectional reaction, as a (reactants, products) pair
    rate constants every reaction's rate constant, evaluated numerically
    RAF activity   the per-species weights of the aggregate RAF kinase activity
    MEK kinetics   the two Michaelis-Menten substrate pools and their Km values

The thermodynamic factors are set to distinct non-unit values (`TEST`) before the
rate constants are compared, so a misplaced energy pattern cannot hide behind a
factor that happens to be 1 in the deposited file.

Three species of the deposited export are absent from the reconstruction by
design and are excluded here: the BRAF V600E allele (zero abundance in both cell
lines, and given no dimerization rules by the export), the `tick()` clock of a
time-dependent RAS-GTP input whose amplitude is zero, and the unused Trametinib
and KSR1-inhibitor species. See the `#@note:` in the `.bngl`.

Usage
-----
    python compare_network_imoto2026.py <generated.net> [<model.sbml>]
"""

from __future__ import annotations

import itertools
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_SBML = HERE.parents[1] / "dev/papers/Imoto2026/mcf7_sbml/ksrmodel_MCF7.sbml"

# Sites present in the deposited SBML but inert at the published parameterization.
DROP_SITES = {"V600", "trametinib", "ksri"}

# Distinct, non-unit test values so that a misplaced energy pattern cannot hide.
TEST = dict(
    p1=49.0, p2=0.004, p3=9.0, p3KSR=7.0,
    pEB=25.0, pEC=17.0, pEK=11.0,
    fKR=1.0e-4, fM=0.1, dKMP=3.0,
    fa=0.005, fb=0.01, g1a=0.044, g1b=0.429,
    g2a=4.346, g2b=102.95, g3a=1.7, g3b=2.3,
    fKa=0.025, fKb=50.0,
)

# Both networks carry the two Michaelis-Menten rate laws as state-dependent
# functions. They are stood in for by markers so the reaction sets still align;
# the functions themselves are compared by `compare_functions`.
FUNC_MARKERS = {
    "Michaelis_Menten_kinetics_uMEK_phosphorylation": 7.0, "MM_uMEK": 7.0,
    "Michaelis_Menten_kinetics_pMEK_phosphorylation": 11.0, "MM_pMEK": 11.0,
}

SAFE = {"ln": math.log, "log": math.log10, "log10": math.log10, "exp": math.exp,
        "sqrt": math.sqrt, "pow": pow, "piecewise": lambda a, c, b: a if c else b}

# The Table S4 activity classes, as (observable name, weight relative to one
# active BRAF monomer).
ACTIVITY_WEIGHTS = {
    "Obs_Act_BRAF_mono": 1, "Obs_Act_BRAF_act_nonact": 5, "Obs_Act_BRAF_homodimer": 10,
    "Obs_Act_RAF1_mono_off": 0.05, "Obs_Act_RAF1_mono_on": 0.2,
    "Obs_Act_RAF1_homodimer_off": 0.5, "Obs_Act_RAF1_homodimer_half": 1.25,
    "Obs_Act_RAF1_homodimer_on": 2, "Obs_Act_heterodimer_partial": 15,
    "Obs_Act_heterodimer_full": 30, "Obs_Act_BRAF_BRAF_inh": 5,
    "Obs_Act_RAF1_RAF1_inh_off": 0.25, "Obs_Act_RAF1_RAF1_inh_on": 1,
    "Obs_Act_RAF1_off_BRAF_inh": 8, "Obs_Act_RAF1_off_inh_BRAF": 8,
    "Obs_Act_RAF1_on_BRAF_inh": 15, "Obs_Act_RAF1_on_inh_BRAF": 15,
    "Obs_Act_BRAF_KSR": 10, "Obs_Act_RAF1_KSR_off": 0.5, "Obs_Act_RAF1_KSR_on": 2,
}


# ------------------------------------------------------------- canonical species


def _split_top(s, sep):
    out, depth, cur = [], 0, ""
    for ch in s:
        depth += (ch == "(") - (ch == ")")
        if ch == sep and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    out.append(cur)
    return out


def parse_bngl_species(text):
    """'@cell::A(x!1,y~u).B(z!1)' -> [(name, {site: (state, bond)}), ...]."""
    text = re.sub(r"^@\w+::", "", text.strip())
    mols = []
    for part in _split_top(text, "."):
        part = re.sub(r"@\w+", "", part.strip()).lstrip("$")
        m = re.match(r"^(\w+)\((.*)\)$", part)
        sites = {}
        for tok in filter(None, (t.strip() for t in m.group(2).split(","))):
            sm = re.match(r"^(\w+)(~\w+)?(![\w+?]+)?$", tok)
            if sm.group(1) in DROP_SITES:
                continue
            sites[sm.group(1)] = (sm.group(2)[1:] if sm.group(2) else None,
                                  sm.group(3)[1:] if sm.group(3) else None)
        mols.append((m.group(1), sites))
    return mols


def parse_pysb_species(text):
    """A PySB species string into the same structure."""
    mols = []
    for part in text.split(" ._br_"):
        m = re.match(r"^(\w+)\((.*)\)$", part.strip())
        sites = {}
        for tok in filter(None, (t.strip() for t in m.group(2).split(","))):
            site, _, val = tok.partition("=")
            site, val = site.strip(), val.strip()
            if site in DROP_SITES:
                continue
            if val == "None":
                sites[site] = (None, None)
            elif val.startswith("'"):
                sites[site] = (val.strip("'"), None)
            else:
                sites[site] = (None, val)
        mols.append((m.group(1), sites))
    return mols


def canonical(mols):
    """A string for a complex, invariant to molecule order and bond numbering."""
    n = len(mols)

    def render(order):
        bondmap, nxt = {}, [1]
        for i in order:
            for site in sorted(mols[i][1]):
                bond = mols[i][1][site][1]
                if bond and bond not in ("+", "?") and bond not in bondmap:
                    bondmap[bond] = nxt[0]
                    nxt[0] += 1
        parts = []
        for i in order:
            name, sites = mols[i]
            toks = []
            for site in sorted(sites):
                state, bond = sites[site]
                t = site + (f"~{state}" if state is not None else "")
                if bond:
                    t += "!" + (str(bondmap[bond]) if bond in bondmap else bond)
                toks.append(t)
            parts.append(f"{name}({','.join(toks)})")
        return ".".join(parts)

    # Only molecules with identical local labels can be permuted; the rest are
    # pinned by the sorted label order, so the search stays tiny.
    local = [(name, tuple(sorted((s, v[0], v[1] is not None) for s, v in sites.items())))
             for name, sites in mols]
    groups = defaultdict(list)
    for i in range(n):
        groups[local[i]].append(i)
    blocks = [groups[k] for k in sorted(groups)]
    best = None
    for combo in itertools.product(*(itertools.permutations(b) for b in blocks)):
        r = render([i for blk in combo for i in blk])
        if best is None or r < best:
            best = r
    return best


# ------------------------------------------------------------------- .net reader


def read_net(path):
    """(species, reactions, params, groups) from a BioNetGen `.net` file."""
    species, reactions, params, groups = [], [], {}, {}
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
            params[f[1] if indexed else f[0]] = f[2] if indexed else f[1]
        elif section == "species":
            species.append((f[1], f[2]) if indexed else (f[0], f[1]))
        elif section == "reactions":
            reactions.append(([int(x) for x in f[1].split(",")],
                              [int(x) for x in f[2].split(",")], f[3]))
        elif section == "groups":
            off = 1 if indexed else 0
            d = {}
            for t in (f[off + 1] if len(f) > off + 1 else "").split(","):
                if not t:
                    continue
                c, i = t.split("*") if "*" in t else ("1", t)
                d[int(i)] = d.get(int(i), 0.0) + float(c)
            groups[f[off]] = d
    return species, reactions, params, groups


# ------------------------------------------------------------------ SBML reading


def sbml_env(model, overrides):
    """Evaluate the SBML constants and assignment rules to numbers."""
    env = dict(SAFE)
    env.update(FUNC_MARKERS)
    for p in model.getListOfParameters():
        if p.getConstant():
            env[p.getId()] = float(p.getValue())
    env.update(overrides)
    import libsbml
    pending = [(r.getVariable(), libsbml.formulaToL3String(r.getMath()))
               for r in model.getListOfRules()]
    pending = [(v, f) for v, f in pending if "__s" not in f]
    for _ in range(10):
        rest = []
        for v, f in pending:
            try:
                env[v] = eval(f.replace("^", "**"), {"__builtins__": {}}, env)
            except Exception:
                rest.append((v, f))
        pending = rest
        if not pending:
            break
    return env


def monomials(law):
    """Split a PySB net-rate expression into (sign, coefficient, [species])."""
    parts, depth, cur, sign, i = [], 0, "", 1, 0
    while i < len(law):
        ch = law[i]
        depth += (ch == "(") - (ch == ")")
        if depth == 0 and ch in "+-" and cur.strip():
            parts.append((sign, cur))
            sign, cur, i = (1 if ch == "+" else -1), "", i + 1
            continue
        if depth == 0 and ch == "-" and not cur.strip():
            sign, i = -sign, i + 1
            continue
        cur += ch
        i += 1
    parts.append((sign, cur))

    out = []
    for sg, expr in parts:
        expr = expr.strip()
        if not expr:
            continue
        while expr.startswith("(") and expr.endswith(")") and _balanced(expr[1:-1]):
            expr = expr[1:-1].strip()
        factors = []
        for f in _split_top(expr, "*"):
            f = f.strip()
            pw = re.fullmatch(r"(__s\d+)\^(\d+)", f)
            factors += [pw.group(1)] * int(pw.group(2)) if pw else [f]
        specs = [f for f in factors if re.fullmatch(r"__s\d+", f)]
        coef = [f for f in factors if not re.fullmatch(r"__s\d+", f)]
        out.append((sg, "*".join(coef), specs))
    return out


def _balanced(s):
    d = 0
    for ch in s:
        d += (ch == "(") - (ch == ")")
        if d < 0:
            return False
    return d == 0


def _stoich(ref):
    v = ref.getStoichiometry()
    return 1 if v is None or v != v else int(v)


# ------------------------------------------------------------------- comparisons


def load(net_path, sbml_path=DEFAULT_SBML, overrides=TEST):
    """Canonicalized state of both networks, ready for the compare_* functions."""
    import libsbml
    model = libsbml.readSBMLFromFile(str(sbml_path)).getModel()
    names = {s.getId(): s.getName() for s in model.getListOfSpecies()}
    scanon = {k: (None if ("V600='E'" in v or v == "tick()") else canonical(parse_pysb_species(v)))
              for k, v in names.items()}
    species, reactions, params, groups = read_net(net_path)
    bcanon = {i + 1: canonical(parse_bngl_species(s)) for i, (s, _) in enumerate(species)}
    return dict(model=model, names=names, scanon=scanon, env=sbml_env(model, overrides),
                species=species, reactions=reactions, params=params, groups=groups,
                bcanon=bcanon, overrides=dict(overrides))


def compare_species(ctx):
    sb = {c for c in ctx["scanon"].values() if c is not None}
    bn = set(ctx["bcanon"].values())
    extra_sbml = {c for c in sb - bn}
    return dict(n_sbml=len(sb), n_bng=len(bn),
                only_sbml=sorted(extra_sbml), only_bng=sorted(bn - sb))


def _sbml_reactions(ctx):
    import libsbml
    out = defaultdict(float)
    for r in ctx["model"].getListOfReactions():
        if r.getName() in ("time", "RAS_signal"):
            continue
        reac, prod = [], []
        for x in r.getListOfReactants():
            reac += [ctx["scanon"][x.getSpecies()]] * _stoich(x)
        for x in r.getListOfProducts():
            prod += [ctx["scanon"][x.getSpecies()]] * _stoich(x)
        if any(c is None for c in reac + prod):
            continue
        law = libsbml.formulaToL3String(r.getKineticLaw().getMath())
        kr, kp = tuple(sorted(reac)), tuple(sorted(prod))
        for sign, coef, specs in monomials(law):
            got = tuple(sorted(ctx["scanon"][s] for s in specs))
            val = eval(coef.replace("^", "**"), {"__builtins__": {}}, ctx["env"]) if coef else 1.0
            if sign > 0 and got == kr:
                out[(kr, kp)] += val
            elif sign < 0 and got == kp:
                out[(kp, kr)] += val
            else:
                raise AssertionError(f"unmatched monomial in {r.getId()}")
    return out


def _bng_reactions(ctx):
    # The `.net` keeps every rate constant symbolic, so overriding the base
    # thermodynamic factors here re-derives all 27457 of them, exactly as the
    # same overrides do on the SBML side. Parameters are in dependency order.
    env = dict(SAFE)
    env.update(FUNC_MARKERS)
    for k, v in ctx["params"].items():
        if k in ctx["overrides"]:
            env[k] = float(ctx["overrides"][k])
            continue
        try:
            env[k] = float(v)
        except ValueError:
            env[k] = eval(v.replace("^", "**"), {"__builtins__": {}}, env)
    out = defaultdict(float)
    for reac, prod, rate in ctx["reactions"]:
        kr = tuple(sorted(ctx["bcanon"][i] for i in reac))
        kp = tuple(sorted(ctx["bcanon"][i] for i in prod))
        out[(kr, kp)] += eval(rate.replace("^", "**"), {"__builtins__": {}}, env)
    return out, env


def compare_reactions(ctx, tol=1e-9):
    sb = _sbml_reactions(ctx)
    bn, _ = _bng_reactions(ctx)
    shared = set(sb) & set(bn)
    worst, worst_key = 0.0, None
    for k in shared:
        a, b = sb[k], bn[k]
        d = abs(a - b) / max(abs(a), abs(b), 1e-300)
        if d > worst:
            worst, worst_key = d, k
    return dict(n_sbml=len(sb), n_bng=len(bn),
                only_sbml=len(set(sb) - set(bn)), only_bng=len(set(bn) - set(sb)),
                n_shared=len(shared), max_rel_err=worst, worst=worst_key,
                n_mismatch=sum(1 for k in shared
                               if abs(sb[k] - bn[k]) / max(abs(sb[k]), abs(bn[k]), 1e-300) > tol))


def rate_constant_pairs(ctx):
    """(SBML, BioNetGen) rate constants for every reaction both networks share.

    The two sequences the aggregate parity panel of the verification figure plots.
    """
    sb = _sbml_reactions(ctx)
    bn, _ = _bng_reactions(ctx)
    keys = sorted(set(sb) & set(bn))
    return ([sb[k] for k in keys], [bn[k] for k in keys])


def compare_activity(ctx):
    """Per-species weights of the aggregate RAF kinase activity."""
    import libsbml
    expr = [libsbml.formulaToL3String(r.getMath()) for r in ctx["model"].getListOfRules()
            if r.getVariable() == "RAF_kinase_activity"][0]
    sb = {}
    for term in expr.split(" + "):
        m = re.fullmatch(r"([\d.]+) \* (__s\d+)", term.strip())
        if m:
            sb[ctx["scanon"][m.group(2)]] = float(m.group(1))
        else:
            sb[ctx["scanon"][re.fullmatch(r"(__s\d+)", term.strip()).group(1)]] = 1.0
    bn = defaultdict(float)
    for g, w in ACTIVITY_WEIGHTS.items():
        for i, c in ctx["groups"][g].items():
            bn[ctx["bcanon"][i]] += w * c
    bn = {k: v for k, v in bn.items() if v}
    diff = [(k, sb.get(k), bn.get(k)) for k in set(sb) | set(bn)
            if abs((sb.get(k) or 0) - (bn.get(k) or 0)) > 1e-9]
    return dict(n_sbml=len(sb), n_bng=len(bn), n_mismatch=len(diff), examples=diff[:5])


def compare_functions(ctx):
    """The two Michaelis-Menten substrate pools and their Michaelis constants."""
    import libsbml
    den = [libsbml.formulaToL3String(r.getMath()) for r in ctx["model"].getListOfRules()
           if r.getVariable() == "Michaelis_Menten_kinetics_uMEK_phosphorylation"][0]
    tail = den[den.rindex(") * ("):]
    pools = re.findall(r"\((__s\d+(?: \+ __s\d+)*)\)", tail)
    out = {}
    for pool, obs in zip(pools, ("Obs_Sub_uMEK", "Obs_Sub_pMEK"), strict=False):
        s = sorted(ctx["scanon"][x] for x in pool.split(" + "))
        b = sorted(ctx["bcanon"][i] for i in ctx["groups"][obs])
        out[obs] = dict(n=len(s), identical=(s == b),
                        unit_weights=set(ctx["groups"][obs].values()) == {1.0})
    _, benv = _bng_reactions(ctx)
    for km in ("Km_RAF_MEK_1", "Km_RAF_MEK_2"):
        out[km] = dict(sbml=ctx["env"][km], bng=benv[km],
                       equal=abs(ctx["env"][km] - benv[km]) < 1e-9)
    return out


def report(net_path, sbml_path=DEFAULT_SBML):
    ctx = load(net_path, sbml_path)
    sp = compare_species(ctx)
    rx = compare_reactions(ctx)
    ac = compare_activity(ctx)
    fn = compare_functions(ctx)
    print(f"species    SBML {sp['n_sbml']:5d}  BNG {sp['n_bng']:5d}  "
          f"only-SBML {len(sp['only_sbml'])}  only-BNG {len(sp['only_bng'])}")
    for c in sp["only_sbml"]:
        print(f"           only in SBML: {c}")
    print(f"reactions  SBML {rx['n_sbml']:5d}  BNG {rx['n_bng']:5d}  "
          f"only-SBML {rx['only_sbml']}  only-BNG {rx['only_bng']}")
    print(f"rate const {rx['n_shared']} shared, {rx['n_mismatch']} mismatched, "
          f"max rel err {rx['max_rel_err']:.3g}")
    print(f"RAF activity  {ac['n_sbml']} weighted species, {ac['n_mismatch']} mismatched")
    for k, v in fn.items():
        print(f"  {k}: {v}")
    return dict(species=sp, reactions=rx, activity=ac, functions=fn)


if __name__ == "__main__":
    report(sys.argv[1], Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_SBML)
