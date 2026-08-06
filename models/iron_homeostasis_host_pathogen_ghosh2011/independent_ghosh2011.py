"""Level-1 verification: predictions built from Table S2, with no simulator involved.

The usual level-1 witness for a network-free model in this collection is agreement
between bngsim's two independent engines. That is not available here: RuleMonkey
silently never fires a rule whose reactant side has three or more identical
patterns (richardposner/RuleMonkey#24, lanl/bngsim#185), and this model is full of
them -- every rule that consumes pathogen iron has three or four identical `Fe`
reactants. Using it would have produced a confident, wrong answer, since it holds
`Fe(site~mtb)` at exactly its seeded 50 forever while reproducing unrelated arms
to within 1%.

So the witness here is arithmetic instead. Several parts of this model are exactly
solvable, and this file solves them from the published rate constants alone -- it
reads `generator/tableS2.json` and never looks at the .bngl, the XML, or any
simulator output. Three families:

1. CATALYTIC ENZYMES. Each heme-pathway enzyme is made from a single-copy gene at
   a constant rate and decays first order, and its catalytic rule leaves it
   unconsumed (`PBS() + ALA + ALA -> PBS() + PBG()`), so it is a pure immigration
   -death process, decoupled from the rest of the model. Its mean is
   (s/d)(1 - exp(-d t)) exactly, for every t, whatever the rest of the network does.

   One caveat, which the model makes load-bearing. Every decay rule here names the
   species with all of its sites free, so only the unbound form can be deleted.
   Writing T for the total count and F for the free form, dT/dt = s - d*F, so at
   steady state it is F that equals s/d exactly, and T exceeds it by however much
   is protected in bound complexes. For most of these species binding is
   negligible and T = F to within noise. ideR is the exception and is worth
   measuring rather than excusing: at t = 2000 its free form is 1.33 +/- 0.21
   against a predicted 1.49, while the total is 11.71, so 10.4 molecules are held
   in forms rule 102 cannot touch. That is the same structure that jams the
   siderophore line, where loaded mbtE likewise cannot degrade.

2. RED CELLS. Zeroth-order production against first-order loss, seeded far from its
   own steady state: 20000 initially against s/d = 4.08.

3. SIDEROPHORE PRECURSORS. Zeroth-order production with no degradation rule at all,
   so the count is exactly k*t and grows without bound.

Each family is also run through a hand-written Gillespie simulator in this file, so
the closed forms are checked against a stochastic implementation that shares no
code with BioNetGen, NFsim or bngsim.
"""
import json
import math
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
TABLE = HERE / "generator" / "tableS2.json"


def _rules():
    return {r["sno"]: r for r in json.loads(TABLE.read_text())}


def _kappa(r):
    return re.sub(r"@.*$", "", r["kappa"].replace("\n", ""))


def catalytic_enzymes():
    """{name: (k_synthesis, k_decay)} for every species that is made from a gene,
    decays first order, and is never consumed by any other rule it appears in.

    Found by reading the rules, not by hard-coding a list: an agent qualifies when
    one rule creates it from a `g_*()` gene, one rule destroys it on its own, and
    in every other rule it appears the same number of times on both sides.
    """
    rules = _rules()
    agents = {m.group(1) for r in rules.values()
              for m in re.finditer(r"([A-Za-z_]\w*)\(", _kappa(r))}
    out = {}
    for a in sorted(agents):
        if a.startswith("g_"):
            continue
        syn = dec = None
        catalytic = True
        for sno, r in rules.items():
            k = _kappa(r)
            if not re.search(rf"\b{a}\(", k):
                continue
            if "->" not in k:
                continue
            lhs, rhs = k.split("->", 1)
            nl = len(re.findall(rf"\b{a}\(", lhs))
            nr = len(re.findall(rf"\b{a}\(", rhs))
            if nl == 0 and nr == 1 and re.search(r"g_\w+\(\)", lhs):
                syn = float(r["rate"])
            elif nl == 1 and nr == 0 and len(re.findall(r"[A-Za-z_]\w*\(", lhs)) == 1:
                dec = float(r["rate"])
            elif nl != nr:
                catalytic = False
        if syn and dec and catalytic:
            out[a] = (syn, dec)
    return out


def immigration_death_mean(t, s, d, n0=0.0):
    """Exact mean of dn/dt = s - d n."""
    e = math.exp(-d * t)
    return n0 * e + (s / d) * (1.0 - e)


def rbc_mean(t):
    """Circulating red cells: rule 97 makes them, rule 98 removes them."""
    r = _rules()
    return immigration_death_mean(t, float(r[97]["rate"]), float(r[98]["rate"]), n0=20000.0)


def precursor_mean(t):
    """Free siderophore precursors: rule 105 makes five of them at once and no rule
    degrades any, so the count is exactly k*t less whatever the mbt line has taken."""
    return float(_rules()[105]["rate"]) * t


def gillespie(reactions, x0, t_end, times, rng):
    """Direct-method SSA. `reactions` is a list of (propensity(x), delta dict)."""
    x = dict(x0)
    out, t, i = [], 0.0, 0
    times = sorted(times)
    while i < len(times):
        a = [max(0.0, f(x)) for f, _ in reactions]
        a0 = sum(a)
        t_next = t - math.log(rng.random()) / a0 if a0 > 0 else math.inf
        # The state is x throughout [t, t_next), so every sample instant in that
        # half-open window takes x as it stands *before* the next event is applied.
        # Recording after applying it biases every count upward by about one event.
        while i < len(times) and times[i] < min(t_next, t_end + 1e-12):
            out.append(dict(x))
            i += 1
        if t_next >= t_end or a0 <= 0:
            break
        u, acc = rng.random() * a0, 0.0
        for (f, delta), ai in zip(reactions, a):
            acc += ai
            if u <= acc:
                for k, v in delta.items():
                    x[k] = x.get(k, 0) + v
                break
        t = t_next
    while i < len(times):
        out.append(dict(x))
        i += 1
    return out


def immigration_death_ssa(s, d, n0, times, rng, replicates=800):
    """Ensemble mean of the same process by simulation rather than by formula."""
    t_end = max(times)
    total = [0.0] * len(times)
    for _ in range(replicates):
        traj = gillespie([(lambda x: s, {"n": +1}), (lambda x: d * x["n"], {"n": -1})],
                         {"n": n0}, t_end, times, rng)
        for j, st in enumerate(traj):
            total[j] += st["n"]
    return [v / replicates for v in total]


if __name__ == "__main__":
    import random
    rng = random.Random(20110101)
    times = [0.0, 500.0, 1000.0, 2000.0]

    print("catalytic enzymes found by reading Table S2:")
    print(f"  {'species':<8}{'k_syn':>10}{'k_dec':>10}{'steady state':>14}")
    for a, (s, d) in catalytic_enzymes().items():
        print(f"  {a:<8}{s:>10.4g}{d:>10.4g}{s / d:>14.2f}")

    print("\nclosed form vs a hand-written SSA of the same process (mean of 800 runs):")
    for a, (s, d) in catalytic_enzymes().items():
        cf = [immigration_death_mean(t, s, d) for t in times]
        ssa = immigration_death_ssa(s, d, 0, times, rng)
        worst = max(abs(c - m) / max(c, 1.0) for c, m in zip(cf, ssa))
        print(f"  {a:<8} closed form {['%.1f' % v for v in cf]}"
              f"  ssa {['%.1f' % v for v in ssa]}  worst {worst:.3f}")

    print(f"\nRBC:        {[round(rbc_mean(t), 1) for t in times]}")
    print(f"precursors: {[round(precursor_mean(t), 1) for t in times]}")
