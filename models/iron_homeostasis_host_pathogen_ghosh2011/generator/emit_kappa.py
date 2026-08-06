"""Write the Ghosh 2011 model as a Kappa file, for TRuML to translate.

The paper ships no .ka file -- the ESI is a PDF, and the rules exist only as the
Kappa text printed in Table S2, recovered into tableS2.json. That text is Kappa 3
(`site!1`, `site!_`, `site~state`, a bare site meaning free); TRuML parses Kappa 4
(`site[1]`, `site[_]`, `site{state}`, `site[.]`). The rewriting here is purely
syntactic. Agent signatures are derived from the rules themselves, the rate
constants are the rate column of tableS2.json, and the initial counts are the
Methods' own numbers in initial_counts.json.
"""
import json
import re
from pathlib import Path

G = Path(__file__).parent
W = G                      # tableS2.json sits beside this script
COUNTS = G / "initial_counts.json"
MANGLE = "st_"    # prefix for internal states whose name starts with a digit

# Rule 163 (sno 168) frees twelve identical Fe agents at once. TRuML computes a
# rule's symmetry factor by brute force over permutations of identically-named
# agents, so 12! candidates make it run forever -- the same wall BioNetGen hits
# reading the rule. It is handled separately in build_ghosh.py.
EXCLUDE = set()


BOND_K4 = {"_": "_", "+": "_", "?": "#", None: "."}   # `!+` and `!?` reach this from the BNGL blocks

# Extracting Table S2 from the PDF displaced some arrow glyphs onto the preceding
# row: rules 5, 7 and 16 come out carrying a second, trailing `->` and rules 6, 9
# and 17 come out with none. The strays are dropped (they sit past the end of the
# products, so they never affected the translation) and the three gaps are filled
# in below. The split points are unambiguous: in 6 and 9 the reactants end where an
# agent name follows a closing parenthesis with no comma, and 17 is a decay rule --
# 16 is the matching synthesis -- whose product side is empty.
ARROW_RECOVERED = {
    6:  ("Fe(site~3ex),Tf(Fe2~a)", "Fe(site~3ex!1),Tf(Fe2~a!1)"),
    9:  ("Tf(Fe1!_,Fe2~a?,tfR~a),TfR(Tf)", "Tf(Fe1!_,Fe2~a?,tfR~a!1),TfR(Tf!1)"),
    17: ("Dcytb(site)", ""),
}

# Rule 31's printed text stops mid-agent, at `fp(site!1` with no closing bracket.
# The completion is forced: bond 1 has to have two ends, and the site it names is
# fp's only one. Every other rule in Table S2 has balanced brackets and an even
# number of ends for every bond, so this is the only rule affected.
TRUNCATED = {
    31: "Fe(site~2cy),fp(site)->Fe(site~fp!1),fp(site!1)",
}


def site_k4(site):
    """One Kappa 3 (or BNGL) site expression -> Kappa 4."""
    site = re.sub(r"~(\w+)\?$", r"~\1!?", site.strip())   # Kappa 3 `Fe2~a?`: state set, bond unspecified
    m = re.fullmatch(r"(\w+)(?:~(\w+))?(?:!(\d+|_|\+|\?))?", site.strip())
    if not m:
        raise ValueError(f"unparsed site: {site!r}")
    name, state, bond = m.groups()
    out = name
    if state:
        out += "{%s}" % state
    return out + "[%s]" % (bond if bond and bond.isdigit() else BOND_K4[bond])


def pattern_k4(text):
    """A Kappa 3 pattern (a comma-separated agent list) -> Kappa 4."""
    agents, consumed = [], []
    for am in re.finditer(r"([A-Za-z_]\w*)\(([^)]*)\)", text):
        sites = [s for s in am.group(2).split(",") if s.strip()]
        agents.append(f"{am.group(1)}({' '.join(site_k4(s) for s in sites)})")
        consumed.append(am.span())
    # An agent the pattern matcher cannot see would be dropped in silence, which is
    # how the truncation in rule 31 first went unnoticed. Anything left over that is
    # not a separator is an error.
    rest, last = "", 0
    for a, b in consumed:
        rest += text[last:a]
        last = b
    rest += text[last:]
    if rest.strip(" ,\t"):
        raise ValueError(f"unconsumed text in pattern: {rest.strip()!r} (from {text!r})")
    return ", ".join(agents)


def signatures(rules):
    """Agent signatures read off the rules themselves.

    Taking these from the prototype BNGL would import that file's own modelling
    choices -- it represents ideR's three iron sites as one counter, which the
    Kappa does not. The rules are the primary source, so the signature is the
    union of every site and internal state each agent is ever written with.
    """
    sites, states = {}, {}
    for r in rules.values():
        for am in re.finditer(r"([A-Za-z_]\w*)\(([^)]*)\)", r["kappa"].replace("\n", "")):
            agent = am.group(1)
            sites.setdefault(agent, dict())     # agents such as ALAS() carry no sites
            for s in am.group(2).split(","):
                m = re.fullmatch(r"(\w+)(?:~(\w+))?(?:!(?:\d+|_|\+|\?))?", s.strip())
                if not m:
                    continue
                sites.setdefault(agent, dict())[m.group(1)] = None
                if m.group(2):
                    states.setdefault((agent, m.group(1)), dict())[m.group(2)] = None
    return [f"%agent: {a}(" + " ".join(
        s + ("{%s}" % " ".join(states[(a, s)]) if (a, s) in states else "")
        for s in sites[a]) + ")" for a in sorted(sites)]


# The eight observables plotted in Fig. 6, in the paper's own terms.
OBSERVABLES = [
    ("Obs_Fe_mtb", "Fe(site{mtb}[.])"),          # free iron in the pathogen
    ("Obs_Fe_3ex", "Fe(site{3ex}[.])"),          # free iron, extracellular
    ("Obs_Fer_Fe", "Fe(site{fer}[_])"),          # iron stored in ferritin
    ("Obs_Fe_myB", "m(site1[_])"),               # iron on mycobactin
    ("Obs_Fe_cmB", "cmex(site1[_])"),            # iron on carboxymycobactin
    ("Obs_act_ideR", "ideR(Fe1[_] Fe2[_] Fe3[_] ideR[_])"),   # loaded ideR dimer
    ("Obs_Tf_TfR", "Tf(tfR{a}[_])"),
    ("Obs_RBC", "RBC(site{a}[.])"),
]


def main():
    rules = {r["sno"]: r for r in json.load(open(W / "tableS2.json"))}
    out = ["# Ghosh et al. 2011, Mol. BioSyst. 7:2750, Table S2 (ESI).",
           "# Kappa 4 transcription of the printed Kappa 3 rules; see emit_kappa.py.",
           ""]
    out += signatures(rules)

    out.append("")
    # The rate constants are the rate column of Table S2 itself; there is no second
    # copy of them to drift out of step.
    for sno in sorted(rules):
        out.append(f"%var: 'k{sno}' {rules[sno]['rate']}")

    out.append("")
    for pat, count in json.loads(COUNTS.read_text())["seed_species"]:
        out.append(f"%init: {count} {pattern_k4(pat)}")

    out.append("")
    for name, pat in OBSERVABLES:
        out.append(f"%obs: '{name}' |{pat}|")

    out.append("")
    for sno in sorted(rules):
        if sno in EXCLUDE:
            continue
        if sno in ARROW_RECOVERED:
            lhs, rhs = ARROW_RECOVERED[sno]
        else:
            k = TRUNCATED.get(sno) or re.sub(r"@.*$", "", rules[sno]["kappa"].replace("\n", ""))
            if k.count("->") == 2:            # a displaced arrow from the next row;
                k = re.sub(r"->\s*$", "", k)      # one trailing arrow is a decay rule
            lhs, rhs = k.split("->", 1)
        lhs, rhs = pattern_k4(lhs), pattern_k4(rhs)
        out.append(f"'R{sno}' {lhs or '.'} -> {rhs or '.'} @ 'k{sno}'"
                   f"  // {rules[sno]['desc']}")

    text = "\n".join(out) + "\n"

    # TRuML's grammar requires an internal state to start with a letter, and this
    # model names iron pools 3ex, 2cy, 2mt and so on. Rename them for the round
    # trip; unmangle_states() puts them back in the BNGL TRuML returns.
    mangled = set()

    def fix_states(m):
        names = m.group(1).split()
        mangled.update(n for n in names if n[0].isdigit())
        return "{%s}" % " ".join(MANGLE + n if n[0].isdigit() else n for n in names)

    text = re.sub(r"\{([^}]*)\}", fix_states, text)     # only inside state braces
    mangled = sorted(mangled)
    (G / "state_mangling.json").write_text(json.dumps(mangled, indent=2))

    (G.parent / "ghosh2011.ka").write_text(text)
    print(f"wrote ghosh2011.ka: {len(rules) - len(EXCLUDE)} rules"
          f" (rule sno {sorted(EXCLUDE)} held back), "
          f"arrow recovered for sno {sorted(ARROW_RECOVERED)}, "
          f"renamed {len(mangled)} digit-initial states {mangled}")


main()
