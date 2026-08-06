"""Build the curated BNGL from the Kappa rules printed in Ghosh et al. (2011), Table S2.

The chain is

    tableS2.json  ->  emit_kappa.py  ->  ghosh2011.ka  ->  TRuML  ->  post-processing here

TRuML (github.com/lanl/TRuML) does the translation. It is Python 2 and its Kappa
parser is Kappa 4 only, so `truml_py3.patch` in this folder ports it and replaces
one algorithm; apply it to a checkout of upstream head. Everything this script does
after TRuML is forced by a BioNetGen or Kappa/BNGL language limit, never by
modelling taste, and each change is commented where it happens.

    python build_ghosh2011.py --truml /path/to/patched/TRuML

`--check` rebuilds the committed files in place and exits non-zero if they drift.
"""
import argparse
import difflib
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FOLDER = HERE.parent
STEM = "iron_homeostasis_host_pathogen_ghosh2011"

KFAST = "1.0e4"

# Rules 137 and 143 move iron from host transferrin and host ferritin onto the
# pathogen's two siderophores, and 138 and 144 complete each transfer. These are the
# four rules the paper calls the host-pathogen interactions, and the non-infected
# model of the Methods is the same rule set with their rates set to zero.
INTERACTION_RATES = ("k142", "k143", "k148", "k149")


def run_truml(truml_dir, python=sys.executable):
    subprocess.run([python, str(HERE / "emit_kappa.py")], check=True)
    work = HERE / "_truml_work"
    work.mkdir(exist_ok=True)
    (work / "ghosh2011.ka").write_text((FOLDER / "ghosh2011.ka").read_text())
    subprocess.run([python, "-m", "truml", "-k", "ghosh2011.ka"], cwd=work, check=True,
                   env={"PYTHONPATH": str(truml_dir), "PATH": "/usr/bin:/bin"},
                   capture_output=True)
    return (work / "ghosh2011.bngl").read_text()


def unmangle(text):
    """Restore state names emit_kappa.py renamed to get them past TRuML's grammar."""
    for s in json.loads((HERE / "state_mangling.json").read_text()):
        text = text.replace("~st_" + s, "~" + s)
    assert "st_" not in text, "a renamed state survived"
    return text


def split_rule_163(text):
    """Rule 163 as published frees twelve identical Fe in one step.

    Both BioNetGen and TRuML size a rule by enumerating permutations of its
    identically-named molecules, so twelve of them is 12! candidates and neither
    finishes reading the rule. 163a performs the published step but leaves each ideR
    holding its iron and flags it; 163b then frees three irons at a time. The flag
    confines the fast step to ideR that has just come off this operator, so no other
    rule can see the intermediate, and the error in the trajectory is first order in
    1/kfast. 163a keeps the published reactant pattern, so TRuML's factor for the
    rule still applies to it unchanged.
    """
    rule = next(l for l in text.split("\n") if "k168" in l and "->" in l)
    lhs, rhs = rule.split("->", 1)
    assert rhs.strip().endswith("k168 * 2592.0"), rhs.strip()[-40:]

    loaded = re.findall(r"ideR\(Fe1!\d+,Fe2!\d+,Fe3!\d+,ideR!\d+,site(?:!\d+)?\)", lhs)
    assert len(loaded) == 4, loaded
    # The published rule breaks the operator bonds and both dimer bonds, so the
    # products are four separate monomers; only the iron bonds carry over.
    products = " + ".join(
        re.sub(r",ideR!\d+,site(!\d+)?\)", ",ideR,site,rel~1)", m) + "." +
        ".".join(f"Fe(site~mtb!{n})" for n in re.findall(r"!(\d+)", m)[:3])
        for m in loaded)
    assert "ideR!" not in products, products
    a = (lhs.replace(",site)", ",site,rel~0)")
         + "-> g_bfr(site1~inc,site2) + " + products + "  k168 * 2592.0")
    a = re.sub(r"(ideR\(Fe1!\d+,Fe2!\d+,Fe3!\d+,ideR!\d+,site!\d+)\)", r"\1,rel~0)", a)
    b = ("\tideR(Fe1!1,Fe2!2,Fe3!3,ideR,site,rel~1).Fe(site~mtb!1).Fe(site~mtb!2).Fe(site~mtb!3)"
         " -> ideR(Fe1,Fe2,Fe3,ideR,site,rel~0) + Fe(site~mtb) + Fe(site~mtb) + Fe(site~mtb)"
         "  kfast")
    text = text.replace(rule, a.rstrip() + "\n" + b)

    head, sep, tail = text.partition("end molecule types")
    text = head.replace("\tideR(Fe1,Fe2,Fe3,ideR,site)\n",
                        "\tideR(Fe1,Fe2,Fe3,ideR,site,rel~0~1)\n") + sep + tail
    # An ideR a rule creates rather than maps must name every component.
    return re.sub(r"(\+ideR\(Fe1,Fe2,Fe3,ideR,site\))(?= k\d)",
                  "+ideR(Fe1,Fe2,Fe3,ideR,site,rel~0)", text)


def name_wildcard_partner(text):
    """Kappa dissolves a bond to an unspecified partner; BNGL cannot.

    Rule 10 dissociates transferrin from its receptor, written in Kappa as
    `tfR{a}[_] -> tfR{a}[.]`, which BioNetGen rejects as illegal wildcard bond
    breaking. TfR is the only agent in the model that ever bonds Tf's tfR site --
    rules 9 and 15 form and hold that bond and nothing else touches it -- so naming
    the partner says exactly what the Kappa says.
    """
    old = "\tTf(Fe1!+,tfR~a!+) -> Tf(Fe1!+,tfR~a) k10"
    new = "\tTf(Fe1!+,tfR~a!1).TfR(Tf!1) -> Tf(Fe1!+,tfR~a)+TfR(Tf) k10"
    assert old in text, "rule 10 is not in the expected form"
    return text.replace(old, new)


def fix_rule_27(text):
    """Ferroportin takes ferritin's iron and the ferritin is destroyed.

    TRuML writes the two products as separate complexes although a bond joins them,
    which BioNetGen rejects as a dangling edge, and omits DeleteMolecules, without
    which BioNetGen destroys the whole complex rather than the one molecule. The
    Kappa is unambiguous:
        Fe(site{fer}[1]), fp(site[.]), Fer(Fe[1]) -> Fe(site{fp}[1]), fp(site[1])
    This is the only rule in the model that deletes a bonded molecule and keeps its
    partner; rules 123 and 124 consume their whole complex and need no flag.
    """
    old = "\tFe(site~fer!1).Fer(Fe!1)+fp(site) -> Fe(site~fp!1)+fp(site!1) k29"
    new = "\tFe(site~fer!1).Fer(Fe!1)+fp(site) -> Fe(site~fp!1).fp(site!1) k29 DeleteMolecules"
    assert old in text, "rule 27 is not in the expected form"
    return text.replace(old, new)


HEADER = """begin model

#@title: Host-pathogen iron homeostasis in tuberculosis (Ghosh et al., 2011){variant_title}

#@description: |
#  A rule-based, network-free stochastic model of the contest for iron between a human host
#  and Mycobacterium tuberculosis, translated from the Kappa model of Ghosh et al. (2011).
#  The host module carries transferrin uptake through the Tf-TfR endocytic cycle, ferritin
#  storage, ferroportin export, the Fe-S cluster and heme biosynthesis pathways, iron
#  regulatory protein control, and red cell turnover. The pathogen module carries siderophore
#  synthesis on the mbt assembly line, iron import through irtAB and esx3, bacterioferritin
#  storage, its own Fe-S cluster pathway, and ideR- and SufR-dependent transcriptional
#  regulation. Four rules connect the two: the siderophores carboxymycobactin and mycobactin
#  strip iron from host transferrin and host ferritin respectively.
#
#  How it works: 194 published rules over 90 molecule types. The species space does not close
#  -- iron and the siderophore precursors are made without bound -- so the model is run
#  network-free, as the authors ran it in Kappa, rather than through generate_network.
#
#  What it does: from one copy of each gene, 5000 extracellular and 50 pathogen iron atoms
#  and 20000 red cells, the system is incubated for 30000 time units, the protocol the
#  Methods use to reach a steady state.{variant_description}

#@keyword: |
#  host-pathogen interaction, iron homeostasis, tuberculosis, Mycobacterium tuberculosis,
#  siderophore, mycobactin, transferrin, ferritin, bacterioferritin, ideR, network-free,
#  stochastic, rule-based, translated from Kappa, Homo sapiens

#@reference: |
#  Ghosh S, Prasad KVSK, Vishveshwara S, Chandra N (2011). Rule-based modelling of iron
#  homeostasis in tuberculosis. Molecular BioSystems 7:2750-2768. doi:10.1039/c1mb05093a

#@note: |
#  SIMULATION INTENT. Network-free and stochastic: every amount is a molecule count, and the
#  rate constants are the per-molecule constants of the published Kappa file. The network is
#  not generable -- rule 100 makes five siderophore precursors from nothing at a constant rate
#  and no rule degrades any of them -- so `generate_network` cannot terminate and the model is
#  run with NFsim. Use bngsim's NFsim core rather than the NFsim binary BNG2.pl ships: seven
#  rules here move a bond from one partner to another, and NFsim v1.14.3 aborts partway
#  through such a run with "sites that are already occupied", at a different time on each
#  seed. bngsim runs all of them. Nothing in this file is shaped around that defect.

#@note: |
#  TRANSLATED FROM KAPPA BY TRuML, NOT BY HAND. The paper publishes no .ka file; the rules
#  exist only as the Kappa text printed in Table S2 of the ESI, recovered into
#  generator/tableS2.json. generator/emit_kappa.py rewrites that Kappa 3 text as the Kappa 4
#  TRuML parses, and TRuML (github.com/lanl/TRuML) performs the translation. This matters
#  quantitatively: Kappa and BNGL count the automorphisms of a reactant pattern differently,
#  and TRuML applies a symmetry factor to 31 of the 194 rules, ranging from 2 to 2592. A hand
#  translation of the same rules applied none of them, i.e. it had those rates wrong by up to
#  three orders of magnitude. generator/truml_py3.patch ports TRuML to Python 3 and replaces
#  its automorphism counter, which enumerates permutations of identically-named molecules and
#  therefore cannot finish on a rule containing twelve identical Fe; the replacement counts
#  the same group with VF2 and reproduces TRuML's own test assertions.

#@note: |
#  RULE 163 IS SPLIT IN TWO, AND IT IS THE ONLY RULE THAT IS. As published it releases four
#  iron-loaded ideR from the bfr operator and frees twelve identical Fe in a single step.
#  Sizing that rule means enumerating 12! permutations, which neither BioNetGen nor TRuML
#  finishes -- BioNetGen never gets past reading the rule. 163a performs the published step
#  with the iron still bound and flags the released ideR; 163b frees three irons at a time at
#  kfast, far above the fastest published rate of 30 /s. The flag is invisible to every other
#  rule, iron is conserved throughout, and the error is first order in 1/kfast.

#@note: |
#  FOUR REPAIRS TO THE PRINTED TABLE. Extracting Table S2 from the PDF displaced some arrow
#  glyphs onto the preceding row, so rules 5, 7 and 16 carry a spurious trailing arrow and
#  rules 6, 9 and 17 have none; the strays sit past the end of the products and the three
#  gaps are filled from the printed text, which is unambiguous. Rule 31's text stops
#  mid-agent at `fp(site!1` with no closing bracket, and its completion is forced because a
#  bond needs two ends. Auditing all 194 rules, rule 31 is the only truncation and every bond
#  has an even number of ends, so the extraction is otherwise sound. emit_kappa.py now fails
#  loudly rather than silently dropping an agent it cannot parse, which is how the rule 31
#  damage went unnoticed at first.

#@note: |
#  WHAT THIS MODEL DOES NOT REPRODUCE, AND WHY THAT IS A PROPERTY OF THE PUBLISHED RULES. The
#  paper reports a steady state with bacterioferritin at 82 molecules, irtAB at 22 and stored
#  Fe-bfr at 41. This model produces none of them: bfr, mycobactin, carboxymycobactin and the
#  assembled siderophore are identically zero for the whole run. The cause is arithmetic, not
#  simulation. Rule 100 synthesises lysine, salicylate, serine and two acyl precursors at
#  k105 = 1.0 per time unit, zeroth order, and no rule degrades any of them, so each grows as
#  exactly k105*t -- 30000 molecules by the end of the incubation. Their only consumers, mbtE
#  and mbtF, come from single-copy genes with steady states of about 10 and 12 molecules, so
#  the assembly line is oversubscribed roughly a thousandfold. Loaded mbtE cannot degrade
#  either, because its decay rule requires every site free, so partially-loaded mbtE
#  accumulates and the assembly never completes. With no siderophore, no iron reaches the
#  pathogen, ideR is never activated, the bfr operator is never switched on. Separately, red
#  cells are seeded at 20000 against an implied steady state of k97/k98 = 4.08 over an
#  incubation lasting 39 decay time constants, so that pool can only collapse. All of this is
#  readable from the rate constants; verify_ghosh2011.ipynb checks each claim in closed form.

#@figure: Fig. 5 and Fig. 6 in Ghosh et al. (2011)
"""


def build(truml_dir):
    text = run_truml(truml_dir)
    text = unmangle(text)
    text = split_rule_163(text)
    text = name_wildcard_partner(text)
    text = fix_rule_27(text)
    text = text.replace("\tk168 0.05",
                        f"\tk168 0.05\n\tkfast {KFAST}  # /time unit, completes split rule 163")
    body = text.split("begin model", 1)[1]

    out = {}
    for variant in ("", "_noninfected"):
        head = HEADER.format(
            variant_title=" -- non-infected control" if variant else "",
            variant_description=(
                "\n#\n#  This is the paper's non-infected control: the same 194 rules with the four\n"
                "#  host-pathogen interaction rates set to zero, so the siderophores cannot reach host\n"
                "#  iron. The Methods build it exactly this way."
                if variant else ""))
        b = body
        if variant:
            for k in INTERACTION_RATES:
                b = re.sub(rf"^\t{k} \S+", f"\t{k} 0.0", b, flags=re.M)
        b += ('\nbegin actions\n\n'
              '# The Methods incubate for 30000 time units to reach a steady state.\n'
              'simulate({method=>"nf",t_end=>30000,n_steps=>300,'
              'complex_bookkeeping=>1,utl=>1000})\n\nend actions\n')
        out[FOLDER / f"{STEM}{variant}.bngl"] = head + b
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truml", required=True, help="checkout of TRuML with truml_py3.patch applied")
    ap.add_argument("--check", action="store_true", help="fail if the committed files would change")
    args = ap.parse_args()

    drift = False
    for path, text in build(args.truml).items():
        old = path.read_text() if path.exists() else ""
        if args.check and old != text:
            drift = True
            print(f"DRIFT {path.name}")
            sys.stdout.writelines(difflib.unified_diff(
                old.splitlines(True), text.splitlines(True), "committed", "rebuilt", n=1))
        path.write_text(text)
        print(f"wrote {path.name}: {len(re.findall(r'->', text.split('begin reaction rules')[1]))} rules")
    sys.exit(1 if drift else 0)


if __name__ == "__main__":
    main()
