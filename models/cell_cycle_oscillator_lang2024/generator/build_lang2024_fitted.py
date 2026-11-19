#!/usr/bin/env python
"""Build `cell_cycle_oscillator_lang2024_fitted.bngl` from the authors' own PEtab tables.

The fitted file is the `_fused` model at the parameter set Lang et al. (2024) obtained by
fitting it: the reaction network is not touched, only numbers change. Every number is read
from the authors' released tables, vendored byte-for-byte under `reference/`:

  * all 205 parameter values (114 kinetic, 16 observation-model, 2 noise, 73 initial
    concentrations) from `results_20220421/parameters_v3.2.0_optimized.tsv` -- the table the
    authors simulated to make `results_20220421/optimisation.tsv`, the v3.2.0 fit figure;
  * the initial state from the PEtab `wt` condition (`experimentalCondition_v3.2.0.tsv`),
    which names, for every one of the 73 species, the `ic_wt_cell__*` parameter that sets it.

What is copied unchanged from `cell_cycle_oscillator_lang2024_fused.bngl` (pinned below by
sha256): the compartments, molecule types, seed species, observables, functions and reaction
rules blocks, byte for byte. What is written here: the header, the parameters block (house
style, with a unit comment per parameter derived from the reaction order of the rules that use
it) and the actions block.

Nothing here runs BioNetGen and nothing is typed by hand: re-running this script must reproduce
the committed `.bngl` exactly (`--check`), and `verify_lang2024_fitted.ipynb` asserts it.

Usage (from anywhere):
    python generator/build_lang2024_fitted.py           # (re)write the .bngl
    python generator/build_lang2024_fitted.py --check   # exit 1 if the committed file differs
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

FOLDER = Path(__file__).resolve().parents[1]
REF = FOLDER / "reference"
FUSED = FOLDER / "cell_cycle_oscillator_lang2024_fused.bngl"
TARGET = FOLDER / "cell_cycle_oscillator_lang2024_fitted.bngl"

# The model body is taken from this exact file; if it ever changes, the fitted file must be
# re-derived and re-verified deliberately rather than drift along with it.
FUSED_SHA256 = "4f4106485535a3ac8ed9e012606ba14e81a40222f46138c96209b131fa33b383"

# Upstream: the authors' repository at the commit that relicensed it under MIT (2024-01-14).
# No file under versions/v3.2.0/ differs between this commit, the earlier d562d3d4 the
# collection first cited, and the repository HEAD 874d13a4 (2024-07-16).
UPSTREAM_REPO = "https://github.com/paulflang/cell_cycle_petab"
UPSTREAM_COMMIT = "362b6667a2096b6c1cab43ddc3a92d4174c5a087"

# vendored file (in reference/) -> (path in the upstream repository, sha256)
PINS = {
    "lang2024_upstream_LICENCE": (
        "LICENCE",
        "d2d58b863fb4b29b94ac30dac13e825911592b031370b494d9d1e2f6a593d4aa",
    ),
    "lang2024_upstream_cell_cycle_v3.2.0.bngl": (
        "versions/v3.2.0/cell_cycle_v3.2.0.bngl",
        "0f23b3e866aff406a0b0ad7c7a404085b01322d2f4a350a6986a57d2589e04ed",
    ),
    "lang2024_upstream_cell_cycle_v3.2.0_sbml.xml": (
        "versions/v3.2.0/cell_cycle_v3.2.0_sbml.xml",
        "5ac5ea95131f7b57ddbed04702d54cce331ec13edc93741a1a0eef50045812db",
    ),
    "lang2024_upstream_cell_cycle_v3.2.0_petab.xml": (
        "versions/v3.2.0/cell_cycle_v3.2.0_petab.xml",
        "f8a15a57f71d334fe187d3905cef7bc2eb5f98a396afc9d2c570318a989da3c7",
    ),
    "lang2024_upstream_parameters_v3.2.0.tsv": (
        "versions/v3.2.0/parameters_v3.2.0.tsv",
        "a0d2cb24f0a16138e8608e500a4fe5518855943c761607ac31c15368c419bf63",
    ),
    "lang2024_upstream_experimentalCondition_v3.2.0.tsv": (
        "versions/v3.2.0/experimentalCondition_v3.2.0.tsv",
        "28c13f5038df7c41deca257b7cc54f295ff8470355821f5dff8f7d5da360942d",
    ),
    "lang2024_upstream_observables_v3.2.0.tsv": (
        "versions/v3.2.0/observables_v3.2.0.tsv",
        "ecd9a4cbd24424eec29c591d7658b474e2a200c73be1ccf033e5c9f3cd4722c1",
    ),
    "lang2024_upstream_results_20220421_parameters_v3.2.0_optimized.tsv": (
        "versions/v3.2.0/results_20220421/parameters_v3.2.0_optimized.tsv",
        "dc8c0d9547880395a8570e85f930e49c8eb24b354dc73134af5be8b59a9d6987",
    ),
    "lang2024_upstream_results_20220421_optimisation.tsv": (
        "versions/v3.2.0/results_20220421/optimisation.tsv",
        "371f0abfde46fa770577ee2a6c85d14925e8799df0be341708dc59f6be2802b6",
    ),
    "lang2024_upstream_results_20220421_v3.2.0_plot.yaml": (
        "versions/v3.2.0/results_20220421/v3.2.0_plot.yaml",
        "e5bf731211fa8cf576e6024713e72a5e9df5c07d1a431cc498fb58c4083c2c35",
    ),
}

OPTIMIZED = "lang2024_upstream_results_20220421_parameters_v3.2.0_optimized.tsv"
CONDITIONS = "lang2024_upstream_experimentalCondition_v3.2.0.tsv"

# The simulation protocol: the horizon of the authors' simulation output (its last sample is
# t = 136471.6126 s), on a uniform grid finer than their smallest sample spacing (164.6 s).
T_END = 136471.6126
N_STEPS = 1200
ATOL, RTOL = 1e-12, 1e-10

# Blocks carried over from _fused unchanged, in canonical order (bngl/skill.md section 2).
BODY_BLOCKS = (
    "compartments",
    "molecule types",
    "seed species",
    "observables",
    "functions",
    "reaction rules",
)

UNITS = {0: "AU/s", 1: "/s", 2: "/AU/s"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_pins() -> None:
    """Every vendored upstream file and the _fused template are the bytes that were pinned."""
    bad = [name for name, (_, digest) in PINS.items() if sha256(REF / name) != digest]
    if sha256(FUSED) != FUSED_SHA256:
        bad.append(FUSED.name)
    if bad:
        raise SystemExit(f"sha256 mismatch, refusing to build: {', '.join(bad)}")


def read_tsv(name: str) -> list[dict[str, str]]:
    with open(REF / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def block_span(text: str, name: str) -> tuple[int, int]:
    """Character span of `begin <name>` .. `end <name>` (both lines included)."""
    m = re.search(rf"^[ \t]*begin {name}[ \t]*\n.*?^[ \t]*end {name}[ \t]*$", text, re.S | re.M)
    if m is None:
        raise ValueError(f"no `{name}` block")
    return m.start(), m.end()


def block_lines(text: str, name: str) -> list[str]:
    start, end = block_span(text, name)
    inner = text[start:end].splitlines()[1:-1]
    return [ln.split("#", 1)[0].strip() for ln in inner if ln.split("#", 1)[0].strip()]


def rate_orders(fused: str) -> dict[str, set[int]]:
    """Reaction order of every rate-constant parameter, read off the reaction rules.

    A rule is `label: A + B -> C rate[, reverse_rate] [modifiers]`. Its forward rate constant
    has order = number of reactant patterns (`0` counts as none), the reverse rate constant
    order = number of product patterns.
    """
    orders: dict[str, set[int]] = {}
    for line in block_lines(fused, "reaction rules"):
        _, _, rule = line.partition(":")
        arrow = " <-> " if " <-> " in rule else " -> "
        lhs, _, rest = rule.partition(arrow)
        tokens = rest.split()
        rhs = [tokens[0]]
        i = 1
        while i < len(tokens) and tokens[i] == "+":
            rhs.append(tokens[i + 1])
            i += 2
        rates = [r.strip() for r in " ".join(tokens[i:]).split(",")]
        rates = [rates[0].split()[0]] + ([rates[1].split()[0]] if arrow == " <-> " else [])
        left = [p.strip() for p in lhs.split(" + ")]
        n_fwd = sum(p != "0" for p in left)
        n_rev = sum(p != "0" for p in rhs)
        for rate, order in zip(rates, (n_fwd, n_rev), strict=False):
            orders.setdefault(rate, set()).add(order)
    return orders


def fmt(value: str) -> str:
    """The table's number as the shortest decimal that round-trips to the same double."""
    return repr(float(value))


HEADER = """\
begin model

#@title: RPE-1 cell cycle model v3.2.0 at its fitted parameters (Lang et al., 2024)

#@description: |
#  The full fused, rule-based human cell-cycle model of Lang et al. (2024), version v3.2.0,
#  at the parameter set the authors estimated by fitting it to single-cell RPE-1 data. It is
#  the network of cell_cycle_oscillator_lang2024_fused.bngl with every number replaced by the
#  authors' fitted value: the Rb/E2F restriction point, cyclin E/A/B-CDK activity, the
#  APC/C:FZR1 (Cdh1), APC/C:CDC20 and FBXO5 (Emi1) degradation switches, the Greatwall
#  (MASTL)-ENSA/ARPP19-PP2A:B55 mitotic phosphatase loop, Wee1/Cdc25 control of CDK1,
#  FOXM1-driven transcription, and the p21 (CDKN1A), p27 (CDKN1B) and Skp2 CDK-inhibitor
#  module. All kinetics are mass action; BioNetGen generates 73 species and 332 reactions.
#
#  From the PEtab wild-type condition the model oscillates without damping. The 0-136471.6 s
#  (37.9 h) window of the fit holds two cell cycles of about 19 h. Within each cycle total
#  cyclin E peaks first (9.5 h into the run), then E2F1 (11.1 h) and Skp2 (16.4 h), then
#  cyclins A and B1 together with phospho-RB (18.5-18.6 h); p21 and p27 peak after mitosis
#  (about 23 h). Longer runs keep oscillating, with the period slowly shortening toward 17 h.
#  The eight functions obs_*__nuc_median() = offset + scale x (species sum) are the authors'
#  observation model for the nuclear median intensities of the Stallaert et al. (2022)
#  single-cell data that the model was fitted to.
#
#  Simulation intent: concentration-based, in arbitrary concentration units (AU) on the unit
#  volume compartment cell, time in seconds. The authors' model has no molar or
#  molecule-count basis, so it is ODE-only; NA and a reference volume are not defined.

#@keyword: |
#  cell cycle, restriction point, G1/S transition, G2/M transition, mitotic exit, E2F,
#  retinoblastoma protein, cyclin E, cyclin A, cyclin B, CDK1, APC/C, Cdh1, Cdc20, Emi1,
#  Greatwall, MASTL, ENSA, PP2A-B55, Wee1, Cdc25, FOXM1, p21, p27, Skp2, oscillator,
#  negative feedback, positive feedback, mass action, rule-based, parameter estimation,
#  PEtab, single-cell data, RPE-1, Homo sapiens

#@reference: |
#  Lang PF, Penas DR, Banga JR, Weindl D, Novak B (2024). Reusable rule-based cell cycle
#  model explains compartment-resolved dynamics of 16 observables in RPE-1 cells. PLoS
#  Comput Biol 20(1):e1011151. doi:10.1371/journal.pcbi.1011151
#
#  Model, PEtab problem and fitted parameters: paulflang/cell_cycle_petab at commit
#  362b6667a2096b6c1cab43ddc3a92d4174c5a087, versions/v3.2.0/ (MIT licence).
#  https://github.com/paulflang/cell_cycle_petab
#
#  Stallaert W, Kedziora KM, Taylor CD, Zikry TM, Ranek JS, Sobon HK, Taylor SR, Young CL,
#  Cook JG, Purvis JE (2022). The structure of the human cell cycle. Cell Syst
#  13(3):230-240.e3. doi:10.1016/j.cels.2021.10.007

#@note: |
#  Provenance. Generated by generator/build_lang2024_fitted.py; do not edit by hand. The
#  compartments, molecule types, seed species, observables, functions and reaction rules
#  blocks are byte-identical to cell_cycle_oscillator_lang2024_fused.bngl (sha256
#  4f4106485535a3ac8ed9e012606ba14e81a40222f46138c96209b131fa33b383), whose body is the
#  authors' cell_cycle_v3.2.0.bngl rule set. Every parameter value is the nominalValue of
#  the authors' versions/v3.2.0/results_20220421/parameters_v3.2.0_optimized.tsv, written as
#  the shortest decimal that reads back as the same double. Each seed species carries the
#  ic_wt_cell__* parameter that the authors' wt condition (experimentalCondition_v3.2.0.tsv)
#  assigns to it, so the initial state is the PEtab wt condition at the fitted values.
#
#  Relation to the _fused file. The network is identical; 187 of the 205 values differ.
#  _fused carries the pre-fit PEtab nominal values, except for 13 fixed promoter species,
#  where it carries the seed amounts of the authors' .bngl (free E2F-target promoters
#  8.448e-5, E2F-bound 1.552e-5; similar for the FOXM1 promoters). The PEtab wt condition,
#  used here, fixes free promoters at 1e-4 and promoter complexes at 0, except CCNB
#  promoter:FOXM1 at 1.9e-7 (the authors' 1.903e-7 rounded).
#
#  Verification (verify_lang2024_fitted.ipynb). Integrated from this file's generated
#  network, the eight observation functions reproduce the authors' own simulation of the
#  fitted model (results_20220421/optimisation.tsv, AMICI; plotted as Fig L in S1 Appendix
#  of the paper, "Supplementary Figure 12" in the repository README) at all 600 time
#  points.
#
#  Units. Rate-constant units follow from the reaction order of the rules that use each one:
#  AU/s (zeroth order), /s (first order), /AU/s (second order); seed amounts are in AU. The
#  compartment cell has volume 1, so BioNetGen's volume scaling of bimolecular rates is the
#  identity.
#
#  Formatting. The six body blocks are kept byte-identical to the authors' rule set, so they
#  keep its tab indentation, long lines and unit-free seed amounts; the header, parameters and
#  actions blocks follow the house style. Parameter lines for the longest ic_wt_cell__* names
#  exceed 100 characters because the names are the PEtab ids verbatim.

"""

ACTIONS = f"""\
begin actions

  generate_network({{overwrite=>1}})  # 73 species, 332 reactions

  # Two cell cycles from the PEtab wild-type condition

  #@protocol: |
  #  ODE time course from the PEtab wt initial state over the horizon of the authors'
  #  simulation output, 0 to 136471.6126 s. The authors sampled 600 non-uniform pseudotime
  #  points; this protocol uses a uniform grid finer than their smallest spacing (164.6 s),
  #  and the verification evaluates the same network on their points directly.
  #  print_functions=>1 writes the eight observation functions.

  #@figure: Fig L in S1 Appendix of Lang et al. (2024)

  simulate({{method=>"ode",suffix=>"ode",t_start=>0,t_end=>{T_END},\\
    n_steps=>{N_STEPS},print_functions=>1,atol=>{ATOL:g},rtol=>{RTOL:g}}})

end actions
"""


def parameters_block(fused: str) -> str:
    table = {r["parameterId"]: r for r in read_tsv(OPTIMIZED)}
    condition = read_tsv(CONDITIONS)
    if len(condition) != 1 or condition[0]["conditionId"] != "wt":
        raise ValueError("expected exactly one PEtab condition, `wt`")
    wt_parameters = {v for k, v in condition[0].items() if k != "conditionId"}

    names = [ln.split()[0] for ln in block_lines(fused, "parameters")]
    if set(names) != set(table) or len(names) != len(table):
        raise ValueError("the _fused parameters and the optimized table name different sets")
    ics = [n for n in names if n.startswith("ic_wt_cell__")]
    if set(ics) != wt_parameters:
        raise ValueError("the ic_wt_cell__* parameters are not the wt condition's parameters")
    seeds = [ln.rsplit(None, 1)[1] for ln in block_lines(fused, "seed species")]
    if sorted(seeds) != sorted(ics):
        raise ValueError("seed species are not each set by exactly one ic_wt_cell__* parameter")

    orders = rate_orders(fused)
    noise = [n for n in names if n.startswith("noiseParameter")]
    observation = [n for n in names if re.match(r"^[os][A-Z]", n)]
    kinetic = [n for n in names if n not in set(ics) | set(noise) | set(observation)]
    for n in kinetic:
        if len(orders.get(n, ())) != 1:
            raise ValueError(f"{n}: used at reaction orders {orders.get(n)}; cannot assign a unit")
    estimated = {n for n in names if table[n]["estimate"] == "1"}
    # The grouping below states which values were fitted; refuse to mislabel a group.
    if not set(kinetic) | set(observation) <= estimated or set(noise) & estimated:
        raise ValueError("estimate flags no longer match the group comments")

    def group(members: list[str], unit) -> list[str]:
        width = max(len(n) for n in members)
        vals = {n: fmt(table[n]["nominalValue"]) for n in members}
        vwidth = max(len(v) for v in vals.values())
        return [f"  {n:<{width}}  {vals[n]:<{vwidth}}  # {unit(n)}" for n in members]

    est_ics = [n for n in ics if n in estimated]
    fixed_ics = [n for n in ics if n not in estimated]
    out = ["begin parameters", ""]
    out += [
        f"  # Kinetic rate constants ({len(kinetic)}, all estimated). Units follow from the",
        "  # reaction order of the rules that use each one; see the header note.",
    ]
    out += group(kinetic, lambda n: UNITS[next(iter(orders[n]))])
    out += [
        "",
        f"  # Observation model ({len(observation)}, all estimated):",
        "  # obs_X__nuc_median() = offset oX + scale sX x (sum of the species measured).",
    ]
    out += group(observation, lambda n: "a.u." if n.startswith("o") else "a.u./AU")
    out += [
        "",
        f"  # Noise model ({len(noise)}, fixed): Gaussian sigma = "
        "noiseParameter1 + noiseParameter2.",
        "  # Not used in simulation; kept so the file carries the whole PEtab parameter table.",
    ]
    out += group(noise, lambda n: "a.u.")
    out += [
        "",
        f"  # Initial concentrations, estimated ({len(est_ics)}): the PEtab wt condition assigns",
        "  # one ic_wt_cell__* parameter to each seed species (names are the PEtab ids).",
    ]
    out += group(est_ics, lambda n: "AU")
    out += [
        "",
        f"  # Initial concentrations, fixed in the PEtab problem ({len(fixed_ics)}): promoters,",
        "  # promoter complexes and inputs, identical in the authors' nominal and optimized",
        "  # tables.",
    ]
    out += group(fixed_ics, lambda n: "AU")
    out += ["", "end parameters"]
    return "\n".join(out)


def build() -> str:
    check_pins()
    fused = FUSED.read_text(encoding="utf-8")
    parts = [HEADER, parameters_block(fused), ""]
    for name in BODY_BLOCKS:
        start, end = block_span(fused, name)
        parts += [fused[start:end], ""]
    parts += ["end model", "", ACTIONS]
    return "\n".join(parts)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true", help="compare with the committed file")
    args = ap.parse_args(argv)
    text = build()
    if args.check:
        same = TARGET.exists() and TARGET.read_text(encoding="utf-8") == text
        print(f"{TARGET.name}: {'matches' if same else 'DIFFERS from'} the generator output")
        return 0 if same else 1
    TARGET.write_text(text, encoding="utf-8")
    print(f"wrote {TARGET.relative_to(FOLDER)} ({sha256(TARGET)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
