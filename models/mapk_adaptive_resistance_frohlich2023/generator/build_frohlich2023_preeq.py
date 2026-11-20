#!/usr/bin/env python
"""Build `mapk_adaptive_resistance_frohlich2023_preeq.bngl` from the primary file and the
authors' own pre-equilibrated state.

The `_preeq` file is the primary model seeded the way the authors' pipeline seeds the Fig. 2A
simulation: with the steady state at 1 uM vemurafenib and no EGF. Nothing about the model
changes, only the seed amounts:

  * every block except the header, the seed species and the actions -- the unit-convention
    comment, parameters, compartments, molecule types, observables, functions, energy
    patterns and reaction rules -- is copied byte for byte from the primary
    `mapk_adaptive_resistance_frohlich2023.bngl` (pinned below by sha256). Before copying,
    its 97 calibrated values are checked against row 0 of the authors' vendored parameter
    table (relative difference <= 1e-12);
  * the seed species block keeps the primary's lines, species and order, and replaces each
    amount (uM; FPKM for the three transcripts) with the authors' steady state, read from
    `reference/frohlich2023_marm2_adapted_state.csv` (column `authors_steady_state_uM`).
    That file is written by `run_frohlich2023.py`, which computes the steady state with the
    authors' own PySB model built from the vendored MARM2.0 release. The two clamped species
    keep the primary's lines: vemurafenib at 1 uM and free EGF at `initEGF`, the EGF step at
    t = 0. Species whose steady state is below SEED_FLOOR stay unseeded, as in the primary.

Values are written to 13 significant digits, the primary's precision. Nothing here runs
BioNetGen or PySB, and nothing is typed by hand: re-running this script must reproduce the
committed `.bngl` exactly (`--check`; `tests/test_frohlich2023_preeq_generator.py`).

Usage (from anywhere):
    python generator/build_frohlich2023_preeq.py           # (re)write the .bngl
    python generator/build_frohlich2023_preeq.py --check   # exit 1 if the committed file differs
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
PRIMARY = FOLDER / "mapk_adaptive_resistance_frohlich2023.bngl"
TARGET = FOLDER / "mapk_adaptive_resistance_frohlich2023_preeq.bngl"
STATE = REF / "frohlich2023_marm2_adapted_state.csv"
TABLE = REF / "frohlich2023_marm2_RTKERK_pRAF_EGF_EGFR_MEKi_PRAFi_RAFi.csv"

# Inputs, pinned: the primary file (whose sha256 is also pinned downstream), the authors'
# parameter table as vendored from labsyspharm/marm2-supplement @ 7e97789 (see
# reference/frohlich2023_marm2_upstream.json), and the steady state run_frohlich2023.py wrote
# from the authors' model. Any change must be re-derived and re-verified deliberately.
PINS = {
    PRIMARY: "51856655ec20a825bb8f457bfa0615a789fad4d01157c839de5a0c4b2ae09452",
    TABLE: "6a25c61d1e6373cf504cee162de561f60fceee92bcc779691e0066ccefa27ee4",
    STATE: "beb91a7b5de8dd66f8441933887d5b3897aab4b463fb8d2d263024855f1d9a91",
}

SEED_FLOOR = 1e-15  # uM, ~6e-10 molecules per 1 pL cell; the primary seeds nothing below it
CLAMPED = ("@CP::$RAFi(raf)", "@CP::$EGF(rtk)")
# Parameters of the primary that are not in the authors' table: the experimental condition
# (doses, constants under the house names) and expressions derived from calibrated values.
CONDITION = {"NA", "V_ref", "h_mem", "m_Da_EGF", "RAFi_0", "EGF_0", "EGFR_crispr"}

BODY_BLOCKS = (
    "parameters",
    "compartments",
    "molecule types",
)
TAIL_BLOCKS = (
    "observables",
    "functions",
    "energy patterns",
    "reaction rules",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_pins() -> None:
    bad = [p.name for p, digest in PINS.items() if sha256(p) != digest]
    if bad:
        raise SystemExit(f"sha256 mismatch, refusing to build: {', '.join(bad)}")


def block_span(text: str, name: str) -> tuple[int, int]:
    """Character span of `begin <name>` .. `end <name>` (both lines included)."""
    m = re.search(rf"^begin {name}[ \t]*\n.*?^end {name}[ \t]*$", text, re.S | re.M)
    if m is None:
        raise ValueError(f"no `{name}` block")
    return m.start(), m.end()


def check_parameters(primary: str) -> int:
    """The primary's calibrated values equal row 0 of the authors' table; returns the count."""
    with open(TABLE, newline="") as fh:
        rows = list(csv.DictReader(fh))
    fvals = [float(r["fval"]) for r in rows]
    if fvals.index(min(fvals)) != 0:
        raise ValueError("row 0 of the authors' table is not the best-ranked set")
    row0 = {k: float(v) for k, v in rows[0].items() if k not in ("", "fval")}
    start, end = block_span(primary, "parameters")
    checked = 0
    for line in primary[start:end].splitlines()[1:-1]:
        fields = line.split("#", 1)[0].split()
        if len(fields) != 2:
            continue
        name, value = fields
        try:
            value = float(value)
        except ValueError:
            continue  # a derived expression
        if name in CONDITION:
            continue
        # MARM.parameters.specialise_par_name(name, prafi=None, rafi="Vemurafenib", meki=None)
        ref = row0[name.replace("RAFi", "Vemurafenib")]
        if abs(value - ref) > 1e-12 * abs(ref):
            raise ValueError(f"{name}: {value!r} is not the authors' row-0 value {ref!r}")
        checked += 1
    if checked != 97:
        raise ValueError(f"expected 97 calibrated values in the primary, found {checked}")
    return checked


def seed_block(primary: str) -> str:
    with open(STATE, newline="") as fh:
        state = {r["species"]: float(r["authors_steady_state_uM"]) for r in csv.DictReader(fh)}
    start, end = block_span(primary, "seed species")
    lines = primary[start:end].splitlines()
    out, seeded = [lines[0]], set()
    for line in lines[1:-1]:
        m = re.fullmatch(r"  (\S+)  (\S+)  # (uM|FPKM)", line)
        if m is None:
            raise ValueError(f"unexpected seed line: {line!r}")
        species = m.group(1)
        seeded.add(species)
        if species in CLAMPED:
            out.append(line)
            continue
        value = state[species]
        if value < SEED_FLOOR:
            raise ValueError(f"{species}: steady state {value:g} uM is below the seed floor")
        out.append(f"  {species}  {value:.13g}  # {m.group(3)}")
    missing = [s for s, v in state.items()
               if s not in seeded and "$" not in s and v >= SEED_FLOOR]
    if missing:
        raise ValueError(f"{len(missing)} species above the floor are not seeded: {missing[:3]}")
    out.append(lines[-1])
    return "\n".join(out)


HEADER = """\
begin model

#@title: Adaptive MAPK signaling from the authors' steady state (Fröhlich et al., 2023)

#@description: |
#  Condition-specific MARM2.0 energy-pattern model of EGFR-RAS-RAF-MEK-ERK
#  signaling in BRAF(V600E) melanoma cells, seeded as the authors seed it for
#  Fig. 2A: at the steady state reached under 1 uM vemurafenib without EGF,
#  their model of the experiment's 24 h drug adaptation. The active ODE
#  protocol applies 100 ng/mL EGF and reproduces the authors' own simulation
#  of the Fig. 2A pERK pulse. Feedback includes ERK-dependent DUSP, SPRY, and
#  EGFR transcription, receptor internalization/recycling, SOS1 and RAF
#  phosphorylation, and allosteric RAF dimer/drug interactions.

#@keyword: |
#  MAPK signaling, adaptive resistance, transient, mass action,
#  energy patterns, transcriptional feedback, steady state,
#  pre-equilibration, Homo sapiens

#@reference: |
#  Fröhlich F, Gerosa L, Muhlich JL, Sorger PK (2023). Mechanistic model of
#  MAPK signaling reveals how allostery and rewiring contribute to drug
#  resistance. Molecular Systems Biology 19:e10988.
#  doi:10.15252/msb.202210988.
#
#  Model, calibrated parameters and protocol: labsyspharm/marm2-supplement at
#  commit 7e977899045e3c82199c827a37d79797f4fcf1a9 (Code EV1 of the article,
#  CC BY 4.0). https://github.com/labsyspharm/marm2-supplement

#@note: |
#  This is a concentration-based deterministic model. NA and V_ref provide
#  the source model's molecules/cell-to-uM conversion. Its compartment sizes
#  are normalized topological factors, while h_mem records the membrane-layer
#  convention needed for population conversion. ODE simulation is the source
#  protocol; discrete simulation is not appropriate for the normalized
#  compartments or the high cellular copy numbers.
#
#  Provenance. Generated by generator/build_frohlich2023_preeq.py; do not
#  edit by hand. Everything but this header, the seed amounts and the actions
#  is copied byte for byte from mapk_adaptive_resistance_frohlich2023.bngl
#  (sha256 51856655ec20a825bb8f457bfa0615a789fad4d01157c839de5a0c4b2ae09452),
#  whose 97 calibrated values the generator checks against row 0 of the
#  authors' parameter table. The seed amounts are the authors' steady state
#  as run_frohlich2023.py computes it from their PySB model
#  (reference/frohlich2023_marm2_RTKERK__pRAF.py at row 0), integrated from
#  the model's own initial values at 1 uM vemurafenib and no EGF: at t = 200 h
#  no species had changed by more than 3.8e-9 (relative, with a 1e-9 uM floor)
#  since t = 100 h, and AMICI's steady-state test with the authors' settings
#  gives a WRMS of dx/dt of 1.5e-5 (converged below 1). They are written to
#  13 significant digits. As in the primary, 333 species are seeded: the 331
#  whose steady state is at or above 1e-15 uM, and the two clamped ones. The
#  other 1,279 are all below 9.5e-16 uM in magnitude (the 252 EGF-bound ones
#  below 2e-27 uM) and stay unseeded. Vemurafenib is clamped at 1 uM and free
#  EGF is seeded at initEGF, which is the step to 100 ng/mL at t = 0.
#
#  Relation to the primary. Same network, parameters, rules and protocol; only
#  the seed amounts differ. The primary's seeds are the model's raw initial
#  values (free proteins only, no transcripts, DUSP, SPRY or EGFR) integrated
#  for 24 h at 1 uM vemurafenib, which is not the authors' protocol: rows 0-9
#  of their processed data table set the *_preeq columns with t_presim = 0,
#  which AMICI runs as a pre-equilibration to steady state, and the paper
#  states that "model simulations were pre-equilibrated to these conditions".
#  The two seed states differ by up to 2.3e-2 (relative; median 8.5e-5; the
#  worst an EGFR dimer), and the two files' run_network pERK_IF_obs
#  trajectories by up to 1.1e-2 (0.0034 IF units).
#
#  Verification. run_frohlich2023.py compares this file with the authors'
#  release under the authors' protocol: parameters, network (identical under
#  a species bijection), equations, seeds (within 4.6e-13 of the authors'
#  state), and the 0-2 h trajectories of pERK_IF_obs, pMEK_IF_obs, all
#  observables and all species. Against the authors' trajectory, the network
#  integrated by independent SciPy code agrees to 2.7e-11 in pERK_IF_obs, and
#  the committed run_network output to 1.8e-6 over pERK_IF_obs, pMEK_IF_obs
#  and the 30 observables (tolerance 1e-5). It draws
#  verify_frohlich2023_preeq.png; every number is in
#  reference/frohlich2023_marm2_comparison.json.
"""


ACTIONS = """\
begin actions
  # EGF pulse from the authors' vemurafenib-adapted steady state; ODE is the
  # source model's simulation mode.

  #@protocol: 100 ng/mL EGF pulse from the steady state at 1 uM vemurafenib

  #@figure: Fig. 2A in Fröhlich et al. (2023)

  generate_network({overwrite=>1})  # 1612 species, 17516 reactions
  simulate({method=>"ode",suffix=>"ode",t_start=>0,t_end=>2,n_steps=>200,\\
    print_functions=>1,atol=>1e-11,rtol=>1e-9})
end actions
"""


def build() -> str:
    check_pins()
    primary = PRIMARY.read_text(encoding="utf-8")
    check_parameters(primary)
    # The primary's unit-convention comment sits between its header and the parameters.
    conventions = re.search(r"^# Concentrations are uM.*?(?=^begin parameters)", primary,
                            re.S | re.M)
    if conventions is None:
        raise ValueError("the primary's unit-convention comment was not found")
    parts = [HEADER, conventions.group(0).rstrip("\n"), ""]
    for name in BODY_BLOCKS:
        start, end = block_span(primary, name)
        parts += [primary[start:end], ""]
    parts += [seed_block(primary), ""]
    for name in TAIL_BLOCKS:
        start, end = block_span(primary, name)
        parts += [primary[start:end], ""]
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
