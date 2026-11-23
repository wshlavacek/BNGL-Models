#!/usr/bin/env python
"""Check cstar_trkb.bngl against the authors' TrkB_S_model.bngl, with BioNetGen alone.

The authors' file is github.com/OleksiiR/cSTAR_Nature @ bfb092f, Trk_AB_models/
TrkB_S_model.bngl, pinned by sha256 in upstream.json and not vendored (the repository grants no
licence). BNG2.pl runs it with its own actions block, which is the published protocol: no BDNF
and no inhibitors for 50000 s, S reset to 2*S0, another 50000 s, saveConcentrations, then for
each arm resetConcentrations, BDNF (TRKl = 10) and at most one inhibitor, 2700 s. The slug's
file runs the same protocol, translated by `cstar_upstream.translate_actions`:

  TRKl(Rec) = 0 / 10   ->  Lig_on = 0 / 1. The slug keeps the ligand at TRKl_conc = 10 and
                           gates every ligand-association rate by Lig_on (the hetero-dimer
                           family included), so during the equilibrations only the
                           free-ligand species differs (0 against 10); that is checked
                           separately ("gated_species").
  any other species    ->  the same setConcentration; an inhibitor dose named by a parameter
                           the slug sets to 0 is written out as the authors' literal.

Writes upstream_comparison.json (every metric, tolerance and verdict) and exits 0 only if
every asserted check passes. About 30 s on one core, nearly all of it BNG2.pl.

Usage (BNGPATH = the BioNetGen folder with BNG2.pl; the run directory is git-ignored):
    uv run python pybnf-jobs/Rukhlenko-2022/cstar_trkb/compare_upstream.py --fetch
    uv run python pybnf-jobs/Rukhlenko-2022/cstar_trkb/compare_upstream.py --upstream DIR
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import cstar_upstream as cu  # noqa: E402

LIGAND = "TRKl(Rec)"
LIGAND_CONC = 10.0  # TRKl_conc in both files
S_LABEL = "S()"  # the phenotype coordinate; sampled over every if() branch


def lig_on(value, text):
    if value == 0.0:
        return {"Lig_on": "0"}
    if value == LIGAND_CONC:
        return {"Lig_on": "1"}
    sys.exit(f"no Lig_on equivalent of TRKl(Rec) = {value}")


def ligand_on(rng, authors, job):
    y = cu.sample_lognormal(rng, authors, spread=(1e-3, 12.0), signed=(S_LABEL,))
    return y, y.copy(), {"Lig_on": 1}


def ligand_off(rng, authors, job):
    y = cu.sample_lognormal(rng, authors, spread=(1e-3, 12.0), signed=(S_LABEL,))
    y_j = y.copy()
    i = cu.species_index(authors, LIGAND)
    y[i], y_j[i] = 0.0, LIGAND_CONC
    return y, y_j, {"Lig_on": 0}


def from_seeds(calls):
    """Pre-equilibrate from the seeds to steady state, with no S reset."""
    calls = cu.set_t_end(calls, "equilibration1", 1e8)
    return cu.drop_calls(
        calls,
        lambda n, a: (
            (n == "setConcentration" and a.startswith('"S"'))
            or (n == "simulate" and cu.simulate_suffix(a) == "equilibration2")
        ),
    )


def after_reset(calls):
    """Both equilibration phases run to steady state; the S reset kept."""
    return cu.set_t_end(cu.set_t_end(calls, "equilibration1", 1e8), "equilibration2", 1e8)


ASSOCIATION = [
    "_rateLaw1",
    "_rateLaw3",
    "_rateLaw5",
    "_rateLaw7",
    "_rateLaw9",
    "_rateLaw19",
    "_rateLaw21",
    "_rateLaw23",
    "_rateLaw25",
    # TrkB only: ligand binding to the TrkB-ErbB hetero-dimer, sqrt(fh*factor)
    "_rateLaw34",
    "_rateLaw36",
    "_rateLaw38",
    "_rateLaw40",
]

SPEC = cu.Spec(
    slug="cstar_trkb",
    job_bngl=HERE / "cstar_trkb.bngl",
    gates={LIGAND: cu.Gate(params=lig_on, keep_species=False)},
    conditions=[
        ("BDNF added: Lig_on = 1, the same state", ligand_on, True),
        (
            "equilibration: authors TRKl(Rec) = 0; slug Lig_on = 0, TRKl(Rec) = 10",
            ligand_off,
            False,
        ),
    ],
    expected={
        "parameters_only_job": ["Lig_on"],
        # The DMSO arm: the slug zeroes the six inhibitor doses, and its ligand-association
        # rate laws carry the Lig_on factor (0 in the file).
        "parameters_differing_value": [
            "I_AKT_conc",
            "I_ERK_conc",
            "I_JNK_conc",
            "I_RSK_conc",
            "I_S6K_conc",
            "I_TRK_conc",
        ]
        + ASSOCIATION,
    },
    probes={
        "steady_state_from_seeds": (from_seeds, [("equilibration2", "equilibration1")]),
        "steady_state_after_S_reset": (after_reset, [("equilibration2", "equilibration2")]),
    },
    watch=["Sval", "FC_pAdRTK", "FC_ppERK", "FC_pAKT"],
    notes={
        "protocol": "The published protocol is the authors' actions block, run verbatim. Its "
        "pre-equilibration is two finite 50000 s phases with S reset to 2*S0 "
        "between them, not a steady state; protocol_sensitivity reports how far "
        "a steady-state pre-equilibration moves the state and the trajectories.",
    },
)

if __name__ == "__main__":
    sys.exit(cu.main(SPEC))
