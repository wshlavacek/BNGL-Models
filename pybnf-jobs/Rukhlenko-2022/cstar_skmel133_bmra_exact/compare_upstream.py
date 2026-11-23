#!/usr/bin/env python
"""Check cstar_skmel133_bmra_exact.bngl against the authors' SKMEL-133-3.bngl, BioNetGen only.

The authors' file is github.com/OleksiiR/cSTAR_Nature @ bfb092f, SKMEL-133_models/
SKMEL-133-3.bngl, pinned by sha256 in upstream.json and not vendored (the repository grants no
licence). BNG2.pl runs it with its own actions block, which is the published protocol:
86400 s with every inhibitor at 0, saveConcentrations, then for each arm resetConcentrations
and the inhibitor doses of that arm: AKTi, IRSi and AKTi + IRSi for 30000 s, and the ERK, AKT,
PKC, SRC, mTOR and CDK inhibitors at dose 1 and dose 2 for 86400 s (15 arms). The slug's file
runs the same protocol, translated by `cstar_upstream.translate_actions`:

  IERK, IAKT, ISRC     ->  the same setConcentration, and before each simulate the parameter
                           I_ERK_conc / I_AKT_conc / I_SRC_conc set to the species' value:
                           the slug's ERK, AKT and SRC phosphorylation rate laws read that
                           parameter where the authors' read the fixed inhibitor species.
  any other species    ->  the same setConcentration; a dose named by a parameter the slug
                           sets to 0 is written out as the authors' literal.

Writes upstream_comparison.json (every metric, tolerance and verdict) and exits 0 only if
every asserted check passes. About 1 minute on one core, nearly all of it BNG2.pl.

Usage (BNGPATH = the BioNetGen folder with BNG2.pl; the run directory is git-ignored):
    uv run python pybnf-jobs/Rukhlenko-2022/cstar_skmel133_bmra_exact/compare_upstream.py --fetch
    uv run python pybnf-jobs/Rukhlenko-2022/cstar_skmel133_bmra_exact/compare_upstream.py \\
        --upstream DIR
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import cstar_upstream as cu  # noqa: E402

S_LABEL = "S()"  # the phenotype coordinate; sampled over every if() branch
# Inhibitor species whose effect the slug routes through a parameter of the same value.
ROUTED = {"IERK(ERKBD)": "I_ERK_conc", "IAKT(AKTBD)": "I_AKT_conc", "ISRC(SRCBD)": "I_SRC_conc"}


def routed(param):
    return lambda value, text: {param: text}


def drugs_as_sampled(rng, authors, job):
    y = cu.sample_lognormal(rng, authors, spread=(1e-3, 12.0), signed=(S_LABEL,))
    params = {p: float(y[cu.species_index(authors, s)]) for s, p in ROUTED.items()}
    return y, y.copy(), params


def steady_state(calls):
    """Pre-equilibrate to steady state instead of for 86400 s."""
    return cu.set_t_end(calls, "equilibration1", 1e8)


SPEC = cu.Spec(
    slug="cstar_skmel133_bmra_exact",
    job_bngl=HERE / "cstar_skmel133_bmra_exact.bngl",
    gates={s: cu.Gate(params=routed(p), keep_species=True) for s, p in ROUTED.items()},
    conditions=[
        (
            "inhibitors as sampled; the slug's I_ERK/AKT/SRC_conc = the authors' species",
            drugs_as_sampled,
            True,
        ),
    ],
    expected={
        # The DMSO baseline: the slug zeroes the eight inhibitor doses.
        "parameters_differing_value": [
            "I_AKT_conc",
            "I_CDK_conc",
            "I_ERK_conc",
            "I_IRS_conc",
            "I_MYC_conc",
            "I_PKC_conc",
            "I_SRC_conc",
            "I_mTOR_conc",
        ],
        # The ERK, AKT and SRC phosphorylation rate laws: 1/(1+I_X_conc) for 1/(1+I_X).
        "functions_differing_text": ["_rateLaw10", "_rateLaw12", "_rateLaw16"],
        # A conserved carrier (total ERK) for the constraint functions; no rate law uses them.
        "observables_only_job": ["totERK"],
    },
    probes={
        "steady_state_pre_equilibration": (steady_state, [("equilibration1", "equilibration1")]),
    },
    watch=["Sval", "tIRS", "FC_pERK", "FC_pAKT"],
    notes={
        "protocol": "The published protocol is the authors' actions block, run verbatim. Its "
        "pre-equilibration is a finite 86400 s phase; protocol_sensitivity reports how "
        "far a steady-state pre-equilibration moves the state and each arm.",
    },
)

if __name__ == "__main__":
    sys.exit(cu.main(SPEC))
