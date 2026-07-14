#!/usr/bin/env python
"""Prepare the SKMEL-133 single-drug .exp for the PAPER-EXACT BMRA job (all SIX inhibitors).

Originals (in ./src/) are the AUTHORS' OWN pyBioNetFit fit targets, verbatim from
github.com/OleksiiR/cSTAR_Nature -> SKMEL-133_preproc/:
  * dose{1,2}_{ERK,AKT,SRC,PKC,mTOR,CDK}inh.exp     -- 9 phospho/total FOLD-CHANGE readouts
  * dose{1,2}_S_{ERK,AKT,SRC,PKC,mTOR,CDK}inh.exp   -- the DPD / state-transition coordinate Sval
Each is the 24 h (86400 s) steady state under one inhibitor at dose1 (published I_<x>_conc) or
dose2 (2x), with a per-point _SD companion -> chi_sq.

Unlike ../cstar_skmel133_bmra (which shipped only the three PARAMETER-driven inhibitors
ERK/AKT/SRC), this EXACT build ships ALL SIX: the mTOR/PKC/CDK inhibitors act by COMPETITIVE
BINDING and are now applied in edition 2 via a quoted-species setConcentration condition
(lanl/PyBNF#474: `perturbations: "IPKC(PKCBD)" = <dose>`), so their single-drug arms join the
training set. The drug-combination arms remain the paper's validation set (not fit here).

TRAINING ARMS (11): dose1 for all six inhibitors + dose2 for ERK/AKT/SRC/PKC/mTOR. The
authors' `dose2_CDKinh.exp` carries NO _SD columns (a single unreplicated measurement), so it
is incompatible with chi_sq and is EXCLUDED (the CDK inhibitor is represented by dose1_CDKinh,
which has per-point SD). The STV (Sval) arms are NOT fit -- the DPD S is in physical distance
units incommensurate with BMRA's standardized space (see VALIDATION.md), and its drug responses
swing across zero (e.g. Sval 10 -> -6.7), so an absolute-Sval fit against the physical betas is
ill-posed; the DPD row is constrained by sign. (The STV originals remain in ./src/ for reference.)

Two mechanical edits per FC file (identical to ../cstar_skmel133_bmra):
  1. Function columns take parentheses in the header (FC_pERK -> FC_pERK()); their _SD
     companions do NOT (PyBNF forms the noise column as <entity>_SD = FC_pERK_SD).
  2. The authors' t=0 row (all 1.0, SD 0.0) is dropped -- SD=0 breaks chi_sq, and
     `normalization = init` (per-observable) supplies the fold-change baseline. A single
     intermediate unmeasured (NaN) sample at 43200 s gives PyBNF >=2 positive grid times.
The STV (Sval) files keep the ABSOLUTE distance-unit value (NO normalization -- the model's
Sval is in the same data units); same t=0-drop + NaN-row treatment. Sval has no parens (it is
a `Species` observable, not a function).

Run from this folder with the authors' originals in ./src/.
"""
import os, re

SRC = "src"
FUNCS = {"FC_tIRS", "FC_IRSI", "FC_pERK", "FC_pAKT", "FC_pSRC",
         "FC_pPKC", "FC_pS6K", "FC_pRB", "FC_MYC"}
INH = ["ERK", "AKT", "SRC", "PKC", "mTOR", "CDK"]
# dose1 for all six + dose2 for all but CDK (authors' dose2_CDKinh.exp lacks _SD -> not chi_sq).
FC = [f"dose1_{x}inh.exp" for x in INH] + [f"dose2_{x}inh.exp" for x in INH if x != "CDK"]


def fixcol(c):
    m = re.match(r"^(FC_[A-Za-z0-9]+)(_SD)?$", c)
    if not (m and m.group(1) in FUNCS):
        return c
    return m.group(1) + ("_SD" if m.group(2) else "()")   # primary: parens; _SD: none


def prep(fn):
    lines = open(os.path.join(SRC, fn)).read().splitlines()
    cols = lines[0].lstrip("#").split()
    header = "# " + "\t".join(fixcol(c) for c in cols)
    data = [l for l in lines[1:] if l.strip() and l.split()[0] not in ("0", "0.0")]
    nan_row = "\t".join(["43200"] + ["NaN"] * (len(cols) - 1))
    open(fn, "w").write("\n".join([header, nan_row] + data) + "\n")
    print(f"wrote {fn}")


def main():
    for fn in FC:
        prep(fn)


if __name__ == "__main__":
    main()
