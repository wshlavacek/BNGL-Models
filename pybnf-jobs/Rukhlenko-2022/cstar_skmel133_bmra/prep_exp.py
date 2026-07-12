#!/usr/bin/env python
"""Prepare the SKMEL-133 single-drug .exp training set for the BMRA-constrained job.

These are the AUTHORS' OWN pyBioNetFit fit targets, verbatim from
github.com/OleksiiR/cSTAR_Nature -> SKMEL-133_preproc/dose{1,2}_{ERK,AKT,SRC}inh.exp.
Each is the 24 h (86400 s) steady-state phospho / total-protein fold change (vs. the
no-drug DMSO baseline) under one inhibitor at dose1 (published I_<x>_conc) or dose2
(2x), with a per-point _SD companion -> chi_sq.

This is the full single-drug TRAINING set for the three PARAMETER-driven inhibitors
(ERK, AKT, SRC), which act purely through the rate-law factor 1/(1+I_<x>_conc) with no
sequestering binding, so an edition-2 Condition (setParameter) applies them exactly.
The mTOR/PKC/CDK inhibitors act by COMPETITIVE BINDING (their dose must re-initialise a
seed-species concentration via setConcentration, which an edition-2 Condition cannot do
-- the documented igf1r/lanl-PyBNF#474 limitation), so those three single-drug arms are
NOT included here; the drug-COMBINATION arms are the paper's validation set either way.

Two mechanical edits adapt each file to PyBNF edition-2 (identical to ../cstar_skmel133):
  1. Function columns take parentheses in the header (FC_pERK -> FC_pERK()); their _SD
     companions do NOT (PyBNF forms the noise column as <entity>_SD = FC_pERK_SD).
  2. The authors' t=0 row (all 1.0, SD 0.0) is dropped -- SD=0 breaks chi_sq, and
     `normalization = init` supplies the baseline. A single intermediate unmeasured
     (NaN) sample point is added so PyBNF has >=2 positive grid times.

Run from this folder with the authors' originals in ./src/ to regenerate the .exp.
"""
import os
import re

SRC = "src"                       # place the authors' original .exp files here
FUNCS = {"FC_tIRS", "FC_IRSI", "FC_pERK", "FC_pAKT", "FC_pSRC",
         "FC_pPKC", "FC_pS6K", "FC_pRB", "FC_MYC"}
USE = ["dose1_ERKinh.exp", "dose2_ERKinh.exp",
       "dose1_AKTinh.exp", "dose2_AKTinh.exp",
       "dose1_SRCinh.exp", "dose2_SRCinh.exp"]


def fixcol(c):
    m = re.match(r"^(FC_[A-Za-z0-9]+)(_SD)?$", c)
    if not (m and m.group(1) in FUNCS):
        return c
    return m.group(1) + ("_SD" if m.group(2) else "()")   # primary: parens; _SD: none


def main():
    for f in USE:
        lines = open(os.path.join(SRC, f)).read().splitlines()
        cols = lines[0].lstrip("#").split()
        header = "# " + "\t".join(fixcol(c) for c in cols)
        data = [l for l in lines[1:] if l.strip() and l.split()[0] not in ("0", "0.0")]
        nan_row = "\t".join(["43200"] + ["NaN"] * (len(cols) - 1))
        open(f, "w").write("\n".join([header, nan_row] + data) + "\n")
        print(f"wrote {f}")


if __name__ == "__main__":
    main()
