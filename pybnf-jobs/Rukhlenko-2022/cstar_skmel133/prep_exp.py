#!/usr/bin/env python
"""Prepare the SKMEL-133 inhibitor .exp files for the edition-2 job.

These are the AUTHORS' OWN pyBioNetFit fit targets, downloaded verbatim from
github.com/OleksiiR/cSTAR_Nature -> SKMEL-133_preproc/{dose1_ERKinh,dose2_ERKinh,
dose1_AKTinh,dose1_SRCinh}.exp. Each is the 24 h (86400 s) steady-state phospho /
total-protein fold change (vs. the no-drug DMSO baseline) under one inhibitor, with a
per-point _SD companion -> chi_sq.

Two mechanical edits adapt them to PyBNF edition-2:
  1. Function columns take parentheses in the header (FC_pERK -> FC_pERK()); their _SD
     companions do NOT (PyBNF forms the noise column as <entity>_SD = FC_pERK_SD).
  2. The authors' t=0 row (all 1.0, SD 0.0) is dropped -- SD=0 breaks chi_sq, and the
     model's own `normalization = init` supplies the baseline. A single intermediate
     unmeasured (NaN) sample point is added so PyBNF has >=2 positive grid times.

Run from this folder with the source files in ./src/ to regenerate the .exp files.
"""
import os
import re

SRC = "src"                       # place the authors' original .exp files here
FUNCS = {"FC_tIRS", "FC_IRSI", "FC_pERK", "FC_pAKT", "FC_pSRC",
         "FC_pPKC", "FC_pS6K", "FC_pRB", "FC_MYC"}
USE = ["dose1_ERKinh.exp", "dose2_ERKinh.exp", "dose1_AKTinh.exp", "dose1_SRCinh.exp"]


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
