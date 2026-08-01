#!/usr/bin/env python3
"""Build the .exp files for the nfkb_tnfa job from the digitized fig. S12 of Korwek 2023.

Source of truth is the committed CSV that
``models/innate_immune_response_korwek2023/digitize_korwek2023.py`` writes; this script only
reshapes it into PyBNF .exp tables, one per blot per replicate.

Four blots, each supplied as two replicates -- the blot reproduced in the figure (``_r1``) and
the additional replicate (``_r2``):

    nfkb_tnfa_wt_nuclear   fig. S12A  nuclear NF-kB (RelA), WT,     9 lanes over 6 h
    nfkb_tnfa_wt_fine      fig. S12B  p-IKK / IkBa / A20, WT,      10 lanes over 3 h
    nfkb_tnfa_wt_long      fig. S12B  p-IKK / IkBa / A20, WT,       8 lanes over 6 h
    nfkb_tnfa_a20ko        fig. S12B  p-IKK / IkBa, A20 KO,         8 lanes over 6 h

Column names are the measurement-model observable ids the conf declares, one per blot per
protein, because each blot carries its own unpublished normalization constant.

A marker flagged in the CSV as sitting on the axis floor or clipped by the axis top is written
as ``NaN``: its value is not readable from the figure, so it is not a measurement. The A20
panel of the A20 KO blot has no data at all (A20 is absent in those cells and the figure draws
a flat placeholder), so it contributes no column.

Per-point sigma
---------------
``<obs>_SD`` is a flat 25% of the measured value, which makes the objective a relative-error
least squares -- the multiplicative-error criterion Korwek et al. use throughout. The 25% is
the paper's own number: table 1 reports an "average multiplicative error between experimental
replicates (for Western blots)" of 1.24, over all 2915 points they fit.

The two replicates digitized here scatter more than that -- multiplicative error 1.55 over the
35 lanes where both are separately readable -- but that figure is biased upward and should not
be used in its place: markers only resolve into two dots when the replicates *disagree*, so
the lanes where the two blots agree are exactly the ones that cannot be measured.
"""
from __future__ import annotations

import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
CSV = (HERE.parents[2] / "models" / "innate_immune_response_korwek2023" / "reference"
       / "korwek2023_figS12_digitized.csv")

REL_SD = 0.25   # table 1 of Korwek et al. (2023): replicate multiplicative error 1.24

# published panel -> (exp stem, {figure observable: conf observable id})
PANELS = {
    "S12A":            ("nfkb_tnfa_wt_nuclear", {"NF-kB nuclear": "NFkBn_nuc"}),
    "S12B_WT_fine":    ("nfkb_tnfa_wt_fine",
                        {"p-IKK": "pIKK_fine", "IkBa": "IkBa_fine", "A20": "A20_fine"}),
    "S12B_WT_long":    ("nfkb_tnfa_wt_long",
                        {"p-IKK": "pIKK_long", "IkBa": "IkBa_long", "A20": "A20_long"}),
    "S12B_A20KO_long": ("nfkb_tnfa_a20ko", {"p-IKK": "pIKK_ko", "IkBa": "IkBa_ko"}),
}
REPLICATES = {"r1": "blot_shown", "r2": "blot_replicate"}


def read_rows():
    with open(CSV) as fh:
        return list(csv.DictReader(fh))


def value(row, series):
    """The measurement, or None when the marker is unreadable."""
    if row[series] == "" or row[f"{series}_on_axis_floor"] == "1" \
            or row[f"{series}_above_axis"] == "1":
        return None
    return float(row[series])


def main():
    rows = read_rows()
    total = 0
    for panel, (stem, obsmap) in PANELS.items():
        sub = [r for r in rows if r["panel"] == panel]
        times = sorted({float(r["time_h"]) for r in sub})
        for tag, series in REPLICATES.items():
            cols, data, n = [], {}, 0
            for fig_obs, obs_id in obsmap.items():
                by_time = {float(r["time_h"]): r for r in sub if r["observable"] == fig_obs}
                vals = [value(by_time[t], series) if t in by_time else None for t in times]
                if all(v is None for v in vals):
                    continue
                cols.append(obs_id)
                data[obs_id] = vals
                n += sum(v is not None for v in vals)
            if not cols:
                continue
            header = ["time"] + [c for obs in cols for c in (obs, f"{obs}_SD")]
            path = HERE / f"{stem}_{tag}.exp"
            with open(path, "w") as fh:
                fh.write("# " + " ".join(header) + "\n")
                for i, t in enumerate(times):
                    fields = [f"{round(t * 3600):7d}"]
                    for obs in cols:
                        v = data[obs][i]
                        fields += ["      NaN", "      NaN"] if v is None else \
                                  [f"{v:9.4f}", f"{v * REL_SD:9.4f}"]
                    fh.write("".join(fields) + "\n")
            print(f"{path.name:34s} {len(times):2d} time points, {n:2d} measurements")
            total += n
    print(f"\ntotal: {total} measurements")


if __name__ == "__main__":
    main()
