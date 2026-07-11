#!/usr/bin/env python
"""Extract the TrkA NGF-stimulation fold-change fit target (Fig. 4A of PMC9644236).

Source: OleksiiR/cSTAR_Nature, RPPA_DA/RPPA_data_trusted.csv (raw RPPA signal, 6
replicates per condition/timepoint) + RPPA_data_Trk_normalized_new.csv (for pTrk,
which is a Western-blot readout not on the RPPA panel).

For each model observable we take the SH-SY5Y/TrkA + NGF, DMSO (no inhibitor) columns
at 0/10/45 min, average replicates (nanmean over r1..r6), and normalize each series to
its own t=0 mean -> fold change (so t=0 == 1 by construction, matching the model's
`normalization = init`). Model time is in seconds: 10 min = 600 s, 45 min = 2700 s.
"""
import csv
import math

SRC = "cstar_src/RPPA_data_trusted.csv"
PTRK_SRC = "cstar_src/RPPA_data_Trk_normalized_new.csv"

# model function (exp header, WITH parens)  ->  antibody substring in the RPPA CSV
MAP = [
    ("FC_pAdRTK()", "ErbB-2/Her2/EGFR P Tyr1248/Tyr1173"),   # ERBB2/HER2 -> AdRTK layer
    ("FC_ppERK()",  "p44/42 MAPK (ERK1/2) P Thr202/Thr185"),  # ERK1/2
    ("FC_pAKT()",   "Akt P Ser473"),
    ("FC_pJNK()",   "SAPK/JNK P Thr183"),
    ("FC_pRSK()",   "p90 S6 kinase (Rsk1-3) P Thr359"),       # p90RSK
    ("FC_pS6K()",   "p70 S6 Kinase P Thr389"),
]
TIMES = [("0", 0), ("10", 600), ("45", 2700)]


def load_rows(path):
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        rows = {r[0].strip().strip('"'): r for r in reader if r}
    return header, rows


def cols_for(header, prefix):
    """Column indices whose name starts with `prefix` (e.g. 'TrkA_DMSO_10_')."""
    return [i for i, h in enumerate(header) if h.startswith(prefix)]


def nanmean(vals):
    xs = [v for v in vals if v is not None and not math.isnan(v)]
    return sum(xs) / len(xs) if xs else float("nan")


def to_f(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return float("nan")


def series(header, rows, antibody, cell="TrkA"):
    # find the row (substring match, quotes/space tolerant)
    key = next((k for k in rows if antibody in k), None)
    if key is None:
        raise KeyError(f"antibody not found: {antibody!r}")
    row = rows[key]
    out = {}
    for label, _sec in TIMES:
        idx = cols_for(header, f"{cell}_DMSO_{label}_")
        out[label] = nanmean([to_f(row[i]) for i in idx])
    t0 = out["0"]
    if not (t0 == t0) or t0 == 0:   # NaN or zero t0 -> fall back to no normalization
        return {k: v for k, v in out.items()}, key, False
    return {k: v / t0 for k, v in out.items()}, key, True


def main():
    header, rows = load_rows(SRC)
    ph, prows = load_rows(PTRK_SRC)

    # pTrk from the normalized Western file (its DMSO_0 columns are already 1)
    ptrk, pkey, pnorm = series(ph, prows, "pTrk")

    fc = {}     # obs -> {timelabel -> fold change}
    prov = []   # provenance table
    for obs, antibody in MAP:
        s, key, normed = series(header, rows, antibody)
        fc[obs] = s
        prov.append((obs, key, normed))
    fc["FC_pTRK()"] = ptrk
    prov.insert(0, ("FC_pTRK()", pkey + "  [normalized_new.csv, Western]", pnorm))

    order = ["FC_pTRK()", "FC_pAdRTK()", "FC_ppERK()", "FC_pAKT()",
             "FC_pJNK()", "FC_pRSK()", "FC_pS6K()"]

    # write the .exp
    lines = ["# time " + " ".join(order)]
    for label, sec in TIMES:
        vals = []
        for obs in order:
            v = fc[obs].get(label, float("nan"))
            vals.append("NaN" if (v != v) else f"{v:.4f}")
        lines.append(f"{sec} " + " ".join(vals))
    exp = "\n".join(lines) + "\n"
    with open("cstar_trka.exp", "w") as fh:
        fh.write(exp)

    print("=== provenance (model obs -> RPPA antibody, init-normalized?) ===")
    for obs, key, normed in prov:
        print(f"  {obs:14s} <- {key}   norm_t0={normed}")
    print("\n=== cstar_trka.exp ===")
    print(exp)


if __name__ == "__main__":
    main()
