#!/usr/bin/env python
"""Extract the WT pp-ATF2(T69/T71) anisomycin time course fit target (PMC7666158).

Source: the paper's Source Data workbook, sheet ``Figure_4`` (Fig. 4a), block
"Anti-ppATF2(T69/T71), treated/untreated" -- the mean western-blot ratio of
anisomycin-stimulated to unstimulated HEK293T cells over 0-41 min, with per-point SD
(error propagation over n=3 replicates). We take the WT column (mean = col C, SD = col C
of the SD sub-block) and write it as the fit target for the model observable ``pT69pT71``.

Download the workbook (5.8 MB) from the article's supplementary materials:
  https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19582-3/MediaObjects/41467_2020_19582_MOESM6_ESM.xlsx
and pass its path as argv[1] (default: ./41467_2020_19582_MOESM6_ESM.xlsx).

"treated/untreated" is a fold change vs the unstimulated control, so the conf compares it
to the model's pT69pT71 under ``normalization = init`` (each simulated series divided by
its own t=0 / basal value). Model time is in seconds: t[s] = t[min] * 60. Only stdlib is
used (the xlsx is unzipped and its sheet XML parsed directly).
"""
import sys
import re
import zipfile
import xml.etree.ElementTree as ET

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
XLSX = sys.argv[1] if len(sys.argv) > 1 else "41467_2020_19582_MOESM6_ESM.xlsx"
OUT = "fig7b_timecourse.exp"

# Figure_4 is the 4th sheet (sheet4.xml); the MEAN block is rows 6-24, the SD block rows
# 29-47; WT is column C (=3), the independent variable t(min) is column B (=2).
SHEET = "xl/worksheets/sheet4.xml"
MEAN_ROW0, SD_ROW0, N = 6, 29, 19
T_COL, WT_COL = 2, 3


def colnum(ref):
    letters = re.match(r"[A-Z]+", ref).group()
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n


def load_grid(xlsx_path):
    with zipfile.ZipFile(xlsx_path) as z:
        shared = [
            "".join(t.text or "" for t in si.iter(NS + "t"))
            for si in ET.fromstring(z.read("xl/sharedStrings.xml"))
        ]
        grid = {}
        for row in ET.fromstring(z.read(SHEET)).find(NS + "sheetData"):
            for c in row:
                v = c.find(NS + "v")
                if v is None:
                    continue
                val = shared[int(v.text)] if c.get("t") == "s" else v.text
                r = int(re.search(r"\d+", c.get("r")).group())
                grid[(r, colnum(c.get("r")))] = val
    return grid


def main():
    grid = load_grid(XLSX)
    rows = []
    for i in range(N):
        t_min = float(grid[(MEAN_ROW0 + i, T_COL)])
        mean = float(grid[(MEAN_ROW0 + i, WT_COL)])
        sd = float(grid[(SD_ROW0 + i, WT_COL)])
        rows.append((round(t_min * 60.0), mean, sd))
    with open(OUT, "w") as fh:
        fh.write("# time pT69pT71 pT69pT71_SD\n")
        for t_s, mean, sd in rows:
            fh.write(f"{t_s} {mean:.4f} {sd:.4f}\n")
    print(f"wrote {OUT} with {len(rows)} points (t = 0..{rows[-1][0]} s)")


if __name__ == "__main__":
    main()
