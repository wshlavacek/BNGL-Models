#!/usr/bin/env python3
"""Digitize the quantification panels of Korwek et al. (2023), Sci. Signal. 16:eabq1173.

The supplementary figures of Korwek et al. (2023) juxtapose quantified Western blot
data with the fitted model's trajectories, but no tabular version of either is
published, so both have to be recovered from the figure. This script does that
reproducibly and writes the CSVs in `reference/`.

Figures digitized
-----------------
fig. S4A  IFN-beta / IFN-beta+poly(I:C) protocol, 11 observables x 10 lanes
fig. S4B  the same experiment's OAS3 blot, 11 lanes
fig. S4C  poly(I:C)-only protocol, 11 observables x 5 lanes + 1 IFN-beta reference lane
fig. S4D  the same experiment's OAS3 blot
fig. S12A TNF-alpha protocol, nuclear NF-kB (RelA), 9 lanes
fig. S12B TNF-alpha protocol, p-IKK / IkBa / A20 for WT (fine and long time courses)
          and for A20 KO cells, 3 observables x 3 panels

Method
------
Every panel is a base-10 log plot whose y axis is ticked at 0.1, 1 and 10 and whose
series are colour-coded: the fitted model is drawn as orange open circles
(RGB ~ 255,75,9) joined by salmon lines, the blot reproduced in the figure as dark
grey filled dots, and the additional replicate as light grey filled dots. The panel
images are pulled out of the PDF at their native resolution (~230 ppi) with
`pdfimages`, so no rasterization step of our own degrades them.

Calibration is taken from the axis itself: the tick marks are short dark strokes
immediately to the left of each panel's vertical spine, and their pixel rows fix the
decade spacing (53.5 px for fig. S12B, 55.0 for fig. S12A, 49.0 for fig. S4). Lane
x-positions are read off the model's open circles, which sit exactly on the lane
centres and are the only orange objects in a panel. A series value is the intensity-
weighted centroid of the fattest contiguous run of that series' colour inside a
5-px-wide strip at the lane centre; open circles are read as the midpoint of the
ring's vertical extent instead.

Accuracy and its limits
-----------------------
Centroids land within about 1 px, i.e. ~4% in value. Two systematic limits are
recorded in the CSVs rather than hidden:

* values at or below ~0.11 sit on the axis floor. Where the published model
  prediction is effectively zero (p-IRF3 in the unstimulated lanes of fig. S4A,
  p-STAT1/2 in the untreated lane, A20 in the A20 KO panel of fig. S12B) the figure
  draws the marker at a small positive display value instead of off-scale, so those
  points carry `on_axis_floor = 1` and must not be read as model predictions.
* values above ~9.5 are clipped by the top of the axis and carry `above_axis = 1`.

Both flags are emitted so downstream users can filter; `verify_korwek2023.ipynb`
excludes flagged points from its quantitative comparison.

Usage
-----
    python digitize_korwek2023.py [--sm PATH_TO_scisignal.abq1173_sm.pdf]

The supplementary PDF is not redistributable and is not part of this repository; it
lives outside the tree (by default under `dev/papers/Korwek-2023/`). The CSVs this
script writes are derived data and *are* committed.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

try:
    from PIL import Image
except ImportError:  # pragma: no cover
    sys.exit("Pillow is required: pip install pillow")

HERE = Path(__file__).resolve().parent
DEFAULT_SM = (HERE.parents[1] / "dev" / "papers" / "Korwek-2023"
              / "scisignal.abq1173_sm.pdf")
OUT = HERE / "reference"

FLOOR, CEIL = 0.115, 9.5          # usable window inside the 0.1 .. 10 axis


# --------------------------------------------------------------------------- #
# image handling
# --------------------------------------------------------------------------- #
def extract_page(pdf: Path, page: int, workdir: Path) -> np.ndarray:
    """Pull the embedded figure bitmap for `page` out of the PDF, unscaled."""
    stem = workdir / f"p{page}"
    subprocess.run(["pdfimages", "-png", "-f", str(page), "-l", str(page),
                    str(pdf), str(stem)], check=True)
    pngs = sorted(workdir.glob(f"p{page}-*.png"))
    if not pngs:
        raise RuntimeError(f"no embedded image found on page {page} of {pdf}")
    return np.array(Image.open(pngs[0]).convert("RGB")).astype(int)


def series_masks(a: np.ndarray):
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    orange = (r - g > 100) & (g - b > 25)                       # model circles
    grey = (np.abs(r - g) < 14) & (np.abs(g - b) < 14)
    return orange, grey & (r < 145), grey & (r >= 145) & (r < 205)


def _runs(idx):
    if len(idx) == 0:
        return []
    out, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] <= 2:
            cur.append(i)
        else:
            out.append(cur)
            cur = [i]
    out.append(cur)
    return out


def lane_centres(orange, x0, x1, ylo, yhi, expected):
    xs = np.where(orange[ylo:yhi, x0:x1].any(axis=0))[0]
    groups = _runs(xs)
    if len(groups) != expected:
        raise RuntimeError(f"found {len(groups)} lanes in x[{x0}:{x1}], "
                           f"expected {expected}")
    return [int(round(float(np.mean(g)))) + x0 for g in groups]


def ring_centre(orange, x, ylo, yhi, hw=3):
    ys = np.where(orange[ylo:yhi, x - hw:x + hw + 1].any(axis=1))[0]
    return None if len(ys) == 0 else ylo + 0.5 * (ys.min() + ys.max())


def blob_centre(mask, x, ylo, yhi, hw=2, min_len=3):
    counts = mask[ylo:yhi, x - hw:x + hw + 1].sum(axis=1)
    idx = np.where(counts >= 3)[0]
    best, weight = None, 0
    for run in _runs(idx):
        if len(run) < min_len:
            continue
        w = counts[run].sum()
        if w > weight:
            best, weight = np.array(run), w
    if best is None:
        return None
    return ylo + float((best * counts[best]).sum() / counts[best].sum())


class LogPanel:
    """A log-scale panel whose top gridline is 10 and whose axis spans 2 decades."""

    def __init__(self, imgs, y_top, px_per_decade, pad=6):
        self.orange, self.dark, self.light = imgs
        self.y_top, self.ppd = float(y_top), float(px_per_decade)
        self.ylo = int(y_top) - pad
        self.yhi = int(y_top + 2 * px_per_decade) + pad

    def _v(self, y):
        return None if y is None else 10.0 * 10.0 ** (-(y - self.y_top) / self.ppd)

    def read(self, x):
        return (self._v(ring_centre(self.orange, x, self.ylo, self.yhi)),
                self._v(blob_centre(self.dark, x, self.ylo, self.yhi)),
                self._v(blob_centre(self.light, x, self.ylo, self.yhi)))


def flags(v):
    if v is None:
        return "", "", ""
    return (f"{v:.4f}", "1" if v <= FLOOR else "0", "1" if v >= CEIL else "0")


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote {path.relative_to(HERE)}  ({len(rows)} rows)")


# --------------------------------------------------------------------------- #
# figure layouts (pixel geometry measured from the embedded bitmaps)
# --------------------------------------------------------------------------- #
S4_OBS = ["p-STAT1", "STAT1", "p-STAT2", "STAT2", "PKR", "p-eIF2a",
          "RIG-I", "p-IRF3", "IkBa", "A20", "RNase L"]
S4_TOPS = [109, 233, 358, 482, 606, 731, 855, 979, 1104, 1228, 1352]
S4_PPD = 49.0
S4_OAS3_TOP, S4_OAS3_PPD = 1595.0, 49.5

# lane -> (IFN-beta exposure in hours or "", poly(I:C) exposure in hours or "")
S4A_LANES = [("", ""), (24, ""), (26, ""), (28, ""), (30, ""), (34, ""),
             (26, 2), (28, 4), (30, 6), (34, 10)]
S4B_LANES = [("", ""), (24, ""), (26, ""), (28, ""), (30, ""), (34, ""), (48, ""),
             (26, 2), (28, 4), (30, 6), (34, 10)]
S4C_LANES = [("", ""), ("", 2), ("", 4), ("", 6), ("", 10), (24, "")]

S12B_ROWS = [("p-IKK", 702.0), ("IkBa", 833.0), ("A20", 964.0)]
S12B_PPD = 53.5
S12A_TOP, S12A_PPD = 226.0, 55.0
S12A_TIMES = [0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6]
S12B_PANELS = [
    # key, genotype, x-search window, lane count, lane times in hours
    ("WT_fine", "WT", 83, 640, 10,
     [0, 5 / 60, 10 / 60, 0.25, 0.5, 0.75, 1, 1.5, 2, 3]),
    ("WT_long", "WT", 846, 1150, 8, [0, 0.25, 0.5, 1, 1.5, 2, 4, 6]),
    ("A20KO_long", "A20 KO", 1192, 1520, 8, [0, 0.25, 0.5, 1, 1.5, 2, 4, 6]),
]


def do_fig_s4(a):
    imgs = series_masks(a)
    orange = imgs[0]
    rows_main, rows_oas3 = [], []
    for panel, x0, x1, lanes in [("S4A", 132, 860, S4A_LANES),
                                 ("S4C", 872, 1290, S4C_LANES)]:
        xs = lane_centres(orange, x0, x1, 100, 1460, len(lanes))
        for obs, top in zip(S4_OBS, S4_TOPS):
            p = LogPanel(imgs, top, S4_PPD)
            for i, (x, (ifnb, pic)) in enumerate(zip(xs, lanes), start=1):
                mo, dk, lt = p.read(x)
                rows_main.append([panel, obs, i, ifnb, pic,
                                  *flags(mo), *flags(dk), *flags(lt)])
    for panel, x0, x1, lanes in [("S4B", 132, 860, S4B_LANES),
                                 ("S4D", 872, 1290, S4C_LANES)]:
        xs = lane_centres(orange, x0, x1, 1590, 1700, len(lanes))
        p = LogPanel(imgs, S4_OAS3_TOP, S4_OAS3_PPD)
        for i, (x, (ifnb, pic)) in enumerate(zip(xs, lanes), start=1):
            mo, dk, lt = p.read(x)
            rows_oas3.append([panel, "OAS3", i, ifnb, pic,
                              *flags(mo), *flags(dk), *flags(lt)])
    header = ["panel", "observable", "lane", "ifnb_h", "polyic_h",
              "model", "model_on_axis_floor", "model_above_axis",
              "blot_shown", "blot_shown_on_axis_floor", "blot_shown_above_axis",
              "blot_replicate", "blot_replicate_on_axis_floor",
              "blot_replicate_above_axis"]
    write_csv(OUT / "korwek2023_figS4_digitized.csv", header,
              rows_main + rows_oas3)


def do_fig_s12(a):
    imgs = series_masks(a)
    orange = imgs[0]
    rows = []
    xs = lane_centres(orange, 569, 1000, 224, 338, len(S12A_TIMES))
    p = LogPanel(imgs, S12A_TOP, S12A_PPD)
    for x, t in zip(xs, S12A_TIMES):
        mo, dk, lt = p.read(x)
        rows.append(["S12A", "NF-kB nuclear", "WT", f"{t:g}",
                     *flags(mo), *flags(dk), *flags(lt)])
    for key, geno, x0, x1, n, times in S12B_PANELS:
        xs = lane_centres(orange, x0, x1, 700, 812, n)
        for obs, top in S12B_ROWS:
            p = LogPanel(imgs, top, S12B_PPD)
            for x, t in zip(xs, times):
                mo, dk, lt = p.read(x)
                rows.append([f"S12B_{key}", obs, geno, f"{t:.5g}",
                             *flags(mo), *flags(dk), *flags(lt)])
    header = ["panel", "observable", "genotype", "time_h",
              "model", "model_on_axis_floor", "model_above_axis",
              "blot_shown", "blot_shown_on_axis_floor", "blot_shown_above_axis",
              "blot_replicate", "blot_replicate_on_axis_floor",
              "blot_replicate_above_axis"]
    write_csv(OUT / "korwek2023_figS12_digitized.csv", header, rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sm", type=Path, default=DEFAULT_SM,
                    help="path to scisignal.abq1173_sm.pdf")
    args = ap.parse_args()
    if not args.sm.exists():
        sys.exit(f"supplementary PDF not found: {args.sm}\n"
                 "It is a primary source and is deliberately not committed; "
                 "pass its location with --sm.")
    OUT.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        do_fig_s4(extract_page(args.sm, 5, tmp))
        do_fig_s12(extract_page(args.sm, 13, tmp))


if __name__ == "__main__":
    main()
