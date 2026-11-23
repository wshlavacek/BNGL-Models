#!/usr/bin/env python3
"""Digitize the nine model curves of Fig. 3 of Faeder et al. (2003).

Fig. 3 of Faeder JR et al., J Immunol 170:3769-3781 (2003), doi:10.4049/jimmunol.170.7.3769,
plots the model's predicted kinetics after RBL cells are given 1.6 nM IgE dimer, 0-70 min:

    a  receptors in aggregates                         % of receptors      1 solid curve
    b  receptors phosphorylated on gamma / on beta     % of receptors      2 solid curves
    c  Lyn free / weakly associated / SH2-bound        % of receptor number  solid, dashed, solid
    d  Syk phosphorylated by Lyn or Syk / by Syk (Auto) / by Lyn
                                                       % of receptor number  solid, dashed, dotted

Panels b and d also carry densitometry points with error bars (one experiment, rescaled by a
common factor; caption). Those are experimental data and are not digitized here: what this
script recovers is the published *model* output, which is what the curated file is checked
against.

Source and route
----------------
The source is the publisher's own raster of the figure, the "Download slide" JPEG of Fig. 3 on
the article page (https://academic.oup.com/jimmunol/article/170/7/3769/8028343; file
im0735521003.jpeg, 630 x 1800 px, sha256 below). No vector route was available: the article
PDF could not be retrieved. The figure is black ink on white, so series are separated by
geometry, not colour (digitization.md §3): vertical order and continuity for the solid curves,
connected components for the dashes and dots, and the edge reading described below where a
dashed curve overlaps a solid one.

Calibration (digitization.md §2)
--------------------------------
Every panel is boxed, with ticks drawn inward from the spines: minor ticks 14-15 px long and
major ticks 30-31 px. The ticks are detected in a band just inside the bottom spine and the
left spine, classed minor or major by length, indexed against the frame, and the pixel-to-data
map is a least-squares fit through the ticks alone (`Axis.from_ticks`); the frame only numbers
them. A tick whose class does not match its index (a curve crossing the tick band) is dropped.
The fit residuals are written into each CSV header: at most 0.06 px in x and 0.55 px in y. x
ticks are 10..60 min, y ticks every 10 % (a), 2.5 % (b), 0.5 % (c) and 1 % (d).

Reading the ink (digitization.md §4)
------------------------------------
A pixel is ink when its grey level is below 128. In one column a stroke is a run of ink; its
top and bottom edges are located to sub-pixel precision where the grey level crosses 128, and
the curve value is the midpoint. The solid strokes are 6-7 px thick, about 2 % of a panel's
height.

* Solid curves are followed column by column from the right edge of the panel, where they are
  well separated, toward t = 0, taking in each column the run nearest the value predicted from
  the columns already read. A run thicker than one stroke (allowing for slope) is merged with
  other ink and is skipped, as is a run too far from the prediction.
* In panels b and d the error bars and markers of the data points (t = 0, 1, 2, 4, 8, 16, 32,
  64 min; caps are about 32 px wide) occlude the curves. Every column within 17 px of a data
  time is excluded outright, and nothing is emitted there: no interpolation across them
  (digitization.md §3). In practice that leaves panels b and d without samples before
  t of about 10 min.
* Panel d's dashed (Auto) and dotted (Lyn) curves are read as connected components: a component
  whose size is that of one dash (about 7-10 x 33 px) or one dot (about 7 x 7 px) gives one
  sample at its ink-weighted centroid. A dash or dot that touches other ink (an error bar, the
  solid curve) is a larger component and is not used.
* Panel c's Free (solid) and Weakly associated (dashed) curves lie about one stroke width apart
  and overlap. In a column between two dashes only the solid stroke is present and its midpoint
  is read as usual (`read = centre`). In a column through a dash the two strokes form one band;
  the solid curve's value is then its top edge plus half a stroke and the dashed curve's is the
  band's bottom edge minus half a stroke (`read = edge`), with the stroke thickness taken as the
  median thickness of the solid-only columns of the same band. That rests on the two strokes
  being drawn equally thick, which the 6-7 px of every solid and dashed stroke in the figure
  supports. Panel c is read from t of about 8 min, where the band is shallow; before that the
  three curves of the panel are steep and cross.

Solid-curve columns are decimated to every fifth column (about 0.72 min), finer than the stroke
width without carrying every correlated raster column (digitization.md §7). Values are written
at %.6g; abscissae are not snapped.

Usage
-----
    python digitize_faeder2003.py [path/to/im0735521003.jpeg]

The default path is dev/papers/Faeder2003/im0735521003.jpeg. dev/ is not committed; the script
checks the image's sha256 and refuses any other rendition, so a re-run regenerates the committed
CSVs byte for byte (digitization.md §8).
"""

from __future__ import annotations

import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "skills/curate-model/scripts"))

from digitize import Axis, index_runs, write_csv  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
IMAGE = REPO / "dev/papers/Faeder2003/im0735521003.jpeg"
IMAGE_SHA256 = "7b7d6d2ac5e8283010cd6b7c01ccc8506cd8ab5dd6d993fee7b877bee61e32c1"
IMAGE_URL = "https://academic.oup.com/jimmunol/article/170/7/3769/8028343"
REFERENCE = HERE / "reference"
CITATION = ("Faeder JR, Hlavacek WS, Reischl I, Blinov ML, Metzger H, Redondo A, Wofsy C, "
            "Goldstein B (2003) J Immunol 170:3769-3781, doi:10.4049/jimmunol.170.7.3769")

INK = 128.0                 # grey level below which a pixel is ink
LEFT_SPINE = (134, 136)     # columns of the shared left spine (all four panels)
RIGHT_FRAME = 621           # column of the shared right frame
STRIDE = 5                  # keep every fifth column of a solid curve
DATA_TIMES = (0, 1, 2, 4, 8, 16, 32, 64)   # min; data points in panels b and d
DATA_HALF_WIDTH = 17        # px; error-bar caps are about 32 px wide


@dataclass(frozen=True)
class Panel:
    key: str
    top: int                # first row of the top frame
    bottom: int             # first row of the bottom spine
    y_step: float           # value between adjacent y ticks (minor and major alternate)
    y_max: float            # value at the top frame, used only to number the ticks
    quantity: str


PANELS = {
    "a": Panel("a", 73, 366, 10.0, 90.0, "receptors in aggregates, % of receptors"),
    "b": Panel("b", 512, 805, 2.5, 20.0, "receptors phosphorylated, % of receptors"),
    "c": Panel("c", 952, 1245, 0.5, 5.0, "Lyn per receptor, % of receptor number"),
    "d": Panel("d", 1391, 1684, 1.0, 10.0, "phosphorylated Syk per receptor, % of receptor number"),
}

MINOR, MAJOR = (13, 16), (29, 32)     # tick lengths in px, measured from the spine


def tick_class(length: int) -> str | None:
    if MINOR[0] <= length <= MINOR[1]:
        return "minor"
    if MAJOR[0] <= length <= MAJOR[1]:
        return "major"
    return None


# --------------------------------------------------------------------------------------------
# Calibration
# --------------------------------------------------------------------------------------------


def x_ticks(grey: np.ndarray, p: Panel) -> list[tuple[float, int]]:
    """(column, length) of each tick rising from the bottom spine of panel `p`."""
    ink = grey < INK
    band = ink[p.bottom - 8: p.bottom - 2, LEFT_SPINE[1] + 3: RIGHT_FRAME - 3]
    cols = np.flatnonzero(band.all(axis=0)) + LEFT_SPINE[1] + 3
    found = []
    for run in index_runs(cols, 1):
        c = int(round(run.mean()))
        length = 0
        while ink[p.bottom - 1 - length, c]:
            length += 1
        span = np.arange(run[0] - 2, run[-1] + 3)
        weight = (255.0 - grey[p.bottom - 8: p.bottom - 2, span]).sum(axis=0)
        found.append((float((span * weight).sum() / weight.sum()), length))
    return found


def y_ticks(grey: np.ndarray, p: Panel) -> list[tuple[float, int]]:
    """(row, length) of each tick extending right from the left spine of panel `p`."""
    ink = grey < INK
    c0 = LEFT_SPINE[1] + 3
    band = ink[p.top + 5: p.bottom - 2, c0: c0 + 6]
    rows = np.flatnonzero(band.all(axis=1)) + p.top + 5
    found = []
    for run in index_runs(rows, 1):
        r = int(round(run.mean()))
        length = 0
        while ink[r, LEFT_SPINE[1] + 1 + length]:
            length += 1
        span = np.arange(run[0] - 2, run[-1] + 3)
        weight = (255.0 - grey[span, c0: c0 + 6]).sum(axis=1)
        found.append((float((span * weight).sum() / weight.sum()), length))
    return found


def calibrate(grey: np.ndarray, p: Panel) -> tuple[Axis, Axis, list[float], list[float]]:
    """Tick-fitted x (min) and y maps for panel `p`, and the tick values used."""
    x0 = float(np.mean(LEFT_SPINE))
    xs, xv = [], []
    for col, length in x_ticks(grey, p):
        kind = tick_class(length)
        value = 10.0 * round((col - x0) / (RIGHT_FRAME - x0) * 7.0)
        want = "major" if value % 20.0 == 0 else "minor"
        if kind == want and 10.0 <= value <= 60.0:
            xs.append(col)
            xv.append(value)
    if xv != [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]:
        raise SystemExit(f"panel {p.key}: x ticks found at {xv}, expected 10..60 min")
    ys, yv = [], []
    y_zero = p.bottom + 1.0
    y_top = p.top + 1.0
    for row, length in y_ticks(grey, p):
        kind = tick_class(length)
        index = round((y_zero - row) / (y_zero - y_top) * p.y_max / p.y_step)
        value = index * p.y_step
        want = "major" if index % 2 == 0 else "minor"
        if kind == want and 0.0 < value < p.y_max:
            ys.append(row)
            yv.append(value)
    if len(yv) < 5 or len(set(yv)) != len(yv):
        raise SystemExit(f"panel {p.key}: y ticks {yv} are too few or ambiguous")
    return Axis.from_ticks(xs, xv), Axis.from_ticks(ys, yv), xv, yv


# --------------------------------------------------------------------------------------------
# Reading strokes
# --------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Stroke:
    top: float
    bottom: float

    @property
    def centre(self) -> float:
        return 0.5 * (self.top + self.bottom)

    @property
    def thickness(self) -> float:
        return self.bottom - self.top


def column_strokes(grey: np.ndarray, col: int, r0: int, r1: int) -> list[Stroke]:
    """Every ink run in rows r0..r1 of one column, with sub-pixel edges at the grey level INK."""
    g = grey[:, col]
    rows = np.flatnonzero(g[r0: r1 + 1] < INK) + r0
    out = []
    for run in index_runs(rows, 1):
        lo, hi = int(run[0]), int(run[-1])
        above, below = g[lo - 1], g[hi + 1]
        top = lo - 1 + (above - INK) / (above - g[lo]) if above > g[lo] else lo - 0.5
        bottom = hi + (INK - g[hi]) / (below - g[hi]) if below > g[hi] else hi + 0.5
        out.append(Stroke(float(top), float(bottom)))
    return out


def excluded(col: int, x_axis: Axis, windows: bool) -> bool:
    """True inside a data point's error-bar window (panels b and d)."""
    if not windows:
        return False
    t_per_px = x_axis(1.0) - x_axis(0.0)
    centres = [(t - x_axis(0.0)) / t_per_px for t in DATA_TIMES]
    return any(abs(col - c) <= DATA_HALF_WIDTH for c in centres)


def follow(grey: np.ndarray, p: Panel, x_axis: Axis, start_rank: int, stop_col: int, *,
           windows: bool, max_thick: float = 9.5, band: bool = False) -> dict[int, Stroke]:
    """Follow one solid curve (or the panel-c band) leftward from the panel's right edge.

    `start_rank` picks the curve in the first column by vertical order (0 = topmost). Returns
    the accepted stroke per column. A run is accepted only if it is within 4 px plus the
    predicted drift of the expected row, and no thicker than one stroke allowing for slope
    (`max_thick` at zero slope; the panel-c band is two strokes, so it gets twice that).
    """
    r0, r1 = p.top + 4, p.bottom - 3
    start = RIGHT_FRAME - 6
    first = column_strokes(grey, start, r0, r1)
    accepted = {start: first[start_rank]}
    history = [(start, first[start_rank].centre)]
    for col in range(start - 1, stop_col - 1, -1):
        if excluded(col, x_axis, windows):
            continue
        recent = history[-12:]
        if len(recent) >= 3:
            slope = np.polyfit([c for c, _ in recent], [r for _, r in recent], 1)[0]
        else:
            slope = 0.0
        last_col, last_row = history[-1]
        predicted = last_row + slope * (col - last_col)
        stretch = float(np.hypot(1.0, slope))
        limit = (2.0 * max_thick if band else max_thick) * stretch + 1.0
        reach = 4.0 + abs(slope) * 1.5 + (limit if band else 0.0)
        candidates = [s for s in column_strokes(grey, col, r0, r1)
                      if abs(s.centre - predicted) <= reach and s.thickness <= limit]
        if not candidates:
            continue
        best = min(candidates, key=lambda s: abs(s.centre - predicted))
        accepted[col] = best
        history.append((col, best.centre))
    return accepted


def components(grey: np.ndarray, p: Panel):
    """Connected ink components inside panel `p`: (row_slice, col_slice, mask, weights)."""
    from scipy import ndimage

    r0, r1 = p.top + 4, p.bottom - 3
    c0, c1 = LEFT_SPINE[1] + 5, RIGHT_FRAME - 3
    sub = grey[r0: r1 + 1, c0: c1 + 1]
    labels, _ = ndimage.label(sub < INK)
    for index, sl in enumerate(ndimage.find_objects(labels), start=1):
        mask = labels[sl] == index
        weight = (255.0 - sub[sl]) * mask
        yield (sl[0].start + r0, sl[0].stop + r0, sl[1].start + c0, sl[1].stop + c0, mask, weight)


def component_samples(grey: np.ndarray, p: Panel, height: tuple[int, int], width: tuple[int, int],
                      area: tuple[int, int]) -> list[tuple[float, float]]:
    """Ink-weighted centroid (col, row) of every component of the given size."""
    out = []
    for ra, rb, ca, cb, mask, weight in components(grey, p):
        h, w, n = rb - ra, cb - ca, int(mask.sum())
        if height[0] <= h <= height[1] and width[0] <= w <= width[1] and area[0] <= n <= area[1]:
            rr, cc = np.mgrid[ra:rb, ca:cb]
            total = weight.sum()
            out.append((float((cc * weight).sum() / total), float((rr * weight).sum() / total)))
    return sorted(out)


def on_stride(cols) -> list[int]:
    x0 = int(np.mean(LEFT_SPINE))
    return [c for c in sorted(cols) if (c - x0) % STRIDE == 0]


# --------------------------------------------------------------------------------------------
# Panels
# --------------------------------------------------------------------------------------------


def solid_rows(grey, p, x_axis, y_axis, series, rank, stop_col, windows):
    strokes = follow(grey, p, x_axis, rank, stop_col, windows=windows)
    return [(series, x_axis(c), y_axis(strokes[c].centre), "centre") for c in on_stride(strokes)]


def panel_c_band(grey, p, x_axis, y_axis, stop_col):
    """Free (solid) and Weakly associated (dashed) from the overlapping band of panel c."""
    strokes = follow(grey, p, x_axis, 1, stop_col, windows=False, band=True)
    single = {c: s for c, s in strokes.items() if s.thickness <= 8.5}
    double = {c: s for c, s in strokes.items() if s.thickness >= 10.0}
    if len(single) < 50 or len(double) < 50:
        raise SystemExit("panel c: could not split the Free / Weakly associated band")
    h = float(np.median([s.thickness for s in single.values()]))
    rows = []
    for c in on_stride(strokes):
        if c in single:
            rows.append(("free", x_axis(c), y_axis(single[c].centre), "centre"))
        elif c in double:
            s = double[c]
            rows.append(("free", x_axis(c), y_axis(s.top + 0.5 * h), "edge"))
            rows.append(("weakly_associated", x_axis(c), y_axis(s.bottom - 0.5 * h), "edge"))
    return rows, h


def main(argv: list[str]) -> None:
    image = Path(argv[1]) if len(argv) > 1 else IMAGE
    if not image.exists():
        raise SystemExit(
            f"source figure not found: {image}\n"
            "dev/papers/ is not committed. Download Fig. 3 ('Download slide', im0735521003.jpeg)\n"
            f"from {IMAGE_URL} and save it there.")
    digest = hashlib.sha256(image.read_bytes()).hexdigest()
    if digest != IMAGE_SHA256:
        raise SystemExit(f"{image}: sha256 {digest} is not the rendition this script reads "
                         f"({IMAGE_SHA256})")
    from PIL import Image

    grey = np.asarray(Image.open(image).convert("L"), dtype=float)
    if grey.shape != (1800, 630):
        raise SystemExit(f"{image}: unexpected size {grey.shape}")

    REFERENCE.mkdir(exist_ok=True)
    common = [
        f"Fig. 3 of {CITATION}.",
        "Model predictions after 1.6 nM IgE dimer (curves only; the densitometry points are "
        "not digitized).",
        f"Source: publisher raster im0735521003.jpeg (630x1800 px, sha256 {IMAGE_SHA256[:16]}...) "
        f"from {IMAGE_URL}.",
        "Written by models/early_fceri_signaling_faeder2003/digitize_faeder2003.py; read = centre "
        "(stroke midpoint), edge (overlapping strokes, see the script), dot or dash (component "
        "centroid).",
    ]
    outputs = {
        "a": ("faeder2003_fig3a_receptor_aggregation_digitized.csv",
              "Panel a: receptors in aggregates (% of receptors)."),
        "b": ("faeder2003_fig3b_receptor_phosphorylation_digitized.csv",
              "Panel b: receptors with phosphorylated gamma / beta ITAMs (% of receptors)."),
        "c": ("faeder2003_fig3c_lyn_redistribution_digitized.csv",
              "Panel c: Lyn free / weakly associated / SH2-bound (% of receptor number)."),
        "d": ("faeder2003_fig3d_syk_phosphorylation_digitized.csv",
              "Panel d: Syk phosphorylated by Lyn or Syk (total) / by Syk (auto) / by Lyn "
              "(% of receptor number)."),
    }
    x0 = int(np.mean(LEFT_SPINE))
    for key, p in PANELS.items():
        x_axis, y_axis, xv, yv = calibrate(grey, p)
        px_per_x = 1.0 / (x_axis(1.0) - x_axis(0.0))
        px_per_y = abs(1.0 / (y_axis(1.0) - y_axis(0.0)))
        note = [f"Calibration: x ticks {xv[0]:g}..{xv[-1]:g} min, residual {x_axis.residual:.3g} "
                f"min ({x_axis.residual * px_per_x:.2f} px); y ticks {len(yv)} of "
                f"{yv[0]:g}..{yv[-1]:g}, residual {y_axis.residual:.3g} "
                f"({y_axis.residual * px_per_y:.2f} px); 1 px = {1 / px_per_y:.4g} in y."]
        if key == "a":
            rows = solid_rows(grey, p, x_axis, y_axis, "aggregated", 0, x0 + 5, False)
        elif key == "b":
            rows = (solid_rows(grey, p, x_axis, y_axis, "gamma", 0, x0 + 5, True)
                    + solid_rows(grey, p, x_axis, y_axis, "beta", 1, x0 + 5, True))
        elif key == "c":
            stop = int(round(x0 + 8.0 * px_per_x))
            rows = solid_rows(grey, p, x_axis, y_axis, "sh2_bound", 0, stop, False)
            band_rows, h = panel_c_band(grey, p, x_axis, y_axis, stop)
            rows += band_rows
            note.append(f"Free / weakly associated read from t >= 8 min; stroke thickness for "
                        f"the edge reading {h:.2f} px.")
        else:
            rows = solid_rows(grey, p, x_axis, y_axis, "total", 0, x0 + 5, True)
            dashes = component_samples(grey, p, (6, 12), (24, 38), (140, 260))
            dots = component_samples(grey, p, (5, 10), (5, 10), (30, 60))
            rows += [("auto", x_axis(c), y_axis(r), "dash") for c, r in dashes]
            rows += [("lyn", x_axis(c), y_axis(r), "dot") for c, r in dots]
        name, caption = outputs[key]
        write_csv(REFERENCE / name, ["series", "t_min", "value_pct", "read"],
                  [(s, float(t), float(v), how) for s, t, v, how in rows],
                  comments=[caption, *common, *note], fmt=".6g")


if __name__ == "__main__":
    main(sys.argv)
