#!/usr/bin/env python3
"""Digitize the plotted curves of Michalski and Loew (2012) from the source PDF.

The published figures are vector graphics, so the plotted data are recovered
exactly rather than estimated from a raster image: every data point is a small
filled or open marker (a circle, diamond or square built from four or eight path
segments), and the axis tick marks are short line segments in a separate path.
Marker centres are read from marker bounding boxes and mapped to data
coordinates with a linear (or log-linear) calibration fitted to the *tick mark*
positions, never to the plot frame, which in this paper is offset from the axis
limits by several points.

Usage:

    python digitize_michalski2012.py [PDF] [OUTDIR]

with defaults PDF = dev/papers/michalski2012/Michalski_2012_Phys._Biol._9_036010.pdf
(relative to the repository root; the dev/ tree is not committed) and
OUTDIR = the reference/ subdirectory next to this script. Re-running the script
regenerates every reference/michalski2012_fig*_digitized.csv byte for byte.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import fitz
import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_PDF = (
    HERE.parents[1] / "dev/papers/michalski2012/Michalski_2012_Phys._Biol._9_036010.pdf"
)

# Series colours used throughout the paper, keyed by holoenzyme size.
COLORS = {
    "dimer": (0.93, 0.12, 0.14),
    "trimer": (0.22, 0.32, 0.64),
    "tetramer": (0.41, 0.74, 0.27),
    "pentamer": (0.96, 0.5, 0.12),
    "hexamer": (0.14, 0.12, 0.13),
}
RED = (0.93, 0.12, 0.14)
BLACK = (0.14, 0.12, 0.13)


def rgb(c):
    return tuple(round(v, 2) for v in c) if c else None


def in_box(r, box):
    x0, y0, x1, y1 = box
    return r.x0 >= x0 - 0.5 and r.x1 <= x1 + 0.5 and r.y0 >= y0 - 0.5 and r.y1 <= y1 + 0.5


def paths(page, box):
    for g in page.get_drawings():
        if in_box(g["rect"], box):
            yield g


def markers(page, box, stroke=None, fill=None, nitems=None, kinds=None, maxsize=8.0):
    """Marker centres inside box, filtered by stroke/fill colour and shape."""
    out = []
    for g in paths(page, box):
        r = g["rect"]
        if r.width > maxsize or r.height > maxsize or r.width == 0:
            continue
        if stroke is not None and rgb(g["color"]) != stroke:
            continue
        if fill is not None and rgb(g["fill"]) != fill:
            continue
        if nitems is not None and len(g["items"]) != nitems:
            continue
        if kinds is not None and "".join(sorted({it[0] for it in g["items"]})) != kinds:
            continue
        out.append(((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2))
    return sorted(set(out))


def polyline(page, box, stroke, nitems):
    """Vertices of a polyline path (BioNetGen-style scatter drawn as segments)."""
    for g in paths(page, box):
        if rgb(g["color"]) == stroke and len(g["items"]) == nitems:
            if all(it[0] == "l" for it in g["items"]):
                pts = [(it[1].x, it[1].y) for it in g["items"]]
                return sorted(pts)
    return []


def ticks(page, box, nitems):
    """Return (xticks, yticks) pixel positions from the tick-mark path."""
    for g in paths(page, box):
        if len(g["items"]) == nitems and all(it[0] == "l" for it in g["items"]):
            xs, ys = [], []
            for it in g["items"]:
                a, b = it[1], it[2]
                if abs(a.x - b.x) < 0.01 and abs(a.y - b.y) > 0.01:
                    xs.append(a.x)
                elif abs(a.y - b.y) < 0.01 and abs(a.x - b.x) > 0.01:
                    ys.append(a.y)
            return sorted(xs), sorted(ys)
    raise LookupError(f"no tick path with {nitems} items in {box}")


def calib(pix, vals, log=False):
    """Least-squares linear map from pixel to data (or log10 data) coordinate."""
    pix = np.asarray(pix, float)
    tgt = np.log10(np.asarray(vals, float)) if log else np.asarray(vals, float)
    if len(pix) != len(tgt):
        raise ValueError(f"{len(pix)} ticks vs {len(tgt)} values")
    slope, icpt = np.polyfit(pix, tgt, 1)
    resid = np.max(np.abs(np.polyval([slope, icpt], pix) - tgt))
    return (lambda p: 10 ** (slope * np.asarray(p) + icpt)) if log else (
        lambda p: slope * np.asarray(p) + icpt
    ), resid


def on_grid(points, fx, tol=0.02):
    """Keep only markers whose abscissa lies on the 1-2-5 decade grid.

    Legend keys are drawn as ordinary markers inside the plot area; they sit at
    arbitrary abscissae and are removed by this filter.
    """
    if not points:
        return []
    vals = np.atleast_1d(fx([p[0] for p in points]))
    snapped = snap_log(vals)
    return [
        p for p, v, s in zip(points, vals, snapped) if abs(np.log10(v) - np.log10(s)) < tol
    ]


def collapse(points, tol=2.0):
    """Average markers that share an abscissa (overlapping copies of a series)."""
    out, cur = [], []
    for x, y in sorted(points):
        if cur and x - cur[-1][0] > tol:
            out.append((float(np.mean([p[0] for p in cur])), float(np.mean([p[1] for p in cur]))))
            cur = []
        cur.append((x, y))
    if cur:
        out.append((float(np.mean([p[0] for p in cur])), float(np.mean([p[1] for p in cur]))))
    return out


def write_csv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="\n") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        for row in rows:
            w.writerow([f"{v:.6g}" for v in row])
    print(f"wrote {path.name}  ({len(rows)} points)")


def snap_log(values, per_decade=(1.0, 2.0, 5.0)):
    """Snap digitized abscissae to the 1-2-5 decade grid the paper plots on."""
    out = []
    for v in values:
        dec = np.floor(np.log10(v) + 1e-9)
        best = min(
            (abs(np.log10(v) - np.log10(m * 10 ** d)), m * 10 ** d)
            for d in (dec - 1, dec, dec + 1)
            for m in per_decade
        )
        out.append(best[1])
    return out


def fig2(doc, out):
    """Fig. 2: FP after 1 s and 6 s versus CaM4, dimer through hexamer."""
    page = doc[6]
    box = (101.0, 70.0, 262.0, 178.0)
    xt, yt = ticks(page, box, 13)
    fx, rx = calib(xt, [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000], log=True)
    fy, ry = calib(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])
    print(f"  fig2 calibration residuals: x {rx:.2e} decades, y {ry:.2e} FP")

    # One colour per holoenzyme size. In this figure the 1 s reaction is drawn
    # with filled circles (16 points, 0.01 to 1000 uM, saturating at 1 - 1/e)
    # and the 6 s reaction with open circles (14 points, 0.001 to 20 uM,
    # saturating near 1); Fig. 7a uses the opposite fill convention, so each
    # assignment here was checked against the saturation value. Legend keys are
    # drawn with four path segments instead of eight and are excluded.
    for name, col in COLORS.items():
        for tag, kw in (
            ("1s", dict(stroke=None, fill=col)),
            ("6s", dict(stroke=col, fill=(1.0, 1.0, 1.0))),
        ):
            pts = markers(page, box, nitems=8, kinds="c", **kw)
            rows = sorted(zip(snap_log(fx([p[0] for p in pts])), fy([p[1] for p in pts])))
            write_csv(
                out / f"michalski2012_fig2_{tag}_{name}_digitized.csv", ["CaM4_uM", "FP"], rows
            )


def fig4(doc, out):
    """Fig. 4: FP after 6 s versus CaM4 in the limit r2 = 0."""
    page = doc[7]
    box = (106.0, 74.0, 261.0, 183.0)
    xt, yt = ticks(page, box, 8)
    fx, rx = calib(xt, [1e-2, 1e-1, 1, 10], log=True)
    fy, ry = calib(yt, [0.6, 0.4, 0.2, 0.0])
    print(f"  fig4 calibration residuals: x {rx:.2e} decades, y {ry:.2e} FP")
    for name, col in COLORS.items():
        pts = polyline(page, box, col, 11)
        rows = sorted(zip(fx([p[0] for p in pts]), fy([p[1] for p in pts])))
        write_csv(out / f"michalski2012_fig4_{name}_digitized.csv", ["CaM4_uM", "FP"], rows)


def fig5a(doc, out):
    """Fig. 5a: equilibrium FP versus CaM4 at 0.1, 1.0 and 10 uM PP1."""
    page = doc[7]
    box = (361.0, 85.0, 524.0, 190.0)
    xt, yt = ticks(page, box, 12)
    fx, rx = calib(xt, [1e-4, 1e-3, 1e-2, 1e-1, 1], log=True)
    fy, ry = calib(yt[:-1], [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])
    print(f"  fig5a calibration residuals: x {rx:.2e} decades, y {ry:.2e} FP")
    # Each PP1 level uses a different marker shape, drawn once per holoenzyme
    # size in the paper's colours. The shape-to-PP1 assignment is not stated in
    # the legend, so it is read off the curve order: more phosphatase means less
    # phosphorylation, so at the right-hand end of the axis the three curves are
    # ordered 0.1 > 1.0 > 10 uM PP1.
    shapes = [("c", 8), ("qu", 1), ("re", 1)]
    by_shape = {}
    for kind, n in shapes:
        series = {}
        for name, col in COLORS.items():
            pts = [p for p in markers(page, box, stroke=col, nitems=n, kinds=kind) if p[0] > 362]
            series[name] = on_grid(pts, fx)
        by_shape[(kind, n)] = series
    order = sorted(by_shape, key=lambda k: -float(fy(max(by_shape[k]["hexamer"])[1])))
    for pp1, key in zip(("0.1", "1.0", "10"), order):
        for name, pts in by_shape[key].items():
            rows = sorted(zip(snap_log(fx([p[0] for p in pts])), fy([p[1] for p in pts])))
            write_csv(
                out / f"michalski2012_fig5a_pp1_{pp1.replace('.', 'p')}uM_{name}_digitized.csv",
                ["CaM4_uM", "FP"],
                rows,
            )


def fig7(doc, out):
    """Fig. 7: ISHA versus hexamer, without (a) and with (b) phosphatase."""
    page = doc[8]

    # Panel (a): 1 s (open) and 6 s (closed) autophosphorylation, no PP1.
    box = (359.0, 70.0, 519.0, 177.0)
    xt, yt = ticks(page, box, 14)
    fx, rx = calib(xt, [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000], log=True)
    fy, ry = calib(yt[:-1], [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])
    print(f"  fig7a calibration residuals: x {rx:.2e} decades, y {ry:.2e} FP")
    # Open symbols are the 1 s reaction and closed symbols the 6 s reaction;
    # the hexamer is black and the ISHA red. The two red legend keys sit inside
    # the plot area at FP ~ 0.71 and are dropped by the abscissa cut.
    series = {
        "hexamer_1s": dict(stroke=BLACK, fill=(1.0, 1.0, 1.0), nitems=8, kinds="c"),
        "hexamer_6s": dict(stroke=None, fill=BLACK, nitems=8, kinds="c"),
        "isha_1s": dict(stroke=RED, nitems=1, kinds="qu"),
        "isha_6s": dict(stroke=None, fill=RED, nitems=4, kinds="l"),
    }
    for name, kw in series.items():
        pts = on_grid(markers(page, box, **kw), fx)
        rows = sorted(zip(snap_log(fx([p[0] for p in pts])), fy([p[1] for p in pts])))
        write_csv(out / f"michalski2012_fig7a_{name}_digitized.csv", ["CaM4_uM", "FP"], rows)

    # Panel (b): equilibrium with 0.1, 1.0 and 10 uM PP1, one marker shape per
    # phosphatase level; the levels are identified from the curve order as in
    # Fig. 5a.
    box = (359.0, 209.0, 518.0, 326.0)
    xt, yt = ticks(page, box, 11)
    fx, rx = calib(xt, [1e-4, 1e-3, 1e-2, 1e-1, 1], log=True)
    fy, ry = calib(yt, [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])
    print(f"  fig7b calibration residuals: x {rx:.2e} decades, y {ry:.2e} FP")
    by_shape = {}
    for kind, n in (("c", 8), ("qu", 1), ("re", 1)):
        series_b = {}
        for name, col in (("hexamer", BLACK), ("isha", RED)):
            series_b[name] = on_grid(markers(page, box, stroke=col, nitems=n, kinds=kind), fx)
        by_shape[(kind, n)] = series_b
    order = sorted(by_shape, key=lambda k: -float(fy(max(by_shape[k]["hexamer"])[1])))
    for pp1, key in zip(("0.1", "1.0", "10"), order):
        for name, pts in by_shape[key].items():
            rows = sorted(zip(snap_log(fx([p[0] for p in pts])), fy([p[1] for p in pts])))
            write_csv(
                out / f"michalski2012_fig7b_pp1_{pp1.replace('.', 'p')}uM_{name}_digitized.csv",
                ["CaM4_uM", "FP"],
                rows,
            )


def fig8(doc, out):
    """Fig. 8: six-state ISHA versus trimer as a function of calcium."""
    page = doc[9]
    panels = {
        "a": ((368.0, 83.0, 524.0, 186.0), 11, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]),
        "b": ((368.0, 185.0, 524.0, 293.0), 11, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]),
    }
    for panel, (box, nt, yvals) in panels.items():
        xt, yt = ticks(page, box, nt)
        fx, rx = calib(xt, [1, 2, 5, 10, 20], log=True)
        fy, ry = calib(yt, list(reversed(yvals)))
        print(f"  fig8{panel} calibration residuals: x {rx:.2e} decades, y {ry:.2e}")
        for name, col in (("trimer", BLACK), ("isha", RED)):
            pts = markers(page, box, stroke=col, nitems=8, kinds="c") + markers(
                page, box, stroke=col, nitems=1, kinds="qu"
            )
            rows = sorted(zip(fx([p[0] for p in pts]), fy([p[1] for p in pts])))
            write_csv(
                out / f"michalski2012_fig8{panel}_{name}_digitized.csv",
                ["Ca_uM", "Frac_pT286"],
                rows,
            )

    # Panels (c) and (d) share a wider calcium axis and a 0 to 0.01 ordinate.
    panels_cd = {
        "c": ((369.0, 331.0, 524.0, 385.0), 6),
        "d": ((369.0, 384.0, 524.0, 440.0), 6),
    }
    for panel, (box, nt) in panels_cd.items():
        try:
            xt, yt = ticks(page, box, nt)
        except LookupError:
            print(f"    note: no tick path found for fig8{panel}")
            continue
        fx, rx = calib(xt, [1, 10, 100], log=True)
        fy, ry = calib(yt, [0.010, 0.005, 0.0][: len(yt)])
        print(f"  fig8{panel} calibration residuals: x {rx:.2e} decades, y {ry:.2e}")
        for name, col in (("trimer", BLACK), ("isha", RED)):
            pts = markers(page, box, stroke=col, nitems=8, kinds="c") + markers(
                page, box, stroke=col, nitems=1, kinds="qu"
            )
            rows = sorted(zip(fx([p[0] for p in pts]), fy([p[1] for p in pts])))
            write_csv(
                out / f"michalski2012_fig8{panel}_{name}_digitized.csv",
                ["Ca_uM", "Frac_pT305"],
                rows,
            )


def main(argv):
    pdf = Path(argv[1]) if len(argv) > 1 else DEFAULT_PDF
    out = Path(argv[2]) if len(argv) > 2 else HERE / "reference"
    if not pdf.exists():
        sys.exit(f"source PDF not found: {pdf}")
    doc = fitz.open(pdf)
    for fn in (fig2, fig4, fig5a, fig7, fig8):
        print(f"{fn.__name__}: {fn.__doc__.splitlines()[0]}")
        fn(doc, out)


if __name__ == "__main__":
    main(sys.argv)
