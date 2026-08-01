#!/usr/bin/env python3
"""Regenerate every ``reference/rohrs2018_*_digitized.csv`` from the source PDF.

Rohrs JA, Zheng D, Graham NA, Wang P, Finley SD (2018), Biophys J 115:1116-1129.
The article and its Supporting Material ship as one PDF, ``mmc8.pdf``, in the
uncommitted ``dev/papers/rohrs2018/`` tree.

Three figures are digitized:

``Fig. S5`` (PDF page 20, **vector**)
    Nine panels of CD3-zeta ITAM site phosphorylation: measured percent
    phosphorylation (filled circles + SD error bars) and the competitive
    inhibition model fit (solid curves), for wild-type CD3-zeta on 10%, 0% and
    45% POPS liposomes and for the six single tyrosine-to-phenylalanine ITAM
    mutants.  Because the panels are vector art, every marker, error bar and
    curve vertex is recovered exactly; nothing is read off a raster.

    Axis calibration uses the *tick marks*, not the plot frame.  Here the two
    happen to coincide (MATLAB draws the box at the axis limits), but the tick
    positions are what the calibration is anchored to.  The x axis is log10
    time in minutes over [-1, 3] and the y axis is percent phosphorylation over
    [0, 100].

    Data-point positions come from the error bars rather than from the marker
    glyphs: the marker paths in this PDF carry a small (~0.15 pt, ~2% in time)
    systematic offset, whereas the error-bar segments land exactly on the axis
    limit for the t = 0.1 min point.  Each error bar is drawn as two vertical
    segments that meet at the mean and two horizontal caps at mean +/- SD, so
    both the mean and the SD are recovered exactly.  Points with no error bar
    (mutated sites, pinned at zero) fall back to the marker centre.

``Fig. 3`` (PDF page 8, **raster**)
    Model-predicted half-maximal time and Hill coefficient for each of the four
    candidate mechanisms (sequential order, random order, phosphate priming,
    competitive inhibition), plus the black dots that mark the corresponding
    values fitted to the experimental data (Fig. 1, D and E).  Bars are read by
    colour segmentation; each panel is calibrated on its own y-axis tick marks.
    The same six data dots are repeated in all four panels, so the spread across
    panels is an empirical estimate of the raster digitization error.

``Fig. 5`` (PDF page 10, **raster**)
    Phosphatase model: normalized percent phosphorylation versus
    log10(LCK/phosphatase) for wild-type and ITAM-mutant CD3-zeta (panel A),
    with the normalized log(EC50) and Hill coefficient of each response
    (panels B and C).

Usage (from the repository root, with the paper folder present)::

    python models/car_cd3zeta_phosphorylation_rohrs2018/digitize_rohrs2018.py
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

try:
    import fitz  # PyMuPDF
except ImportError:  # pragma: no cover
    sys.exit("PyMuPDF is required: pip install pymupdf")
try:
    from PIL import Image
except ImportError:  # pragma: no cover
    sys.exit("Pillow is required: pip install pillow")

HERE = Path(__file__).resolve().parent
OUT = HERE / "reference"
PDF = HERE.parents[1] / "dev" / "papers" / "rohrs2018" / "mmc8.pdf"

SITES = ["A1", "A2", "B1", "B2", "C1", "C2"]

# ---------------------------------------------------------------------------
# Fig. S5 -- vector
# ---------------------------------------------------------------------------

FRAME_GREY = (0.14900435507297516,) * 3
# MATLAB default colour order, in the order the Fig. S5 legend lists the sites.
S5_COLOR = {
    "A1": (0.0, 0.44699999690055847, 0.7409999966621399),
    "A2": (0.8510000109672546, 0.32499998807907104, 0.09799999743700027),
    "B1": (0.9290000200271606, 0.6940000057220459, 0.125),
    "B2": (0.49399998784065247, 0.18400000035762787, 0.5569999814033508),
    "C1": (0.46700000762939453, 0.675000011920929, 0.18799999356269836),
    "C2": (0.3019999861717224, 0.7450000047683716, 0.9330000281333923),
}
# Panel letter -> (file slug, human-readable condition).  Fig. S5 skips "G".
S5_PANELS = [
    ("A", "wildtype_10pct_pops_rep2", "wild-type CD3z, 10% POPS (biological replicate 2)"),
    ("B", "wildtype_0pct_pops", "wild-type CD3z, 0% POPS"),
    ("C", "wildtype_45pct_pops", "wild-type CD3z, 45% POPS"),
    ("D", "mutant_A1", "CD3z A1 tyrosine-to-phenylalanine mutant, 10% POPS"),
    ("E", "mutant_B1", "CD3z B1 tyrosine-to-phenylalanine mutant, 10% POPS"),
    ("F", "mutant_C1", "CD3z C1 tyrosine-to-phenylalanine mutant, 10% POPS"),
    ("H", "mutant_A2", "CD3z A2 tyrosine-to-phenylalanine mutant, 10% POPS"),
    ("I", "mutant_B2", "CD3z B2 tyrosine-to-phenylalanine mutant, 10% POPS"),
    ("J", "mutant_C2", "CD3z C2 tyrosine-to-phenylalanine mutant, 10% POPS"),
]


def _s5_panel_boxes(drawings):
    """Reassemble the nine plot frames from their four stroked edges."""
    edges = [d for d in drawings
             if d.get("color") == FRAME_GREY
             and (d["rect"].width > 50 or d["rect"].height > 50)]
    vertical = sorted({(round(d["rect"].x0, 2), round(d["rect"].y0, 2),
                        round(d["rect"].y1, 2))
                       for d in edges if d["rect"].width < 1})
    boxes = []
    for i, left in enumerate(vertical):
        for right in vertical[i + 1:]:
            if abs(right[1] - left[1]) > 1 or abs(right[2] - left[2]) > 1:
                continue
            if 60 < right[0] - left[0] < 75:      # one panel wide
                boxes.append((left[0], right[0], left[1], left[2]))
    boxes.sort(key=lambda b: (round(b[2]), b[0]))   # row-major
    return boxes


def _s5_ticks(drawings, box, axis):
    """Major + minor tick positions along one edge of a panel."""
    x0, x1, y0, y1 = box
    found = []
    for d in drawings:
        if d.get("color") != FRAME_GREY or d["type"] != "s":
            continue
        r = d["rect"]
        if axis == "x":
            if r.width > 0.5 or not (y1 - 1 < r.y0 < y1 + 1):
                continue
            if not 0.3 < r.height < 1.0:
                continue
            if x0 - 1 <= r.x0 <= x1 + 1:
                found.append(round(r.x0, 3))
        else:
            if r.height > 0.5 or not (x0 - 1 < r.x0 < x0 + 1):
                continue
            if not 0.3 < r.width < 1.0:
                continue
            if y0 - 1 <= r.y0 <= y1 + 1:
                found.append(round(r.y0, 3))
    return sorted(set(found))


def _s5_major(ticks, n):
    """Keep the n evenly spaced major ticks out of a log-minor tick ladder."""
    span = (ticks[-1] - ticks[0]) / (n - 1)
    major = [ticks[0]]
    for want in range(1, n):
        target = ticks[0] + want * span
        major.append(min(ticks, key=lambda t: abs(t - target)))
    return major


def digitize_fig_s5(doc):
    page = doc[19]                       # 0-based; "Figure S5" page
    drawings = page.get_drawings()
    boxes = _s5_panel_boxes(drawings)
    if len(boxes) != len(S5_PANELS):
        raise RuntimeError(f"expected 9 Fig. S5 panels, found {len(boxes)}")

    for (letter, slug, caption), box in zip(S5_PANELS, boxes):
        x0, x1, y0, y1 = box
        xt = _s5_major(_s5_ticks(drawings, box, "x"), 5)   # log10 t = -1..3
        yt = _s5_major(_s5_ticks(drawings, box, "y"), 5)   # 100..0 %
        px_per_decade = (xt[-1] - xt[0]) / 4.0
        px_per_pct = (yt[-1] - yt[0]) / 100.0
        to_logt = lambda px: -1.0 + (px - xt[0]) / px_per_decade  # noqa: E731
        to_pct = lambda py: (yt[-1] - py) / px_per_pct            # noqa: E731

        def in_panel(rect, pad=3.0):
            return (rect.x0 >= x0 - pad and rect.x1 <= x1 + pad
                    and rect.y0 >= y0 - pad and rect.y1 <= y1 + pad)

        markers, errbars, curves = {}, {}, {}
        for site in SITES:
            colour = S5_COLOR[site]
            circles, verticals, curve = [], [], None
            for d in drawings:
                if not in_panel(d["rect"]):
                    continue
                if d.get("color") != colour and d.get("fill") != colour:
                    continue
                items = d["items"]
                if len(items) == 5 and d["type"] in ("f", "fs"):
                    quad = [seg[3] for seg in items if seg[0] == "c"]
                    if len(quad) == 4:
                        circles.append((
                            (min(p.x for p in quad) + max(p.x for p in quad)) / 2,
                            (min(p.y for p in quad) + max(p.y for p in quad)) / 2))
                elif len(items) > 40:
                    curve = [seg[1] for seg in items if seg[0] == "l"]
                    curve.append(items[-1][2])
                elif len(items) == 2 and d["type"] == "s" and d["rect"].width < 0.05:
                    verticals.append(d["rect"])
            by_x = defaultdict(list)
            for rect in verticals:
                by_x[round(rect.x0, 2)] += [rect.y0, rect.y1]
            bars = {}
            for xk, ys in by_x.items():
                ys = sorted(ys)
                if len(ys) == 4:                       # two segments meeting at the mean
                    bars[xk] = ((ys[1] + ys[2]) / 2, (ys[3] - ys[0]) / 2)
            markers[site] = sorted(circles)
            errbars[site] = bars
            curves[site] = curve

        # Time grid: error-bar x positions are exact, and every site in a panel
        # shares the same sampling times, so pool them across sites.
        grid = sorted({xk for site in SITES for xk in errbars[site]})
        if not grid:
            grid = sorted({round(mx, 1) for mx in
                           (m[0] for site in SITES for m in markers[site])})

        rows = []
        for xk in grid:
            row = {"time_min": 10.0 ** to_logt(xk)}
            for site in SITES:
                hit = [k for k in errbars[site] if abs(k - xk) < 0.5]
                if hit:
                    mean, half = errbars[site][hit[0]]
                    row[site] = to_pct(mean)
                    row[site + "_SD"] = half / px_per_pct
                else:
                    near = [m for m in markers[site] if abs(m[0] - xk) < 1.2]
                    row[site] = to_pct(near[0][1]) if near else float("nan")
                    row[site + "_SD"] = float("nan")
            rows.append(row)

        header = ["time_min"] + [c for s in SITES for c in (s, s + "_SD")]
        _write(OUT / f"rohrs2018_figS5{letter.lower()}_{slug}_data_digitized.csv",
               header, rows, f"Fig. S5{letter}: measured % phosphorylation, {caption}")

        # Model curves share one x grid; resample onto the A1 vertices.
        base = [to_logt(p.x) for p in curves["A1"]]
        mrows = []
        for i, lt in enumerate(base):
            row = {"time_min": 10.0 ** lt}
            for site in SITES:
                pts = curves[site]
                row[site] = to_pct(pts[i].y) if pts and i < len(pts) else float("nan")
            mrows.append(row)
        _write(OUT / f"rohrs2018_figS5{letter.lower()}_{slug}_model_digitized.csv",
               ["time_min"] + SITES, mrows,
               f"Fig. S5{letter}: competitive inhibition model fit, {caption}")


# ---------------------------------------------------------------------------
# raster helpers, shared by Fig. 3 and Fig. 5
# ---------------------------------------------------------------------------

def _raster(doc, page_index, name):
    page = doc[page_index]
    images = page.get_images(full=True)
    if len(images) != 1:
        raise RuntimeError(f"expected one image on page {page_index} ({name})")
    pix = fitz.Pixmap(doc, images[0][0])
    if pix.n > 4:
        pix = fitz.Pixmap(fitz.csRGB, pix)
    return np.frombuffer(pix.samples, dtype=np.uint8).reshape(
        pix.height, pix.width, pix.n)[:, :, :3].astype(int)


def _longest_run(mask_1d):
    best = (0, 0, 0)
    i = 0
    while i < len(mask_1d):
        if mask_1d[i]:
            j = i
            while j < len(mask_1d) and mask_1d[j]:
                j += 1
            if j - i > best[0]:
                best = (j - i, i, j - 1)
            i = j
        else:
            i += 1
    return best


def _runs(mask_1d):
    out, i = [], 0
    while i < len(mask_1d):
        if mask_1d[i]:
            j = i
            while j < len(mask_1d) and mask_1d[j]:
                j += 1
            out.append((i + j - 1) / 2.0)
            i = j
        else:
            i += 1
    return out


def _tick_positions(sub, axis_col, n_ticks):
    """Evenly spaced tick centres in the narrow strip left of a y axis."""
    for lo, hi in ((6, 1), (8, 2), (10, 1), (12, 2)):
        ticks = _runs(sub[:, max(axis_col - lo, 0):axis_col - hi + 1].any(axis=1))
        if len(ticks) != n_ticks:
            continue
        gaps = np.diff(ticks)
        if gaps.size and np.ptp(gaps) < 0.25 * np.mean(gaps):
            return ticks
    raise RuntimeError(f"could not find {n_ticks} evenly spaced ticks")


def _bar_top(panel, colour, tol=35, min_height=6):
    """Topmost row of a solid bar of the given fill colour, or None."""
    mask = np.abs(panel - np.array(colour)).sum(axis=2) < tol
    keep = mask.sum(axis=0) > min_height
    if keep.sum() < 5:
        return None
    mask[:, ~keep] = False
    return int(np.where(mask.any(axis=1))[0].min())


def _round_blobs(mask, lo=60, hi=400, size=(8, 22)):
    """Connected components that are roughly circular -- the black data dots."""
    seen = np.zeros_like(mask)
    dots = []
    for r in range(mask.shape[0]):
        for c in range(mask.shape[1]):
            if not mask[r, c] or seen[r, c]:
                continue
            stack, pix = [(r, c)], []
            seen[r, c] = True
            while stack:
                a, b = stack.pop()
                pix.append((a, b))
                for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    na, nb = a + da, b + db
                    if (0 <= na < mask.shape[0] and 0 <= nb < mask.shape[1]
                            and mask[na, nb] and not seen[na, nb]):
                        seen[na, nb] = True
                        stack.append((na, nb))
            rows = [p[0] for p in pix]
            cols = [p[1] for p in pix]
            h, w = max(rows) - min(rows) + 1, max(cols) - min(cols) + 1
            if (lo < len(pix) < hi and size[0] < h < size[1] and size[0] < w < size[1]
                    and abs(h - w) <= 3 and len(pix) > 0.6 * h * w):
                dots.append(((min(cols) + max(cols)) / 2, (min(rows) + max(rows)) / 2))
    dots.sort()
    return dots


# ---------------------------------------------------------------------------
# Fig. 3 -- raster bar charts
# ---------------------------------------------------------------------------

# Column crops sit inside the vertical rules that separate the four mechanisms.
FIG3_COLUMNS = [(0, 518), (530, 1022), (1034, 1528), (1540, 2018)]
FIG3_MECHANISMS = ["sequential_order", "random_order",
                   "phosphate_priming", "competitive_inhibition"]
FIG3_BAR = {"A1": (57, 90, 171), "A2": (157, 76, 49), "B1": (200, 166, 69),
            "B2": (86, 39, 119), "C1": (123, 154, 60), "C2": (124, 170, 222)}
FIG3_BANDS = {
    "half_max_time_min": (620, 1010, 250.0, 6),
    "hill_coefficient": (1010, 1428, 3.0, 4),
}


def digitize_fig3(doc):
    img = _raster(doc, 7, "Fig. 3")
    dark = img.max(axis=2) < 90
    rows, dot_rows = [], defaultdict(list)

    for quantity, (y0, y1, vmax, n_ticks) in FIG3_BANDS.items():
        for mech, (cx0, cx1) in zip(FIG3_MECHANISMS, FIG3_COLUMNS):
            sub = dark[y0:y1, cx0:cx1]
            panel = img[y0:y1, cx0:cx1]
            runs = sorted(((_longest_run(sub[:, c]), c) for c in range(sub.shape[1])),
                          reverse=True)
            axis_col = min(c for (ln, _, _), c in runs if ln == runs[0][0][0])
            ticks = _tick_positions(sub, axis_col, n_ticks)
            top, bottom = ticks[0], ticks[-1]
            scale = vmax / (bottom - top)

            row = {"quantity": quantity, "source": mech}
            for site in SITES:
                t = _bar_top(panel, FIG3_BAR[site])
                row[site] = float("nan") if t is None else (bottom - t) * scale
            rows.append(row)

            blobs = panel.max(axis=2) < 60
            blobs[:, :axis_col + 6] = False
            blobs[int(bottom) - 4:, :] = False
            blobs[:int(top) - 4, :] = False
            dots = _round_blobs(blobs)
            if len(dots) == len(SITES):
                for site, (_, dy) in zip(SITES, dots):
                    dot_rows[quantity].append({site: (bottom - dy) * scale})

    # the same experimental dots are repeated in all four panels -- average them
    for quantity in FIG3_BANDS:
        pooled = defaultdict(list)
        for entry in dot_rows[quantity]:
            for site, value in entry.items():
                pooled[site].append(value)
        row = {"quantity": quantity, "source": "experiment_sigmoid_fit"}
        row.update({s: float(np.mean(pooled[s])) for s in SITES})
        rows.append(row)
        spread = {s: float(np.ptp(pooled[s])) for s in SITES}
        print(f"  Fig. 3 {quantity}: data-dot spread across the four panels = "
              + ", ".join(f"{s} {spread[s]:.2f}" for s in SITES))

    rows.sort(key=lambda r: (r["quantity"], r["source"]))
    _write(OUT / "rohrs2018_fig3_sigmoid_summary_digitized.csv",
           ["quantity", "source"] + SITES, rows,
           "Fig. 3, middle and bottom rows: half-maximal time (min) and Hill "
           "coefficient of each mechanism's predicted response, plus the "
           "sigmoid fit to the experimental data (Fig. 1, D and E)")


# ---------------------------------------------------------------------------
# Fig. 5 -- raster, phosphatase model
# ---------------------------------------------------------------------------

FIG5_VARIANTS = ["ABC", "xBC", "AxC", "ABx", "Axx", "xBx", "xxC"]
FIG5_COLOR = {"ABC": (57, 90, 171), "xBC": (170, 71, 32), "AxC": (207, 165, 45),
              "ABx": (91, 31, 121), "Axx": (117, 157, 45), "xBx": (112, 171, 229),
              "xxC": (115, 23, 38)}


def digitize_fig5(doc):
    img = _raster(doc, 9, "Fig. 5")
    dark = img.max(axis=2) < 90

    # --- panel A: boxed axes, x = log10(LCK/phosphatase) in [-2, 2], y = 0..100
    left_half = dark[:, :550]
    col_runs = sorted(((_longest_run(left_half[:, c]), c)
                       for c in range(left_half.shape[1])), reverse=True)
    axis_x = min(c for (ln, _, _), c in col_runs if ln == col_runs[0][0][0])
    (_, ytop, ybot) = col_runs[0][0]
    right_x = _longest_run(left_half[ybot, :])[2]
    # Panel A is the one place where the frame, not the ticks, sets the scale:
    # MATLAB draws these ticks *inside* the axes, where the plotted curves cover
    # them.  The frame is a safe substitute here -- in the vector Fig. S5 panels,
    # where both are measurable, the outermost ticks sit on the frame to better
    # than 0.01 pt, i.e. MATLAB boxes the axes at their limits.
    to_x = lambda px: -2.0 + 4.0 * (px - axis_x) / (right_x - axis_x)  # noqa: E731
    to_y = lambda py: 100.0 * (ybot - py) / (ybot - ytop)              # noqa: E731

    # Only the wild-type (ABC) trace is recovered.  The six mutant traces are
    # drawn as three overlapping dashed pairs (xBC/AxC/ABx and Axx/xBx/xxC) at
    # 974 px panel width, so a colour scan cannot separate the branches where
    # they cross; the mutant comparison uses the Fig. 5, B and C bar charts
    # instead, which report the same responses as log(EC50) and Hill.
    # Walk the trace right to left, following the run nearest the previous one.
    # The inset legend carries a blue sample line of its own in the upper left,
    # so a column can hold two blue runs; continuity picks the curve.
    rows, previous = [], None
    for px in range(right_x - 1, axis_x, -1):
        column = img[ytop:ybot + 1, px]
        d = np.abs(column - np.array(FIG5_COLOR["ABC"])).sum(axis=1)
        centres = _runs(d < 60)
        if not centres:
            continue
        centre = (centres[0] if previous is None
                  else min(centres, key=lambda c: abs(c - previous)))
        if previous is not None and abs(centre - previous) > 12:
            continue                      # a jump this large is not the curve
        previous = centre
        rows.append({"log10_lck_over_phosphatase": to_x(px),
                     "ABC": to_y(ytop + centre)})
    rows.reverse()
    _write(OUT / "rohrs2018_fig5a_wildtype_response_digitized.csv",
           ["log10_lck_over_phosphatase", "ABC"], rows,
           "Fig. 5A: normalized % phosphorylation of wild-type CD3z (ABC trace only) "
           "vs log10(LCK/phosphatase)")

    # --- panels B and C: bar charts stacked in the right-hand column
    bar_rows = []
    for quantity, (y0, y1, vmax, n_ticks) in {
            "log_ec50_normalized": (0, 196, 0.20, 5),
            "hill_coefficient_normalized": (196, 447, 1.2, 5)}.items():
        sub = dark[y0:y1, 560:]
        panel = img[y0:y1, 560:]
        runs = sorted(((_longest_run(sub[:, c]), c) for c in range(sub.shape[1])),
                      reverse=True)
        axis_col = min(c for (ln, _, _), c in runs if ln == runs[0][0][0])
        ticks = _tick_positions(sub, axis_col, n_ticks)
        top, bottom = ticks[0], ticks[-1]
        scale = vmax / (bottom - top)
        row = {"quantity": quantity}
        for variant in FIG5_VARIANTS:
            t = _bar_top(panel, FIG5_COLOR[variant])
            row[variant] = 0.0 if t is None else (bottom - t) * scale
        bar_rows.append(row)
    _write(OUT / "rohrs2018_fig5bc_phosphatase_summary_digitized.csv",
           ["quantity"] + FIG5_VARIANTS, bar_rows,
           "Fig. 5, B and C: normalized log(EC50) and Hill coefficient of the "
           "phosphatase-model response for each ITAM mutant")


# ---------------------------------------------------------------------------

def _write(path, header, rows, caption):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        fh.write(f"# {caption}\n")
        fh.write("# Digitized from Rohrs et al. (2018) Biophys J 115:1116-1129 by "
                 "digitize_rohrs2018.py\n")
        # lineterminator: csv defaults to CRLF, which the repo's mixed-line-ending
        # pre-commit hook rejects.
        writer = csv.DictWriter(fh, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if isinstance(row.get(k), float)
                                 and np.isnan(row[k]) else
                                 (f"{row[k]:.6g}" if isinstance(row.get(k), float)
                                  else row.get(k)))
                             for k in header})
    print(f"  wrote {path.relative_to(HERE)} ({len(rows)} rows)")


def main():
    if not PDF.exists():
        sys.exit(f"source PDF not found: {PDF}\n"
                 "dev/papers/ is not committed; obtain mmc8.pdf from the publisher.")
    doc = fitz.open(PDF)
    print("Fig. S5 (vector)")
    digitize_fig_s5(doc)
    print("Fig. 3 (raster)")
    digitize_fig3(doc)
    print("Fig. 5 (raster)")
    digitize_fig5(doc)


if __name__ == "__main__":
    main()
