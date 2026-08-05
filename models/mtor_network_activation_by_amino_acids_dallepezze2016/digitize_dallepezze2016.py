"""Recover the simulated curves of Fig. 2b, 2c and 2d from Dalle Pezze et al. (2016).

The paper publishes no table and no source-data file for its simulations: the fitted
trajectories exist only as the red curves of Fig. 2. This script digitizes them so the
curation can compare against numbers rather than against an impression.

Extraction route (digitization.md §1): **rasterized page, 600 dpi**, not the embedded
bitmaps. Page 5 does carry the three panel grids as JPEGs at 220 ppi, and
``embedded_raster`` would hand back the publisher's own pixels -- but each JPEG is cropped
flush to its y-axis, so the y tick marks, which point outward, are not in it. Only the
composed page has them. Calibrating instead on the top of the y-axis line is not a
substitute: the axis overshoots the largest tick by ten pixels in the Fig. 2b AMPK panel
and by two in the Insulin panel, so that route silently rescales several readout panels by
10-25% -- the Fig. 2b p70-S6K-pT389 plateau reads 28.2 that way against a true 35.7.
Rasterizing costs about a pixel of positional noise and buys a calibration anchored where
digitization.md §2 requires.

Calibration differs by axis, because the two are legible in different ways.

**y: the tick marks.** Each panel's three-to-five y ticks are found as notches on its
y-axis, matched to the tick label values from the page's text layer, and fitted by least
squares. The residual is reported per panel and is the dominant uncertainty these CSVs
carry -- the artwork places its ticks imprecisely, and the residual reaches 3-4% of full
scale on a few panels (12.7 units on the 0-400 Akt-pT308 panel of Fig. 2b).

**x: the printed tick labels.** The x tick marks point inward, so the red curve and the
blue error bars paint over them and several panels expose only one -- a per-panel
detection then locks onto a seven-window shifted by a whole tick. The labels 0, 20 ... 120
are drawn once per figure and centred on their ticks, and their horizontal centres are
exactly evenly spaced (9.250 pt in Fig. 2b, 9.167 in 2c, 9.333 in 2d, identical across all
six gaps to three decimals). That grid is placed on each panel by an offset measured from
a per-COLUMN median of the x-axis start; single panels disagree by several points where a
y-axis run breaks up, but the columns are aligned so the median is stable. Every panel is
then checked against whatever tick marks it does expose, and all three figures agree to
better than a point.

The y calibration is checked against landmarks whose answer is known exactly
(digitization.md §4): the two input panels of every figure plot Insulin and Amino_Acids at
values the model fixes. The value-1 landmarks sit in open panel space and are gated on. A
value-0 landmark lies on top of the x-axis line, so its visible ink is displaced upward by
the axis stroke and it measures the axis rather than the calibration; it is reported and
not gated.

Series separation is by hue: the simulated curve is red, the experimental mean +- s.e.m.
error bars are blue, and the two input traces are green. Each is classified by channel
dominance rather than distance to a nominal colour, because JPEG ringing spreads the red
of a thin stroke over (210,72,72)...(255,188,189). Where a blue error bar crosses the red
curve the red stays visible on both sides and the reader follows continuity.

Usage:
    python digitize_dallepezze2016.py [path/to/NatComm.pdf]

Re-running must leave `git diff` empty (digitization.md §8).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "skills/curate-model/scripts"))
from digitize import Axis, index_runs, page_raster, write_csv  # noqa: E402

HERE = Path(__file__).resolve().parent
DEFAULT_PDF = HERE.parents[1] / "dev" / "papers" / "DallePezze2016" / "NatComm.pdf"
DPI = 600
PAGE = 5

CITATION = "Dalle Pezze P et al. (2016) Nature Communications 7:13254. doi:10.1038/ncomms13254"

# Each panel grid, located by its page-point bounding box, with the readout of every cell
# in reading order. The boxes come from the placement rectangles of the three JPEGs
# (PyMuPDF `page.get_image_rects`), padded left and below to take in the tick marks and
# labels that sit outside the image crop. Panel identity is the title printed above each
# cell; re-check with `pdftotext -f 5 -l 5 -bbox`.
FIGURES = [
    dict(
        fig="2b", condition="aa_insulin",
        box=(60.0, 386.0, 300.0, 580.0),
        caption="Fig. 2b - simulated time courses, amino acids + insulin",
        panels=[
            "Insulin", "Amino_Acids", "IR_pY1146",
            "AMPK_pT172", "Akt_pT308", "Akt_pS473",
            "TSC_pS1387", "mTOR_pS2448", "mTOR_pS2481",
            "IRS_pS636", "PRAS40_pS183", "PRAS40_pT246",
            "S6K_pT389", "S6K_pT229",
        ],
    ),
    dict(
        fig="2c", condition="aa",
        box=(305.0, 386.0, 545.0, 580.0),
        caption="Fig. 2c - simulated time courses, amino acids only",
        panels=[
            "Insulin", "Amino_Acids", "Akt_pS473",
            "AMPK_pT172", "Akt_pT308", "mTOR_pS2481",
            "TSC_pS1387", "mTOR_pS2448", "PRAS40_pT246",
            "IRS_pS636", "PRAS40_pS183",
            "S6K_pT389", "S6K_pT229",
        ],
    ),
    dict(
        fig="2d", condition="aa_wortmannin",
        box=(58.0, 596.0, 215.0, 682.0),
        caption="Fig. 2d - simulated time courses, amino acids + wortmannin",
        panels=["Insulin", "Amino_Acids", "S6K_pT389", "S6K_pT229"],
    ),
]

INPUT_TRUTH = {
    "aa": {"Insulin": 0.0, "Amino_Acids": 1.0},
    "aa_insulin": {"Insulin": 1.0, "Amino_Acids": 1.0},
    "aa_wortmannin": {"Insulin": 0.0, "Amino_Acids": 1.0},
}

X_TICKS = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0])  # minutes


def masks(img: np.ndarray):
    r = img[:, :, 0].astype(int)
    g = img[:, :, 1].astype(int)
    b = img[:, :, 2].astype(int)
    return (
        (r > g + 40) & (r > b + 40),        # red: the simulated curve
        (g > r + 30) & (g > b + 30),        # green: the fixed inputs
        # structure: axes and ticks. Dark AND achromatic -- 175 rather than a true
        # black because the Fig. 2b row-4 axes render grey (max ~170) and are missed
        # entirely at 120, and achromatic because at that threshold the green input
        # lines (46,134,48) would otherwise be read as axes. Text is rejected by the
        # run-length filter in find_panels, not by the threshold.
        (img.max(axis=2) < 175) & (img.max(axis=2).astype(int)
                                   - img.min(axis=2).astype(int) < 40),
    )


def numeric_labels(page):
    """Every numeric word on the page as (value, x0, x1, box-centre y), in page points."""
    return [
        (float(w[4]), w[0], w[2], (w[1] + w[3]) / 2.0)
        for w in page.get_text("words")
        if re.fullmatch(r"\d+", w[4])
    ]


def find_panels(black: np.ndarray, box, s: float):
    """Every panel inside `box` as (x0, x1, xaxis_row, yaxis_col, ytop), raster pixels.

    A panel is an L: a long horizontal black run is the x-axis, and the tall vertical run
    rising from its left end is the y-axis.
    """
    X0, Y0, X1, Y1 = (int(v * s) for v in box)
    sub = black[Y0:Y1, X0:X1]
    _h, w = sub.shape
    # An axis row is one carrying a contiguous horizontal run at least 18% of the box
    # wide; that is longer than any glyph and shorter than the shortest panel axis.
    minrun = int(0.18 * w)
    axis_rows = []
    for y in range(sub.shape[0]):
        idx = np.where(sub[y])[0]
        if len(idx) and any(len(r) >= minrun for r in index_runs(idx, 3)):
            axis_rows.append(y)
    panels = []
    for band in index_runs(np.array(axis_rows), 3):
        yb = int(np.mean(band))
        strip = sub[band[0]: band[-1] + 1].any(axis=0)
        for run in index_runs(np.where(strip)[0], int(6 * s)):
            if len(run) < minrun:
                continue
            x0, x1 = int(run[0]), int(run[-1])
            # The y-axis is the tallest vertical run that reaches the x-axis, searched a
            # few points either side of the run's left end -- the x-axis line overhangs
            # the corner by a point or two, so x0 is not exactly the axis column. A run
            # must be at least 8 pt tall to qualify, which rejects the tick stubs and the
            # corner itself.
            best_col, best_top, best_h = None, None, 0
            lo = max(0, yb - int(40 * s))
            for x in range(x0, min(x0 + int(6 * s), sub.shape[1]) + 1):
                col = sub[lo:yb, x]
                if not col.any():
                    continue
                tail = index_runs(np.where(col)[0], 2)[-1]
                if tail[-1] < len(col) - 3:          # must reach the x-axis
                    continue
                if len(tail) > best_h:
                    best_col, best_top, best_h = x, lo + int(tail[0]), len(tail)
            if best_col is None or best_h < int(8 * s):
                # No tall vertical run reached the x-axis: the panel is still real (its
                # x-axis is unambiguous), the y-axis is just broken up in the render.
                # ytop only bounds the search windows below -- both axis calibrations come
                # from tick marks and printed labels -- so fall back to the nominal height
                # rather than dropping a panel.
                best_col = best_col if best_col is not None else x0
                best_top = max(0, yb - int(26 * s))
            panels.append((X0 + x0, X0 + x1, Y0 + yb, X0 + best_col, Y0 + best_top))
    return panels


def ticks_along(black: np.ndarray, *, axis: str, line: int, lo: int, hi: int, s: float,
                want: int | None = None):
    """Tick-mark centres on one axis line.

    Ticks may be drawn inward or outward and the figure does both, so both sides are
    searched. `want` is the number expected: the side that produces exactly that many is
    taken, which is what keeps a stray curve or axis stub on the other side from winning
    simply by contributing more detections. Without `want`, the busier side wins.
    """
    reach = max(3, int(2.4 * s))
    gap = max(1, int(0.5 * s))
    best: list[float] = []
    for sign in (-1, +1):
        a, b = line + sign * (gap + reach), line + sign * gap
        lo_i, hi_i = min(a, b), max(a, b)
        if axis == "x":
            band = black[lo_i:hi_i, lo:hi + 1].any(axis=0)
        else:
            band = black[lo:hi + 1, lo_i:hi_i].any(axis=1)
        centres = [float(np.mean(t)) + lo for t in index_runs(np.where(band)[0], 2)]
        if want is not None and len(centres) == want:
            return centres
        if len(centres) > len(best):
            best = centres
    return best


def x_tick_labels(labels, box, axis_row: float) -> list[float]:
    """Centres of the printed x tick labels 0, 20 ... 120, in page points.

    They are printed once per figure, in a single row immediately under its top-left
    readout panel, and are centred on their ticks. Their horizontal centres are exactly
    evenly spaced, so they give the tick pitch to better than a hundredth of a point.

    `axis_row` is the page y of that panel's x-axis: the search band is the eight points
    below it. Without that bound the y-axis "0" of every panel in the figure also matches
    the value filter, and the evenly-spaced-subset step then locks onto the wrong seven.
    """
    X0, _Y0, X1, _Y1 = box
    xs = sorted(
        lab[1] + (lab[2] - lab[1]) / 2.0
        for lab in labels
        if X0 < lab[1] < X1 and axis_row + 0.5 < lab[3] < axis_row + 9.0
        and lab[0] in (0, 20, 40, 60, 80, 100, 120)
    )
    xs = [x for i, x in enumerate(xs) if i == 0 or x - xs[i - 1] > 2.0]
    if len(xs) < 7:
        raise RuntimeError(f"found {len(xs)} x tick labels under y={axis_row:.1f}, want 7")
    return regular_subset(xs, 7)


def grid_from_pitch(candidates, pitch: float, want: int, span, tol: float):
    """`want` ticks at the known `pitch`, positioned by whichever offset the marks support.

    Each candidate is tried as each tick index; the offset agreeing with the most other
    candidates wins, and the whole grid must fit inside the panel's own axis line. Three
    agreeing marks are enough because the pitch is already known exactly -- which is what
    lets the reader survive a panel whose middle ticks are painted over by an error bar.
    """
    if not candidates:
        return None
    lo, hi = span
    best, score = None, 0
    for c in candidates:
        for j in range(want):
            t0 = c - j * pitch
            if t0 < lo - tol or t0 + (want - 1) * pitch > hi + tol:
                continue
            hits = sum(
                1 for i in range(want)
                if min(abs(t0 + i * pitch - q) for q in candidates) < tol
            )
            if hits > score:
                best, score = t0, hits
    if best is None or score < 3:
        return None
    return [best + i * pitch for i in range(want)]


def regular_subset(centres: list[float], want: int) -> list[float]:
    """The `want` consecutive candidates that are most evenly spaced.

    Tick detection on an axis line also picks up its two end caps, so a panel can offer
    nine candidates for seven ticks. Real ticks are equally spaced and the end caps are
    not, so the window with the smallest spacing variance is the tick set.
    """
    if len(centres) <= want:
        return centres
    c = sorted(centres)
    best, score = None, None
    for i in range(len(c) - want + 1):
        win = c[i: i + want]
        gaps = np.diff(win)
        v = float(np.std(gaps) / max(np.mean(gaps), 1e-9))
        if score is None or v < score:
            best, score = win, v
    return best


def digitize(pdf: Path):
    import fitz

    img = page_raster(pdf, PAGE, dpi=DPI)
    red, green, black = masks(img)
    s = DPI / 72.0
    page = fitz.open(pdf)[PAGE - 1]
    labels = numeric_labels(page)

    written = []
    for spec in FIGURES:
        panels = find_panels(black, spec["box"], s)
        if len(panels) != len(spec["panels"]):
            raise RuntimeError(
                f"{spec['fig']}: found {len(panels)} panels, expected {len(spec['panels'])}"
            )
        # The x grid is taken from the printed tick labels, not from the tick marks.
        # The labels 0, 20 ... 120 are drawn once per figure under its top-left readout
        # panel and are centred on their ticks; their horizontal centres are exactly
        # evenly spaced (9.250 pt in Fig. 2b, 9.167 in 2c, 9.333 in 2d, identical across
        # all six gaps to three decimals). The tick MARKS are not usable per panel --
        # they point inward, so the red curve and the blue error bars paint over them,
        # and several panels expose only one. Every panel in a figure is the same grid
        # translated horizontally, so the label grid is placed on each panel by the
        # offset from its own y-axis line, which is a long unambiguous stroke. Panels
        # that do expose enough marks are then checked against that placement.
        top_row = min(p[2] for p in panels) / s
        xlab = x_tick_labels(labels, spec["box"], top_row)
        pitch = float(np.mean(np.diff(xlab)))
        # Place the grid off a per-COLUMN median of the x-axis start, not off any single
        # panel. Individual panels disagree by several points -- a broken y-axis run or an
        # antialiased corner moves the detected start -- but the columns are aligned, so
        # the median over the four or five panels in a column is stable to a fraction of a
        # point.
        cols: dict[int, list[float]] = {}
        for q in panels:
            cols.setdefault(round(q[0] / s / 40), []).append(q[0] / s)
        col_x0 = {k: float(np.median(v)) for k, v in cols.items()}
        anchor_key = min(col_x0, key=lambda k: abs(col_x0[k] - xlab[0]))
        offset = xlab[0] - col_x0[anchor_key]

        rows, checks, xres, yres, xcheck, skipped = [], [], [], [], [], []
        for (x0, x1, yb, xax, ytop), name in zip(panels, spec["panels"], strict=True):
            base = col_x0[round(x0 / s / 40)] + offset
            xt = [(base + i * pitch) * s for i in range(len(X_TICKS))]
            xaxis = Axis.from_ticks(xt, X_TICKS)
            xres.append(xaxis.residual)
            # Check the placement wherever the panel exposes enough tick marks.
            cand = ticks_along(black, axis="x", line=yb, lo=xax + int(0.3 * s), hi=x1, s=s)
            agree = [min(abs(c - q) for q in xt) / s for c in cand]
            close = [a for a in agree if a < 1.0]
            if len(close) >= 3:
                xcheck.append(max(close))

            near = [
                lab for lab in labels
                if xax / s - 16 < lab[2] < xax / s + 0.6 and ytop / s - 2 < lab[3] < yb / s + 3
            ]
            yt = ticks_along(black, axis="y", line=xax, lo=ytop, hi=yb, s=s,
                             want=len(near) or None)
            px, vals = [], []
            for value, _lx0, _lx1, ly in sorted(near, key=lambda z: z[3]):
                if not yt:
                    break
                p = min(yt, key=lambda q: abs(q / s - ly))
                if abs(p / s - ly) < 3.0:
                    px.append(p)
                    vals.append(value)
            if len(px) < 2:
                if name in ("Insulin", "Amino_Acids"):
                    # An input panel is a calibration landmark, not data. Losing one
                    # costs a check, not a curve, so say so and carry on; losing a
                    # readout panel would mean emitting uncalibrated numbers.
                    skipped.append(name)
                    continue
                raise RuntimeError(
                    f"{spec['fig']} {name}: matched {len(px)} y ticks to labels, want >= 2"
                )
            yaxis = Axis.from_ticks(px, vals)
            yres.append(yaxis.residual)

            if name in ("Insulin", "Amino_Acids"):
                ys, _ = np.nonzero(green[ytop - 2: yb + 3, x0: x1 + 1])
                if len(ys) == 0:
                    raise RuntimeError(f"{spec['fig']} {name}: no green trace")
                checks.append((name, float(yaxis(float(np.median(ys)) + ytop - 2)),
                               INPUT_TRUTH[spec["condition"]][name]))
                continue

            top = max(0, ytop - 4)
            prev = None
            for x in range(x0, x1 + 1):
                col = np.where(red[top: yb + 3, x])[0]
                if len(col) == 0:
                    continue
                runs = index_runs(col, 2)
                run = (max(runs, key=len) if prev is None
                       else min(runs, key=lambda r: abs(float(np.mean(r)) + top - prev)))
                prev = float(np.mean(run)) + top
                t = float(xaxis(float(x)))
                if -0.5 <= t <= 120.5:
                    # Where the curve is steep, one column of ink spans many rows and its
                    # centre is not a value -- any ordinate in the run is equally
                    # consistent with the artwork. Emit the run's extent alongside the
                    # centre so a comparison can score against the interval instead of
                    # against a point that the figure never claimed.
                    hi = float(yaxis(float(run[0]) + top))
                    lo = float(yaxis(float(run[-1]) + top))
                    rows.append((name, float(np.clip(t, 0.0, 120.0)), float(yaxis(prev)),
                                 min(lo, hi), max(lo, hi)))

        for name, got, want in checks:
            # 10%: the landmark is a line one stroke thick on a panel whose y-axis
            # spans 25 pt in Figs. 2b/2c and 26 pt in 2d, so a single point of tick
            # placement error is already 4-8% of it. A miss larger than this is a
            # calibration fault rather than the artwork's own precision.
            if want != 0.0 and abs(got - want) / want > 0.10:
                raise RuntimeError(
                    f"{spec['fig']} {name}: calibration check failed, read {got:.4f} "
                    f"where the model fixes {want}"
                )
        print(
            f"  {spec['fig']}: {len(panels)} panels, x residual <= {max(xres):.3f} min, "
            f"y residual <= {max(yres):.4f} units, "
            + ", ".join(
                f"{n}={g:.3f} (exact {w}" + ("" if w else ", on-axis, not gated") + ")"
                for n, g, w in checks
            )
            + (f"; landmark not calibrated: {', '.join(skipped)}" if skipped else "")
            + (f"; x grid confirmed by marks to {max(xcheck):.2f} pt" if xcheck else "")
        )

        path = HERE / "reference" / (
            f"dallepezze2016_fig{spec['fig'][1:]}_{spec['condition']}_digitized.csv"
        )
        write_csv(
            path, ["readout", "time_min", "value_au", "ink_lo_au", "ink_hi_au"], rows,
            comments=[
                spec["caption"],
                CITATION,
                "Digitized by digitize_dallepezze2016.py from a 600 dpi render of page 5.",
                "Both axes calibrated on their tick marks; the per-panel residuals printed "
                "by the script bound the accuracy of these values.",
                "value_au is the centre of the red ink in that pixel column; ink_lo_au "
                "and ink_hi_au are its extent, which is wide wherever the curve is steep "
                "and the centre is therefore not a value the figure claims.",
                "Values are the paper's arbitrary relative units and are not comparable "
                "between readouts.",
            ],
            fmt=".6g",
        )
        written.append(path)
    return written


def main() -> int:
    pdf = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_PDF
    if not pdf.exists():
        print(
            f"{pdf} not found.\nThe source PDF is not committed (dev/ is gitignored). "
            "Place the Nature Communications article there, or pass its path.",
            file=sys.stderr,
        )
        return 1
    print(f"digitizing {pdf}")
    for p in digitize(pdf):
        print(f"  wrote {p.relative_to(HERE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
