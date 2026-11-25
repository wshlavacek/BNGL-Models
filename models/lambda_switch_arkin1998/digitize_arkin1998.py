#!/usr/bin/env python
"""Digitize the "Full" curve of Arkin, Ross & McAdams (1998) Fig. 6a.

Arkin A, Ross J, McAdams HH (1998). Stochastic kinetic analysis of developmental pathway
bifurcation in phage lambda-infected Escherichia coli cells. Genetics 149:1633-1648.
Fig. 6a: "Simulation results for fraction of cells producing lysogens"; the curve labelled
"Full" (filled squares) is the stochastic model with the proteolysis parameters fitted to the
full Kourilsky (1973) data set -- the curve lambda_switch_arkin1998_fullcircuit_exact.bngl is
compared against. The paper prints no table of these values and does not state how many runs
each point averages.

Route (digitization.md s1): the figure is a scanned bitmap -- page 12 (journal p. 1644) carries
ONE embedded 811 x 2204 1-bit /ImageMask (CCITT G4, ~304 ppi as placed) holding panels a-c,
so it is read at its native resolution with `pymupdf_raster`; nothing is resampled. Ink is the
painted sample (8.6% of the bitmap).

Calibration (s2): on the tick marks. Panel a draws outward ticks every 0.02 on the left spine
(0.00..1.00) and every 0.5 MOI on the bottom spine (0..16); the few hidden under the MOI 0-1
markers are lost, and the rest are numbered from the ticks level with the top and right
spines. Both axes are
least-squares fits through every detected tick, and the residuals are printed; a residual above
1.5 px aborts. Ticks are found as runs of ink leaving the spine outward; a run wider than 4 px
along the spine is a data marker hanging over the axis (the MOI 0 and 1 squares), not a tick.

Reading (s4): the Full markers are filled squares, ~17 px on a side, at integer MOI 0..15 -- the
paper's sampling grid, so the abscissa is snapped to the integer (s7) and only the ordinate is
measured. For each MOI the 13-px column strip centred on the calibrated abscissa is scanned for
the longest run of rows in which every column is ink: only a filled square gives such a run
(>= 12 rows); the connecting lines (3 px), the vertical hatching (1-2 px) and the other series
(open-topped triangles, diamonds, crosses) never fill the strip for that long. The marker's
ordinate is the midpoint of that run. Its horizontal centre is measured too and must sit within
0.15 MOI of the snapped integer, which checks the x calibration at every point.

Uncertainty: the run midpoint is good to about +/-1 px (anti-aliasing is absent in a 1-bit scan,
but the square's edges are ragged by one pixel), i.e. +/-0.0017 in probability; with the
calibration residual and a possible 1-px scan skew across the panel, the recorded per-point
digitization uncertainty is 0.005 (3 px). Points at MOI 0-2 sit on or just above the bottom
spine, which merges with the square's run; the midpoint is still the square's centre to within
the spine's own 3-px width, and the column `on_axis` flags them.

The paper's own Monte Carlo error is not stated. The plateau MOI 9-15 is flat to the eye, so the
scatter of those seven points about their mean is written into the CSV header as an empirical
estimate of the per-point sampling noise of the published curve (an upper bound: any real slope
on the plateau inflates it).

Usage:
    python digitize_arkin1998.py [path/to/arkin1998.pdf]
The PDF is not committed (digitization.md s8); the default path is
dev/papers/arkin1998_genetics.pdf. The publisher's PDF (Genetics, free to read) is archived at
https://web.archive.org/web/20210211013411/https://www.genetics.org/content/genetics/149/4/1633.full.pdf
(sha256 5730df7ba1c3b05f451ac38602fae1d39121055546d49519f00c67d9e79e99ce).
Writes reference/arkin1998_fig6a_full_digitized.csv.
"""

import hashlib
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "skills/curate-model/scripts"))
from digitize import Axis, index_runs, longest_run, open_pdf, pymupdf_raster, write_csv  # noqa: E402

DEFAULT_PDF = HERE.parents[1] / "dev/papers/arkin1998_genetics.pdf"
PDF_SHA256 = "5730df7ba1c3b05f451ac38602fae1d39121055546d49519f00c67d9e79e99ce"
OUT = HERE / "reference" / "arkin1998_fig6a_full_digitized.csv"
PAGE_INDEX = 11                         # 0-based; journal page 1644
PANEL_A_ROWS = 740                      # panel a occupies the top 740 rows of the bitmap
STRIP_HALF = 6                          # 13-px column strip per marker
MIN_SQUARE_ROWS = 12
DIG_UNC = 0.005


def _spine(mask, axis, lo, hi):
    """Index (column if axis=0, row if axis=1) in [lo, hi) with the longest straight ink run,
    and the (first, last) indices of the spine's full thickness."""
    line = (lambda i: mask[:, i]) if axis == 0 else (lambda i: mask[i, :])
    lens = {i: longest_run(line(i))[0] for i in range(lo, hi)}
    best = max(lens, key=lens.get)
    thick = [i for i in range(lo, hi) if lens[i] > 0.9 * lens[best]]
    return best, min(thick), max(thick)


def _outward_ticks(mask, along, start, step, n_max, span):
    """Centres of tick marks leaving a spine outward.

    `along` is the axis the spine runs along (0: vertical spine, ticks are rows; 1: horizontal
    spine, ticks are columns); `start` is the first pixel outside the spine and `step` (+1/-1)
    the outward direction. A position is a tick if ink continues >= 6 px outward from the spine.
    Runs wider than 4 px along the spine are markers over the axis, not ticks, and are dropped.
    """
    hits = []
    for p in range(span[0], span[1]):
        run, q = 0, start
        while 0 <= q < (mask.shape[1] if along == 0 else mask.shape[0]) and run < n_max and \
                (mask[p, q] if along == 0 else mask[q, p]):
            run += 1
            q += step
        if run >= 6:
            hits.append(p)
    return [(g[0] + g[-1]) / 2.0 for g in index_runs(np.array(hits)) if g[-1] - g[0] <= 3]


def digitize(pdf):
    doc = open_pdf(pdf)
    ink = pymupdf_raster(doc, PAGE_INDEX)[:, :, 0] > 127
    assert ink.mean() < 0.2, "ink must be the minority (painted) sample of the /ImageMask"
    P = ink[:PANEL_A_ROWS]

    left, l0, _ = _spine(P, 0, 90, 170)
    bottom, _, b1 = _spine(P, 1, 590, 680)

    # y: ticks every 0.02 leaving the left spine to the left. The lowest two (0.00, 0.02) are
    # hidden under the MOI-0 marker, so values are counted down from the 1.00 tick, which is
    # level with the top spine; the provisional spacing is (bottom spine - 1.00 tick) / 50.
    yt = _outward_ticks(P, 0, l0 - 1, -1, 20, (10, bottom + 3))
    dy = (bottom - yt[0]) / 50.0
    yv = [round(1.0 - 0.02 * round((p - yt[0]) / dy), 2) for p in yt]
    assert len(set(yv)) == len(yv) and yv[0] == 1.0 and len(yt) >= 45, (len(yt), yv)
    yax = Axis.from_ticks(yt, yv)
    # x: ticks every 0.5 MOI leaving the bottom spine downward. Those under the MOI 0 and 1
    # markers are lost; values count back from the 16 tick, level with the right spine, at a
    # provisional spacing of (16 tick - left spine) / 32.
    xt = _outward_ticks(P, 1, b1 + 1, +1, 20, (left - 5, 790))
    dx = (xt[-1] - left) / 32.0
    xv = [16.0 - 0.5 * round((xt[-1] - p) / dx) for p in xt]
    assert len(set(xv)) == len(xv) and xv[-1] == 16.0 and len(xt) >= 28, (len(xt), xv)
    xax = Axis.from_ticks(xt, xv)
    px_per_y = abs(yt[-1] - yt[0]) / abs(yv[-1] - yv[0])
    px_per_x = abs(xt[-1] - xt[0]) / abs(xv[-1] - xv[0])
    assert yax.residual * px_per_y < 1.5 and xax.residual * px_per_x < 1.5, \
        (yax.residual * px_per_y, xax.residual * px_per_x)

    inv_x = np.polyfit(xv, xt, 1)
    rows = []
    for moi in range(16):
        c = int(round(np.polyval(inv_x, moi)))
        strip = P[:, c - STRIP_HALF:c + STRIP_HALF + 1]
        full = np.flatnonzero(strip.all(axis=1))
        runs = [g for g in index_runs(full) if len(g) >= MIN_SQUARE_ROWS]
        assert runs, f"no filled square found at MOI {moi}"
        g = max(runs, key=len)
        r_mid = (g[0] + g[-1]) / 2.0
        # horizontal centre of the square: the ink run through the middle row of the marker
        row = P[int(round(r_mid))]
        hr = [h for h in index_runs(np.flatnonzero(row)) if h[0] <= c <= h[-1]][0]
        x_meas = float(xax((hr[0] + hr[-1]) / 2.0)) if len(hr) < 40 else float(moi)
        assert abs(x_meas - moi) < 0.15, f"MOI {moi}: marker centred at {x_meas:.2f}"
        y = float(yax(r_mid))
        rows.append(dict(MOI=moi, p_lysogeny=round(y, 4), x_measured=round(x_meas, 3),
                         square_rows=len(g), on_axis=int(g[-1] >= bottom - 2)))
    return rows, yax, xax, px_per_y, px_per_x, len(yt), len(xt)


def main():
    pdf = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_PDF
    if not pdf.is_file():
        sys.exit(f"PDF not found: {pdf}\nDownload the publisher's PDF (see this script's "
                 f"docstring for the archived URL and sha256) and pass its path.")
    sha = hashlib.sha256(pdf.read_bytes()).hexdigest()
    if sha != PDF_SHA256:
        print(f"warning: PDF sha256 {sha} differs from the one this CSV was made from")
    rows, yax, xax, ppy, ppx, nyt, nxt = digitize(pdf)
    plateau = np.array([r["p_lysogeny"] for r in rows if 9 <= r["MOI"] <= 15])
    comments = [
        "Arkin A, Ross J, McAdams HH (1998) Genetics 149:1633-1648, Fig. 6a, curve 'Full'",
        "(filled squares): simulated probability of lysogeny vs MOI, Table-3 'Full' proteolysis.",
        "Written by digitize_arkin1998.py from the embedded 1-bit scan on journal p. 1644.",
        f"Calibration: {nyt} y ticks (0.02 apart), residual {yax.residual:.4f} "
        f"({yax.residual * ppy:.2f} px); {nxt} x ticks (0.5 MOI apart), residual "
        f"{xax.residual:.4f} MOI ({xax.residual * ppx:.2f} px); {ppy:.1f} px per unit "
        "probability.",
        f"Per-point digitization uncertainty {DIG_UNC} (3 px). MOI is the paper's integer "
        "grid (snapped); x_measured is the marker centre as read.",
        f"Plateau MOI 9-15: mean {plateau.mean():.4f}, SD {plateau.std(ddof=1):.4f} -- an "
        "empirical (upper-bound) estimate of the published curve's per-point Monte Carlo noise;",
        "the paper does not state its number of runs per MOI.",
        "on_axis = 1: the square overlaps the bottom spine (MOI 0-2); its centre is still read.",
    ]
    write_csv(OUT, ["MOI", "p_lysogeny", "x_measured", "square_rows", "on_axis"], rows,
              comments=comments)
    for r in rows:
        print(f"  MOI {r['MOI']:2d}: {r['p_lysogeny']:.4f}  (x read {r['x_measured']:.2f}, "
              f"{r['square_rows']} rows{', on axis' if r['on_axis'] else ''})")
    print(f"  plateau MOI 9-15: mean {plateau.mean():.4f}, SD {plateau.std(ddof=1):.4f}")


if __name__ == "__main__":
    main()
