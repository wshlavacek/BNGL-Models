#!/usr/bin/env python3
"""Digitize Figure 15.4 of Mu et al. (2010), Handbook of Chemoinformatics Algorithms.

Figure 15.4 is the only quantitative result reported for the four-flux example, and
no tabular version of it is published, so the four plotted P01 trajectories have to
be recovered from the figure itself.

Method
------
The panel is line art, not a raster image: the axis box, the tick marks and all four
trajectories are stroked vector paths in the PDF. `pdftocairo -svg` (poppler) rewrites
that page as SVG without resampling anything, so the cubic Bezier control points of
each trajectory survive exactly as the typesetter emitted them. This script therefore
reads geometry, not pixels, and the only approximation left is the Bezier fit that the
figure's own producer used to represent the simulated polyline.

Calibration comes from the axis box of the same SVG, so any global offset or scale
that poppler introduces cancels: the box spans 0 to 20 in Time and 0 to 0.05 in P01.
The four interior x tick marks (5, 10, 15) and eight interior y tick marks (0.01 to
0.04, one set on each spine) are then used as an independent check and must land on
their nominal values to within `TICK_TOL_*`. This is the frame-then-check calibration
of `skills/curate-model/references/digitization.md` §2: the ticks remain the authority,
they are just used as an assertion rather than as a fit.

Series identification
---------------------
The four trajectories are distinguished by stroke dash pattern, which the caption of
Figure 15.4 maps onto simulation conditions:

    solid          pool size of M2 = 1   (v1 = v2 = 0.5, v3 = v4 = 0.4)
    dashed         pool size of M2 = 2   (same fluxes)
    dash-dot       pool size of M2 = 0.5 (same fluxes)
    dotted         all pool sizes = 1, fluxes v1 = 0.6, v2 = 0.4, v3 = 0.5, v4 = 0.3

Output
------
`reference/mu2010_fig15_4_digitized.csv` with columns

    series, condition, time, P01

`series` is the dash pattern name, `condition` the suffix of the matching simulation
protocol in `four_flux_network_isotopomer_labeling_mu2010.bngl`.

The SVG parsing, axis calibration and CSV writing come from
`skills/curate-model/scripts/digitize.py`; what is left here is what is particular to
this figure. Requires the poppler `pdftocairo` executable on PATH; no third-party
Python packages beyond numpy.

Usage
-----
    python digitize_mu2010.py [path/to/HandbookChemoinformatics10.pdf]
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "skills/curate-model/scripts"))

from digitize import (  # noqa: E402
    Axis,
    check_ticks,
    dedupe_by_x,
    split_tick_segments,
    svg_paths,
    write_csv,
)

# Page 414 of the chapter, i.e. the 16th page of the PDF, carries Figure 15.4.
FIGURE_PAGE = 16

# Axis limits of the panel, read off the printed tick labels.
X_MIN, X_MAX = 0.0, 20.0
Y_MIN, Y_MAX = 0.0, 0.05

# Interior tick values that the calibration is checked against.
X_TICKS = (5.0, 10.0, 15.0)
Y_TICKS = (0.01, 0.02, 0.03, 0.04)
TICK_TOL_X = 0.05  # time units
TICK_TOL_Y = 5e-5  # P01 units

# Dash patterns as emitted by poppler, keyed by the name used in the caption.
DASH_TO_SERIES = {
    "": "solid",
    "3 2": "dashed",
    "3 2 1 2": "dash-dot",
    "0.3 1 0.3 1": "dotted",
}

# Simulation protocol each plotted curve corresponds to.
SERIES_TO_CONDITION = {
    "solid": "ode",
    "dashed": "ode_m2high",
    "dash-dot": "ode_m2low",
    "dotted": "ode_altflux",
}

# Samples per cubic Bezier segment. The curves have 3-4 segments each, so this
# yields a few hundred points per trajectory, well beyond the resolution of the
# underlying simulation (41 output times).
SAMPLES_PER_SEGMENT = 200


def axis_paths(paths):
    """Return the undashed 0.5 pt paths, i.e. the axis box and the tick marks."""
    return [p for p in paths if not p.dash and abs(p.width - 0.5) < 1e-9]


def find_frame(paths) -> tuple[float, float, float, float]:
    """Return (x_left, x_right, y_bottom, y_top) of the axis box, in page units.

    The box and the tick marks share a bounding box, so they are told apart by
    segment count: the box is four segments (eight endpoints), the ticks fourteen.
    """
    cands = []
    for p in axis_paths(paths):
        x_min, x_max, y_min, y_max = p.bbox
        if len(p.points) == 8 and x_max - x_min > 100 and y_max - y_min > 100:
            cands.append((x_min, x_max, y_max, y_min))
    if len(cands) != 1:
        raise SystemExit(f"expected exactly one axis box, found {len(cands)}")
    return cands[0]


def check_calibration(paths, to_x: Axis, to_y: Axis) -> None:
    """Assert that the interior tick marks land on their nominal data values."""
    # The tick marks are the other undashed 0.5 pt path of the panel, drawn as a run
    # of two-point strokes: a vertical stroke is an x tick and a horizontal one a y tick.
    ticks = max(axis_paths(paths), key=lambda p: len(p.points))
    segments = [ticks.points[i:i + 2] for i in range(0, len(ticks.points), 2)]
    xs, ys = split_tick_segments(segments, tol=1e-6)
    check_ticks(to_x, xs, X_TICKS, TICK_TOL_X, label="x")
    check_ticks(to_y, ys, Y_TICKS, TICK_TOL_Y, label="y")
    print(f"calibration check: {len(X_TICKS)} x ticks and {len(Y_TICKS)} y ticks within tolerance")


def main(argv: list[str]) -> int:
    here = Path(__file__).resolve().parent
    default_pdf = (
        here.parents[1] / "dev" / "papers" / "Mu2010" / "HandbookChemoinformatics10.pdf"
    )
    pdf = Path(argv[1]) if len(argv) > 1 else default_pdf
    if not pdf.exists():
        raise SystemExit(f"source PDF not found: {pdf}")

    paths = svg_paths(pdf, FIGURE_PAGE, samples=SAMPLES_PER_SEGMENT)
    xl, xr, yb, yt = find_frame(paths)
    to_x = Axis.from_limits(xl, xr, X_MIN, X_MAX)
    to_y = Axis.from_limits(yb, yt, Y_MIN, Y_MAX)
    check_calibration(paths, to_x, to_y)

    rows = []
    for dash, series in DASH_TO_SERIES.items():
        # The trajectories are the 0.75 pt strokes; the solid one is the only
        # undashed 0.75 pt path (the axis box and ticks are drawn at 0.5 pt).
        cands = [p for p in paths if p.dash == dash and abs(p.width - 0.75) < 1e-9]
        if len(cands) != 1:
            raise SystemExit(f"expected exactly one {series} trajectory, found {len(cands)}")
        pts = sorted((to_x(px), to_y(py)) for px, py in cands[0].points)
        # Consecutive Bezier segments share a knot; drop the duplicated endpoint.
        dedup = dedupe_by_x(pts)
        for t, v in dedup:
            rows.append((series, SERIES_TO_CONDITION[series], f"{t:.6f}", f"{v:.8f}"))
        print(f"{series:9s} -> {SERIES_TO_CONDITION[series]:12s} "
              f"{len(dedup):4d} points, t in [{dedup[0][0]:.2f}, {dedup[-1][0]:.2f}], "
              f"P01 in [{min(v for _, v in dedup):.5f}, {max(v for _, v in dedup):.5f}]")

    # lineterminator: this file was first written with csv.writer's CRLF default and is
    # committed that way. It is consistently CRLF, so the mixed-line-ending hook accepts
    # it; preserved here so re-running does not rewrite every line. New work uses LF.
    write_csv(here / "reference" / "mu2010_fig15_4_digitized.csv",
              ["series", "condition", "time", "P01"], rows, lineterminator="\r\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
