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
their nominal values to within `TICK_TOL_*`.

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

Requires the poppler `pdftocairo` executable on PATH; no third-party Python packages.

Usage
-----
    python digitize_mu2010.py [path/to/HandbookChemoinformatics10.pdf]
"""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

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

NUM = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def page_to_svg(pdf: Path, page: int) -> str:
    """Render one page of `pdf` to SVG with poppler, preserving vector geometry."""
    if shutil.which("pdftocairo") is None:
        raise SystemExit("pdftocairo (poppler) is required but was not found on PATH")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "page.svg"
        subprocess.run(
            ["pdftocairo", "-svg", "-f", str(page), "-l", str(page), str(pdf), str(out)],
            check=True,
            capture_output=True,
        )
        return out.read_text()


def parse_matrix(transform: str | None) -> tuple[float, float, float, float, float, float]:
    """Return the affine `matrix(a b c d e f)` of an SVG transform, identity if absent."""
    if not transform:
        return (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
    m = re.search(r"matrix\(([^)]*)\)", transform)
    if m is None:
        raise ValueError(f"unsupported transform: {transform}")
    a, b, c, d, e, f = (float(v) for v in NUM.findall(m.group(1)))
    return (a, b, c, d, e, f)


def apply(mat: tuple[float, ...], x: float, y: float) -> tuple[float, float]:
    a, b, c, d, e, f = mat
    return (a * x + c * y + e, b * x + d * y + f)


def path_points(d_attr: str, mat: tuple[float, ...]) -> list[tuple[float, float]]:
    """Flatten an SVG path of M/L/C commands into page-space points.

    Only the absolute commands poppler emits for stroked line art are handled; a
    cubic segment is sampled at `SAMPLES_PER_SEGMENT` parameter values.
    """
    tokens = re.findall(r"[MLCZmlcz]|" + NUM.pattern, d_attr)
    pts: list[tuple[float, float]] = []
    cur: tuple[float, float] | None = None
    i = 0
    while i < len(tokens):
        op = tokens[i]
        if op in "Zz":
            i += 1
            continue
        if op not in "MLCmlc":
            raise ValueError(f"unexpected path token {op!r}")
        n = {"M": 2, "L": 2, "C": 6}[op.upper()]
        i += 1
        coords = [float(tokens[i + k]) for k in range(n)]
        i += n
        if op.upper() == "M":
            cur = (coords[0], coords[1])
            pts.append(apply(mat, *cur))
        elif op.upper() == "L":
            cur = (coords[0], coords[1])
            pts.append(apply(mat, *cur))
        else:
            assert cur is not None
            p0 = cur
            p1, p2, p3 = (coords[0], coords[1]), (coords[2], coords[3]), (coords[4], coords[5])
            for k in range(1, SAMPLES_PER_SEGMENT + 1):
                t = k / SAMPLES_PER_SEGMENT
                u = 1.0 - t
                x = u**3 * p0[0] + 3 * u**2 * t * p1[0] + 3 * u * t**2 * p2[0] + t**3 * p3[0]
                y = u**3 * p0[1] + 3 * u**2 * t * p1[1] + 3 * u * t**2 * p2[1] + t**3 * p3[1]
                pts.append(apply(mat, x, y))
            cur = p3
    return pts


def collect_paths(svg: str) -> list[dict]:
    """Return every stroked <path> of the page as page-space points plus its dash style."""
    root = ET.fromstring(svg)
    ns = "{http://www.w3.org/2000/svg}"
    out = []
    for el in root.iter(f"{ns}path"):
        d_attr = el.get("d")
        if not d_attr or el.get("stroke", "none") == "none":
            continue
        mat = parse_matrix(el.get("transform"))
        out.append(
            {
                "dash": (el.get("stroke-dasharray") or "").strip(),
                "width": float(el.get("stroke-width") or 0.0),
                "points": path_points(d_attr, mat),
            }
        )
    return out


def axis_paths(paths: list[dict]) -> list[dict]:
    """Return the undashed 0.5 pt paths, i.e. the axis box and the tick marks."""
    return [p for p in paths if not p["dash"] and abs(p["width"] - 0.5) < 1e-9]


def find_frame(paths: list[dict]) -> tuple[float, float, float, float]:
    """Return (x_left, x_right, y_bottom, y_top) of the axis box, in page units.

    The box and the tick marks share a bounding box, so they are told apart by
    segment count: the box is four segments (eight endpoints), the ticks fourteen.
    """
    cands = []
    for p in axis_paths(paths):
        xs = [q[0] for q in p["points"]]
        ys = [q[1] for q in p["points"]]
        if len(p["points"]) == 8 and max(xs) - min(xs) > 100 and max(ys) - min(ys) > 100:
            cands.append((min(xs), max(xs), max(ys), min(ys)))
    if len(cands) != 1:
        raise SystemExit(f"expected exactly one axis box, found {len(cands)}")
    return cands[0]


def check_ticks(paths: list[dict], frame: tuple[float, float, float, float]) -> None:
    """Assert that the interior tick marks land on their nominal data values."""
    xl, xr, yb, yt = frame
    to_x = lambda px: X_MIN + (px - xl) / (xr - xl) * (X_MAX - X_MIN)  # noqa: E731
    to_y = lambda py: Y_MIN + (py - yb) / (yt - yb) * (Y_MAX - Y_MIN)  # noqa: E731

    # The tick marks are the other undashed 0.5 pt path of the panel.
    ticks = max(axis_paths(paths), key=lambda p: len(p["points"]))
    seg = [ticks["points"][i : i + 2] for i in range(0, len(ticks["points"]), 2)]
    xs, ys = [], []
    for (x0, y0), (x1, y1) in seg:
        if abs(x0 - x1) < 1e-6:  # vertical stroke -> x-axis tick
            xs.append(to_x(x0))
        elif abs(y0 - y1) < 1e-6:  # horizontal stroke -> y-axis tick
            ys.append(to_y(y0))
    for want in X_TICKS:
        err = min(abs(v - want) for v in xs)
        if err > TICK_TOL_X:
            raise SystemExit(f"x tick {want} off by {err:.4f} > {TICK_TOL_X}")
    for want in Y_TICKS:
        err = min(abs(v - want) for v in ys)
        if err > TICK_TOL_Y:
            raise SystemExit(f"y tick {want} off by {err:.2e} > {TICK_TOL_Y:.0e}")
    print(f"calibration check: {len(X_TICKS)} x ticks and {len(Y_TICKS)} y ticks within tolerance")


def main(argv: list[str]) -> int:
    here = Path(__file__).resolve().parent
    default_pdf = (
        here.parents[1] / "dev" / "papers" / "Mu2010" / "HandbookChemoinformatics10.pdf"
    )
    pdf = Path(argv[1]) if len(argv) > 1 else default_pdf
    if not pdf.exists():
        raise SystemExit(f"source PDF not found: {pdf}")

    paths = collect_paths(page_to_svg(pdf, FIGURE_PAGE))
    frame = find_frame(paths)
    check_ticks(paths, frame)
    xl, xr, yb, yt = frame

    rows = []
    for dash, series in DASH_TO_SERIES.items():
        # The trajectories are the 0.75 pt strokes; the solid one is the only
        # undashed 0.75 pt path (the axis box and ticks are drawn at 0.5 pt).
        cands = [p for p in paths if p["dash"] == dash and abs(p["width"] - 0.75) < 1e-9]
        if len(cands) != 1:
            raise SystemExit(f"expected exactly one {series} trajectory, found {len(cands)}")
        pts = []
        for px, py in cands[0]["points"]:
            t = X_MIN + (px - xl) / (xr - xl) * (X_MAX - X_MIN)
            v = Y_MIN + (py - yb) / (yt - yb) * (Y_MAX - Y_MIN)
            pts.append((t, v))
        pts.sort()
        # Drop duplicated segment endpoints.
        dedup = [pts[0]]
        for t, v in pts[1:]:
            if t - dedup[-1][0] > 1e-9:
                dedup.append((t, v))
        for t, v in dedup:
            rows.append((series, SERIES_TO_CONDITION[series], f"{t:.6f}", f"{v:.8f}"))
        print(f"{series:9s} -> {SERIES_TO_CONDITION[series]:12s} "
              f"{len(dedup):4d} points, t in [{dedup[0][0]:.2f}, {dedup[-1][0]:.2f}], "
              f"P01 in [{min(v for _, v in dedup):.5f}, {max(v for _, v in dedup):.5f}]")

    out = here / "reference" / "mu2010_fig15_4_digitized.csv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["series", "condition", "time", "P01"])
        w.writerows(rows)
    print(f"wrote {out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
