#!/usr/bin/env python3
"""Digitize the two red curves in Fig. 3b of Malleshaiah et al. (2010).

Fig. 3b plots mean Ste5 phosphorylation against alpha-factor for one and four
phosphorylation sites, and is published only as a figure. The panel is vector art but
is deliberately rasterized (`pdftoppm`, 300 dpi, page 3): both curves are the same red,
so they are separated by vertical position within a column rather than by any path
attribute, and a pixel column is the natural place to do that.

Calibration is two landmarks per axis, read off the printed tick labels: x is
logarithmic over 10^-3 to 10^1 uM, y linear over 0 to 4 sites.

The inset legend overlays x = 367..430 px, where a column can hold legend rule pixels
as well as curve pixels. Over that interval the branch logic below keeps only the
lower-panel components, and any column left unresolved is filled by interpolation
across a known, bounded obstruction -- see
`skills/curate-model/references/digitization.md` §3 and §4.

The raster route, run grouping, gap filling, decimation, axis calibration and CSV
writing come from `skills/curate-model/scripts/digitize.py`.

Usage
-----
    python digitize_malleshaiah2010.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "skills/curate-model/scripts"))

from digitize import Axis, fill_gaps, index_runs, page_raster, stride_index, write_csv  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PDF = REPO / "dev/papers/Malleshaiah2010/nature08946.pdf"
REFERENCE = HERE / "reference"

# The source PDF is rendered at 300 dpi. Fig. 3b is then cropped from page 3.
PAGE, DPI = 3, 300
CROP = (1360, 1770, 2010, 2290)

# Plot calibration in cropped-image pixels. The x axis is logarithmic.
X_LEFT = 101.0       # 10^-3 uM
X_RIGHT = 626.0      # 10^1 uM
Y_ZERO = 428.0       # 0 sites
Y_FOUR = 92.0        # 4 sites

TO_ALPHA = Axis.from_limits(X_LEFT, X_RIGHT, 1e-3, 10.0, log=True)
TO_SITES = Axis.from_limits(Y_ZERO, Y_FOUR, 0.0, 4.0)


def column_centres(red: np.ndarray, x: int) -> list[float]:
    """Median row of each run of red pixels in one column of the panel.

    Runs are cut at gaps wider than three pixels, which bridges the anti-aliased
    dropouts inside a stroke without merging the two curves; a run of one pixel is
    noise and is dropped.
    """
    ys = np.flatnonzero(red[75:435, x]) + 75
    return [float(np.median(run)) for run in index_runs(ys, max_gap=3) if len(run) >= 2]


def write_curve(path: Path, x_px: np.ndarray, y_px: np.ndarray) -> None:
    write_csv(path, ["alpha_factor_uM", "mean_ste5_phosphorylation"],
              zip(TO_ALPHA(x_px), TO_SITES(y_px), strict=True))


def main() -> None:
    if not PDF.exists():
        raise SystemExit(f"source PDF not found: {PDF}\n"
                         "dev/papers/ is not committed; obtain nature08946.pdf from the publisher.")
    REFERENCE.mkdir(exist_ok=True)
    rgb = page_raster(PDF, PAGE, dpi=DPI, crop=CROP)

    red = (rgb[:, :, 0] > 175) & (rgb[:, :, 1] < 100) & (rgb[:, :, 2] < 120)
    x_all = np.arange(int(X_LEFT), int(X_RIGHT) + 1)
    solid = np.full(len(x_all), np.nan)
    dashed = np.full(len(x_all), np.nan)

    for index, x in enumerate(x_all):
        centers = column_centres(red, x)
        # The legend overlays x=367..430 and y~160/195. For the continuous
        # solid curve, omit that interval and interpolate across it. Outside
        # the overlay the solid curve is the upper (smaller-y) component.
        if centers and x < 340:
            candidates = [center for center in centers if center < 250]
            if candidates:
                solid[index] = min(candidates)
        elif centers and x < 367:
            candidates = [center for center in centers if center < 360]
            if candidates:
                solid[index] = min(candidates)
        elif centers and x <= 430:
            candidates = [center for center in centers if center > 300]
            if len(candidates) >= 2:
                solid[index] = min(candidates)
        elif centers and x > 430:
            lower_panel = [center for center in centers if center > 300]
            if lower_panel:
                solid[index] = min(lower_panel)

        # Retain dashed-line pixels only where their identity is unambiguous.
        # Missing dash gaps are interpolated after the full panel is scanned.
        if x < 340:
            candidates = [center for center in centers if center > 250]
            if candidates:
                dashed[index] = max(candidates)
        elif x < 367:
            candidates = [center for center in centers if center > 300]
            if candidates:
                dashed[index] = max(candidates)
        else:
            candidates = [center for center in centers if center > 300]
            if len(candidates) >= 2:
                dashed[index] = max(candidates)

    solid = fill_gaps(solid)
    dashed = fill_gaps(dashed)

    # Five-pixel spacing is finer than the plotted line width while avoiding
    # pseudo-replication of every raster column. Include the right endpoint.
    keep = stride_index(len(x_all), 5)
    write_curve(REFERENCE / "malleshaiah2010_fig3b_4ps_digitized.csv", x_all[keep], solid[keep])
    write_curve(REFERENCE / "malleshaiah2010_fig3b_1ps_digitized.csv", x_all[keep], dashed[keep])
    print(f"Digitized {len(keep)} points per curve from Fig. 3b")


if __name__ == "__main__":
    main()
