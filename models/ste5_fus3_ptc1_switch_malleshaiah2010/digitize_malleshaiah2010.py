#!/usr/bin/env python3
"""Digitize the two red curves in Fig. 3b of Malleshaiah et al. (2010)."""

from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PDF = REPO / "dev/papers/Malleshaiah2010/nature08946.pdf"
REFERENCE = HERE / "reference"

# The source PDF is rendered at 300 dpi. Fig. 3b is then cropped from page 3.
CROP = (1360, 1770, 2010, 2290)

# Plot calibration in cropped-image pixels. The x axis is logarithmic.
X_LEFT = 101.0       # 10^-3 uM
X_RIGHT = 626.0      # 10^1 uM
Y_ZERO = 428.0       # 0 sites
Y_FOUR = 92.0        # 4 sites


def groups(values: np.ndarray) -> list[np.ndarray]:
    """Split sorted integer coordinates at gaps larger than three pixels."""
    if not len(values):
        return []
    cuts = np.flatnonzero(np.diff(values) > 3) + 1
    return list(np.split(values, cuts))


def interpolate_missing(values: np.ndarray) -> np.ndarray:
    good = np.isfinite(values)
    return np.interp(np.arange(len(values)), np.flatnonzero(good), values[good])


def write_curve(path: Path, x_px: np.ndarray, y_px: np.ndarray) -> None:
    alpha_uM = 10 ** (-3.0 + 4.0 * (x_px - X_LEFT) / (X_RIGHT - X_LEFT))
    mean_sites = 4.0 * (Y_ZERO - y_px) / (Y_ZERO - Y_FOUR)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["alpha_factor_uM", "mean_ste5_phosphorylation"])
        writer.writerows(zip(alpha_uM, mean_sites, strict=True))


def main() -> None:
    REFERENCE.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="malleshaiah2010_") as temporary:
        prefix = Path(temporary) / "page"
        subprocess.run(
            [
                "pdftoppm", "-f", "3", "-l", "3", "-r", "300", "-png",
                str(PDF), str(prefix),
            ],
            check=True,
        )
        image = Image.open(f"{prefix}-3.png").convert("RGB").crop(CROP)

    rgb = np.asarray(image)
    red = (rgb[:, :, 0] > 175) & (rgb[:, :, 1] < 100) & (rgb[:, :, 2] < 120)
    x_all = np.arange(int(X_LEFT), int(X_RIGHT) + 1)
    solid = np.full(len(x_all), np.nan)
    dashed = np.full(len(x_all), np.nan)

    for index, x in enumerate(x_all):
        ys = np.flatnonzero(red[75:435, x]) + 75
        centers = [float(np.median(group)) for group in groups(ys) if len(group) >= 2]
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

    solid = interpolate_missing(solid)
    dashed = interpolate_missing(dashed)

    # Five-pixel spacing is finer than the plotted line width while avoiding
    # pseudo-replication of every raster column. Include the right endpoint.
    keep = np.unique(np.append(np.arange(0, len(x_all), 5), len(x_all) - 1))
    write_curve(
        REFERENCE / "malleshaiah2010_fig3b_4ps_digitized.csv",
        x_all[keep],
        solid[keep],
    )
    write_curve(
        REFERENCE / "malleshaiah2010_fig3b_1ps_digitized.csv",
        x_all[keep],
        dashed[keep],
    )
    print(f"Digitized {len(keep)} points per curve from Fig. 3b")


if __name__ == "__main__":
    main()
