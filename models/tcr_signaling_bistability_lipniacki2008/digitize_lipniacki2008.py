#!/usr/bin/env python3
"""Digitize the d = 0.05/s steady-state branches of Fig. 7A of Lipniacki et al. (2008).

Fig. 7A plots the stable steady states of ppERK against the number N of agonist peptide-MHC
per cell, log-log, for three agonist dissociation rate constants (d = 0.1, 0.05 and 0.0333/s;
"The unstable states connecting saddle-node bifurcation points are not shown"). This model's
d1 is 0.05/s, so only the blue curves are read: the upper (ERK-high) branch, which begins near
N = 10, and the lower (ERK-low) branch, which ends near N = 700.

Source. The only copy obtainable is the bitmap of the NIH author manuscript, PMC2577002
(NIHMS68204), Fig. 7, served as nihms68204f7.jpg (800 x 660 px, JPEG):

    https://cdn.ncbi.nlm.nih.gov/pmc/blobs/6110/2577002/b97d85a38719/nihms68204f7.jpg

The image is not committed (it is the publisher's figure); its sha256 is pinned below. There is
no vector source, so this is the embedded-bitmap route of digitization.md §1, at ~54 px per
decade of N and ~85 px per decade of ppERK.

Calibration (digitization.md §2). Both axes are logarithmic and carry dotted grid lines at the
labelled ticks. Each interior grid line is located as the centroid of a light-grey column
(row) profile inside the axes, with the legend box masked, and the pixel-to-log10 map is a
least-squares fit through them: N = 10, 100, 1000, 10000, and for ppERK the four lines
labelled 3000, 10000, 30000 and 100000. The lines labelled "3000" and "30000" are drawn at half
decades, 10^3.5 and 10^4.5, and the labels are rounded: with the printed values the four y
lines do not fit one log scale (residual 0.0135 decade, 1.1 px) and the frame's bottom edge
maps to ppERK = 948, while with half decades the residual is 0.0043 decade (0.36 px) and the
bottom edge maps to 999.7. Both fits are computed and the script stops unless the half-decade
one is at least twice as good. The frame is used only as a check: its left and right edges
must map to N = 1 and 1e5, and its bottom and top edges to ppERK = 1000 and 300000, the printed
axis limits, each within one pixel.

Extraction. Pixels are assigned to red, blue or green by ink direction (`classify_ink`,
rana2020), with the three legend-key colours as references. In every pixel column the blue
runs above ppERK ~ 4e4 belong to the upper branch and those below ~ 1.5e4 to the lower one
(nothing blue lies between, once the legend box is masked). A run's centre is the curve's
value in that column. Two flags are written per point, so the comparison can drop what the
raster cannot measure (digitization.md §6):

* `occluded` = 1 where red or green ink lies within two pixels of the blue run, or the run
  is thinner than the stroke less one pixel: another curve is drawn over part of it, so its
  centre is biased (the three upper branches converge above N ~ 150);
* `steep` = 1 where the run is taller than 1.5 x the stroke width, i.e. the curve is close
  to vertical there and a column centre is not a point on it;
* `end_cap` = 1 within half a stroke width plus one pixel of the first or last inked column
  of the branch, where the column holds the end of the stroke rather than the curve.

Each output row is one pixel column of one branch. The stroke width (the median run height
where no other ink is near) is printed and recorded in the CSV header: where a branch has a
vertical tangent, at a fold, its centreline ends half a stroke width inside the inked end.

Usage
-----
    python digitize_lipniacki2008.py [path/to/nihms68204f7.jpg]
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "skills/curate-model/scripts"))

from digitize import Axis, classify_ink, index_runs, write_csv  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SOURCE = REPO / "dev/papers/Lipniacki2008/nihms68204f7.jpg"
SOURCE_URL = "https://cdn.ncbi.nlm.nih.gov/pmc/blobs/6110/2577002/b97d85a38719/nihms68204f7.jpg"
SOURCE_SHA256 = "b7bc01934079a88989ebf4c2cda02b35bbd762ac9375e44b997335cbc598cdab"
OUT = HERE / "reference" / "lipniacki2008_fig7a_d0p05_ppERK_steady_states_digitized.csv"

# Panel A axes frame (image pixels), asserted against the image below.
FRAME_LEFT, FRAME_RIGHT, FRAME_TOP, FRAME_BOTTOM = 291, 561, 59, 269
LEGEND_BOX = (423, 109, 553, 178)                       # x0, y0, x1, y1, masked everywhere
X_TICKS = {10.0: 345, 100.0: 399, 1000.0: 453, 10000.0: 507}        # nominal pixel of each
Y_TICKS = {10**3.5: 227, 1e4: 184, 10**4.5: 142, 1e5: 100}  # grid line, refined below
Y_LABELS = (3000.0, 10000.0, 30000.0, 100000.0)    # what the four lines are labelled
# Legend-key colours, read off the key strokes (row means across the key).
COLOURS = {"red": (220, 37, 35), "blue": (29, 29, 206), "green": (29, 205, 75)}
UPPER_MAX_ROW, LOWER_MIN_ROW = 128, 175   # ppERK ~ 4e4 and ~ 1.5e4


def load(path: Path) -> np.ndarray:
    from PIL import Image

    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != SOURCE_SHA256:
        raise SystemExit(f"{path}: sha256 {digest} is not the pinned {SOURCE_SHA256}")
    return np.asarray(Image.open(path).convert("RGB"), dtype=float)


def frame_edges(img: np.ndarray) -> tuple[int, int, int, int]:
    dark = (img.mean(-1) < 170) & (np.ptp(img, axis=-1) < 60)
    cols = dark[FRAME_TOP:FRAME_BOTTOM + 1].sum(0)
    rows = dark[:, FRAME_LEFT:FRAME_RIGHT + 1].sum(1)
    left = int(np.argmax(cols[280:300])) + 280
    right = int(np.argmax(cols[550:570])) + 550
    top = int(np.argmax(rows[50:70])) + 50
    bottom = int(np.argmax(rows[260:280])) + 260
    return left, right, top, bottom


def grid_lines(img: np.ndarray) -> tuple[dict, dict]:
    """Centroid of each dotted grid line's light-grey profile, near its nominal pixel."""
    grey = (np.ptp(img, axis=-1) < 45) & (img.mean(-1) < 225)
    x0, y0, x1, y1 = LEGEND_BOX
    grey[y0:y1, x0:x1] = False
    cols = grey[FRAME_TOP + 3:FRAME_BOTTOM - 2].sum(0).astype(float)
    rows = grey[:, FRAME_LEFT + 3:FRAME_RIGHT - 2].sum(1).astype(float)

    def centroid(profile, nominal):
        idx = np.arange(nominal - 2, nominal + 3)
        w = profile[idx]
        if w.max() < 15:
            raise SystemExit(f"no grid line near pixel {nominal}")
        return float((idx * w).sum() / w.sum())

    return ({v: centroid(cols, p) for v, p in X_TICKS.items()},
            {v: centroid(rows, p) for v, p in Y_TICKS.items()})


def main() -> None:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else SOURCE
    if not path.exists():
        raise SystemExit(f"source figure not found: {path}\n"
                         f"dev/papers/ is not committed; download Fig. 7 of PMC2577002 from\n"
                         f"{SOURCE_URL}")
    img = load(path)

    edges = frame_edges(img)
    assert edges == (FRAME_LEFT, FRAME_RIGHT, FRAME_TOP, FRAME_BOTTOM), edges
    xt, yt = grid_lines(img)
    to_n = Axis.from_ticks(list(xt.values()), list(xt.keys()), log=True)
    to_erk = Axis.from_ticks(list(yt.values()), list(yt.keys()), log=True)
    print("x grid lines (px):", {f"{k:g}": round(v, 2) for k, v in xt.items()},
          f"fit residual {to_n.residual:.4f} decade")
    print("y grid lines (px):", {f"{k:g}": round(v, 2) for k, v in yt.items()},
          f"fit residual {to_erk.residual:.4f} decade")
    as_labelled = Axis.from_ticks(list(yt.values()), Y_LABELS, log=True)
    print(f"  (with the printed labels instead: residual {as_labelled.residual:.4f} decade, "
          f"bottom edge -> {as_labelled(FRAME_BOTTOM):.1f})")
    if not to_erk.residual < 0.5 * as_labelled.residual:
        raise SystemExit("y grid lines are not at half decades; recheck the calibration")
    # The frame must land on the printed axis limits (1, 1e5; 1000), within one pixel.
    x_frame = max(abs(np.log10(to_n(FRAME_LEFT)) - 0.0), abs(np.log10(to_n(FRAME_RIGHT)) - 5.0))
    y_frame = max(abs(np.log10(to_erk(FRAME_BOTTOM)) - 3.0),
                  abs(np.log10(to_erk(FRAME_TOP)) - np.log10(3e5)))
    if x_frame > 1 / 54 or y_frame > 1 / 85:
        raise SystemExit(f"frame misses the axis limits by {x_frame:.4f} / {y_frame:.4f} decade")
    print(f"frame vs axis limits: x {x_frame:.4f} decade, y {y_frame:.4f} decade")

    x0, y0, x1, y1 = LEGEND_BOX
    window = img[FRAME_TOP + 1:FRAME_BOTTOM, FRAME_LEFT + 1:FRAME_RIGHT]
    label, keep = classify_ink(window, list(COLOURS.values()))
    names = list(COLOURS)
    masks = {}
    for i, name in enumerate(names):
        m = np.zeros(img.shape[:2], bool)
        m[FRAME_TOP + 1:FRAME_BOTTOM, FRAME_LEFT + 1:FRAME_RIGHT] = keep & (label == i)
        m[y0:y1, x0:x1] = False
        masks[name] = m
    blue, other = masks["blue"], masks["red"] | masks["green"]

    points = []
    for x in range(FRAME_LEFT + 1, FRAME_RIGHT):
        for run in index_runs(np.flatnonzero(blue[:, x]), max_gap=1):
            if len(run) < 2:
                continue
            top, bottom = int(run[0]), int(run[-1])
            if bottom <= UPPER_MAX_ROW:
                branch = "high"
            elif top >= LOWER_MIN_ROW:
                branch = "low"
            else:
                raise SystemExit(f"blue run at x={x}, rows {top}..{bottom} is on neither branch")
            near = bool(other[max(top - 2, 0):bottom + 3, max(x - 1, 0):x + 2].any())
            points.append((branch, x, top, bottom, near))

    # Stroke width: the median run height where no other ink is near.
    width = float(np.median([bt - tp + 1 for _, _, tp, bt, near in points if not near]))
    ends = {b: (min(p[1] for p in points if p[0] == b), max(p[1] for p in points if p[0] == b))
            for b in ("high", "low")}
    cap = int(np.ceil(width / 2)) + 1          # columns of an end cap: half a stroke + 1 px
    rows = []
    for branch, x, top, bottom, near in points:
        height = bottom - top + 1
        occluded = int(near or height < width - 1)
        steep = int(height > 1.5 * width)
        end_cap = int(x - ends[branch][0] < cap or ends[branch][1] - x < cap)
        rows.append([branch, x, float(to_n(x)), float(to_erk(0.5 * (top + bottom))), height,
                     occluded, steep, end_cap])
    rows.sort(key=lambda r: (r[0], r[1]))
    for b in ("high", "low"):
        sel = [r for r in rows if r[0] == b]
        clean = [r for r in sel if not (r[5] or r[6] or r[7])]
        print(f"{b}: {len(sel)} columns, ink N = {sel[0][2]:.4g} .. {sel[-1][2]:.4g}; "
              f"{len(clean)} unflagged, N = {clean[0][2]:.4g} .. {clean[-1][2]:.4g}")
    print(f"stroke width: {width:g} px (median run height of columns with no other ink near)")

    comments = [
        "Lipniacki T, Hat B, Faeder JR, Hlavacek WS (2008) J Theor Biol 254:110-122, Fig. 7A: "
        "stable steady states of ppERK vs N (agonist pMHC per cell), blue curves (d = 0.05/s)",
        f"source: PMC2577002 (NIHMS68204) figure bitmap nihms68204f7.jpg, sha256 {SOURCE_SHA256}",
        f"calibration: log-log, least squares through the dotted grid lines (x residual "
        f"{to_n.residual:.4f}, y residual {to_erk.residual:.4f} decade; the y lines labelled "
        f"3000 and 30000 are at 10^3.5 and 10^4.5); frame edges map to N = 1, 1e5 and "
        f"ppERK = 1000, 3e5 within {max(x_frame, y_frame):.4f} decade",
        f"one row per pixel column per branch; ppERK at the centre of the blue run; stroke width "
        f"{width:g} px; occluded = red/green ink within 2 px or run thinner than stroke - 1; "
        f"steep = run taller than 1.5 strokes; end_cap = within {cap} px of the branch's inked end",
        "written by models/tcr_signaling_bistability_lipniacki2008/digitize_lipniacki2008.py",
    ]
    write_csv(OUT, ["branch", "x_px", "N_agonist", "ppERK_molecules", "run_height_px",
                    "occluded", "steep", "end_cap"], rows, comments=comments, fmt=".6g")


if __name__ == "__main__":
    main()
