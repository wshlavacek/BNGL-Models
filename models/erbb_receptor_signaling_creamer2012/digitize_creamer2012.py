#!/usr/bin/env python3
"""Digitize the Fig. 4 heat map of Creamer et al. (2012) from the authors' original figure file.

Source: "Authors' original file for figure 4", the electronic supplementary material
12918_2012_954_MOESM6_ESM.pdf of doi:10.1186/1752-0509-6-107 (BMC Systems Biology). The heat
map is ONE embedded JPEG bitmap (2432 x 1115 px, the GenePattern export, dendrogram included);
site labels, time labels and the colorbar labels are vector text laid over it. The bitmap is
read at its native resolution (digitization.md section 1: embedded-bitmap route), not a
resampled page render.

Grid. GenePattern draws every matrix cell as an equal block. The heat-map area of the bitmap
is found as the block of non-white pixels (x = 150..2311, y = 15..1114), and the cell grid is
measured from the bitmap itself: colour changes between neighbouring pixel columns fall only
on every second column, so cells are 2 px wide, giving 2162 / 2 = 1081 time columns; rows are
1100 / 55 = 20 px. 1081 columns is exactly a 0-10800 s response sampled every 10 s, so column
k is t = 10 k s. The time labels are NOT used for calibration (digitization.md section 2
allows this with evidence): there are no tick marks, and the hand-placed labels "0" and
"180" sit 2.3 min inside the first column and 6.2 min inside the last one, so a label
calibration would start the heat map 2.4 min BEFORE ligand addition. The label positions are
checked against the grid and the offsets printed. Each of the 55 site labels (vector text)
is asserted to sit within half a row of its grid row, and mapped to the model observable by
molecule and residue (Shc = p52Shc1) using the committed .bngl's observables block.

Value. The figure's colorbar is labelled green = 0, black = 0.14, red = 1 (per-site
max-normalized level). The stops are measured on the heat map: saturated green (51,254,50),
black (0,0,0), red (254,0,0). A cell's colour is the median over its interior 10 x 2 px
(rows 5-14 of the 20, away from JPEG ringing at row edges) and is mapped by projection onto
the green-black segment (value 0-0.14) when G > R, else onto the black-red segment
(0.14-1); this assumes GenePattern blends linearly in RGB between stops, which the figure
cannot confirm. What the figure does fix without that assumption is which side of black a
cell is on (red- or green-dominant, i.e. above or below 0.14 of the site's maximum), and the
verification uses that. JPEG noise: between-cell SD of a saturated-green cell is 9.7 in G,
i.e. 0.005 in value on the green side and up to 0.03 on the red side; written to the header.

Output: reference/creamer2012_fig4_digitized.csv, one row per site (Fig. 4 order, top to
bottom), 1081 columns t = 0, 1/6, ..., 180 min. This replaces an earlier 60-column CSV with no
committed script, which mapped the whole bitmap width (dendrogram margin included) onto
0-180 min: its first four columns (0-9.2 min) were the white margin, read as 0, and every
receptor peak sat near 15 min instead of 1-5 min.

Usage:

    python digitize_creamer2012.py [PDF]

default PDF = dev/papers/creamer2012/12918_2012_954_MOESM6_ESM.pdf (relative to the
repository root; dev/ is not committed). Re-running regenerates the CSV byte for byte (the
JPEG is decoded by Pillow; a different libjpeg could move a cell by a count or so, below the
recorded noise). Needs PyMuPDF and Pillow (`uv sync --group digitize`).
"""

from __future__ import annotations

import io
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "skills/curate-model/scripts"))
from digitize import open_pdf, write_csv  # noqa: E402

DEFAULT_PDF = HERE.parents[1] / "dev/papers/creamer2012/12918_2012_954_MOESM6_ESM.pdf"
MODEL = HERE / "erbb_receptor_signaling_creamer2012.bngl"
OUT = HERE / "reference/creamer2012_fig4_digitized.csv"
N_ROWS, N_COLS, DT_S = 55, 1081, 10.0
BLACK_AT = 0.14


def heatmap_bitmap(doc):
    """The largest embedded image of page 1 (the heat map), decoded, and its placement."""
    from PIL import Image

    page = doc[0]
    imgs = sorted(page.get_images(full=True), key=lambda im: -im[2] * im[3])
    xref = imgs[0][0]
    raw = doc.extract_image(xref)
    arr = np.asarray(Image.open(io.BytesIO(raw["image"])).convert("RGB")).astype(float)
    rect = page.get_image_rects(xref)[0]
    return arr, rect


def spans(page):
    out = []
    for b in page.get_text("dict")["blocks"]:
        for ln in b.get("lines", []):
            for s in ln["spans"]:
                t = s["text"].strip().rstrip("!").strip()
                if t:
                    out.append((t, s["bbox"]))
    return out


def label_to_observable():
    """'ErbB2 Y1139' -> 'ErbB2_1139', from the .bngl's `Molecules NAME MOL(SITE~P!?)` lines."""
    alias = {"Shc": "p52Shc1"}
    obs = {}
    for name, mol, site in re.findall(r"(?m)^\s*Molecules\s+(\S+)\s+(\w+)\((\w+)~P!\?\)",
                                      MODEL.read_text()):
        obs[(mol, site)] = name
    return lambda label: obs[(alias.get(label.split()[0], label.split()[0]), label.split()[1])]


def main(argv):
    pdf = Path(argv[1]) if len(argv) > 1 else DEFAULT_PDF
    if not pdf.exists():
        sys.exit(f"{pdf} not found: download 'Authors' original file for figure 4' "
                 "(12918_2012_954_MOESM6_ESM.pdf) from doi:10.1186/1752-0509-6-107 "
                 "(dev/papers/ is not committed).")
    doc = open_pdf(pdf)
    img, rect = heatmap_bitmap(doc)
    H, W, _ = img.shape
    sx, sy = W / rect.width, H / rect.height  # px per pt

    # Heat-map block: columns / rows that are non-white in almost every row / column.
    nonwhite = img.min(axis=2) < 200
    cols = np.where(nonwhite.mean(axis=0) > 0.95)[0]
    x0, x1 = int(cols.min()), int(cols.max()) + 1
    rows = np.where(nonwhite[:, x0:x1].mean(axis=1) > 0.95)[0]
    y0, y1 = int(rows.min()), int(rows.max()) + 1
    assert np.all(np.diff(cols) == 1) and np.all(np.diff(rows) == 1), "heat map not one block"
    hm = img[y0:y1, x0:x1]
    # Cell width from where the colour jumps (|dR|+|dG|+|dB| > 100 between neighbouring
    # pixel columns): only between columns 2k+1 and 2k+2, so cells are 2 px wide.
    jump = np.abs(np.diff(hm, axis=1)).sum(axis=2) > 100
    odd, even = int(jump[:, 1::2].sum()), int(jump[:, 0::2].sum())
    width = 2 if even < 0.01 * odd else None
    assert width == 2, f"column edges not on a 2-px grid ({even} even vs {odd} odd jumps)"
    print(f"colour jumps between pixel columns: {odd} at odd offsets, {even} at even")
    assert (x1 - x0) == width * N_COLS, f"heat map {x1 - x0} px wide, expected {width * N_COLS}"
    assert (y1 - y0) % N_ROWS == 0, f"heat map {y1 - y0} px tall, not {N_ROWS} equal rows"
    rh = (y1 - y0) // N_ROWS
    print(f"heat map block x {x0}..{x1 - 1}, y {y0}..{y1 - 1}: {N_COLS} columns x {width} px, "
          f"{N_ROWS} rows x {rh} px -> t = 0..{(N_COLS - 1) * DT_S:.0f} s every {DT_S:.0f} s")

    # Labels: sites (right of the map) and times (below it); checked against the grid.
    page = doc[0]
    sp = spans(page)
    to_obs = label_to_observable()
    site_labels = [(t, bb) for t, bb in sp if re.fullmatch(r"\w+ [STY]\d+", t)]
    if len(site_labels) != N_ROWS:
        sys.exit(f"found {len(site_labels)} site labels, expected {N_ROWS}")
    site_labels.sort(key=lambda s: s[1][1])
    sites = []
    for j, (t, bb) in enumerate(site_labels):
        cy = ((bb[1] + bb[3]) / 2 - rect.y0) * sy
        r = (cy - y0) / rh - 0.5
        assert abs(r - j) < 0.5, f"label {t} at row {r:.2f}, expected {j}"
        sites.append((t, to_obs(t)))
    assert len({o for _, o in sites}) == N_ROWS
    for t, bb in sp:
        if t in ("0", "30", "60", "90", "120", "150", "180") and rect.y1 - 1 < bb[1] < rect.y1 + 10:
            cx = ((bb[0] + bb[2]) / 2 - rect.x0) * sx
            tl = ((cx - x0) / width - 0.5) * DT_S / 60.0
            if t in ("0", "180"):
                print(f"  time label {t:>3s} min sits over the column at t = {tl:6.2f} min")

    # Stops, measured on the map: saturated green and red.
    cells = hm.reshape(N_ROWS, rh, N_COLS, width, 3)[:, 5:15]
    med = np.median(cells.transpose(0, 2, 1, 3, 4).reshape(N_ROWS, N_COLS, -1, 3), axis=2)
    green = np.median(med[med[..., 1] > 200], axis=0)
    red = np.array([254.0, 0.0, 0.0])
    noise_g = float(med[med[..., 1] > 200][:, 1].std())
    print(f"  green stop {green.round(0)}, red stop {red}; green-cell noise sd {noise_g:.1f}")
    greenside = med[..., 1] > med[..., 0]
    ug = np.clip(med @ green / (green @ green), 0.0, 1.0)
    ur = np.clip(med @ red / (red @ red), 0.0, 1.0)
    val = np.where(greenside, BLACK_AT * (1.0 - ug), BLACK_AT + (1.0 - BLACK_AT) * ur)
    print(f"  per-site max value: min {val.max(axis=1).min():.3f}, "
          f"median {np.median(val.max(axis=1)):.3f} (1 = the site's own maximum)")

    tmin = np.arange(N_COLS) * DT_S / 60.0
    write_csv(OUT, ["site", *[f"{t:.6g}" for t in tmin]],
              [(o, *[float(v) for v in row]) for (_, o), row in zip(sites, val)],
              comments=[
                  "Creamer MS et al. (2012) BMC Syst Biol 6:107, doi:10.1186/1752-0509-6-107, "
                  "Fig. 4 (authors' original file, ESM 12918_2012_954_MOESM6_ESM.pdf)",
                  "per-site max-normalized phosphorylation, one NFsim run at compartment "
                  "scale 0.2; rows in Fig. 4 order; columns t in min (grid: 1081 x 10 s)",
                  "colour map green=0, black=0.14, red=1 (linear RGB between stops assumed); "
                  f"JPEG noise ~{BLACK_AT * noise_g / green[1]:.3f} (green side) to "
                  f"~{(1 - BLACK_AT) * noise_g / 254:.3f} (red side)",
                  "written by digitize_creamer2012.py"], fmt=".3g")


if __name__ == "__main__":
    main(sys.argv)
