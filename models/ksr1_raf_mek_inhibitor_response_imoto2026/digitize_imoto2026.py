#!/usr/bin/env python
"""Recover the six model-prediction panels of Fig. 3 in Imoto et al. (2026).

Panel geometry
--------------
Fig. 3 occupies page 9 of `dev/papers/Imoto2026/npjSystBiolAppl.pdf` as a **single
embedded bitmap** (980 x 648 px, ~146 ppi) plus a soft mask. `pdfimages` therefore
gives the publisher's own pixels, which is strictly better than re-rendering the
page (references/digitization.md §1); the RGB image is composited over white
through its soft mask before anything else, because the mask is what makes the
page background white and the raw RGB carries black there.

Each of the three figure rows holds, for MCF7 (left) and PSN1 (right), an
experimental dose response and — the panels digitized here — a model-predicted
dose response at five KSR1 abundances:

    A / B   Type II RAF inhibitor  (x = [RAFi]/Kd, 0-50)
    C / D   Type I-1/2 RAF inhibitor (x = [RAFi]/Kd, 0-50)
    E / F   Cobimetinib             (x = [MEKi]/Kd, 0-10 for MCF7, 0-20 for PSN1)

Every curve starts at exactly 1.0 at zero dose, so the ordinate is ppERK
normalized to its own no-drug steady state at the same KSR1 abundance. That is
what makes these panels usable as reported data without knowing the absolute
ppERK scale.

Calibration
-----------
The pixel-to-data map is fitted by least squares to the **tick marks** (§2),
never to the plot frame: x ticks are read from the row just below the abscissa
and y ticks from the column just left of the ordinate. The fit residual is
printed per panel and is a free unit test of the calibration.

Series separation
-----------------
The five curves are separated by **ink direction** (§3): anti-aliasing and JPEG
ringing composite ink C over white as `a*C + (1-a)*255`, so `255 - observed` is
parallel to `255 - C` whatever the coverage `a`. Reference inks are read from the
legend key strokes of each panel rather than assumed. The key strokes themselves
are masked out before tracing, since they are the only plotted ink inside the
legend (its text is grey and never classifies as a series).

Uncertainty
-----------
The panels are ~176 px wide for 50 Kd units and ~105-140 px tall for the full
ordinate, so one pixel is ~0.28 Kd in x and ~0.005-0.019 in normalized ppERK.
Curve strokes are ~2 px, and the source is JPEG-compressed. Treat the recovered
values as good to about +/-0.02 in normalized ppERK, degrading where curves
converge; where two series overlap within a stroke width the point is dropped
rather than guessed, which is why some columns are blank in the CSVs.

Usage
-----
    python digitize_imoto2026.py [path/to/npjSystBiolAppl.pdf]

Writes `reference/imoto2026_fig3<panel>_<cell>_<drug>_digitized.csv`. Re-running
must leave `git diff` empty.
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "skills/curate-model/scripts"))
from digitize import Axis, write_csv  # noqa: E402

HERE = Path(__file__).resolve().parent
DEFAULT_PDF = HERE.parents[1] / "dev/papers/Imoto2026/npjSystBiolAppl.pdf"
PAGE = 9

KSR_LEVELS = (0, 15, 30, 50, 100)

# Panel boxes and axis calibration. `xticks`/`yticks` are the nominal tick values
# the detected tick pixels are fitted against.
PANELS = {
    "A": dict(cell="MCF7", drug="type_II_RAFi", box=(287, 100, 530, 245),
              xticks=(0, 10, 20, 30, 40, 50), yticks=(0.0, 0.5, 1.0, 1.5)),
    "B": dict(cell="PSN1", drug="type_II_RAFi", box=(790, 100, 980, 250),
              xticks=(0, 10, 20, 30, 40, 50), yticks=(0.0, 0.5, 1.0, 1.5, 2.0)),
    "C": dict(cell="MCF7", drug="type_I_half_RAFi", box=(288, 285, 530, 428),
              xticks=(0, 10, 20, 30, 40, 50), yticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
    "D": dict(cell="PSN1", drug="type_I_half_RAFi", box=(790, 285, 980, 432),
              xticks=(0, 10, 20, 30, 40, 50), yticks=(0.0, 0.25, 0.5, 0.75, 1.0, 1.25)),
    "E": dict(cell="MCF7", drug="cobimetinib", box=(286, 470, 530, 613),
              xticks=(0, 2, 4, 6, 8, 10), yticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
    "F": dict(cell="PSN1", drug="cobimetinib", box=(790, 470, 980, 619),
              xticks=(0, 5, 10, 15, 20), yticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)),
}


# --------------------------------------------------------------------------- input


def load_figure(pdf: Path) -> np.ndarray:
    """Embedded RGB bitmap of page 9, composited over white through its soft mask."""
    from PIL import Image

    with tempfile.TemporaryDirectory() as tmp:
        stem = Path(tmp) / "fig"
        subprocess.run(
            ["pdfimages", "-f", str(PAGE), "-l", str(PAGE), "-png", str(pdf), str(stem)],
            check=True,
        )
        rgb = np.asarray(Image.open(f"{stem}-000.png").convert("RGB"), dtype=float)
        alpha = np.asarray(Image.open(f"{stem}-001.png").convert("L"), dtype=float) / 255.0
    return rgb * alpha[..., None] + 255.0 * (1.0 - alpha[..., None])


# --------------------------------------------------------------------- calibration


def _dark(img: np.ndarray) -> np.ndarray:
    hi, lo = img.max(axis=2), img.min(axis=2)
    return (hi < 160) & (hi - lo < 50)


def _runs(idx: np.ndarray) -> list[float]:
    out: list[list[int]] = []
    for i in idx:
        if out and i == out[-1][-1] + 1:
            out[-1].append(int(i))
        else:
            out.append([int(i)])
    return [(r[0] + r[-1]) / 2 for r in out]


def calibrate(img: np.ndarray, panel: dict) -> tuple[int, int, Axis, Axis, float, float]:
    x0, y0, x1, y1 = panel["box"]
    dark = _dark(img)
    sub = dark[y0:y1, x0:x1]
    xaxis_col = int(np.argmax(sub.sum(axis=0))) + x0
    yaxis_row = int(np.argmax(sub.sum(axis=1))) + y0

    # Tick strokes sit just outside the axis line; their exact offset varies by a
    # pixel or two between panels, so take whichever offset resolves the most ticks.
    xpix = max(
        ([p for p in _runs(np.where(dark[yaxis_row + off, x0:x1])[0] + x0) if p >= xaxis_col - 1]
         for off in (2, 3, 4, 5)), key=len)
    ypix = max(
        ([p for p in _runs(np.where(dark[y0:y1, xaxis_col - off])[0] + y0) if p <= yaxis_row + 1]
         for off in (2, 3, 4, 5, 6)), key=len)
    xpix, xval = _match(xpix, panel["xticks"])
    ypix, yval = _match(ypix, panel["yticks"], descending=True)

    ax = Axis.from_ticks(xpix, xval)
    ay = Axis.from_ticks(ypix, yval)
    return xaxis_col, yaxis_row, ax, ay, ax.residual, ay.residual, xpix


def _match(pix, nominal, descending=False):
    """Keep the detected tick pixels that form the expected evenly-spaced ladder.

    A label glyph occasionally survives the dark mask and a tick occasionally sits
    under a curve, so the ladder is rebuilt from the two extreme ticks and each
    nominal value is matched to the nearest detected pixel within half a step.
    """
    pix = np.asarray(sorted(pix, reverse=descending), dtype=float)
    nominal = np.asarray(nominal, dtype=float)
    span = (pix[-1] - pix[0]) / (nominal[-1] - nominal[0])
    ideal = pix[0] + (nominal - nominal[0]) * span
    keep_p, keep_v = [], []
    for p_ideal, v in zip(ideal, nominal, strict=False):
        d = np.abs(pix - p_ideal)
        if d.min() <= abs(span) * (nominal[1] - nominal[0]) / 2:
            keep_p.append(float(pix[int(np.argmin(d))]))
            keep_v.append(float(v))
    return np.asarray(keep_p), np.asarray(keep_v)


# ---------------------------------------------------------------- series handling


def legend_inks(img: np.ndarray, panel: dict) -> tuple[np.ndarray, list[tuple[int, int, int]]]:
    """Reference ink colours and the pixel rectangles of the legend key strokes."""
    x0, y0, x1, y1 = panel["box"]
    box = img[y0:y1, x0:x1]
    # A thin, heavily anti-aliased key stroke can fall below a saturation of 55
    # (the dark-blue key of the Type I-1/2 PSN1 panel sits at 35), and missing one
    # key is worse than admitting a few extra candidate rows: the even-spacing
    # test below discards those.
    sat = box.max(axis=2) - box.min(axis=2)
    mask = (sat > 30) & (box.max(axis=2) > 80)
    keys = []
    for r in range(mask.shape[0]):
        cols = np.where(mask[r])[0]
        if len(cols) < 8:
            continue
        seg: list[list[int]] = []
        for c in cols:
            if seg and c == seg[-1][-1] + 1:
                seg[-1].append(int(c))
            else:
                seg.append([int(c)])
        s = max(seg, key=len)
        if len(s) >= 9:
            keys.append((r + y0, s[0] + x0, s[-1] + x0,
                         np.median(box[r, s[0]:s[-1] + 1], axis=0), len(s)))

    # Collapse adjacent rows of the same stroke, keeping the longest run.
    merged: list[tuple] = []
    for k in keys:
        if merged and k[0] - merged[-1][0] < 6:
            if k[4] > merged[-1][4]:
                merged[-1] = k
        else:
            merged.append(k)

    # The legend keys are five *evenly spaced* strokes. Requiring even spacing is
    # what keeps a nearly horizontal stretch of a curve — which also reads as a
    # long single-colour run — from being taken for a key, as happened in the
    # Type I-1/2 PSN1 panel where one true key falls just under the run-length
    # threshold and the next candidate row is 43 px further down.
    n = len(KSR_LEVELS)

    # The legend keys are five evenly spaced strokes of visibly different colours.
    # Both conditions are needed: a nearly horizontal stretch of a curve also reads
    # as a long single-colour run, so candidate rows can interleave with the real
    # keys (they do in the MCF7 Cobimetinib panel) and one real key can fall below
    # the run-length threshold (it does in the PSN1 Type I-1/2 panel). So match a
    # ladder rather than five consecutive candidates.
    rows = np.array([m[0] for m in merged], dtype=float)

    def distinct(window):
        d = 255.0 - np.array([w[3] for w in window])
        u = d / np.linalg.norm(d, axis=1)[:, None]
        return np.max(u @ u.T - np.eye(n)) <= 0.995

    picked = None
    for i in range(len(merged)):
        steps = sorted({rows[j] - rows[i] for j in range(i + 1, len(merged))
                        if 6 <= rows[j] - rows[i] <= 20})
        for step in steps:
            idx = []
            for k in range(n):
                want = rows[i] + k * step
                j = int(np.argmin(np.abs(rows - want)))
                if abs(rows[j] - want) > 2.0:
                    idx = []
                    break
                idx.append(j)
            if len(idx) == n and len(set(idx)) == n:
                window = [merged[j] for j in idx]
                if distinct(window):
                    picked = window
                    break
        if picked is not None:
            break
    if picked is None:
        raise RuntimeError(f"no evenly spaced run of {n} distinctly coloured legend keys")

    # Take each key's colour from the most strongly inked row of its stroke.
    # Blending with white preserves ink *direction*, which is what the classifier
    # uses, but a saturated sample is the more honest thing to record.
    inks = []
    for row, c0, c1, colour, _ in picked:
        best = colour
        for r in range(row - 2, row + 3):
            if y0 <= r < y1:
                cand = np.median(img[r, c0:c1 + 1], axis=0)
                if np.linalg.norm(255.0 - cand) > np.linalg.norm(255.0 - best):
                    best = cand
        inks.append(best)
    inks = np.array(inks)
    rects = [(p[0], p[1], p[2]) for p in picked]
    return inks, rects


def trace(img, panel, inks, rects, xaxis_col, yaxis_row, ax, ay, xpix):
    """Weighted ink centre of each series in every pixel column of the panel."""
    x0, y0, x1, y1 = panel["box"]
    work = img.copy()
    for row, c0, c1 in rects:  # blank the legend key strokes
        work[row - 3:row + 4, c0 - 2:c1 + 3] = 255.0

    directions = 255.0 - inks
    directions /= np.linalg.norm(directions, axis=1)[:, None]

    x_hi = min(x1, int(round(max(xpix))) + 1)
    top = y0
    rows = slice(top, yaxis_row)

    xs, series = [], [[] for _ in KSR_LEVELS]
    for col in range(xaxis_col, x_hi):
        strip = work[rows, col]
        d = 255.0 - strip
        mag = np.linalg.norm(d, axis=1)
        ink = mag > 55
        if not ink.any():
            xs.append(ax(col))
            for s in series:
                s.append(np.nan)
            continue
        unit = np.zeros_like(d)
        unit[ink] = d[ink] / mag[ink][:, None]
        cos = unit @ directions.T                     # (pixel, series)
        best = np.argmax(cos, axis=1)
        good = ink & (cos.max(axis=1) > 0.985)
        xs.append(ax(col))
        for k in range(len(KSR_LEVELS)):
            sel = good & (best == k)
            if not sel.any():
                series[k].append(np.nan)
                continue
            w = mag[sel]
            y = np.arange(top, yaxis_row)[sel]
            series[k].append(ay(float((w * y).sum() / w.sum())))
    return np.asarray(xs), [np.asarray(s) for s in series]


# --------------------------------------------------------------------------- main


def main(pdf: Path) -> None:
    img = load_figure(pdf)
    out_dir = HERE / "reference"
    out_dir.mkdir(exist_ok=True)

    for key, panel in PANELS.items():
        xaxis_col, yaxis_row, ax, ay, xres, yres, xpix = calibrate(img, panel)
        inks, rects = legend_inks(img, panel)
        xs, series = trace(img, panel, inks, rects, xaxis_col, yaxis_row, ax, ay, xpix)

        # keep only columns inside the plotted x range
        inside = (xs >= panel["xticks"][0] - 1e-9) & (xs <= panel["xticks"][-1] + 1e-9)
        xs = xs[inside]
        series = [s[inside] for s in series]

        name = f"imoto2026_fig3{key}_{panel['cell']}_{panel['drug']}_digitized.csv"
        header = ["dose_over_Kd"] + [f"ppERK_norm_KSR_{k}nM" for k in KSR_LEVELS]
        rows = [[f"{x:.4f}"] + ["" if np.isnan(s[i]) else f"{s[i]:.4f}" for s in series]
                for i, x in enumerate(xs)]
        write_csv(
            out_dir / name, header, rows,
            comments=[
                f"Fig. 3{key} (model panel), {panel['cell']} cells, {panel['drug']},",
                "from Imoto H, Rauch N, Nagasato AI, Okada M, Kolch W, Rukhlenko OS,",
                "Kholodenko BN (2026), npj Syst Biol Appl, doi:10.1038/s41540-026-00710-6.",
                "",
                "Source: embedded bitmap of page 9 (pdfimages, 980x648 px, ~146 ppi),",
                "composited over white through its soft mask.",
                "Calibration: least squares on the tick marks;",
                f"  max residual x {xres:.4f} Kd, y {yres:.4f} normalized ppERK.",
                "Series separation: ink-direction classification against these legend",
                "  key colours, read from the panel itself:",
                *[f"    KSR = {k:3d} nM   RGB {tuple(int(v) for v in c)}"
                  for k, c in zip(KSR_LEVELS, inks, strict=False)],
                "",
                "Ordinate is ppERK normalized to the no-drug steady state at the same",
                "KSR1 abundance (every curve starts at 1.0). Blank cells are columns",
                "where the series could not be separated from an overlapping one.",
                "Uncertainty ~0.02 in normalized ppERK; see digitize_imoto2026.py.",
            ],
        )
        counts = ", ".join(f"KSR={k}:{int(np.isfinite(s).sum())}"
                           for k, s in zip(KSR_LEVELS, series, strict=False))
        print(f"Fig. 3{key} {panel['cell']:5s} {panel['drug']:17s} "
              f"tick residual x {xres:.4f} y {yres:.4f}  points {counts}")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_PDF)
