"""Digitize the ThT time courses in Fig. 1 of Rana et al. (2020).

Rana P, Bose P, Vaidya A, Rangachari V, Ghosh P (2020). Global fitting and parameter
identifiability for amyloid-beta aggregation with competing pathways. 2020 IEEE 20th
International Conference on BioInformatics and BioEngineering (BIBE), pp. 73-78.
doi:10.1109/BIBE50027.2020.00020

Fig. 1 has no tabulated counterpart -- the ThT traces exist only as plotted curves -- so the
four panels are recovered from the publisher PDF by colour separation.

Method
------
1. Page 4 of the PDF (article p. 76) is rasterized at 600 dpi with ``pdftoppm``.
2. Each panel's axis box is located as the outermost pair of long horizontal and vertical
   grey runs; the box edges coincide with the axis limits (verified against the printed tick
   labels of panel (a): predicted tick pixels 709.0/842.1/975.2/1108.3/1241.4/1374.5 for
   0/10/20/30/40/50 h, observed label centres within 4 px).
3. Ink-vector classification. Anti-aliasing composites ink C over a white page, so an observed
   pixel is alpha*C + (1-alpha)*255 and the vector 255-obs stays parallel to 255-C whatever the
   coverage alpha. Pixels are therefore assigned to the curve whose 255-C direction is nearest,
   which is stable against the strong anti-aliasing in a 600-dpi rasterization of vector art.
4. For each curve and each pixel column the median row is converted to (time, ThT) through the
   axis calibration, then binned to a uniform grid (``BIN_H`` hours, median within the bin).

Caveats recorded with the output
--------------------------------
* Data markers and the overlaid EKS model line share a colour, so the per-column median tracks
  the centre of the marker cloud with the model line buried inside it. The half-height of the
  marker cloud is reported per point as ``spread_au`` and is the dominant uncertainty
  (typically 0.02-0.05 ThT a.u.).
* Where two curves overlap, pixels go to the nearer colour, so the occluded curve simply has no
  samples over that interval (e.g. panel (b) red is only resolvable after the 24 h dilution,
  before which it runs underneath blue).
* The published caption labels panel (b) "micelle addition event" and (c) "micelle removal
  event", but the body text cites Fig. 1c for the micelle-addition fit ("the three
  experiments", SSE 4.12, Sec. IV-B) and Fig. 1b for the micelle-removal fit (SSE 1.22,
  Sec. IV-C). The kinetics agree with the body text: every trace in panel (b) rises without a
  lag, as expected when fatty acid is present from t=0 and later diluted, while panel (c)
  contains two traces that stay flat until 3 h and 24 h respectively, as expected when fatty
  acid is added to a lag-phase Abeta42 sample. Curves are named after the body text.

Usage
-----
The publisher PDF is a primary source and is not redistributed; point ``pdftoppm`` at your
own copy. Page 4 of the PDF is article page 76, which carries Fig. 1.

    pdftoppm -r 600 -f 4 -l 4 -png <Rana2020>.pdf page
    python digitize_rana2020.py page-4.png reference/

The committed ``reference/rana2020_fig1*_digitized.csv`` files are the output of exactly
that, and are what the verification notebook and the PyBNF jobs read.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from PIL import Image

BIN_H = 0.5  # hours per output bin

# Pure MATLAB primaries (panels a-c) and the default MATLAB colour order (panel d).
RED, GREEN, BLUE = (255, 0, 0), (0, 255, 0), (0, 0, 255)
MBLUE, MYELLOW, MGREEN, MDARKRED = (0, 114, 189), (237, 177, 32), (119, 172, 48), (162, 20, 47)

# frame = (left, right, top, bottom) in pixels of the 600-dpi rasterization of page 4.
PANELS = {
    "fig1a": dict(
        frame=(709.0, 1374.5, 729.5, 1229.5), xlim=(0.0, 50.0), ylim=(0.0, 1.2),
        curves={"on_pathway": RED},
    ),
    "fig1b": dict(
        frame=(1574.5, 2290.0, 714.5, 1253.0), xlim=(0.0, 50.0), ylim=(0.0, 1.5),
        curves={"fa_control": RED, "removal_5h": GREEN, "removal_24h": BLUE},
    ),
    "fig1c": dict(
        frame=(684.5, 1387.0, 1403.5, 1935.0), xlim=(0.0, 80.0), ylim=(0.0, 1.2),
        curves={"fa_control": RED, "addition_3h": GREEN, "addition_24h": BLUE},
    ),
    "fig1d": dict(
        frame=(1586.5, 2261.0, 1426.5, 1935.0), xlim=(0.0, 60.0), ylim=(0.0, 1.5),
        curves={"removal_5h": MBLUE, "addition_3h": MYELLOW, "removal_24h": MGREEN,
                "fa_control": RED, "on_pathway": MDARKRED},
    ),
}


def ink_direction(rgb):
    """Unit vector along 255 - RGB (invariant to anti-aliasing coverage)."""
    inv = 255.0 - np.asarray(rgb, dtype=float)
    return inv / np.maximum(np.linalg.norm(inv, axis=-1, keepdims=True), 1e-9)


def digitize(page_png, tol=0.16, min_ink=90.0, min_sat=25.0):
    page = np.asarray(Image.open(page_png).convert("RGB")).astype(float)
    inked = (np.linalg.norm(255.0 - page, axis=-1) > min_ink) & (page.max(2) - page.min(2) > min_sat)

    rows = []
    for panel, spec in PANELS.items():
        left, right, top, bottom = spec["frame"]
        x0, x1 = spec["xlim"]
        y0, y1 = spec["ylim"]
        # inset by a few pixels so the axis frame itself is never sampled
        li, ri, ti, bi = int(left) + 3, int(right) - 2, int(top) + 3, int(bottom) - 2

        names = list(spec["curves"])
        ref = ink_direction([spec["curves"][n] for n in names])
        window = page[ti:bi, li:ri]
        dist = np.linalg.norm(ink_direction(window)[:, :, None, :] - ref[None, None, :, :], axis=-1)
        label = dist.argmin(2)
        keep = inked[ti:bi, li:ri] & (dist.min(2) < tol)

        for idx, name in enumerate(names):
            mask = keep & (label == idx)
            samples = []
            for col in np.nonzero(mask.any(0))[0]:
                pix = np.nonzero(mask[:, col])[0]
                tht = y0 + (bottom - (pix + ti)) / (bottom - top) * (y1 - y0)
                time = x0 + (col + li - left) / (right - left) * (x1 - x0)
                samples.append((time, float(np.median(tht)),
                                float(tht.max() - tht.min()) / 2.0))
            if not samples:
                continue
            samples = np.array(samples)
            edges = np.arange(x0, x1 + BIN_H, BIN_H)
            which = np.digitize(samples[:, 0], edges) - 1
            for b in np.unique(which):
                sel = samples[which == b]
                rows.append((panel, name, round(float(np.median(sel[:, 0])), 4),
                             round(float(np.median(sel[:, 1])), 5),
                             round(float(np.median(sel[:, 2])), 5), len(sel)))
    return rows


def main():
    page_png, outdir = Path(sys.argv[1]), Path(sys.argv[2])
    outdir.mkdir(parents=True, exist_ok=True)
    rows = digitize(page_png)
    for panel in PANELS:
        sel = [r for r in rows if r[0] == panel]
        path = outdir / f"rana2020_{panel}_digitized.csv"
        with open(path, "w") as fh:
            fh.write("# Digitized from Fig. 1 of Rana et al. (2020), IEEE BIBE 2020, "
                     "doi:10.1109/BIBE50027.2020.00020\n")
            fh.write(f"# 600 dpi rasterization of PDF page 4; colour separation; "
                     f"{BIN_H} h bins. See digitize_rana2020.py for method and caveats.\n")
            fh.write("# spread_au is the half-height of the marker cloud in the bin "
                     "(digitization uncertainty); n_px_cols is the number of pixel columns "
                     "merged into the bin.\n")
            fh.write("curve,time_h,ThT_au,spread_au,n_px_cols\n")
            for _, curve, t, v, s, n in sorted(sel, key=lambda r: (r[1], r[2])):
                fh.write(f"{curve},{t},{v},{s},{n}\n")
        print(f"{path}  {len(sel)} rows  curves={sorted({r[1] for r in sel})}")


if __name__ == "__main__":
    main()
