#!/usr/bin/env python
"""Digitize the experimental points of Cortes et al. (2017) Fig. S3, with their error bars.

Cortes MG, Trinh JT, Zeng L, Balazsi G (2017). Late-Arriving Signals Contribute Less to
Cell-Fate Decisions. Biophys J 113:2110-2120. Supporting Material (Document S1, mmc1.pdf),
Figure S3, physical page 21: "Probabilities of lysogeny from the detailed model versus
experimental data." Filled circles with error bars are the single-cell lysogenization
frequencies of Zeng et al. (2010, Cell 141:682), for initially small (red), average (blue) and
large (green) cells at MOI 1-5; solid and dashed lines are the paper's detailed-model curves
(criteria #1 and #2). The caption does not say what the error bars are; they are recorded here
as drawn (half-lengths), and the verification states how it uses them.

Route (digitization.md s1): the figure is not vector art. Page 21 carries it as SEVEN embedded
bitmaps stacked into horizontal strips (six 1951 x 133 px and one 3902 x 263 px, all spanning
the same 468.2 pt width, i.e. 300 ppi and 600 ppi), so no single embedded image holds the panel.
The page is therefore rendered with PyMuPDF at 300 dpi over the strips' union -- the native
resolution of six of the seven strips -- which reassembles them without resampling those six.

Calibration (s2): on the tick marks, which matplotlib draws inward, short and grey on all four
spines: six y ticks (0..100 %) on the left spine and seven x ticks (MOI 0..6) on the bottom one,
each a least-squares fit; the top and right spines' ticks must land on the same pixels (no
skew), and a residual above 1 px aborts.

Separation (s3): ink-direction classification (`classify_ink`), which keeps the thin
anti-aliased error-bar strokes that a near-colour test loses, against the three series colours.

Reading (s4): per series and MOI, the marker is the one filled-circle blob (17 px across) left
in a +/-0.2-MOI column window after a 3-pixel erosion has removed the curves (5 px) and the
error bar (1-2 px); its bounding-box centre is the point (erosion is symmetric, so the centre is
preserved). The error bar is then traced up and down from the marker along its centre column
(+/-1 px) to its last inked row, and each end must carry a cap -- a row whose ink spans >= 11 px
centred on the bar -- or be the plot frame (a bar clipped at 100 %, flagged `clipped_hi`).
The abscissa is snapped to the integer MOI (the experimental grid); the measured centre is kept
as `x_measured` and must be within 0.1 MOI of it.

Uncertainty: centres and bar ends are read to about +/-1 px = +/-0.13 % lysogeny; allowing for
the curve ink that touches some markers, the recorded per-point digitization uncertainty is
0.5 % lysogeny (4 px).

Usage:
    python digitize_cortes2017.py [path/to/mmc1.pdf]
The PDF is not committed (digitization.md s8); the default path is
dev/papers/cortes2017_bpj_mmc1.pdf. The publisher's file is
https://ars.els-cdn.com/content/image/1-s2.0-S000634951731024X-mmc1.pdf
(sha256 a2dda3584e21126ac24de4e804891f7fd3b0316600da17665c86c946af9d28c2).
Writes reference/cortes2017_figS3_lysogeny_vs_moi_digitized.csv (the data points) and
reference/cortes2017_figS3_model_crit1_digitized.csv (the paper's solid criteria-#1 curves at MOI
1-5, traced as described in `_trace_solid`), so a model can be scored against the data on the
same footing as the paper's own model.
"""

import hashlib
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "skills/curate-model/scripts"))
from digitize import Axis, classify_ink, index_runs, longest_run, open_pdf, write_csv  # noqa: E402

DEFAULT_PDF = HERE.parents[1] / "dev/papers/cortes2017_bpj_mmc1.pdf"
PDF_SHA256 = "a2dda3584e21126ac24de4e804891f7fd3b0316600da17665c86c946af9d28c2"
OUT = HERE / "reference" / "cortes2017_figS3_lysogeny_vs_moi_digitized.csv"
OUT_CURVES = HERE / "reference" / "cortes2017_figS3_model_crit1_digitized.csv"
PAGE_INDEX = 20                                  # 0-based; physical page 21 of mmc1.pdf
DPI = 300
SERIES = (("small", 0.5, (255, 0, 0)), ("average", 1.0, (0, 0, 255)),
          ("large", 1.5, (0, 128, 0)))
DIG_UNC = 0.5                                    # % lysogeny


def render(pdf):
    """The page region covered by the seven figure strips, rendered at 300 dpi, as RGB."""
    import pymupdf
    doc = open_pdf(pdf)
    page = doc[PAGE_INDEX]
    rects = [r for im in page.get_images(full=True) for r in page.get_image_rects(im[0])]
    assert len(rects) == 7, f"expected the 7 figure strips on page 21, found {len(rects)}"
    clip = pymupdf.Rect(min(r.x0 for r in rects), min(r.y0 for r in rects),
                        max(r.x1 for r in rects), max(r.y1 for r in rects))
    pix = page.get_pixmap(dpi=DPI, clip=clip)
    img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    return img[:, :, :3].astype(int)


def _ticks(grey, along, start, step, span):
    """Centres of the short grey ticks leaving a spine inward (>= 4 px of ink)."""
    hits = []
    for p in range(*span):
        run, q = 0, start
        while run < 12 and (grey[p, q] if along == 0 else grey[q, p]):
            run += 1
            q += step
        if run >= 4:
            hits.append(p)
    return [(g[0] + g[-1]) / 2.0 for g in index_runs(np.array(hits))]


def _erode(mask, k):
    for _ in range(k):
        e = mask.copy()
        for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1), (1, -1), (-1, 1)):
            e &= np.roll(np.roll(mask, dr, 0), dc, 1)
        mask = e
    return mask


def _components(mask):
    """Connected components of a boolean mask (8-connected): list of (rows, cols) arrays."""
    seen = np.zeros_like(mask)
    out = []
    for r0, c0 in zip(*np.nonzero(mask), strict=True):
        if seen[r0, c0]:
            continue
        stack, pix = [(r0, c0)], []
        seen[r0, c0] = True
        while stack:
            r, c = stack.pop()
            pix.append((r, c))
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    rr, cc = r + dr, c + dc
                    if 0 <= rr < mask.shape[0] and 0 <= cc < mask.shape[1] and mask[rr, cc] \
                            and not seen[rr, cc]:
                        seen[rr, cc] = True
                        stack.append((rr, cc))
        p = np.array(pix)
        out.append((p[:, 0], p[:, 1]))
    return out


def _bar_end(ink, anyink, cy, xc, r_marker, top, bottom, tag, direction):
    """Walk an error bar outward from its marker (centre row cy, column xc), row by row, through
    rows where the bar column is inked by the series OR covered by another series' line (an
    occluder), stopping at three blank rows. Returns (end row, cap found, clipped): the end is the
    outermost cap-like row -- a series-colour run 13-24 px wide centred on the bar to 2 px -- or,
    if the cap itself is hidden under another line, the last series-colour bar row; clipped if
    the walk reaches the frame."""
    r = int(round(cy + direction * r_marker))
    last, cap, blank = None, None, 0
    while top + 1 < r < bottom and blank < 3:
        row = np.flatnonzero(ink[r, xc - 14:xc + 15]) + xc - 14
        run = [g for g in index_runs(row) if g[0] - 1 <= xc <= g[-1] + 1]
        if run:
            g = run[0]
            last, blank = r, 0
            if 13 <= g[-1] - g[0] + 1 <= 24 and abs((g[0] + g[-1]) / 2.0 - xc) <= 2:
                cap = r
        elif anyink[r, xc - 1:xc + 2].any():
            blank = 0
        else:
            blank += 1
        r += direction
    assert last is not None, f"{tag}: no error bar"
    clipped = (last - top) <= 4 if direction < 0 else (bottom - last) <= 4
    return (cap if cap is not None else last), cap is not None, clipped


def digitize(pdf):
    img = render(pdf)
    grey = (img.max(2) < 215) & (img.max(2) - img.min(2) < 30)
    H, W = grey.shape
    lens_c = {c: longest_run(grey[:, c])[0] for c in range(W)}
    lens_r = {r: longest_run(grey[r, :])[0] for r in range(H)}
    left = min(c for c in range(W // 2) if lens_c[c] > 0.7 * H)
    right = max(c for c in range(W // 2, W) if lens_c[c] > 0.7 * H)
    top = min(r for r in range(H // 2) if lens_r[r] > 0.7 * W)
    bottom = max(r for r in range(H // 2, H) if lens_r[r] > 0.7 * W)

    yt = [t for t in _ticks(grey, 0, left + 2, +1, (top - 2, bottom + 3))]
    yt_r = [t for t in _ticks(grey, 0, right - 1, -1, (top - 2, bottom + 3))]
    xt = [t for t in _ticks(grey, 1, bottom - 1, -1, (left - 2, right + 3))]
    xt_t = [t for t in _ticks(grey, 1, top + 2, +1, (left - 2, right + 3))]
    assert len(yt) == 6 and len(xt) == 7, (yt, xt)
    assert np.allclose(yt, yt_r, atol=1.0) and np.allclose(xt, xt_t, atol=1.0), "skewed frame"
    yax = Axis.from_ticks(yt, [100, 80, 60, 40, 20, 0])
    xax = Axis.from_ticks(xt, [0, 1, 2, 3, 4, 5, 6])
    ppy, ppx = (yt[-1] - yt[0]) / 100.0, (xt[-1] - xt[0]) / 6.0
    assert yax.residual * ppy < 1.0 and xax.residual * ppx < 1.0, \
        (yax.residual * ppy, xax.residual * ppx)
    x_of = np.polyfit([0, 1, 2, 3, 4, 5, 6], xt, 1)

    label, keep = classify_ink(img, [c for _, _, c in SERIES], min_ink=40.0)
    anyink = img.max(2) < 200
    rows, marker_cols = [], {}
    for k, (name, v0, _) in enumerate(SERIES):
        ink = keep & (label == k)
        core = _erode(ink, 3)
        for moi in range(1, 6):
            c = np.polyval(x_of, moi)
            c0, c1 = int(c - 0.2 * ppx), int(c + 0.2 * ppx)
            win = np.zeros_like(core)
            win[top + 1:bottom, c0:c1] = core[top + 1:bottom, c0:c1]
            blobs = [b for b in _components(win) if len(b[0]) >= 40]
            assert len(blobs) == 1, f"{name} MOI {moi}: {len(blobs)} marker blobs"
            br, bc = blobs[0]
            cy = (br.min() + br.max()) / 2.0
            cx = (bc.min() + bc.max()) / 2.0
            xc = int(round(cx))
            r_marker = (br.max() - br.min()) / 2.0 + 3 + 2      # eroded radius + erosion + edge

            bar = (ink, anyink, cy, xc, r_marker, top, bottom, f"{name} MOI {moi}")
            r_hi, cap_hi, clip_hi = _bar_end(*bar, -1)
            r_lo, cap_lo, clip_lo = _bar_end(*bar, +1)
            assert not clip_lo, (name, moi)
            y, y_hi, y_lo = float(yax(cy)), float(yax(r_hi)), float(yax(r_lo))
            x_meas = float(xax(cx))
            assert abs(x_meas - moi) < 0.1, f"{name} MOI {moi}: marker at x = {x_meas:.3f}"
            half_lo, half_hi = y - y_lo, y_hi - y
            # matplotlib draws symmetric bars, so sigma_bar is read from the capped, unclipped
            # halves (their mean; their half-spread is its read error, > 1 % aborts). Where no
            # end shows a clean cap (both covered by other ink, or clipped), the longer
            # unclipped half is used and its read error is set to 0.5 % (4 px, the width of a
            # covering line).
            ends = [(half_lo, cap_lo, clip_lo), (half_hi, cap_hi, clip_hi)]
            capped = [h for h, cap, clip in ends if cap and not clip]
            if capped:
                sigma = float(np.mean(capped))
                sigma_unc = 0.5 * (max(capped) - min(capped))
            else:
                sigma = max(h for h, _, clip in ends if not clip)
                sigma_unc = 0.5
            assert sigma_unc <= 1.0, f"{name} MOI {moi}: bar halves {half_lo:.2f} / {half_hi:.2f}"
            marker_cols.setdefault(k, []).append(cx)
            rows.append(dict(cell_class=name, v0=v0, MOI=moi, pct_lysogeny=round(y, 2),
                             err_lo=round(y_lo, 2), err_hi=round(y_hi, 2),
                             sigma_bar=round(sigma, 2), sigma_bar_unc=round(sigma_unc, 2),
                             cap_lo=int(cap_lo), cap_hi=int(cap_hi), clipped_hi=int(clip_hi),
                             x_measured=round(x_meas, 3)))
    curves = []
    for k, (name, v0, _) in enumerate(SERIES):
        ink = keep & (label == k)
        trace = _trace_solid(ink, marker_cols[k], np.polyval(x_of, 1), np.polyval(x_of, 5),
                             np.polyval(x_of, 2.5), top, bottom)
        xs = np.array([float(xax(c)) for c, _ in trace])
        ys = np.array([float(yax(r)) for _, r in trace])
        for moi in range(1, 6):
            sel = np.abs(xs - moi) <= 0.15
            assert sel.sum() >= 15, f"{name} solid curve: {sel.sum()} samples near MOI {moi}"
            slope, icpt = np.polyfit(xs[sel], ys[sel], 1)
            resid = float(np.abs(np.polyval([slope, icpt], xs[sel]) - ys[sel]).max())
            curves.append(dict(cell_class=name, v0=v0, MOI=moi,
                               pct_lysogeny_crit1=round(float(np.polyval([slope, icpt], moi)), 2),
                               n_samples=int(sel.sum()), fit_resid=round(resid, 2)))
    return rows, yax, xax, ppy, ppx, curves


def _trace_solid(ink, marker_cols, c_lo, c_hi, c_start, top, bottom, skip=14, jump=4):
    """Follow one series' SOLID line (criteria #1) column by column from c_start out to c_lo and
    c_hi. Columns within `skip` px of one of the series' markers are skipped (bar, caps and
    marker cover the line there). At c_start the solid line is told from the dashed one
    (criteria #2) by continuity: each candidate run is followed 80 px each way accepting only
    runs within `jump` px of the previous position, and the solid line is the candidate found
    in the most columns (a dashed line breaks off at every gap)."""
    def runs(c):
        return [(g[0] + g[-1]) / 2.0 for g in index_runs(np.flatnonzero(ink[top + 1:bottom, c]))
                if len(g) <= 12]

    def follow(c0, r0, c_end, step):
        """Accept, column by column, the run within `jump` px of the position predicted from
        the local slope (least squares over the last 12 accepted samples), so the trace carries
        straight across a skipped marker zone."""
        out = []
        for c in range(int(c0), int(c_end), step):
            if min(abs(c - m) for m in marker_cols) <= skip:
                continue
            if len(out) >= 4:
                cc, rr = np.array(out[-12:]).T
                pred = np.polyval(np.polyfit(cc, rr, 1), c)
            else:
                pred = out[-1][1] if out else r0
            near = [r for r in runs(c) if abs(r - pred) <= jump]
            if near:
                out.append((c, min(near, key=lambda r: abs(r - pred))))
        return out

    c0 = int(round(c_start))
    cands = runs(c0)
    assert len(cands) >= 2, f"expected solid and dashed lines at column {c0}, found {cands}"
    score = {r: len(follow(c0, r, c0 + 80, 1)) + len(follow(c0, r, c0 - 80, -1)) for r in cands}
    r0 = max(score, key=score.get)
    best, second = sorted(score.values())[-1], sorted(score.values())[-2]
    assert best >= 0.9 * 150 and second < 0.85 * best, f"solid line not distinguishable {score}"
    trace = follow(c0, r0, c_hi + 1, 1) + follow(c0 - 1, r0, c_lo - 1, -1)
    return sorted((c, top + 1 + r) for c, r in trace)


def main():
    pdf = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_PDF
    if not pdf.is_file():
        sys.exit(f"PDF not found: {pdf}\nDownload the publisher's Supporting Material (URL and "
                 f"sha256 in this script's docstring) and pass its path.")
    sha = hashlib.sha256(pdf.read_bytes()).hexdigest()
    if sha != PDF_SHA256:
        print(f"warning: PDF sha256 {sha} differs from the one this CSV was made from")
    rows, yax, xax, ppy, ppx, curves = digitize(pdf)
    comments = [
        "Experimental single-cell %lysogeny vs MOI (Zeng et al. 2010) as plotted in Cortes et al.",
        "(2017) Biophys J 113:2110, Supporting Material Fig. S3 (mmc1.pdf p. 21): red = small",
        "(v0 ~ 0.5), blue = average (v0 ~ 1), green = large (v0 ~ 1.5) initial cell volume.",
        "Written by digitize_cortes2017.py from a 300-dpi render of the seven embedded strips.",
        f"Calibration: 6 y ticks, residual {yax.residual:.3f} % ({yax.residual * ppy:.2f} px); "
        f"7 x ticks, residual {xax.residual:.4f} MOI ({xax.residual * ppx:.2f} px); "
        f"{ppy:.3f} px per % lysogeny.",
        f"Per-point digitization uncertainty {DIG_UNC} % lysogeny. err_lo / err_hi are the "
        "error-bar ends as traced; cap_lo / cap_hi = 1 where that end carries a clean cap.",
        "sigma_bar is the bar half-length, from the capped ends (the mean if both); an end is",
        "uncapped where the 100 % frame clips it (clipped_hi = 1) or another series' line covers",
        "the cap. sigma_bar_unc is its read error (half-spread of two capped halves, or 0.5 % when",
        "no end is capped and the longer unclipped half is used). The caption",
        "does not say whether the bars are SE, SD or a confidence interval. MOI is snapped to the",
        "integer grid; x_measured is the marker centre as read.",
    ]
    write_csv(OUT, ["cell_class", "v0", "MOI", "pct_lysogeny", "err_lo", "err_hi", "sigma_bar",
                    "sigma_bar_unc", "cap_lo", "cap_hi", "clipped_hi", "x_measured"], rows,
              comments=comments)
    write_csv(OUT_CURVES, ["cell_class", "v0", "MOI", "pct_lysogeny_crit1", "n_samples",
                           "fit_resid"], curves, comments=[
        "The paper's own detailed-model prediction, decision criteria #1 (Fig. S3 SOLID lines),",
        "read at MOI 1-5 by digitize_cortes2017.py from the same render as the data points.",
        "Each value is a straight-line fit through the traced line within 0.15 MOI of the",
        "integer (marker zones excluded; one-sided at MOI 1 and 5, where the lines end);",
        "fit_resid is the fit's largest residual in % lysogeny. Read error ~0.5 % lysogeny;",
        "at MOI 1 and 5 the blue solid and dashed lines nearly meet, so the blue ends are",
        "uncertain by the lines' separation there (< 1 % lysogeny)."])
    for r in rows:
        print(f"  {r['cell_class']:8s} MOI {r['MOI']}: {r['pct_lysogeny']:6.2f} "
              f"[{r['err_lo']:6.2f}, {r['err_hi']:6.2f}]  sigma {r['sigma_bar']:5.2f}"
              f"{'  (clipped at 100)' if r['clipped_hi'] else ''}  x {r['x_measured']:.3f}")


if __name__ == "__main__":
    main()
