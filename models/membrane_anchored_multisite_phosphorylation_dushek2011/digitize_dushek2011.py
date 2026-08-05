"""Recover the 80 simulated curves of Fig. 2 in Dushek et al. (2011).

The paper reports no tabulated output: its result is four panels of twenty dose-response
curves each, total phosphorylation against the kinase/phosphatase balance for substrates with
N = 1 ... 20 sites. This script extracts all eighty so the curation can compare against
numbers.

Extraction route (digitization.md §1): **vector paths**, via PyMuPDF `page.get_drawings()`.
Fig. 2 is line art -- 378 drawing objects on the page, no embedded images -- and every curve
is a separate stroked path of ~50 segments, so the typesetter's own coordinates come out
exactly and no colour separation is needed. That matters here: MATLAB's default colour order
repeats every seven curves, so N = 1, 8 and 15 are all the same blue, and a raster route
keyed on colour could not tell them apart.

Series identity comes from **drawing order within a panel**, which is N ascending, and is
cross-checked two ways: the seven-colour MATLAB cycle must line up with it, and the curves
must be nested (the N = 1 curve is the shallowest and each larger N is strictly steeper).
Both are asserted, so a reordered figure would fail rather than silently mislabel.

Calibration (digitization.md §2): the axes box is the plot's own white background rectangle,
which in this export coincides with the data limits, x in [-2, 2] and y in [0, 1]. It is not
trusted on its own -- it is validated against a landmark whose answer is known exactly. The
kinase and phosphatase carry identical rate constants, so the model has the exact symmetry
<S>(r) + <S>(1/r) = 1, which puts *every* curve through 0.5 at log10([E]/[F]) = 0. The figure
shows them crossing at a single point there. The extracted crossing is measured and reported,
and the script aborts if it is off by more than 0.01 of full scale.

Usage:
    python digitize_dushek2011.py [path/to/PIIS0006349511001457.pdf]

Re-running must leave `git diff` empty (digitization.md §8).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "skills/curate-model/scripts"))
from digitize import write_csv  # noqa: E402

HERE = Path(__file__).resolve().parent
DEFAULT_PDF = HERE.parents[1] / "dev" / "papers" / "Dushek2011" / "PIIS0006349511001457.pdf"
PAGE = 4                      # journal p. 1192, the Fig. 2 page (1-based)
N_MAX = 20

CITATION = (
    "Dushek O, van der Merwe PA, Shahrezaei V (2011) Biophysical Journal 100:1189-1197. "
    "doi:10.1016/j.bpj.2011.01.060"
)

# MATLAB's default colour order, which the figure cycles through. Used only as a
# cross-check on the drawing-order identification of N.
MATLAB_ORDER = [
    (0.22, 0.325, 0.639), (0.035, 0.502, 0.251), (0.925, 0.118, 0.141),
    (0.102, 0.733, 0.737), (0.643, 0.243, 0.588), (0.741, 0.745, 0.196),
    (0.25, 0.25, 0.25),
]

# The four panels of Fig. 2, keyed by the position of their axes box on the page.
# `regime`/`refractory` name the model that produces them.
PANELS = [
    dict(key="2A", regime="reaction",  refractory=False, col=0, row=0),
    dict(key="2B", regime="diffusion", refractory=False, col=1, row=0),
    dict(key="2C", regime="reaction",  refractory=True,  col=0, row=1),
    dict(key="2D", regime="diffusion", refractory=True,  col=1, row=1),
]


def _close(a, b, tol=0.02):
    return a is not None and max(abs(x - y) for x, y in zip(a, b, strict=False)) < tol


def axes_boxes(drawings):
    """The four main axes rectangles, as (x0, y0, x1, y1) in page points.

    Each panel is drawn as a white filled rectangle (the axes background). The four largest
    such rectangles are the panels; the smaller ones are the Hill-number insets.
    """
    rects = []
    for d in drawings:
        if d.get("fill") and _close(d["fill"], (1.0, 1.0, 1.0)) and \
                any(i[0] == "re" for i in d["items"]):
            r = d["rect"]
            if r.width > 100 and r.height > 70:
                rects.append((r.x0, r.y0, r.x1, r.y1))
    # de-duplicate rectangles drawn twice (fill then stroke)
    out = []
    for r in sorted(rects, key=lambda t: (round(t[1]), t[0])):
        if not any(abs(r[0] - s[0]) < 2 and abs(r[1] - s[1]) < 2 for s in out):
            out.append(r)
    if len(out) != 4:
        raise RuntimeError(f"expected 4 panel axes boxes, found {len(out)}")
    return out


def inset_boxes(drawings):
    """The Hill-number inset axes, which are drawn INSIDE their panel and must be excluded."""
    out = []
    for d in drawings:
        if d.get("fill") and _close(d["fill"], (1.0, 1.0, 1.0)) and \
                any(i[0] == "re" for i in d["items"]):
            r = d["rect"]
            if 20 < r.width <= 100 and 15 < r.height <= 70:
                if not any(abs(r.x0 - s[0]) < 2 and abs(r.y0 - s[1]) < 2 for s in out):
                    out.append((r.x0, r.y0, r.x1, r.y1))
    return out


def curve_paths(drawings, box, insets=()):
    """Every plotted curve inside `box`, in drawing order, as (colour, [(x, y), ...])."""
    x0, y0, x1, y1 = box
    out = []
    for d in drawings:
        segs = [i for i in d["items"] if i[0] == "l"]
        if len(segs) < 12:
            continue
        col = d.get("color")
        if col is None or _close(col, (1.0, 1.0, 1.0)):
            continue
        r = d["rect"]
        if not (r.x0 >= x0 - 1 and r.x1 <= x1 + 1 and r.y0 >= y0 - 1 and r.y1 <= y1 + 1):
            continue
        if any(r.x0 >= i[0] - 1 and r.x1 <= i[2] + 1 and r.y0 >= i[1] - 1 and r.y1 <= i[3] + 1
               for i in insets):
            continue                      # a Hill-number inset line, not a dose-response
        pts = [(segs[0][1].x, segs[0][1].y)] + [(s[2].x, s[2].y) for s in segs]
        out.append((tuple(round(c, 3) for c in col), pts))
    return out


def digitize(pdf: Path):
    import fitz

    page = fitz.open(pdf)[PAGE - 1]
    drawings = page.get_drawings()
    if page.get_images(full=True):
        raise RuntimeError("page carries raster images; the vector route assumes line art")
    boxes = axes_boxes(drawings)
    insets = inset_boxes(drawings)

    # order the four boxes into a 2x2 grid: rows top->bottom, columns left->right
    ys = sorted({round(b[1]) for b in boxes})
    xs = sorted({round(b[0]) for b in boxes})
    grid = {}
    for b in boxes:
        grid[(ys.index(round(b[1])), xs.index(round(b[0])))] = b

    written = []
    for spec in PANELS:
        box = grid[(spec["row"], spec["col"])]
        bx0, by0, bx1, by1 = box
        curves = curve_paths(drawings, box, insets)
        if len(curves) != N_MAX:
            raise RuntimeError(f"{spec['key']}: found {len(curves)} curves, expected {N_MAX}")

        def to_data(pts, bx0=bx0, bx1=bx1, by0=by0, by1=by1):
            xs_ = np.array([q[0] for q in pts])
            ys_ = np.array([q[1] for q in pts])
            lx = -2.0 + 4.0 * (xs_ - bx0) / (bx1 - bx0)
            v = (by1 - ys_) / (by1 - by0)
            o = np.argsort(lx)
            return lx[o], v[o]

        # Identity: drawing order is N ascending. Cross-check it against the MATLAB
        # colour cycle, and against the requirement that the curves be strictly nested.
        for i, (col, _) in enumerate(curves):
            want = MATLAB_ORDER[i % len(MATLAB_ORDER)]
            if not _close(col, want, 0.03):
                raise RuntimeError(
                    f"{spec['key']}: curve {i + 1} is {col}, not the MATLAB colour "
                    f"{want} its drawing position implies -- N identification is unsafe")
        # Nesting: for N ascending the curve at the left edge must fall monotonically.
        # The tolerance scales with the panel's own curve separation, because in the
        # diffusion-limited non-refractory panel the whole point is that the curves
        # COLLAPSE onto the single-site curve -- there the spread is a few thousandths and
        # ordering carries no information, so the colour cycle above is the only identity
        # check that means anything.
        left = [to_data(p)[1][0] for _, p in curves]
        spread = max(left) - min(left)
        nested = "collapsed"
        if spread > 0.01:
            if not all(left[i] > left[i + 1] - 1e-9 for i in range(len(left) - 1)):
                raise RuntimeError(
                    f"{spec['key']}: curves are not nested in drawing order "
                    f"(left-edge spread {spread:.4f})")
            nested = "nested"

        # Calibration check: every curve passes through the point (0, 0.5).
        # Measured as the distance from that point to the curve, in units of panel size,
        # because neither interpolation alone is well conditioned across the family: the
        # panel-D curves are near-vertical at the crossing (Hill ~5.4), so interpolating the
        # ordinate is ill conditioned, while the panel-A N=1 curve is nearly flat (it spans
        # 0.47 to 0.53), so interpolating the abscissa is. A distance is conditioned on
        # neither.
        def miss(pts, to_data=to_data):
            lx, v = to_data(pts)
            s = np.linspace(0, 1, 4000)
            xi = np.interp(s, np.linspace(0, 1, len(lx)), lx) / 4.0
            vi = np.interp(s, np.linspace(0, 1, len(v)), v)
            return float(np.min(np.hypot(xi, vi - 0.5)))
        cross = [miss(p) for _, p in curves]
        err = max(cross)
        if err > 0.01:
            raise RuntimeError(
                f"{spec['key']}: a curve misses the point (0, 0.5) by {err:.4f} of panel "
                f"size -- axes calibration is wrong")

        rows = []
        for n, (_, pts) in enumerate(curves, start=1):
            lx, v = to_data(pts)
            for a, b in zip(lx, v, strict=True):
                rows.append((n, float(a), float(b)))
        print(f"  {spec['key']} ({spec['regime']}-limited, "
              f"{'refractory' if spec['refractory'] else 'non-refractory'}): "
              f"{len(curves)} curves, {len(rows)} points; "
              f"worst miss of the exact (0, 0.5) crossing {err:.4f} of panel size; "
              f"left-edge spread {spread:.4f} ({nested})")

        path = HERE / "reference" / f"dushek2011_fig{spec['key'][1:].lower()}_digitized.csv"
        write_csv(
            path, ["n_sites", "log10_E_over_F", "total_phosphorylation"], rows,
            comments=[
                f"Fig. {spec['key']} - total phosphorylation (Eq. 2) vs log10([E]/[F]), "
                f"N = 1..{N_MAX} sites",
                f"{spec['regime']}-limited regime, "
                f"{'with' if spec['refractory'] else 'without'} enzyme refractory state",
                CITATION,
                "Digitized by digitize_dushek2011.py from the vector paths of page 4.",
                "Series identity is drawing order, cross-checked against the MATLAB colour "
                "cycle and against curve nesting; the axes calibration is validated on the "
                "exact 0.5 crossing at log10([E]/[F]) = 0.",
            ],
            fmt=".6g",
        )
        written.append(path)
    return written


def main() -> int:
    pdf = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_PDF
    if not pdf.exists():
        print(
            f"{pdf} not found.\nThe source PDF is not committed (dev/ is gitignored). "
            "Place the Biophysical Journal article there, or pass its path.",
            file=sys.stderr,
        )
        return 1
    print(f"digitizing {pdf}")
    for p in digitize(pdf):
        print(f"  wrote {p.relative_to(HERE)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
