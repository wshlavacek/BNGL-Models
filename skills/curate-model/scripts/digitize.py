#!/usr/bin/env python
"""Shared machinery for recovering plotted data from a published figure.

The method, and the judgment calls this module deliberately does not make for you, are in
`../references/digitization.md`. Read that first; this is the toolbox it describes.

Everything here was factored out of the seven committed digitizers, and each helper names the
script it came from so there is a worked precedent for every one of them:

    §1 sources      page_raster, embedded_raster, pymupdf_raster, svg_paths, drawings
    §2 calibration  Axis (from_limits / from_ticks), check_ticks, split_tick_segments
    §3 separation   near_colour, ink_direction, classify_ink
    §4 extraction   index_runs, run_centres, longest_run, trace_curve, blobs,
                    ring_centre, weighted_centre, bar_top, fill_gaps
    §5 post         snap_log, on_grid, collapse, bin_median, dedupe_by_x, stride_index
    §6 output       write_csv

Use from a digitizer in a model folder:

    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "skills/curate-model/scripts"))
    from digitize import Axis, index_runs, page_raster, write_csv

numpy is the only hard dependency. Pillow is needed for the raster routes and PyMuPDF for
`pymupdf_raster`/`drawings`; both are imported lazily, so a stdlib-plus-numpy digitizer such as
mu2010's pays for neither.

**This module must never change a committed digitized CSV.** After editing it, re-run the ported
digitizers and confirm `git diff` comes back empty:

    python models/four_flux_network_isotopomer_labeling_mu2010/digitize_mu2010.py \\
        dev/papers/archived_reference_materials/Mu2010/HandbookChemoinformatics10.pdf
    python models/ste5_fus3_ptc1_switch_malleshaiah2010/digitize_malleshaiah2010.py
    git diff --stat models/

Run this file directly for a self-test that exercises every helper on synthetic input.
"""

from __future__ import annotations

import csv
import math
import re
import shutil
import subprocess
import tempfile
import xml.etree.ElementTree as ET
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

__all__ = [
    "Axis",
    "SvgPath",
    "bar_top",
    "bbox_centres",
    "bin_median",
    "blobs",
    "check_ticks",
    "classify_ink",
    "collapse",
    "dedupe_by_x",
    "drawing_segments",
    "drawings",
    "embedded_raster",
    "fill_gaps",
    "in_box",
    "index_runs",
    "ink_direction",
    "longest_run",
    "near_colour",
    "on_grid",
    "open_pdf",
    "page_raster",
    "pymupdf_raster",
    "rgb",
    "ring_centre",
    "run_centres",
    "select_drawings",
    "snap_log",
    "split_tick_segments",
    "stride_index",
    "svg_paths",
    "trace_curve",
    "weighted_centre",
    "write_csv",
]


def _require(module: str, pip_name: str):
    """Import an optional dependency, or fail with the install line instead of a traceback."""
    try:
        return __import__(module)
    except ImportError as exc:  # pragma: no cover - exercised only without the dependency
        raise SystemExit(f"{pip_name} is required for this route: pip install {pip_name}") from exc


def _tool(name: str) -> str:
    """Locate a poppler executable, or fail with a message a reader can act on."""
    path = shutil.which(name)
    if path is None:
        raise SystemExit(f"{name} (poppler) is required for this route but is not on PATH")
    return path


# --------------------------------------------------------------------------------------------
# §1  Sources -- getting pixels or geometry out of the PDF
# --------------------------------------------------------------------------------------------


def page_raster(pdf: Path | str, page: int, dpi: int = 600, crop=None) -> np.ndarray:
    """Rasterize one 1-based PDF page with `pdftoppm`, as an (H, W, 3) uint8 RGB array.

    Use for vector art you have *chosen* to rasterize because colour separation is simpler than
    path surgery -- and say why in the digitizer's docstring (digitization.md §1). When the panel
    is genuinely a bitmap, `embedded_raster` is strictly better.

    `crop` is a PIL box `(left, upper, right, lower)` in pixels of the rendered page.
    From rana2020 (600 dpi) and malleshaiah2010 (300 dpi, cropped).
    """
    image_mod = _require("PIL.Image", "pillow")
    with tempfile.TemporaryDirectory(prefix="digitize_") as tmp:
        prefix = Path(tmp) / "page"
        subprocess.run(
            [_tool("pdftoppm"), "-f", str(page), "-l", str(page), "-r", str(dpi), "-png",
             str(pdf), str(prefix)],
            check=True,
        )
        # pdftoppm zero-pads the page number according to the document length, so glob for it.
        rendered = sorted(Path(tmp).glob("page-*.png"))
        if len(rendered) != 1:
            raise RuntimeError(f"expected 1 rendered page, got {len(rendered)} from {pdf} p{page}")
        img = image_mod.Image.open(rendered[0]).convert("RGB")
        if crop is not None:
            img = img.crop(crop)
        return np.asarray(img)


def embedded_raster(pdf: Path | str, page: int, index: int = 0) -> np.ndarray:
    """Pull an embedded bitmap off a 1-based PDF page with `pdfimages`, at its native resolution.

    Preferred over `page_raster` whenever the panel *is* an image: it hands back the publisher's
    own pixels rather than a resampling of them. From korwek2023 (~230 ppi figure bitmaps).
    """
    image_mod = _require("PIL.Image", "pillow")
    with tempfile.TemporaryDirectory(prefix="digitize_") as tmp:
        stem = Path(tmp) / f"p{page}"
        subprocess.run(
            [_tool("pdfimages"), "-png", "-f", str(page), "-l", str(page), str(pdf), str(stem)],
            check=True,
        )
        pngs = sorted(Path(tmp).glob(f"p{page}-*.png"))
        if not pngs:
            raise RuntimeError(f"no embedded image found on page {page} of {pdf}")
        if index >= len(pngs):
            raise RuntimeError(f"page {page} of {pdf} has {len(pngs)} images; asked for #{index}")
        return np.asarray(image_mod.Image.open(pngs[index]).convert("RGB"))


def open_pdf(pdf: Path | str):
    """Open a PDF with PyMuPDF. Use when one document supplies both vector and raster panels."""
    fitz = _require("fitz", "pymupdf")
    return fitz.open(str(pdf))


def pymupdf_raster(doc, page_index: int) -> np.ndarray:
    """The single embedded bitmap of a 0-based page of an open PyMuPDF document, as uint8 RGB.

    The `embedded_raster` equivalent for when the document is already open for vector work, which
    is why rohrs2018 uses it: Figs. 3 and 5 are bitmaps in the same PDF as the vector Fig. S5.
    """
    fitz = _require("fitz", "pymupdf")
    page = doc[page_index]
    images = page.get_images(full=True)
    if len(images) != 1:
        raise RuntimeError(f"expected one image on page {page_index}, found {len(images)}")
    pix = fitz.Pixmap(doc, images[0][0])
    if pix.n > 4:
        pix = fitz.Pixmap(fitz.csRGB, pix)
    raw = np.frombuffer(pix.samples, dtype=np.uint8)
    return raw.reshape(pix.height, pix.width, pix.n)[:, :, :3]


_NUM = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")
_SVG_NS = "{http://www.w3.org/2000/svg}"


@dataclass(frozen=True)
class SvgPath:
    """One stroked `<path>` of a poppler-rendered page, flattened to page-space points.

    `dash` is the raw `stroke-dasharray` string, which is a series key in its own right -- mu2010
    identifies its four simulation conditions from nothing else.
    """

    points: tuple[tuple[float, float], ...]
    stroke: str
    fill: str
    width: float
    dash: str

    @property
    def bbox(self) -> tuple[float, float, float, float]:
        """(x_min, x_max, y_min, y_max) in page units."""
        xs = [p[0] for p in self.points]
        ys = [p[1] for p in self.points]
        return (min(xs), max(xs), min(ys), max(ys))


def _parse_matrix(transform: str | None) -> tuple[float, ...]:
    if not transform:
        return (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
    m = re.search(r"matrix\(([^)]*)\)", transform)
    if m is None:
        raise ValueError(f"unsupported transform: {transform}")
    return tuple(float(v) for v in _NUM.findall(m.group(1)))


def _apply(mat: tuple[float, ...], x: float, y: float) -> tuple[float, float]:
    a, b, c, d, e, f = mat
    return (a * x + c * y + e, b * x + d * y + f)


def _path_points(d_attr: str, mat: tuple[float, ...], samples: int) -> list[tuple[float, float]]:
    """Flatten an SVG path of absolute M/L/C commands, sampling each cubic at `samples` points."""
    tokens = re.findall(r"[MLCZmlcz]|" + _NUM.pattern, d_attr)
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
        if op.upper() in "ML":
            cur = (coords[0], coords[1])
            pts.append(_apply(mat, *cur))
        else:
            assert cur is not None
            p0 = cur
            p1, p2, p3 = (coords[0], coords[1]), (coords[2], coords[3]), (coords[4], coords[5])
            for k in range(1, samples + 1):
                t = k / samples
                u = 1.0 - t
                x = u**3 * p0[0] + 3 * u**2 * t * p1[0] + 3 * u * t**2 * p2[0] + t**3 * p3[0]
                y = u**3 * p0[1] + 3 * u**2 * t * p1[1] + 3 * u * t**2 * p2[1] + t**3 * p3[1]
                pts.append(_apply(mat, x, y))
            cur = p3
    return pts


def svg_paths(pdf: Path | str, page: int, samples: int = 200) -> list[SvgPath]:
    """Every stroked path of a 1-based PDF page, via `pdftocairo -svg`. Stdlib plus numpy only.

    poppler rewrites the page without resampling, so Bezier control points survive exactly as the
    typesetter emitted them; the only approximation left is the producer's own curve fit. From
    mu2010.
    """
    with tempfile.TemporaryDirectory(prefix="digitize_") as tmp:
        out = Path(tmp) / "page.svg"
        subprocess.run(
            [_tool("pdftocairo"), "-svg", "-f", str(page), "-l", str(page), str(pdf), str(out)],
            check=True,
            capture_output=True,
        )
        svg = out.read_text()

    root = ET.fromstring(svg)
    paths = []
    for el in root.iter(f"{_SVG_NS}path"):
        d_attr = el.get("d")
        if not d_attr or el.get("stroke", "none") == "none":
            continue
        paths.append(
            SvgPath(
                points=tuple(_path_points(d_attr, _parse_matrix(el.get("transform")), samples)),
                stroke=el.get("stroke") or "",
                fill=el.get("fill") or "none",
                width=float(el.get("stroke-width") or 0.0),
                dash=(el.get("stroke-dasharray") or "").strip(),
            )
        )
    return paths


def drawings(doc, page_index: int) -> list[dict]:
    """`page.get_drawings()` for a 0-based page: PyMuPDF's own dicts, unwrapped.

    Kept in PyMuPDF's currency on purpose. Filtering these is what `select_drawings` is for, but
    marker *shape* lives in the raw `items`, and abstracting that away would make the richest
    vector key in the corpus unreachable. From michalski2012 and rohrs2018.
    """
    return list(doc[page_index].get_drawings())


# --------------------------------------------------------------------------------------------
# §2  Calibration -- pixel to data coordinate
# --------------------------------------------------------------------------------------------


class Axis:
    """A pixel-to-data map for one plot axis, linear or log10.

    Build it from the **tick marks** wherever they are measurable (digitization.md §2); the plot
    frame is drawn wherever the plotting library likes, and in michalski2012 it is offset from the
    axis limits by several points.

        x = Axis.from_ticks(tick_pixels, [1e-3, 1e-2, 1e-1, 1.0], log=True)
        print(x.residual)          # calibration disagreeing with itself: a free unit test
        times = x(marker_pixels)

    `from_limits` is the two-landmark form, for the case where the frame is the only measurable
    landmark -- pair it with `check_ticks` so the ticks remain the authority.
    """

    def __init__(self, to_t: Callable[[np.ndarray], np.ndarray], *, log: bool = False,
                 residual: float = 0.0):
        self._to_t = to_t
        self.log = log
        self.residual = residual

    def __call__(self, px):
        """Map pixel coordinate(s) to data coordinate(s).

        Accepts a scalar, list or array, and gives back a plain `float` for a scalar so the result
        drops straight into an f-string format spec.
        """
        arr = np.asarray(px, dtype=float)
        t = self._to_t(arr)
        value = 10.0**t if self.log else t
        return float(value) if arr.ndim == 0 else value

    @classmethod
    def from_limits(cls, px_lo: float, px_hi: float, v_lo: float, v_hi: float,
                    log: bool = False) -> Axis:
        """Two-landmark map: pixel `px_lo` is value `v_lo`, pixel `px_hi` is value `v_hi`."""
        t_lo = math.log10(v_lo) if log else float(v_lo)
        t_hi = math.log10(v_hi) if log else float(v_hi)
        lo, span, width = float(px_lo), t_hi - t_lo, float(px_hi) - float(px_lo)
        if width == 0.0:
            raise ValueError("px_lo and px_hi are the same pixel")
        return cls(lambda px: t_lo + span * (px - lo) / width, log=log)

    @classmethod
    def from_ticks(cls, px: Sequence[float], values: Sequence[float], log: bool = False) -> Axis:
        """Least-squares map through every (tick pixel, tick value) pair.

        `residual` is the largest absolute deviation of a tick from the fit, in data units (or in
        decades for a log axis). On a vector route it should be near machine noise; anything else
        means the ticks were misidentified. From michalski2012, which prints it per panel.
        """
        pix = np.asarray(px, dtype=float)
        tgt = np.log10(np.asarray(values, dtype=float)) if log else np.asarray(values, dtype=float)
        if len(pix) != len(tgt):
            raise ValueError(f"{len(pix)} ticks vs {len(tgt)} values")
        slope, icpt = np.polyfit(pix, tgt, 1)
        residual = float(np.max(np.abs(np.polyval([slope, icpt], pix) - tgt)))
        return cls(lambda p: slope * p + icpt, log=log, residual=residual)


def split_tick_segments(segments: Iterable, tol: float = 0.01) -> tuple[list[float], list[float]]:
    """Sort straight segments into (x of the vertical ones, y of the horizontal ones).

    Tick marks are usually one path of short strokes: a vertical stroke marks an x tick and a
    horizontal one marks a y tick. Segments are `((x0, y0), (x1, y1))` pairs, from
    `drawing_segments` or from an `SvgPath`'s points taken two at a time.
    From michalski2012 (`tol=0.01`) and mu2010 (`tol=1e-6`).
    """
    xs, ys = [], []
    for (x0, y0), (x1, y1) in segments:
        if abs(x0 - x1) < tol and abs(y0 - y1) > tol:
            xs.append(x0)
        elif abs(y0 - y1) < tol and abs(x0 - x1) > tol:
            ys.append(y0)
    return sorted(xs), sorted(ys)


def check_ticks(axis: Axis, tick_px: Sequence[float], nominal: Sequence[float], tol: float,
                label: str = "axis") -> float:
    """Assert every nominal tick value is matched by a detected tick, and return the worst error.

    The assertion form of §2: calibrate on the frame, then make the ticks prove the calibration.
    Raises `SystemExit` -- a digitizer that has mis-located its axis should stop, not emit numbers.
    From mu2010.
    """
    got = [float(axis(p)) for p in tick_px]
    if not got:
        raise SystemExit(f"{label}: no ticks detected")
    worst = 0.0
    for want in nominal:
        err = min(abs(v - want) for v in got)
        if err > tol:
            raise SystemExit(f"{label}: tick {want:g} off by {err:.4g} > {tol:g}")
        worst = max(worst, err)
    return worst


# --------------------------------------------------------------------------------------------
# §3  Separation -- which ink belongs to which series
# --------------------------------------------------------------------------------------------


def near_colour(img: np.ndarray, colour, tol: float = 60.0) -> np.ndarray:
    """Boolean mask of pixels within L1 distance `tol` of `colour`. From rohrs2018.

    Adequate for saturated colours on white. On thin anti-aliased strokes, where most pixels are
    partial coverage, prefer `classify_ink`.
    """
    return np.abs(np.asarray(img, dtype=float) - np.asarray(colour, dtype=float)).sum(-1) < tol


def ink_direction(rgb_values) -> np.ndarray:
    """Unit vector along 255 - RGB, which is invariant to anti-aliasing coverage.

    Anti-aliasing composites ink C over a white page, so an observed pixel is a*C + (1-a)*255 for
    unknown coverage a; 255 - observed therefore stays parallel to 255 - C whatever a is.
    From rana2020.
    """
    inv = 255.0 - np.asarray(rgb_values, dtype=float)
    return inv / np.maximum(np.linalg.norm(inv, axis=-1, keepdims=True), 1e-9)


def classify_ink(window: np.ndarray, colours: Sequence, tol: float = 0.16, min_ink: float = 90.0,
                 min_sat: float = 25.0) -> tuple[np.ndarray, np.ndarray]:
    """Assign each pixel of an RGB window to the nearest of `colours` by ink direction.

    Returns `(label, keep)`: `label[r, c]` indexes `colours`, and `keep[r, c]` is True where the
    pixel is inked enough (`min_ink`), saturated enough (`min_sat`) and close enough (`tol`) to be
    trusted. Mask a series with `keep & (label == i)`. From rana2020.
    """
    win = np.asarray(window, dtype=float)
    inked = (np.linalg.norm(255.0 - win, axis=-1) > min_ink) & (
        win.max(-1) - win.min(-1) > min_sat
    )
    ref = ink_direction(list(colours))
    dist = np.linalg.norm(ink_direction(win)[:, :, None, :] - ref[None, None, :, :], axis=-1)
    return dist.argmin(2), inked & (dist.min(2) < tol)


def rgb(colour, ndigits: int = 2):
    """Round a PDF colour tuple so it compares equal to a literal. From michalski2012.

    PDF colour components are floats: the blue you read out of a legend is 0.4469999969, not 0.447.
    """
    return tuple(round(v, ndigits) for v in colour) if colour else None


def in_box(rect, box, pad: float = 0.5) -> bool:
    """Whether a PyMuPDF rect lies inside `box = (x0, y0, x1, y1)`, with `pad` slack."""
    x0, y0, x1, y1 = box
    return (rect.x0 >= x0 - pad and rect.x1 <= x1 + pad
            and rect.y0 >= y0 - pad and rect.y1 <= y1 + pad)


def select_drawings(paths: Iterable[dict], *, box=None, pad: float = 0.5, stroke=None, fill=None,
                    n_items: int | None = None, kinds: str | None = None, type: str | None = None,
                    max_size: float | None = None, where: Callable[[dict], bool] | None = None):
    """Filter PyMuPDF drawings by position, colour and path structure.

    `stroke`/`fill` are compared after `rgb()` rounding. `n_items` is the segment count and `kinds`
    the sorted distinct segment letters ("c" for a circle of four curves, "l" for a polygon) --
    together they are the marker-shape key that lets michalski2012 pull three unlabelled
    phosphatase levels out of one colour. `where` takes any further predicate.
    """
    out = []
    for g in paths:
        r = g["rect"]
        if box is not None and not in_box(r, box, pad):
            continue
        if max_size is not None and (r.width > max_size or r.height > max_size or r.width == 0):
            continue
        if stroke is not None and rgb(g.get("color")) != stroke:
            continue
        if fill is not None and rgb(g.get("fill")) != fill:
            continue
        if type is not None and g.get("type") != type:
            continue
        if n_items is not None and len(g["items"]) != n_items:
            continue
        if kinds is not None and "".join(sorted({it[0] for it in g["items"]})) != kinds:
            continue
        if where is not None and not where(g):
            continue
        out.append(g)
    return out


def bbox_centres(paths: Iterable[dict]) -> list[tuple[float, float]]:
    """Bounding-box centres of PyMuPDF drawings, sorted and de-duplicated -- i.e. marker centres.

    Cross-check these against a landmark whose answer you know before trusting them: the Fig. S5
    marker paths of rohrs2018 carry a systematic ~0.15 pt offset, so that script reads its points
    off the error bars instead. From michalski2012.
    """
    out = {((g["rect"].x0 + g["rect"].x1) / 2, (g["rect"].y0 + g["rect"].y1) / 2) for g in paths}
    return sorted(out)


def drawing_segments(path: dict) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """The straight `("l", p0, p1)` segments of one PyMuPDF drawing, as coordinate pairs."""
    return [((it[1].x, it[1].y), (it[2].x, it[2].y)) for it in path["items"] if it[0] == "l"]


# --------------------------------------------------------------------------------------------
# §4  Extraction -- ink to pixel coordinates
# --------------------------------------------------------------------------------------------


def index_runs(values, max_gap: int = 1) -> list[np.ndarray]:
    """Split a sorted array of indices into runs, cutting wherever the gap exceeds `max_gap`.

    The primitive under every raster reader here: `max_gap=1` is strict contiguity (rohrs2018),
    and a larger gap tolerates the dropouts of a dashed or anti-aliased stroke (korwek2023 uses 2,
    malleshaiah2010 uses 3).
    """
    v = np.asarray(values)
    if v.size == 0:
        return []
    return list(np.split(v, np.flatnonzero(np.diff(v) > max_gap) + 1))


def run_centres(mask, max_gap: int = 1) -> list[float]:
    """Midpoints of the runs of True in a 1-D boolean mask. From rohrs2018 and korwek2023."""
    return [(g[0] + g[-1]) / 2.0 for g in index_runs(np.flatnonzero(np.asarray(mask)), max_gap)]


def longest_run(mask) -> tuple[int, int, int]:
    """`(length, first, last)` of the longest run of True in a 1-D mask, or `(0, 0, 0)`.

    How rohrs2018 finds an axis spine: it is the column with the longest dark run.
    """
    runs = index_runs(np.flatnonzero(np.asarray(mask)), 1)
    if not runs:
        return (0, 0, 0)
    best = max(runs, key=len)
    return (len(best), int(best[0]), int(best[-1]))


def trace_curve(img: np.ndarray, colour, x_from: int, x_to: int, y_lo: int, y_hi: int,
                tol: float = 60.0, max_jump: float = 12.0) -> list[tuple[int, float]]:
    """Follow a curve column by column, preferring the run nearest the previous column's.

    Returns `(column, row)` pairs in image coordinates, walking from `x_from` towards `x_to`
    (either direction). Continuity is what lets the trace survive a legend sample line of the same
    colour sitting elsewhere in the column, and `max_jump` rejects a run too far away to be the
    curve. From rohrs2018 (Fig. 5A).
    """
    step = 1 if x_to >= x_from else -1
    found, previous = [], None
    for px in range(x_from, x_to, step):
        centres = run_centres(near_colour(img[y_lo:y_hi + 1, px], colour, tol))
        if not centres:
            continue
        centre = centres[0] if previous is None else min(centres, key=lambda c: abs(c - previous))
        if previous is not None and abs(centre - previous) > max_jump:
            continue
        previous = centre
        found.append((px, y_lo + centre))
    return found


def blobs(mask: np.ndarray, area=(60, 400), size=(8, 22), max_aspect: int = 3,
          min_fill: float = 0.6) -> list[tuple[float, float]]:
    """Centres of the roughly circular connected components of a boolean mask, sorted by x.

    The "find the data dots" reader: components are kept only if their pixel count is in `area`,
    both bounding-box dimensions are in `size`, the dimensions differ by at most `max_aspect`, and
    they fill at least `min_fill` of their box. From rohrs2018.
    """
    mask = np.asarray(mask, dtype=bool)
    seen = np.zeros_like(mask)
    out = []
    for r in range(mask.shape[0]):
        for c in range(mask.shape[1]):
            if not mask[r, c] or seen[r, c]:
                continue
            stack, pix = [(r, c)], []
            seen[r, c] = True
            while stack:
                a, b = stack.pop()
                pix.append((a, b))
                for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    na, nb = a + da, b + db
                    if (0 <= na < mask.shape[0] and 0 <= nb < mask.shape[1]
                            and mask[na, nb] and not seen[na, nb]):
                        seen[na, nb] = True
                        stack.append((na, nb))
            rows = [p[0] for p in pix]
            cols = [p[1] for p in pix]
            h, w = max(rows) - min(rows) + 1, max(cols) - min(cols) + 1
            if (area[0] < len(pix) < area[1] and size[0] < h < size[1] and size[0] < w < size[1]
                    and abs(h - w) <= max_aspect and len(pix) > min_fill * h * w):
                out.append(((min(cols) + max(cols)) / 2, (min(rows) + max(rows)) / 2))
    return sorted(out)


def ring_centre(mask: np.ndarray, x: int, y_lo: int, y_hi: int, half_width: int = 3):
    """Centre of an **open** marker: the midpoint of the ring's vertical extent, not its centroid.

    A ring's ink centroid is meaningless -- the middle of the glyph is empty. From korwek2023.
    """
    ys = np.where(mask[y_lo:y_hi, x - half_width:x + half_width + 1].any(axis=1))[0]
    return None if len(ys) == 0 else y_lo + 0.5 * (ys.min() + ys.max())


def weighted_centre(mask: np.ndarray, x: int, y_lo: int, y_hi: int, half_width: int = 2,
                    min_len: int = 3, min_count: int = 3):
    """Intensity-weighted centre of the fattest run of a **filled** marker in a column strip.

    "Fattest" is by summed pixel count, so a solid dot wins over a stray stroke crossing the strip.
    From korwek2023.
    """
    counts = mask[y_lo:y_hi, x - half_width:x + half_width + 1].sum(axis=1)
    best, weight = None, 0
    for run in index_runs(np.where(counts >= min_count)[0], 2):
        if len(run) < min_len:
            continue
        w = counts[run].sum()
        if w > weight:
            best, weight = run, w
    if best is None:
        return None
    return y_lo + float((best * counts[best]).sum() / counts[best].sum())


def bar_top(panel: np.ndarray, colour, tol: float = 35.0, min_height: int = 6,
            min_columns: int = 5):
    """Topmost row of a solid bar of `colour`, or None if the bar is absent.

    Columns with fewer than `min_height` matching pixels are discarded first, so a legend swatch
    or an axis label in the same colour cannot set the bar top. From rohrs2018.
    """
    mask = near_colour(panel, colour, tol)
    keep = mask.sum(axis=0) > min_height
    if keep.sum() < min_columns:
        return None
    mask = mask.copy()
    mask[:, ~keep] = False
    return int(np.where(mask.any(axis=1))[0].min())


def fill_gaps(values: np.ndarray) -> np.ndarray:
    """Linearly interpolate over the NaNs of a 1-D array, against the index.

    Only for a *known, bounded* obstruction -- malleshaiah2010 uses it across an inset legend, and
    says so in its docstring. Interpolating across an interval where you simply could not separate
    two series is inventing data (digitization.md §3).
    """
    values = np.asarray(values, dtype=float)
    good = np.isfinite(values)
    return np.interp(np.arange(len(values)), np.flatnonzero(good), values[good])


# --------------------------------------------------------------------------------------------
# §5  Post-processing
# --------------------------------------------------------------------------------------------


def snap_log(values, per_decade=(1.0, 2.0, 5.0)) -> list[float]:
    """Snap values to the nearest point of a repeating per-decade grid.

    Legitimate on an **abscissa** whose sampling grid the paper states -- the true value is known
    exactly, so this removes digitization noise. Never apply it to the ordinate: that is the
    quantity being measured (digitization.md §7). From michalski2012.
    """
    out = []
    for v in values:
        dec = np.floor(np.log10(v) + 1e-9)
        out.append(min(
            (abs(np.log10(v) - np.log10(m * 10**d)), m * 10**d)
            for d in (dec - 1, dec, dec + 1)
            for m in per_decade
        )[1])
    return out


def on_grid(points, axis: Callable, tol: float = 0.02, per_decade=(1.0, 2.0, 5.0)):
    """Keep only points whose abscissa lies on the per-decade sampling grid.

    An in-panel legend draws its keys with the same path objects as the data, so every colour and
    shape filter picks them up; the keys sit at arbitrary abscissae and this drops them
    (digitization.md §5). From michalski2012.
    """
    if not points:
        return []
    vals = np.atleast_1d(axis([p[0] for p in points]))
    return [p for p, v, s in zip(points, vals, snap_log(vals), strict=True)
            if abs(np.log10(v) - np.log10(s)) < tol]


def collapse(points, tol: float = 2.0):
    """Average points sharing an abscissa to within `tol` -- overlapping copies of one series."""
    out, cur = [], []
    for x, y in sorted(points):
        if cur and x - cur[-1][0] > tol:
            out.append((float(np.mean([p[0] for p in cur])), float(np.mean([p[1] for p in cur]))))
            cur = []
        cur.append((x, y))
    if cur:
        out.append((float(np.mean([p[0] for p in cur])), float(np.mean([p[1] for p in cur]))))
    return out


def bin_median(samples, lo: float, hi: float, width: float) -> list[tuple[np.ndarray, int]]:
    """Bin rows of `samples` by their first column and take the column-wise median of each bin.

    Returns `(median_row, n_rows)` per non-empty bin. Reduces a per-pixel-column sampling that was
    never independent to begin with; it does not change a value. From rana2020.
    """
    samples = np.asarray(samples, dtype=float)
    edges = np.arange(lo, hi + width, width)
    which = np.digitize(samples[:, 0], edges) - 1
    out = []
    for b in np.unique(which):
        sel = samples[which == b]
        out.append((np.median(sel, axis=0), len(sel)))
    return out


def dedupe_by_x(points, tol: float = 1e-9):
    """Drop points whose abscissa repeats the previous one -- duplicated segment endpoints.

    Expects `points` sorted by abscissa. From mu2010, where consecutive Bezier segments share a
    knot and would otherwise emit it twice.
    """
    points = list(points)
    if not points:
        return []
    out = [points[0]]
    for p in points[1:]:
        if p[0] - out[-1][0] > tol:
            out.append(p)
    return out


def stride_index(n: int, step: int, include_last: bool = True) -> np.ndarray:
    """Every `step`-th index of `n`, plus the last one -- a decimation that keeps the endpoint.

    Choose `step` finer than the plotted line width but coarse enough to avoid pseudo-replicating
    every raster column. From malleshaiah2010 (5 px at 300 dpi).
    """
    idx = np.arange(0, n, step)
    if include_last and n > 0:
        idx = np.append(idx, n - 1)
    return np.unique(idx)


# --------------------------------------------------------------------------------------------
# §6  Output
# --------------------------------------------------------------------------------------------


def write_csv(path: Path | str, header: Sequence[str], rows: Iterable, *, comments: Sequence[str] =
              (), fmt: str | None = None, lineterminator: str = "\n", nan: str = "",
              quiet: bool = False) -> int:
    """Write a digitized CSV with its provenance header, and return the row count.

    `comments` are emitted as `#` lines above the header: the caption, the citation, and the script
    that wrote the file, so the CSV is usable by someone who cannot see the figure.

    `fmt` is a format spec applied to floats -- `".6g"` is the house choice, well past the
    precision of any raster route. `None` writes values unchanged, which is what a caller that has
    already formatted them wants. NaNs become `nan` (empty by default), because a blank cell is a
    missing measurement and `nan` in a CSV is a value some readers will happily parse.

    `lineterminator` defaults to LF. `csv.writer`'s own default is CRLF, which mixes with any
    `fh.write` header and gets the file rejected by the `mixed-line-ending` pre-commit hook
    (digitization.md §8). Pass `"\\r\\n"` only to preserve an existing committed file.

    Rows may be sequences or dicts keyed by `header`.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    def cell(value):
        if isinstance(value, float) and math.isnan(value):
            return nan
        if fmt is not None and isinstance(value, float):
            return format(value, fmt)
        return value

    count = 0
    with path.open("w", newline="", encoding="utf-8") as fh:
        for line in comments:
            fh.write(f"# {line}{lineterminator}")
        writer = csv.writer(fh, lineterminator=lineterminator)
        writer.writerow(header)
        for row in rows:
            values = [row.get(k) for k in header] if isinstance(row, dict) else list(row)
            writer.writerow([cell(v) for v in values])
            count += 1
    if not quiet:
        print(f"wrote {path} ({count} rows)")
    return count


# --------------------------------------------------------------------------------------------
# Self-test
# --------------------------------------------------------------------------------------------


def _self_test() -> None:
    """Exercise every helper on synthetic input. No PDF, no optional dependency."""
    import tempfile as _tf

    # §2 calibration -----------------------------------------------------------------------
    lin = Axis.from_limits(100.0, 600.0, 0.0, 50.0)
    assert abs(float(lin(350.0)) - 25.0) < 1e-12
    log = Axis.from_limits(101.0, 626.0, 1e-3, 10.0, log=True)
    assert abs(float(log(101.0)) - 1e-3) < 1e-15
    assert abs(float(log(626.0)) - 10.0) < 1e-12
    fitted = Axis.from_ticks([10.0, 20.0, 30.0, 40.0], [1.0, 10.0, 100.0, 1000.0], log=True)
    assert fitted.residual < 1e-12
    assert abs(float(fitted(25.0)) - 10**1.5) < 1e-9
    assert np.allclose(lin([100.0, 600.0]), [0.0, 50.0])
    segs = [((5.0, 0.0), (5.0, 3.0)), ((0.0, 7.0), (2.0, 7.0)), ((1.0, 1.0), (1.0, 1.0))]
    assert split_tick_segments(segs, tol=1e-6) == ([5.0], [7.0])
    assert check_ticks(lin, [200.0, 350.0], [10.0, 25.0], tol=1e-9) < 1e-9
    try:
        check_ticks(lin, [200.0], [11.0], tol=0.1, label="x")
    except SystemExit:
        pass
    else:  # pragma: no cover
        raise AssertionError("check_ticks accepted a tick that is off by 1.0")

    # §3 separation ------------------------------------------------------------------------
    img = np.full((9, 9, 3), 255, dtype=np.uint8)
    img[4, 2:7] = (255, 0, 0)
    assert near_colour(img, (255, 0, 0), tol=30).sum() == 5
    half = np.array([255, 128, 128], dtype=float)          # 50 % coverage of pure red on white
    assert np.allclose(ink_direction(half), ink_direction([255, 0, 0]), atol=1e-9)
    img[6, 2:7] = (0, 0, 255)
    label, keep = classify_ink(img, [(255, 0, 0), (0, 0, 255)])
    assert (keep & (label == 0)).sum() == 5 and (keep & (label == 1)).sum() == 5
    assert rgb((0.4469999969, 0.0, 1.0)) == (0.45, 0.0, 1.0) and rgb(None) is None

    # §4 extraction ------------------------------------------------------------------------
    assert [list(r) for r in index_runs([0, 1, 2, 7, 8], 1)] == [[0, 1, 2], [7, 8]]
    assert [list(r) for r in index_runs([0, 1, 2, 4], 3)] == [[0, 1, 2, 4]]
    assert index_runs([]) == []
    assert run_centres([True, True, False, False, True]) == [0.5, 4.0]
    assert longest_run([False, True, True, True, False, True]) == (3, 1, 3)
    assert longest_run([False, False]) == (0, 0, 0)
    ramp = np.full((20, 20, 3), 255, dtype=np.uint8)
    for c in range(3, 17):
        ramp[c, c] = (255, 0, 0)
    traced = trace_curve(ramp, (255, 0, 0), 3, 17, 0, 19)
    assert traced == [(c, float(c)) for c in range(3, 17)]
    dots = np.zeros((30, 30), dtype=bool)
    dots[4:14, 4:14] = True
    assert blobs(dots, area=(60, 400), size=(8, 22)) == [(8.5, 8.5)]
    ring = np.zeros((20, 20), dtype=bool)
    ring[5, 8:12] = ring[11, 8:12] = True
    ring[5:12, 8] = ring[5:12, 11] = True
    assert ring_centre(ring, 10, 0, 20) == 8.0
    blob = np.zeros((20, 20), dtype=bool)
    blob[8:12, 8:13] = True
    assert weighted_centre(blob, 10, 0, 20) == 9.5
    bars = np.full((20, 20, 3), 255, dtype=np.uint8)
    bars[7:20, 5:15] = (57, 90, 171)
    assert bar_top(bars, (57, 90, 171)) == 7
    assert bar_top(bars, (200, 20, 20)) is None
    assert np.allclose(fill_gaps([1.0, np.nan, 3.0]), [1.0, 2.0, 3.0])

    # §5 post-processing -------------------------------------------------------------------
    assert snap_log([0.0104, 1.93, 480.0]) == [0.01, 2.0, 500.0]
    identity = Axis.from_limits(0.0, 1.0, 1.0, 10.0, log=True)
    kept = on_grid([(0.0, 9.9), (0.30103, 9.8), (0.5, 9.7)], identity)
    assert [round(p[0], 5) for p in kept] == [0.0, 0.30103]   # 1 and 2 survive, 10**0.5 does not
    assert collapse([(1.0, 10.0), (2.0, 20.0), (9.0, 90.0)], tol=2.0) == [(1.5, 15.0), (9.0, 90.0)]
    binned = bin_median([[0.1, 5.0], [0.3, 7.0], [1.2, 9.0]], 0.0, 2.0, 1.0)
    assert [(list(m), n) for m, n in binned] == [([0.2, 6.0], 2), ([1.2, 9.0], 1)]
    assert dedupe_by_x([(0.0, 1.0), (0.0, 2.0), (1.0, 3.0)]) == [(0.0, 1.0), (1.0, 3.0)]
    assert dedupe_by_x([]) == []
    assert list(stride_index(11, 5)) == [0, 5, 10]
    assert list(stride_index(12, 5)) == [0, 5, 10, 11]

    # §6 output ----------------------------------------------------------------------------
    with _tf.TemporaryDirectory() as tmp:
        out = Path(tmp) / "sub" / "t.csv"
        write_csv(out, ["a", "b"], [(1.5, float("nan")), {"a": 2.5, "b": 3.5}],
                  comments=["caption", "citation"], fmt=".6g", quiet=True)
        assert out.read_bytes() == b"# caption\n# citation\na,b\n1.5,\n2.5,3.5\n"
        write_csv(out, ["a"], [("0.100000",)], lineterminator="\r\n", quiet=True)
        assert out.read_bytes() == b"a\r\n0.100000\r\n"

    print("digitize.py self-test: all helpers OK")


if __name__ == "__main__":
    _self_test()
