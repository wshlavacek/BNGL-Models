#!/usr/bin/env python3
"""Digitize Figs. 2, 5 and 7 of Vilar et al. (2002) into reference/*.csv.

Route
-----
The figures are *embedded bitmaps*, not vector art: PyMuPDF's `get_drawings()`
returns five paths on the Fig. 2 page and none on the Fig. 7 page, while
`pdfimages -list` reports one grayscale JPEG per figure at the publisher's native
150 ppi. Per `skills/curate-model/references/digitization.md` §1 the right route
for a panel that *is* a bitmap is `pdfimages`, which returns those pixels
unresampled, rather than `pdftoppm`, which would rasterize the page a second time.

150 ppi is modest, and that -- not curve tracing -- sets the accuracy. Each Fig. 2
panel is 151 px tall for 2000 molecules, so one pixel is 13.2 molecules. The
images are JPEG, so ink is soft-edged: frame detection and tick finding use a
loose threshold (< 210 of 255) because one panel's top frame is compressed to a
minimum of 186, while curve tracing uses a strict one (< 160) so that JPEG
ringing around the axes is not read as data.

Frame detection
---------------
A panel frame is a row whose longest *contiguous* ink run matches the run shared
by all the other frame rows. Requiring contiguity and a shared extent is what
separates a frame from a flat curve: in Fig. 5a the repressor decays to its fixed
point and lies flat across three quarters of the panel, which any ink-count test
reads as a horizontal line. Its run starts at column 166; every real frame in
that figure runs 91 to 491. The modal span is found, not assumed, so the test
does not depend on where a particular figure puts its axes.

Calibration
-----------
Anchored on tick marks, never on the frame alone (§2). The left and right frame
columns of these panels coincide with the first and last x tick, which the script
*verifies* rather than assumes: it takes t = 0 and t = t_max from those columns,
infers the tick interval, and then asserts that every interior tick mark hanging
inward from a top frame lands within one pixel of the implied grid. If a frame
were offset from its axis limits -- michalski2012's failure mode, the reason for
the rule -- the interior ticks would not land on the grid and the assertion
fires. The same check runs on the y axis against the printed axis maximum.

What is extracted
-----------------
These are relaxation oscillators: within one pixel column the curve can span most
of the panel, so a single y per x is not defined. For every column the script
records the topmost and bottommost ink row, giving an upper and a lower envelope.
The upper envelope is the quantity to compare against, and the comparison in
verify_vilar2002.ipynb reduces the BioNetGen trajectory the same way -- column
binning, then max and min within the bin -- so that the two sides are the same
statistic rather than two different ones that happen to be plotted alike.

Frame rows are excluded from the trace, so the lower envelope saturates one pixel
above the axis (13 molecules in Fig. 2) rather than reaching 0. Treat any lower
envelope at that floor as "at or below the floor", which is how the notebook
treats it.

Usage
-----
    uv run python digitize_vilar2002.py

Writes reference/vilar2002_fig<N>_<panel>_digitized.csv. Re-running must leave
`git diff` empty (§8).
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
PDF = (
    HERE.parents[1]
    / "dev/papers/Vilar2002"
    / "vilar-et-al-2002-mechanisms-of-noise-resistance-in-genetic-oscillators.pdf"
)
OUT = HERE / "reference"

INK_CURVE = 160  # strict: curve ink only
INK_LINE = 210  # loose: JPEG-washed frames and tick marks

# Figure layout. `image` indexes the pdfimages extraction of pages 2-5. `t_max`
# and each panel's y_max are read off the printed axis labels; y_min is 0 in
# every panel of all three figures.
FIGURES = {
    2: dict(
        image=1,
        t_max=400.0,
        panels=[
            ("a", 2000.0, "deterministic free activator A"),
            ("b", 2000.0, "deterministic free repressor R"),
            ("c", 2000.0, "stochastic free activator A"),
            ("d", 2000.0, "stochastic free repressor R"),
        ],
    ),
    5: dict(
        image=4,
        t_max=400.0,
        panels=[
            ("a", 3000.0, "deterministic free repressor R at delta_R = 0.05 /h"),
            ("b", 3000.0, "stochastic free repressor R at delta_R = 0.05 /h"),
        ],
    ),
    7: dict(
        image=6,
        t_max=400.0,
        panels=[
            ("a", 1500.0, "stochastic free activator A, rescaled mRNA rates"),
            ("b", 3.0, "stochastic activator transcript M_A, rescaled mRNA rates"),
            ("c", 1500.0, "stochastic free repressor R, rescaled mRNA rates"),
            ("d", 4.0, "stochastic repressor transcript M_R, rescaled mRNA rates"),
        ],
    ),
}


def extract_images(pdf: Path, tmp: Path) -> list[Path]:
    subprocess.run(
        ["pdfimages", "-png", "-f", "2", "-l", "5", str(pdf), str(tmp / "img")],
        check=True,
    )
    return sorted(tmp.glob("img-*.png"))


def longest_run(row: np.ndarray) -> tuple[int, int]:
    """(start, end) of the longest contiguous True run; (0, -1) if none."""
    best = (0, -1)
    cur: list[int] | None = None
    for i, v in enumerate(row):
        if v:
            cur = [i, i] if cur is None else [cur[0], i]
        else:
            if cur is not None and cur[1] - cur[0] > best[1] - best[0]:
                best = (cur[0], cur[1])
            cur = None
    if cur is not None and cur[1] - cur[0] > best[1] - best[0]:
        best = (cur[0], cur[1])
    return best


def group(values, gap: int = 2) -> list[float]:
    out: list[list[int]] = []
    for v in sorted(int(x) for x in values):
        if out and v - out[-1][-1] <= gap:
            out[-1].append(v)
        else:
            out.append([v])
    return [float(np.mean(g)) for g in out]


def frame_rows(line_ink: np.ndarray) -> tuple[list[int], tuple[int, int]]:
    """Panel frame rows, and the column span they share."""
    width = line_ink.shape[1]
    spans = {}
    for r in range(line_ink.shape[0]):
        s, e = longest_run(line_ink[r])
        if e - s > 0.5 * width:
            spans[r] = (s, e)
    if not spans:
        raise AssertionError("no candidate frame rows found")
    mode = Counter(spans.values()).most_common(1)[0][0]
    keep = [r for r, s in spans.items() if abs(s[0] - mode[0]) <= 3 and abs(s[1] - mode[1]) <= 3]
    return [int(round(v)) for v in group(keep)], mode


def axis_columns(line_ink: np.ndarray, top: int, bottom: int):
    """Locate the two vertical axes of a panel.

    Returns the centre of each axis, which is what the time calibration is
    anchored on, and the first and last interior column, which is where tracing
    may begin and end. The two are not the same: these axes are drawn one or two
    pixels wide, and treating an axis column as data makes that column read as a
    full-height excursion.
    """
    counts = line_ink[top : bottom + 1, :].sum(axis=0)
    full = np.where(counts > 0.9 * (bottom - top))[0]
    if full.size < 2:
        raise AssertionError(f"panel {top}-{bottom}: found {full.size} vertical axis columns, need 2")
    groups: list[list[int]] = []
    for c in full:
        if groups and c - groups[-1][-1] <= 2:
            groups[-1].append(int(c))
        else:
            groups.append([int(c)])
    if len(groups) < 2:
        raise AssertionError(f"panel {top}-{bottom}: vertical axes did not separate: {groups}")
    lo, hi = groups[0], groups[-1]
    return float(np.mean(lo)), float(np.mean(hi)), max(lo) + 1, min(hi) - 1


def inward_ticks(line_ink: np.ndarray, top: int, lo: int, hi: int) -> list[float]:
    """Tick marks hanging down from a panel's top frame, within [lo, hi]."""
    seg = line_ink[top + 1 : top + 6, :].sum(axis=0)
    cand = np.where(seg >= 3)[0]
    cand = cand[(cand >= lo) & (cand <= hi)]
    return group(cand)


def side_ticks(line_ink: np.ndarray, top: int, bottom: int, axis_col: int) -> list[float]:
    """Tick marks pointing inward from the left axis, within the panel."""
    seg = line_ink[:, axis_col + 1 : axis_col + 8].sum(axis=1)
    rows = np.arange(seg.size)
    sel = rows[(seg >= 5) & (rows >= top) & (rows <= bottom)]
    return group(sel)


def verify_grid(found, lo: float, hi: float, n_units: float, what: str, tol: float = 1.0):
    """Assert detected ticks land on a regular grid spanning [lo, hi].

    Tries each plausible tick interval and keeps the one matching the most
    detected marks. Candidates that match nothing are curve crossings picked up
    by the tick window and are reported, not silently dropped.
    """
    best: tuple | None = None
    for n_int in (2, 3, 4, 5, 6, 8, 10, 12, 16, 20):
        grid = lo + (hi - lo) * np.arange(n_int + 1) / n_int
        matched, err = 0, 0.0
        for g in grid:
            near = min(found, key=lambda c: abs(c - g)) if found else None
            if near is not None and abs(near - g) <= tol:
                matched += 1
                err = max(err, abs(near - g))
        if best is None or matched > best[1]:
            best = (n_int, matched, err, grid)
    assert best is not None
    n_int, matched, err, grid = best
    if matched < 3:
        raise AssertionError(
            f"{what}: only {matched} of {len(found)} detected marks land on any "
            f"regular grid spanning the axis; detected {np.round(found, 1).tolist()}"
        )
    ignored = [c for c in found if min(abs(c - g) for g in grid) > tol]
    return dict(
        n_intervals=n_int,
        step_units=n_units / n_int,
        matched=matched,
        max_err_px=err,
        ignored=[round(c, 1) for c in ignored],
        grid=grid,
    )


GUARD = 5  # columns beside each vertical axis that the y ticks reach into


def mask_x_ticks(curve_ink, top, bottom, x_marks, reach=6, pad=1):
    """Blank the tick marks hanging inward from the top and bottom frames.

    Left unmasked they are read as data: in Fig. 2 they put a spurious
    1,987-molecule point in every panel, one per gridline. The masked positions
    are the marks actually *detected* on the image rather than the fitted major
    grid, because these panels carry minor ticks between the labelled ones and a
    mask built from the labelled grid alone leaves the minor ones behind.

    Only each mark's own few pixels are cleared, so a curve crossing the same
    column further down the panel is untouched.

    The y ticks, which reach inward from the two vertical axes, are not masked;
    the GUARD columns beside each axis are dropped from the trace instead. They
    cannot be masked safely because a tick and the curve can coincide exactly:
    the first activator peak of Fig. 2a is 1,745 molecules at t = 2.4 h, which
    lands on both the minor tick at 1,750 and the column beside the axis. The
    cost is the first and last few hours of each panel, which is why the
    comparison in the notebook is made on the established cycle.
    """
    out = curve_ink.copy()
    for m in x_marks:
        c0, c1 = int(round(m)) - pad, int(round(m)) + pad + 1
        out[top : top + reach, c0:c1] = False
        out[bottom - reach + 1 : bottom + 1, c0:c1] = False
    return out


def envelope(curve_ink, top: int, bottom: int, c0: int, c1: int, y_max: float):
    """Per-column ink extent in data units, over the panel interior.

    Frame rows are excluded, as are the vertical axes and the GUARD columns
    beside them. A column with no interior ink is left NaN, which means the
    curve coincides with the bottom axis -- at or below the one-pixel floor --
    not that the data are missing.
    """
    cols = np.arange(c0, c1 + 1)
    hi = np.full(cols.size, np.nan)
    lo = np.full(cols.size, np.nan)
    span = bottom - top
    for i, c in enumerate(cols):
        rows = np.where(curve_ink[top + 1 : bottom, c])[0]
        if rows.size == 0:
            continue
        hi[i] = y_max * (span - (rows.min() + 1)) / span
        lo[i] = y_max * (span - (rows.max() + 1)) / span
    return cols, hi, lo


def write_csv(path: Path, header, rows, comments) -> None:
    lines = [f"# {c}" for c in comments]
    lines.append(",".join(header))
    for row in rows:
        lines.append(",".join("" if np.isnan(v) else f"{v:.4g}" for v in row))
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    if not PDF.exists():
        print(f"source PDF not found: {PDF}", file=sys.stderr)
        return 1
    OUT.mkdir(exist_ok=True)

    with tempfile.TemporaryDirectory() as td:
        images = extract_images(PDF, Path(td))

        for fig, spec in FIGURES.items():
            grey = np.array(Image.open(images[spec["image"]]).convert("L"))
            line_ink = grey < INK_LINE
            curve_ink = grey < INK_CURVE

            rows, _ = frame_rows(line_ink)
            panels = spec["panels"]
            if len(rows) != 2 * len(panels):
                raise AssertionError(
                    f"Fig. {fig}: found {len(rows)} frame rows, "
                    f"expected {2 * len(panels)}: {rows}"
                )
            pairs = [(rows[2 * i], rows[2 * i + 1]) for i in range(len(panels))]

            left, right, in_lo, in_hi = axis_columns(line_ink, *pairs[0])
            px_per_h = (right - left) / spec["t_max"]
            trace_lo = max(in_lo, int(round(left)) + GUARD)
            trace_hi = min(in_hi, int(round(right)) - GUARD)

            ticks: list[float] = []
            for top, _ in pairs:
                ticks.extend(inward_ticks(line_ink, top, int(round(left)), int(round(right))))
            x_marks = group(ticks)
            x_grid = verify_grid(x_marks, left, right, spec["t_max"], f"Fig. {fig} x axis")

            for (label, y_max, quantity), (top, bottom) in zip(panels, pairs):
                span = bottom - top
                y_grid = verify_grid(
                    side_ticks(line_ink, top, bottom, int(round(left))),
                    top,
                    bottom,
                    y_max,
                    f"Fig. {fig}{label} y axis",
                )

                clean = mask_x_ticks(curve_ink, top, bottom, x_marks)
                cols, hi, lo = envelope(clean, top, bottom, trace_lo, trace_hi, y_max)
                t = (cols - left) / px_per_h
                floor = y_max / span
                write_csv(
                    OUT / f"vilar2002_fig{fig}_{label}_digitized.csv",
                    ["time_h", "upper", "lower"],
                    np.column_stack([t, hi, lo]),
                    [
                        f"Vilar et al. (2002) Fig. {fig}, panel {label}: {quantity}",
                        "Source: embedded 150 ppi grayscale JPEG extracted with pdfimages "
                        "from the article PDF; no resampling",
                        f"Panel: frame rows {top}-{bottom} ({span} px for 0-{y_max:g}), "
                        f"axis columns {left:g}-{right:g} (0-{spec['t_max']:g} h)",
                        f"Traced columns {trace_lo}-{trace_hi}: the {GUARD} columns beside each "
                        "vertical axis are dropped because the y tick marks reach into them and "
                        "can coincide with the curve, so the first and last few hours are absent",
                        f"x ticks: {x_grid['matched']} marks on a "
                        f"{x_grid['step_units']:g} h grid, max deviation "
                        f"{x_grid['max_err_px']:.2f} px, ignored {x_grid['ignored']}",
                        f"y ticks: {y_grid['matched']} marks on a "
                        f"{y_grid['step_units']:g} grid, max deviation "
                        f"{y_grid['max_err_px']:.2f} px, ignored {y_grid['ignored']}",
                        "upper/lower are the topmost and bottommost curve ink in each pixel "
                        "column; a relaxation oscillator is not single valued within a column",
                        f"Frame rows and the two axis columns are excluded, so 'lower' "
                        f"saturates at {floor:.3g}. An empty row is a column with no interior "
                        "ink, meaning the curve lies on the bottom axis: read it as "
                        "'at or below the floor', not as missing data",
                        f"Uncertainty: 1 px = {y_max / span:.3g} vertically, "
                        f"{1 / px_per_h:.2f} h horizontally",
                        "Regenerated by digitize_vilar2002.py",
                    ],
                )
                print(
                    f"Fig. {fig}{label}: {np.isfinite(hi).sum()}/{cols.size} columns traced, "
                    f"upper max {np.nanmax(hi):7.1f}, "
                    f"x-grid {x_grid['step_units']:g} h (err {x_grid['max_err_px']:.2f} px), "
                    f"y-grid {y_grid['step_units']:g} (err {y_grid['max_err_px']:.2f} px)"
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
