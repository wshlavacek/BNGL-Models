# Digitizing a published figure

How to recover numbers from a plotted panel when the paper publishes no table, and how to record
what you recovered so that the result is a measurement rather than an impression. A digitized CSV
in `reference/` is *derived data*: something downstream compares against it and states a tolerance,
so it needs a stated method, a stated uncertainty, and a script that regenerates it.

28 of the 54 model folders carry at least one `*_digitized.csv`, 130 files in all — **and only six
of them carry the script that produced it.** The other 22 folders hold numbers nobody can
regenerate, check or correct. That is the gap this document exists to stop widening: §8 is the
part that matters most, and it is the cheapest.

The protocol below is generalized from the seven committed digitizers — `rohrs2018` (532 lines),
`michalski2012` (356), `korwek2023` (282), `mu2010` (178), `rana2020` (158), `malleshaiah2010`
(134) and `Kirsch-2020/ppatf2_phospho` (40) — whose shared machinery is factored into
`scripts/digitize.py`, described in §9.

The companion rule for fitting data is `validate-pybnf-job/references/data-provenance.md`, which
uses this method to trace an `.exp` back to its primary source.

## Contents
1. Choose the extraction route
2. Calibrate on the tick marks
3. Separating series
4. Reading a value out of ink
5. Legend keys are inside the plot area
6. What to record alongside the numbers
7. Precision you can claim, and snapping
8. The reproducibility contract
9. The `digitize.py` library
10. Worked examples

---

## 1. Choose the extraction route

**Look for vector art before you rasterize anything.** Most journal PDFs carry line plots as
stroked paths, and a path gives you the typesetter's own coordinates — every marker centre, curve
vertex and Bézier control point, exactly. A raster gives you an estimate of them. Five of the seven
figures here turned out to be vector art: three are read as geometry, and two more are vector but
rasterized on purpose (see below). Only korwek2023's panels are genuinely bitmaps.

| route | tool | use when | scripts |
|---|---|---|---|
| **SVG paths** | `pdftocairo -svg`, parse with `xml.etree` | vector art; you want dash patterns and stroke widths as series keys; stdlib only | mu2010 |
| **PDF drawings** | PyMuPDF `page.get_drawings()` | vector art; you want per-path fill/stroke colour and segment structure (marker shape) | michalski2012, rohrs2018 (Fig. S5) |
| **Embedded bitmap** | `pdfimages -png` | the panel *is* a bitmap (a blot, a screenshot, a figure pasted as an image) | korwek2023, rohrs2018 (Figs. 3, 5) |
| **Rasterized page** | `pdftoppm -r <dpi> -png` | vector art you have *chosen* to rasterize — see below | rana2020 (600 dpi), malleshaiah2010 (300 dpi) |

To decide: open the page's drawings and count stroked paths inside the panel. Dozens of small paths
means vector markers. One large image object means a bitmap, and `pdfimages` is then strictly
better than `pdftoppm` because it hands you the publisher's own pixels at their native resolution
(~230 ppi in korwek2023) instead of resampling them.

**Rasterizing vector art is allowed, but say why.** rana2020 renders at 600 dpi even though the
figure is vector, because five curves overlap and separating them by colour is far simpler than
disambiguating five overlapping path objects. That is a good reason. "I did not check" is not.
Whichever way you go, the choice belongs in the script docstring, because it sets the achievable
accuracy: a vector route has no digitization error in the geometry at all, and a raster route has
about a pixel.

## 2. Calibrate on the tick marks

**Fit the pixel-to-data map to the tick marks, never to the plot frame.** This is the single most
repeated rule in the corpus — five of seven scripts state it — and the reason is that a plot frame
is drawn wherever the plotting library likes. michalski2012's frame is offset from the axis limits
by several points; calibrating on it would have put every point in the paper wrong by a fixed
fraction of a decade, which looks like a real disagreement and is not.

Two acceptable shapes:

- **Fit the map to the ticks.** Least-squares a line through (tick pixel, tick value), in `log10`
  space for a log axis, and **report the residual**. michalski2012 prints a residual per panel in
  decades and in ordinate units. That residual is a free unit test: it is the calibration
  disagreeing with itself, and it should be near machine noise on a vector route.
- **Fit the map to the frame, then check it against the ticks.** mu2010 calibrates on the axis box
  and then requires the three interior x ticks and four interior y ticks to land on their nominal
  values within a hard tolerance, aborting if not. The ticks are still the authority; they are just
  being used as an assertion instead of as a fit.

The rule may be broken with evidence, and rohrs2018 shows how. In Fig. 5A MATLAB draws the ticks
*inside* the axes, where the plotted curves cover them, so the frame is the only measurable
landmark. The script uses the frame — and justifies it by measuring, in the vector Fig. S5 panels
of the same paper, that the outermost ticks sit on the frame to better than 0.01 pt. An exception
backed by a measurement in the same document is fine. An exception backed by "they look the same"
is how you get a systematic error.

**Log axes are where digitization error concentrates.** Calibrate in log space, plot in log space,
and check a decade landmark by hand. A 2 % error in pixel position on a linear axis is a 2 % error;
on a four-decade axis it is a factor of 1.2.

## 3. Separating series

**Vector.** Filter paths by any of: stroke colour, fill colour, stroke width, dash pattern, segment
count, segment kind. Round the colours before comparing — PDF colour components are floats and will
not compare equal to the value you read out of a legend (michalski2012 rounds to two decimals). The
richest keys in practice:

- *fill vs. stroke* separates filled from open markers of the same colour, which is how
  michalski2012 tells the 1 s series from the 6 s series.
- *segment count and kind* separates marker shapes — a circle is four curve segments, a square four
  lines — which is how michalski2012 recovers three phosphatase levels drawn in one colour.
- *dash pattern* is a series key in its own right: mu2010 identifies its four simulation conditions
  entirely from `stroke-dasharray`, because the caption maps dash style onto condition.
- *stroke width* separates data from decoration: in mu2010 the axis box and ticks are 0.5 pt and
  every trajectory is 0.75 pt.

**Raster.** Two techniques, and the second is better:

- *Near-colour distance.* `|pixel − colour|₁ < tol`. Simple, and adequate for saturated primaries
  on white.
- *Ink-direction classification* (rana2020). Anti-aliasing composites ink `C` over a white page, so
  an observed pixel is `α·C + (1−α)·255` for unknown coverage `α`. The vector `255 − observed` is
  therefore parallel to `255 − C` **whatever the coverage**. Normalize both and assign each pixel to
  the nearest direction. This is stable exactly where near-colour distance fails — on thin strokes,
  where most pixels are partial coverage — and it is the right default for a rasterized vector
  figure.

**When two series of the same colour overlap, there is no technique that separates them.** The
corpus contains three honest responses and you should pick one, not invent points:

- **Drop the series and use a different panel.** rohrs2018 abandons six mutant traces in Fig. 5A
  (drawn as three overlapping dashed pairs at 974 px panel width) and takes the same quantities from
  the Fig. 5B/C bar charts instead.
- **Emit nothing over the occluded interval.** rana2020 assigns each pixel to the nearer colour, so
  the hidden curve simply has no samples there, and the docstring names the interval.
- **Interpolate across a *known, bounded* obstruction and say so.** malleshaiah2010 interpolates the
  solid curve across the inset legend that overlays x = 367…430 px.

## 4. Reading a value out of ink

Once a mask isolates one series, a value is a summary of the ink in one column (or the position of
one glyph). Which summary depends on what is drawn:

| what is drawn | read it as | script |
|---|---|---|
| a thin curve | centre of the run of coloured pixels in the column | malleshaiah2010, rohrs2018 |
| a curve crossing other art | the run *nearest the previous column's* — follow continuity | rohrs2018 (Fig. 5A) |
| a filled marker | intensity-weighted centroid of the fattest run | korwek2023 |
| an open marker (a ring) | midpoint of the ring's vertical extent, not its centroid | korwek2023 |
| a bar | topmost row of the fill colour, over columns with enough height to be the bar | rohrs2018 |
| a scatter of round dots | connected components, filtered on area, aspect ratio and fill fraction | rohrs2018 |
| a vector marker of any shape | centre of the path's bounding box | michalski2012, rohrs2018 |
| a marker cloud with a fitted line inside it | per-column median, and report the cloud's half-height as the uncertainty | rana2020 |

**Cross-check the glyph against a landmark whose answer you know.** rohrs2018 found that the marker
paths in Fig. S5 carry a systematic offset of about 0.15 pt — roughly 2 % in time on that log axis —
while the error-bar segments land exactly on the axis limit for the t = 0.1 min point. It therefore
reads the data-point positions off the **error bars** (two vertical segments meeting at the mean,
two caps at mean ± SD, so the mean and the SD both come out exactly) and falls back to the marker
centre only where no error bar is drawn. A landmark you can predict — an axis limit, a point pinned
at zero, a value stated in the text — is worth more than a careful look at the marker.

## 5. Legend keys are inside the plot area

An in-panel legend draws its keys with the *same* path objects as the data, so every colour and
shape filter picks them up. They bite every route:

- michalski2012 drops them with `on_grid`: the paper samples on a 1-2-5 decade grid, legend keys sit
  at arbitrary abscissae, so any marker whose abscissa is not within tolerance of the grid is a key.
- rohrs2018's Fig. 5A curve walker has to survive a blue legend sample line in the upper left of the
  panel, which puts two blue runs in the same column; continuity picks the curve.
- malleshaiah2010 excludes the legend's pixel window outright and interpolates across it.

If a series comes out with one or two more points than the figure shows, suspect the legend before
suspecting the extraction.

## 6. What to record alongside the numbers

The CSV is read by someone who cannot see the figure. Give them what they need to use it correctly.

- **A provenance header.** Two or three `#` lines: the caption, the citation, and the script that
  wrote the file. rohrs2018 and rana2020 both do this.
- **A per-point uncertainty when the figure supports one.** rana2020 emits `spread_au`, the
  half-height of the marker cloud in each bin, which is the dominant error term and is typically
  0.02–0.05 ThT a.u. — far more useful than a single global claim.
- **Range flags, not silent clipping.** korwek2023 emits `on_axis_floor` and `above_axis` per value.
  The reason matters: on a log axis a figure cannot draw zero, so where the published model
  prediction is effectively zero the figure parks the marker at a small positive display value.
  **Those points are not model predictions and must not be compared as if they were** — the flag is
  what lets the verification notebook exclude them instead of fitting to an artifact of the plot.
- **The caveat that is currently in your head.** If the assignment of a marker shape to a condition
  was inferred rather than read from the legend, say how it was inferred: michalski2012 orders three
  unlabelled phosphatase levels by the physical argument that more phosphatase means less
  phosphorylation. If the caption and the body text disagree about which panel is which, record the
  contradiction and which one you followed, as rana2020 does at length for Fig. 1b/1c. That is a
  finding — see `when-the-paper-is-wrong.md`.

The digitization uncertainty you record here is what later justifies the tolerance in the
verification figure. `stochastic-verification.md` §4 puts a digitized curve at RMSE ≤ 0.06 of the
panel height for a theory comparison and ≤ 0.08 for a sampled one; those numbers are affordable
*because* line thickness and digitization error are of that size. A tolerance with no stated source
of error is not a tolerance.

## 7. Precision you can claim, and snapping

Do not emit more digits than the measurement supports. `%.6g` is the house format; it is well past
the precision of any raster route and costs nothing.

**Snapping an abscissa to a known grid is legitimate; snapping an ordinate is not.** When the paper
states its sampling grid — michalski2012 samples decades at 1, 2 and 5 — the true abscissa is known
exactly and the digitized one is known to be that value plus noise, so snapping removes error
rather than inventing precision. The ordinate is the quantity you are measuring, so it must keep
whatever noise the route gives it. Applying `snap_log` to a measured value would be fabricating
agreement.

The same asymmetry governs binning and decimation. rana2020 bins to a uniform 0.5 h grid and takes
the median within the bin; malleshaiah2010 keeps every fifth pixel column, "finer than the plotted
line width while avoiding pseudo-replication of every raster column". Both reduce a per-pixel
sampling that was never independent to begin with. Neither changes a value.

## 8. The reproducibility contract

- **The script is committed; the source PDF is not.** `dev/papers/` is outside the tree. So the
  script takes the PDF path as an argument with a default under `dev/papers/`, and when the file is
  missing it exits with a message that says where to get it — not a traceback. Every digitizer in
  the corpus does this.
- **Re-running regenerates the committed CSVs byte for byte.** This is the actual test of a
  digitizer, and it is worth running before you commit: `git diff --stat reference/` must come back
  empty. It catches the whole class of scripts that were run once, hand-edited, and can never be run
  again.
- **Name it `digitize_<author><year>.py`, next to the model, with `role: verification` in
  `metadata.yaml`.** See SKILL.md, "Committed Scripts".
- **Write every line ending as `\n`.** `csv.writer` defaults to CRLF. If you write a `#` provenance
  header with `fh.write` and the rows with a default `csv.writer`, the file gets LF headers and CRLF
  rows, and the `mixed-line-ending` pre-commit hook rejects it. Pass `lineterminator="\n"`. (Two
  committed CSVs — mu2010's and korwek2023's — are wholly CRLF for this reason and are left alone;
  they are consistent, so the hook accepts them. New work should be LF.)

## 9. The `digitize.py` library

`scripts/digitize.py` holds the machinery that was duplicated across the seven scripts: the four
source routes, tick-and-frame axis calibration with residuals, colour and ink-direction separation,
the run/blob/centroid readers of §4, the post-processing of §7, and a CSV writer that emits the
provenance header and the line endings of §8. Import it from a digitizer in a model folder:

```python
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "skills/curate-model/scripts"))
from digitize import Axis, page_raster, index_runs, write_csv
```

Run the module directly for a self-test that exercises every helper on synthetic input. numpy is
the only hard dependency: the SVG route needs nothing but the poppler CLI, and Pillow and PyMuPDF
are imported lazily for the raster and `get_drawings()` routes — `uv sync --group digitize`.

**The library must never change a committed digitized CSV.** It exists to stop the eighth digitizer
from re-deriving `index_runs` for the eighth time, not to reprocess the seven that already agree
with their figures. `tests/test_digitize.py` enforces that: it runs the self-test and replays both
ported digitizers (§10) against their committed bytes, skipping a replay when the uncommitted PDF
or an optional dependency is absent. Run it after any edit here, and if the PDFs are not on the
machine you are working on, re-run the digitizers by hand and confirm `git diff` comes back empty.

Reach for the library first, but do not contort a panel to fit it. Every figure has something
particular about it, and the particular part belongs in the model's own `digitize_*.py`, spelled
out, with the library supplying the parts that are not particular.

## 10. Worked examples

- **`four_flux_network_isotopomer_labeling_mu2010`** — the vector template, and stdlib-only.
  `pdftocairo -svg`, four trajectories keyed by dash pattern, frame calibration checked against
  seven interior ticks with a hard abort. **Ported to the library.**
- **`ste5_fus3_ptc1_switch_malleshaiah2010`** — the compact raster template. 300 dpi page render,
  red-channel threshold, two red curves separated by vertical position with a legend overlay
  interpolated across, log x and linear y from two landmarks each. **Ported to the library.**
- **`car_cd3zeta_phosphorylation_rohrs2018`** — the widest-range example: one paper digitized by
  three routes (vector Fig. S5, embedded bitmaps for Figs. 3 and 5). Read it for the error-bar
  trick of §4, the documented frame-instead-of-ticks exception of §2, the abandoned overlapping
  series of §3, and the empirical error estimate it gets for free — the same six data dots are
  repeated in all four Fig. 3 panels, so their spread across panels *is* the raster digitization
  error.
- **`camkii_holoenzyme_activation_michalski2012`** — vector marker recovery at its most demanding:
  five series × three unlabelled phosphatase levels, told apart by fill, segment count and segment
  kind, with calibration residuals printed per panel and legend keys filtered off the sampling grid.
- **`innate_immune_response_korwek2023`** — the embedded-bitmap route, and the reference for honest
  range flags: two-decade log panels, open model rings vs. filled blot dots, `on_axis_floor` and
  `above_axis` emitted per value so the notebook can exclude what the figure could not draw.
- **`amyloid_beta_competing_aggregation_pathways_rana2020`** — ink-direction classification, a
  per-point uncertainty column, and a caption-vs-body-text contradiction resolved from the kinetics
  and recorded in full.
- **`pybnf-jobs/Kirsch-2020/ppatf2_phospho`** — the honest floor. Four points per condition read by
  eye against the gridlines, transcribed as literals in a committed script, with the tolerance
  (±0.03 µM, the marker half-width) and the consequence (unweighted `sos`, since the panel draws no
  error bars) both stated. Reading by eye is acceptable when the panel is small and clean; hiding
  that you did is not.
