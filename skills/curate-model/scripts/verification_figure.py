#!/usr/bin/env python
"""Helpers for building `verify_<author><year>.png` to the contract in SKILL.md.

The verification figure has to answer one question — is the BioNetGen output correct? — from
the figure alone. That means four things, and this module supplies one helper for each:

    overlay_axis()    BioNetGen against a reference, in the house style (bngl/skill.md §2)
    residual_axis()   the disagreement, with the tolerance drawn as a band
    parity_axis()     many points collapsed onto an identity line, with the pooled metric
    verdict_footer()  the metrics and the PASS/FAIL, printed on the figure

`two_level_figure()` assembles the default layout: overlays on row one, residuals on row two,
metrics and verdict in the footer. Reach for it first; drop to the individual helpers when a
model needs a layout the default cannot express.

Nothing here decides whether a fit is good. You pass the tolerance and its justification; the
helpers make the comparison legible and label it honestly.

Use from a verification notebook in a model folder:

    import sys; sys.path.insert(0, "../../skills/curate-model/scripts")
    from verification_figure import Comparison, two_level_figure

    fig, checks = two_level_figure(
        t, bng_values,
        comparisons=[
            Comparison("independent SciPy ODE", scipy_values, tol=1e-4),
            Comparison("Fig. 2A (digitized)", paper_values, tol=0.15, x=paper_t),
        ],
        title="...", ylabel="...", xlabel="time (h)",
    )
    fig.savefig("verify_author2026.png", dpi=200, bbox_inches="tight")

Run this file directly for a synthetic self-test that exercises every helper.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# House overlay convention (bngl/skill.md §2): BNG solid, reference open markers.
BNG_KW = {"lw": 1.5, "zorder": 3}
REF_KW = {"mfc": "none", "ms": 5.5, "ls": "none", "zorder": 4}
MARKERS = ("o", "s", "^", "D", "v", "P")
CONVENTION = "Lines = BioNetGen; open markers = reference"


@dataclass
class Comparison:
    """One reference series to check BioNetGen against.

    `tol` is the acceptance threshold on the metric named by `kind`, and `why` is the
    justification that belongs in the figure or the notebook — data precision, figure
    resolution, digitization uncertainty. `x` defaults to the BioNetGen abscissa; supply it
    when the reference is sampled at different points (a digitized curve usually is).
    """

    label: str
    values: np.ndarray
    tol: float
    kind: str = "rel"          # "rel" -> relative error, "abs" -> absolute
    x: np.ndarray | None = None
    why: str = ""
    floor: float | None = None  # denominator floor for relative error
    stat: str = "max"          # "max" or "median" — which statistic the tol applies to


@dataclass
class Metrics:
    """Summary of one comparison. `passed` is measured against the Comparison's tol."""

    n: int
    max_rel: float
    median_rel: float
    max_abs: float
    rmse: float
    stat: str
    tol: float
    kind: str
    passed: bool
    label: str = ""
    detail: list[str] = field(default_factory=list)

    def line(self) -> str:
        value = self.median_rel if self.stat == "median" else self.max_rel
        if self.kind == "abs":
            value = self.max_abs
            shown = f"{self.stat} abs err {value:.3g}"
            tol = f"{self.tol:.3g}"
        else:
            shown = f"{self.stat} rel err {value:.3g}"
            tol = f"{self.tol:.3g}"
        return f"{self.label}: {shown} (tol {tol}) {'PASS' if self.passed else 'FAIL'}"


def _align(x_model, model, comparison):
    """Interpolate the model onto the reference abscissa when they differ."""
    model = np.asarray(model, dtype=float)
    ref = np.asarray(comparison.values, dtype=float)
    if comparison.x is None:
        if model.shape != ref.shape:
            raise ValueError(
                f"{comparison.label!r}: model has {model.shape} points and the reference "
                f"{ref.shape}; pass Comparison(x=...) so they can be aligned"
            )
        return np.asarray(x_model, dtype=float), model, ref
    xr = np.asarray(comparison.x, dtype=float)
    return xr, np.interp(xr, np.asarray(x_model, dtype=float), model), ref


def summarize(x_model, model, comparison: Comparison) -> Metrics:
    """Metrics for one comparison, on the points where the reference is finite.

    Relative error uses a denominator floor so a near-zero reference does not manufacture a
    huge error. The default floor is 1e-6 x max|reference|, matching the reproducibility rule
    in SKILL.md; pass `floor` to override.
    """
    _, m, r = _align(x_model, model, comparison)
    good = np.isfinite(m) & np.isfinite(r)
    m, r = m[good], r[good]
    if m.size == 0:
        raise ValueError(f"{comparison.label!r}: no finite points to compare")
    scale = comparison.floor
    if scale is None:
        scale = 1e-6 * float(np.max(np.abs(r))) if np.any(r) else 1.0
    denom = np.maximum(np.abs(r), scale)
    rel = np.abs(m - r) / denom
    abs_err = np.abs(m - r)
    stat_value = float(np.median(rel)) if comparison.stat == "median" else float(np.max(rel))
    if comparison.kind == "abs":
        stat_value = float(np.median(abs_err)) if comparison.stat == "median" else float(
            np.max(abs_err)
        )
    return Metrics(
        n=int(m.size),
        max_rel=float(np.max(rel)),
        median_rel=float(np.median(rel)),
        max_abs=float(np.max(abs_err)),
        rmse=float(np.sqrt(np.mean((m - r) ** 2))),
        stat=comparison.stat,
        tol=comparison.tol,
        kind=comparison.kind,
        passed=bool(stat_value <= comparison.tol),
        label=comparison.label,
    )


def overlay_axis(ax, x, model, comparisons, *, model_label="curated BioNetGen",
                 every=None, ylabel="", xlabel="", title="", logy=False):
    """BioNetGen as a solid line, each reference as open markers (bngl/skill.md §2).

    `every` subsamples the reference markers; the house default is about every 15th point, so
    that perfect agreement reads as markers sitting on the line rather than one curve hiding
    another.
    """
    x = np.asarray(x, dtype=float)
    ax.plot(x, model, "-", label=model_label, color="tab:blue", **BNG_KW)
    for i, comp in enumerate(comparisons):
        xr, _, r = _align(x, model, comp)
        step = every if every is not None else max(1, len(xr) // 15)
        ax.plot(xr[::step], r[::step], MARKERS[i % len(MARKERS)],
                label=comp.label, mec=f"C{i + 1}", **REF_KW)
    if logy:
        ax.set_yscale("log")
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, framealpha=0.9)
    return ax


def residual_axis(ax, x, model, comparison: Comparison, *, xlabel="", title=None):
    """The disagreement, with the tolerance as a shaded band.

    This is the panel that makes an overlay meaningful: two curves drawn on top of each other
    look identical whether the error is 1e-6 or 5%.
    """
    xr, m, r = _align(x, model, comparison)
    met = summarize(x, model, comparison)
    if comparison.kind == "abs":
        resid, ylabel = m - r, "absolute residual"
    else:
        scale = comparison.floor
        if scale is None:
            scale = 1e-6 * float(np.max(np.abs(r))) if np.any(r) else 1.0
        resid = (m - r) / np.maximum(np.abs(r), scale)
        ylabel = "relative residual"
    ax.axhspan(-comparison.tol, comparison.tol, color="tab:green", alpha=0.15,
               label=f"tolerance ±{comparison.tol:.3g}")
    ax.axhline(0.0, color="0.4", lw=0.8)
    ax.plot(xr, resid, ".", ms=4, color="tab:red", zorder=3)
    ax.set(xlabel=xlabel, ylabel=ylabel,
           title=title if title is not None else f"residual vs {comparison.label}")
    # Headroom above the data for the metric banner, so it cannot sit on the points; the
    # legend goes low, where a residual scatter centred on zero leaves room.
    span = max(float(np.max(np.abs(resid))), comparison.tol)
    ax.set_ylim(-1.45 * span, 2.1 * span)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="lower right", framealpha=0.9)
    ax.text(0.02, 0.97, met.line(), transform=ax.transAxes, fontsize=8, va="top",
            color="darkgreen" if met.passed else "darkred",
            bbox={"fc": "white", "ec": "darkgreen" if met.passed else "darkred",
                  "alpha": 0.95, "pad": 2.5})
    return ax


def parity_axis(ax, model, reference, *, groups=None, tol=None, label="", logscale=False,
                xlabel="reference", ylabel="BioNetGen"):
    """Every model value against every reference value, on the identity line.

    The aggregate panel required when a paper has many panels: it collapses what would be a
    survey of eyeball judgements into one readable claim. `groups` colours the points by
    panel/condition.
    """
    model = np.asarray(model, dtype=float).ravel()
    reference = np.asarray(reference, dtype=float).ravel()
    good = np.isfinite(model) & np.isfinite(reference)
    model, reference = model[good], reference[good]
    if groups is not None:
        groups = np.asarray(groups).ravel()[good]
        for i, g in enumerate(dict.fromkeys(groups.tolist())):
            sel = groups == g
            ax.plot(reference[sel], model[sel], MARKERS[i % len(MARKERS)], ms=4,
                    mfc="none", mec=f"C{i}", ls="none", label=str(g))
        ax.legend(fontsize=7, ncol=2, framealpha=0.9)
    else:
        ax.plot(reference, model, "o", ms=4, mfc="none", mec="tab:blue", ls="none")
    lo = float(min(reference.min(), model.min()))
    hi = float(max(reference.max(), model.max()))
    pad = 0.05 * (hi - lo or 1.0)
    line = np.array([lo - pad, hi + pad])
    ax.plot(line, line, "-", color="0.3", lw=1.0, zorder=1)
    if tol:
        ax.fill_between(line, line * (1 - tol), line * (1 + tol), color="tab:green",
                        alpha=0.13, zorder=0)
    if logscale:
        ax.set_xscale("log")
        ax.set_yscale("log")
    scale = 1e-6 * float(np.max(np.abs(reference))) if np.any(reference) else 1.0
    rel = np.abs(model - reference) / np.maximum(np.abs(reference), scale)
    note = f"n={model.size}  median rel err {np.median(rel):.3g}  max {np.max(rel):.3g}"
    ax.set(xlabel=xlabel, ylabel=ylabel, title=label or "aggregate agreement")
    ax.text(0.03, 0.96, note, transform=ax.transAxes, fontsize=8, va="top",
            bbox={"fc": "white", "ec": "0.7", "alpha": 0.9, "pad": 2.5})
    ax.grid(alpha=0.25)
    return ax


def verdict_footer(fig, metrics, *, extra=""):
    """Print each comparison's metric and PASS/FAIL along the bottom of the figure.

    A reader should never have to open the notebook to learn how well it agreed.
    """
    lines = [m.line() if isinstance(m, Metrics) else str(m) for m in metrics]
    ok = all(m.passed for m in metrics if isinstance(m, Metrics))
    text = "   ·   ".join(lines)
    if extra:
        text += f"\n{extra}"
    fig.text(0.5, 0.005, text, ha="center", va="bottom", fontsize=9,
             color="darkgreen" if ok else "darkred",
             bbox={"fc": "white", "ec": "darkgreen" if ok else "darkred", "pad": 4})
    return ok


def two_level_figure(x, model, comparisons, *, title="", ylabel="", xlabel="",
                     logy=False, convention=CONVENTION, figsize=None):
    """The default layout: overlays on row 1, residuals on row 2, verdict in the footer.

    Returns `(fig, metrics)`. Save with `fig.savefig(..., dpi=200, bbox_inches="tight")`.
    """
    import matplotlib.pyplot as plt

    n = len(comparisons)
    if n == 0:
        raise ValueError("a verification figure needs at least one comparison")
    figsize = figsize or (5.6 * n, 7.4)
    fig, axes = plt.subplots(2, n, figsize=figsize, squeeze=False,
                             gridspec_kw={"height_ratios": [2.0, 1.0]})
    metrics = []
    for j, comp in enumerate(comparisons):
        overlay_axis(axes[0][j], x, model, [comp], ylabel=ylabel, xlabel=xlabel,
                     title=f"BioNetGen vs {comp.label}", logy=logy)
        residual_axis(axes[1][j], x, model, comp, xlabel=xlabel)
        metrics.append(summarize(x, model, comp))
    header = title + (f"\n{convention}" if convention else "")
    fig.suptitle(header, fontsize=11)
    fig.tight_layout(rect=(0, 0.045, 1, 0.97))
    verdict_footer(fig, metrics)
    return fig, metrics


def _selftest() -> int:
    """Synthetic exercise of every helper; writes nothing outside the temp dir."""
    import tempfile
    from pathlib import Path

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = np.linspace(0, 10, 400)
    bng = 2.0 * np.exp(-0.4 * t) * np.sin(1.3 * t) + 1.0
    scipy_like = bng * (1 + 3e-7)                       # agrees to ~3e-7
    paper_t = np.linspace(0.2, 9.5, 14)
    paper = np.interp(paper_t, t, bng) * (1 + 0.06 * np.cos(paper_t))

    comps = [
        Comparison("independent SciPy ODE", scipy_like, tol=1e-4),
        Comparison("Fig. 2A (digitized)", paper, tol=0.15, x=paper_t, stat="median"),
    ]
    fig, mets = two_level_figure(t, bng, comps, title="self-test",
                                 ylabel="value (a.u.)", xlabel="time (h)")
    assert len(mets) == 2
    assert mets[0].passed and mets[0].max_rel < 1e-5, mets[0]
    assert mets[1].passed, mets[1]
    # identity: a series compared to itself is exactly zero error
    zero = summarize(t, bng, Comparison("self", bng, tol=0.0))
    assert zero.max_rel == 0.0 and zero.passed, zero
    # a comparison outside tolerance must FAIL, not quietly pass
    bad = summarize(t, bng, Comparison("bad", bng * 1.5, tol=0.01))
    assert not bad.passed and "FAIL" in bad.line(), bad
    # interpolation path: mismatched lengths without x= is an error, not a silent truncation
    try:
        summarize(t, bng, Comparison("mismatch", paper, tol=0.1))
    except ValueError as exc:
        assert "aligned" in str(exc)
    else:
        raise AssertionError("mismatched lengths should raise")

    fig2, ax = plt.subplots(figsize=(5, 5))
    parity_axis(ax, np.interp(paper_t, t, bng), paper, tol=0.15,
                groups=np.array(["S4A"] * 7 + ["S12B"] * 7), label="aggregate")
    with tempfile.TemporaryDirectory() as d:
        for i, f in enumerate((fig, fig2)):
            p = Path(d) / f"selftest_{i}.png"
            f.savefig(p, dpi=110, bbox_inches="tight")
            assert p.stat().st_size > 8000, p
    plt.close("all")
    print("verification_figure self-test: OK")
    print("  " + mets[0].line())
    print("  " + mets[1].line())
    return 0


if __name__ == "__main__":
    raise SystemExit(_selftest())
