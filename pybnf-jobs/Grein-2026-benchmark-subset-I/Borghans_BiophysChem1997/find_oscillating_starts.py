#!/usr/bin/env python3
"""Find prior-box draws that actually oscillate -- the starts a search discards.

The point of this script is the class of start it selects, and why that class is invisible
to every arm of #563's acceptance benchmark. On Borghans a correctly-shaped oscillator whose
period is wrong scores ~25 NLL units WORSE than a flat line, so a global search ranking by
the objective walks away from oscillating points and toward the no-dynamics region. The
retained trajectories of the completed baseline runs contain no oscillating point at all.

Multiple shooting's whole mechanism is that over one short segment a period error cannot
accumulate, so period information moves out of a saturated residual term and into continuity
defects, which carry a direction. That mechanism has content only where there IS a period to
get wrong. So the honest test of it needs starts of exactly the kind no search would hand it.

Draws are scored for oscillation on the observable's own trajectory (peaks + relative
amplitude), never against the data, so nothing here uses the answer.

Usage:  python find_oscillating_starts.py --n 400 --keep 12 --seed 20260814
"""

from __future__ import annotations

import argparse
import json
import signal
import os
import sys
import time
from pathlib import Path

import numpy as np

# Run with the PyBNF environment's interpreter, per this collection's convention.
# Set PYBNF_SRC to prepend a source checkout instead (that is how these runs were made).
_pybnf_src = os.environ.get("PYBNF_SRC")
if _pybnf_src:
    sys.path.insert(0, _pybnf_src)

HERE = Path(__file__).resolve().parent
CONF = HERE / "Borghans_bench_a3_ms.conf"

#: The benchmark box for the loguniform observation parameters.
BOX_LO, BOX_HI = 1e-3, 1e5


class _DrawTimeout(Exception):
    pass


def _alarm(signum, frame):
    raise _DrawTimeout()


def build():
    """The real Configuration + models the benchmark arms use, so a start drawn here is a
    start those arms could have drawn."""
    from pybnf.parse import load_config
    return load_config(str(CONF))


def observation_fit(pred, obs):
    """Least-squares ``(scale, offset)`` putting ``Ca = Z_state*scale + offset`` on the
    data's own scale, plus the resulting relative residual.

    Selecting draws on the *state* alone was not enough, and that was the first version's
    mistake: a draw can oscillate beautifully while ``scale`` and ``offset`` are eight
    decades out, and then the objective is dominated by the measurement model rather than by
    the dynamics. Start 1 of that screen scored **+726** where a flat line scores -165.98 and
    a wrong-period oscillator scores ~-141 -- so every method's fastest descent was to fix
    the observation parameters, and flattening the trajectory came along for free. Nothing
    about the transcription was being tested.

    This is a *starting value* only -- both stay free parameters of the fit -- and it uses
    the data's own scale, which any modeller has, rather than the fitted solution.
    """
    pred = np.asarray(pred, dtype=float)
    obs = np.asarray(obs, dtype=float)
    design = np.column_stack([pred, np.ones_like(pred)])
    try:
        (scale, offset), *_ = np.linalg.lstsq(design, obs, rcond=None)
    except np.linalg.LinAlgError:
        return None
    if not (np.isfinite(scale) and np.isfinite(offset)):
        return None
    # Both are loguniform_var over the benchmark's own [1e-3, 1e5] box, so a fit outside it
    # is not a start those arms could have taken. Silently keeping it would be worse than
    # dropping it: the conf carries the value as a truncated-normal mean, and a mean below
    # the lower bound puts the actual start at the BOUND rather than where it was intended --
    # a different point than the one that was screened, with nothing to show it changed.
    if not (BOX_LO <= scale <= BOX_HI and BOX_LO <= offset <= BOX_HI):
        return None
    residual = design @ [scale, offset] - obs
    return float(scale), float(offset), float(np.sqrt(np.mean(residual ** 2))
                                             / max(np.sqrt(np.mean(obs ** 2)), 1e-30))


def load_observations(path):
    """``(times, Ca)`` from the experiment file."""
    rows = [line.split() for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    arr = np.array([[float(x) for x in r[:2]] for r in rows], dtype=float)
    return arr[:, 0], arr[:, 1]


def oscillation(times, values):
    """``(n_peaks, relative amplitude)`` of one trajectory.

    A peak is a strict interior local maximum that rises at least 5 % of the trajectory's
    range above the lower of its two flanking minima -- enough to reject numerical ripple on
    a monotone or flat curve without tuning to any particular period.
    """
    values = np.asarray(values, dtype=float)
    if values.size < 5 or not np.all(np.isfinite(values)):
        return 0, 0.0
    span = float(np.max(values) - np.min(values))
    scale = max(abs(float(np.mean(values))), 1e-30)
    relative = span / scale
    if span <= 0.0:
        return 0, 0.0
    interior = np.flatnonzero((values[1:-1] > values[:-2]) & (values[1:-1] > values[2:])) + 1
    peaks = [i for i in interior
             if values[i] - min(np.min(values[:i]), np.min(values[i:])) > 0.05 * span]
    return len(peaks), relative


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=400, help="box draws to screen")
    parser.add_argument("--keep", type=int, default=12)
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--min-peaks", type=int, default=2)
    parser.add_argument("--min-amplitude", type=float, default=0.2)
    parser.add_argument("--out", default="oscillating_starts.json")
    parser.add_argument("--report-every", type=int, default=25)
    parser.add_argument("--draw-timeout", type=float, default=5.0,
                        help="seconds before a single draw's integration is abandoned")
    parser.add_argument("--max-relative-residual", type=float, default=0.9,
                        help="reject a draw whose best (scale, offset) still cannot get the "
                             "observable within this relative RMS of the data")
    args = parser.parse_args()

    config = build()
    variables = config.variables
    rng = np.random.default_rng(args.seed)
    model = list(config.models.values())[0]

    from pybnf.pset import PSet

    obs_times, obs_values = load_observations(HERE / "experiment1.exp")

    out_path = HERE / args.out
    kept, screened, integrated, stalled = [], 0, 0, 0
    started = time.time()

    def report():
        rate = f"{integrated / screened:.1%}" if screened else "-"
        hit = f"1 in {screened / len(kept):.0f}" if kept else f"0 in {screened}"
        print(f"  {screened:5d} draws  {integrated:5d} integrated ({rate})  "
              f"{stalled:3d} timed out  {len(kept):2d} oscillating ({hit})  "
              f"{time.time() - started:5.0f}s", flush=True)
        # Written as we go: a screen that has to be stopped early is still evidence, and
        # the *rate* is a finding in its own right (if oscillating points are 1-in-1000 of
        # this box, no arm of #563's benchmark can reach the regime the method operates in).
        out_path.write_text(json.dumps(kept, indent=2) + "\n")

    for _ in range(args.n):
        if len(kept) >= args.keep:
            break
        screened += 1
        if screened % args.report_every == 0:
            report()
        pset = PSet([v.set_value(v.value_from_quantile(rng.random()).value)
                     for v in variables])
        # The model class is BngsimSbmlModelNoTimeout -- it has no wall-clock guard, and a
        # pathological draw can sit inside CVODE indefinitely. The first version of this
        # screen burned 41 minutes on a single such point while reporting nothing. SIGALRM
        # works here because bngsim releases the GIL during integration, so the Python
        # handler actually gets to run.
        signal.signal(signal.SIGALRM, _alarm)
        signal.setitimer(signal.ITIMER_REAL, args.draw_timeout)
        try:
            model.param_set = pset
            # execute()'s third argument IS the per-simulation timeout, and it is only
            # applied when > 0. Passing 0 here (the first version did) makes every draw
            # unbounded, and roughly 1 in 100 Borghans draws then sits in CVODE forever --
            # which cost 41 minutes and was very nearly misfiled as a PyBNF defect. The
            # real fit path passes wall_time_sim; so does this now.
            data = model.execute(str(HERE / "_screen"), "screen", args.draw_timeout)
        except _DrawTimeout:
            stalled += 1
            continue
        except Exception:
            continue
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
        integrated += 1
        series = list(data.values())[0]
        column = "Ca" if "Ca" in series.cols else "Z_state"
        if column not in series.cols:
            continue
        arr = np.asarray(series.data, dtype=float)
        times = arr[:, series.cols[series.indvar]]
        values = arr[:, series.cols[column]]
        peaks, amplitude = oscillation(times, values)
        if peaks < args.min_peaks or amplitude < args.min_amplitude:
            continue
        # Put the observable on the data's scale before keeping the draw; see
        # observation_fit's docstring for why selecting on the state alone was not enough.
        fitted = observation_fit(np.interp(obs_times, times, values), obs_values)
        if fitted is None:
            continue
        scale, offset, relative = fitted
        if relative > args.max_relative_residual:
            continue
        params = {v.name: pset[v.name] for v in variables}
        params["scale"], params["offset"] = scale, offset
        kept.append({"params": params, "peaks": peaks, "amplitude": amplitude,
                     "relative_residual": relative})

    report()
    print(f"\nFINAL: screened {screened} draws; {integrated} integrated "
          f"({integrated / max(screened, 1):.1%}); {len(kept)} oscillate "
          f"(>= {args.min_peaks} peaks, relative amplitude >= {args.min_amplitude})")
    if kept:
        print(f"       oscillating draws are ~1 in {screened / len(kept):.0f} of this box")
    else:
        print(f"       NO oscillating draw in {screened} -- the regime multiple shooting "
              f"operates in is not reachable by uniform sampling of this prior box")
    for i, item in enumerate(kept, 1):
        print(f"  start {i:2d}: {item['peaks']} peaks, amplitude "
              f"{item['amplitude']:.2f}, relative residual "
              f"{item['relative_residual']:.2f}, scale {item['params']['scale']:.3g}, "
              f"offset {item['params']['offset']:.3g}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
