#!/usr/bin/env python3
"""Paired ms-vs-gntr from starts that oscillate -- the starts a search discards.

#563's acceptance benchmark asks four arms to fit Borghans from uninformed starts. On this
problem that measures the *global search*, which nothing in the multiple-shooting feature
was built to address: a correctly-shaped oscillator with the wrong period scores ~25 NLL
units WORSE than a flat line, so a search ranking by the objective is correctly pushed into
the no-dynamics region and away from the only region the transcription can convert. The
completed baseline runs' retained trajectories bear that out -- 5000 points each, every one
in the no-dynamics band.

So "refine the best fit" hands multiple shooting a flat trajectory, which is trivially
continuous: no period, therefore no continuity defects, therefore nothing for the
transcription to work with. This script tests the claim the issue's motivation actually
makes -- "multiple shooting can enlarge the useful convergence region" -- on the class of
start where it has content.

Both methods get the SAME start, the SAME box (the benchmark's own bounds, carried by a
bounded normal prior whose median IS the start point), the SAME budget and the SAME
objective. The only difference is the transcription.

Usage:  python run_paired_oscillating.py --budget 300
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
import time
import shutil
import os
from pathlib import Path

HERE = Path(__file__).resolve().parent
# Set PYBNF_BIN to the pybnf entry point; otherwise it is resolved on PATH, which is
# what activating the PyBNF environment gives you.
PYBNF = Path(os.environ.get("PYBNF_BIN") or shutil.which("pybnf") or "pybnf")
STARTS = HERE / "oscillating_starts.json"
TAG = "osc"

# The benchmark's own box, per parameter family.
LOG_BOUNDS = (0.001, 100000)
LIN_BOUNDS = (0.0, 1.0)
LINEAR = {"init_A_state", "init_Y_state", "init_Z_state"}

PROBLEM = """
edition = 2
model: model_Borghans_BiophysChem1997.xml
noise_model = lognormal, sigma = fit sigma
observable: Ca, formula: Z_state*scale + offset
experiment: experiment1, method: ode, data: experiment1.exp
"""

HEADER = """\
# Paired ms-vs-gntr from an OSCILLATING box draw (#563).
# The start is the prior median, so both methods begin at exactly the same point, inside
# exactly the benchmark's box. sd is irrelevant at population_size = 1 (no scatter is drawn).
"""


def parameter_lines(params: dict) -> list[str]:
    out = []
    for name, value in sorted(params.items()):
        # sd is deliberately TINY, and that is load-bearing rather than cosmetic. The start
        # is the box centre, which is the prior's MEDIAN -- and for a normal truncated to
        # [lower, upper] the median is not the mean whenever the truncation is asymmetric.
        # At sd = 0.2 on init_Z_state (mean 0.088, bounds [0, 1]) the left tail is cut at
        # -0.44 sd and the start lands at 0.173 -- a factor of 1.97 from the point that was
        # intended, silently. A negligible sd makes the truncation negligible, so median ==
        # mean and the fit starts where it was told. Nothing else reads sd here:
        # population_size = 1 draws no scatter, and priors do not enter an optimizer's
        # objective.
        if name in LINEAR:
            lo, hi = LIN_BOUNDS
            out.append(f"parameter: {name}, prior: normal, mean: {value!r}, sd: 1e-6, "
                       f"lower: {lo}, upper: {hi}")
        else:
            lo, hi = LOG_BOUNDS
            out.append(f"parameter: {name}, prior: normal, parameter_scale: log10, "
                       f"mean: {math.log10(value)!r}, sd: 1e-6, lower: {lo}, upper: {hi}")
    return out


def make_conf(index: int, method: str, params: dict, budget: int,
              ms_segments: int = 4) -> tuple[Path, Path]:
    stem = f"{TAG}{index:02d}_{method}"
    output = HERE / f"output_{stem}"
    conf = HERE / f"Borghans_{stem}.conf"
    body = [f"job_type = {method}", "population_size = 1", "parallel_count = 3",
            "noise_profiling = 1",
            f"wall_time_fit = {budget}", "wall_time_sim = 10", "verbosity = 1",
            "sbml_backend = bngsim", "random_seed = 1", f"output_dir = {output.name}"]
    if method == "ms":
        # ms_segments is the FINEST rung. The shipped default is 4, but the #563 prototype's
        # one and only solve of this problem came off an 8-4-2-1 ladder
        # (`m=8: -132.32  m=4: -166.52  m=2: -221.91  m=1: -248.07`), so replicating that
        # solve means replicating its ladder rather than the default.
        body += [f"max_iterations = {ms_segments * 6}", f"ms_segments = {ms_segments}",
                 "ms_coarsening = 2", "ms_max_iterations = 25",
                 "ms_inner_iterations = 50"]
    else:
        body += ["max_iterations = 2000", "gntr_max_iterations = 2000"]
    conf.write_text(HEADER + "\n".join(body) + "\n" + PROBLEM
                    + "\n".join(parameter_lines(params)) + "\n")
    return conf, output


def verify_start(conf: Path, params: dict, tolerance: float = 1e-4) -> float:
    """The largest relative gap between where the fit will start and where it was told to.

    A start is expressed as a bounded prior whose MEDIAN the box-start path takes, and a
    truncated normal's median is not its mean. That displaced an intended start by a factor
    of 1.97 once already, silently, and produced a whole experiment's worth of runs that did
    not begin where their table said they did. Checked every time now, and loudly, because
    the failure has no symptom of its own -- the fits run perfectly happily from the wrong
    point.
    """
    from pybnf.parse import load_config
    config = load_config(conf.name)
    worst, culprit = 0.0, None
    for v in config.variables:
        # Measured in SAMPLING space and normalised by the parameter's own box width. A
        # plain relative error is wrong here in both directions: it divides by zero for a
        # parameter legitimately clipped to a bound (a perturbation can put init_Z_state at
        # exactly 0, and 6.7e-7 / 0 "displaces" by 1e23), and it is scale-blind for a log
        # parameter, where what matters is the distance in decades rather than in value.
        # lower_bound / upper_bound, NOT p1 / p2: for a normal prior p1 and p2 are the
        # mean and the sd, so using them here made the "box width" 1e-6 and turned a 6.7e-7
        # displacement into 0.674.
        lo = v.to_sampling_space(v.lower_bound)
        hi = v.to_sampling_space(v.upper_bound)
        width = abs(hi - lo) or 1.0
        intended = v.to_sampling_space(params[v.name])
        actual = v.to_sampling_space(v.value_from_quantile(0.5).value)
        gap = abs(actual - intended) / width
        if gap > worst:
            worst, culprit = gap, v.name
    if worst > tolerance:
        raise SystemExit(
            f"start displaced by {worst:.3g} of the box width in {culprit} "
            f"(> {tolerance:g}) in {conf.name} -- the fit would not start where it was told")
    return worst


def score(output: Path) -> tuple[str, str]:
    """``(reduced objective, OG)`` through the job's own score.py."""
    proc = subprocess.run([sys.executable, str(HERE / "score.py"), str(output)],
                          cwd=HERE, capture_output=True, text=True)
    if proc.returncode != 0:
        return "missing", "missing"
    reduced = og = "missing"
    for line in proc.stdout.splitlines():
        # score.py appends a parenthetical gloss to the reduced-objective line, so take the
        # first token after the '=' rather than the rest of the line.
        if line.startswith("PyBNF reduced objective"):
            reduced = line.split("=")[-1].split()[0]
        if line.startswith("OPTIMALITY GAP"):
            og = line.split("=")[-1].split()[0]
    return reduced, og


def continuity(output: Path) -> str:
    path = output / "Results" / "continuity_defects.txt"
    if not path.exists():
        return "-"
    for line in path.read_text().splitlines():
        if line.startswith("scaled_defect_norm_inf\t"):
            return line.split("\t", 1)[1].strip()
    return "-"


def stage_trace(log: Path) -> str:
    """The last homotopy stage trace in a run's stdout, or ``-`` for a non-``ms`` run.

    The single most informative artifact a multiple-shooting run produces (ADR-0109): it
    shows whether the *coarsening* is converting the segmented stages, which is the mechanism
    the method rests on. On the prototype's solving run it reads
    ``m=8: -132.32   m=4: -166.52   m=2: -221.91   m=1: -248.07`` -- every segmented stage
    worse than a flat line (-165.98), and the coarsening cashing them in. A run whose whole
    trace sits inside the flat-line band never converted anything, which the final objective
    alone does not distinguish from a run that converted a little.
    """
    if not log.exists():
        return "-"
    found = "-"
    for line in log.read_text(errors="replace").splitlines():
        if "stage trace:" in line:
            found = line.split("stage trace:", 1)[1].strip()
    return found


def simulations(output: Path) -> str:
    path = output / "Results" / "method_chain.json"
    if not path.exists():
        return "?"
    try:
        doc = json.loads(path.read_text())
    except Exception:
        return "?"
    return str(sum(int(p.get("simulations") or 0) for p in doc.get("phases", [])))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--budget", type=int, default=300)
    parser.add_argument("--limit", type=int, default=0, help="only the first N starts")
    parser.add_argument("--starts", default=str(STARTS))
    parser.add_argument("--tag", default="osc", help="prefix for conf/output names")
    parser.add_argument("--ms-segments", type=int, default=4,
                        help="finest rung of the ms ladder (the prototype's solve used 8)")
    parser.add_argument("--methods", default="gntr,ms",
                        help="comma-separated methods to run for each start")
    args = parser.parse_args()

    global TAG
    TAG = args.tag
    starts = json.loads(Path(args.starts).read_text())
    if args.limit:
        starts = starts[: args.limit]
    if not starts:
        print("no oscillating starts found -- run find_oscillating_starts.py first")
        return 1

    status = HERE / f"paired_{TAG}_status.tsv"
    new = not status.exists() or status.stat().st_size == 0
    with status.open("a", buffering=1) as fh:
        writer = csv.writer(fh, delimiter="\t")
        if new:
            writer.writerow(("start", "peaks", "amplitude", "method", "seconds",
                             "reduced", "OG", "simulations", "continuity_norm",
                             "stage_trace", "output"))
        for i, item in enumerate(starts, 1):
            for method in [m.strip() for m in args.methods.split(",") if m.strip()]:
                conf, output = make_conf(i, method, item["params"], args.budget,
                                         ms_segments=args.ms_segments)
                verify_start(conf, item["params"])
                prefix = conf.with_suffix("").name.replace("Borghans_", "bnf_")
                started = time.time()
                with (HERE / f"{prefix}.out").open("w") as out:
                    subprocess.run([str(PYBNF), "-c", conf.name, "-o", "-l", prefix,
                                    "-L", "critical"],
                                   cwd=HERE, stdout=out, stderr=subprocess.STDOUT,
                                   check=False)
                reduced, og = score(output)
                row = (i, item["peaks"], f"{item['amplitude']:.3f}", method,
                       int(time.time() - started), reduced, og, simulations(output),
                       continuity(output), stage_trace(HERE / f"{prefix}.out"),
                       output.name)
                writer.writerow(row)
                print("  ".join(str(x) for x in row), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
