#!/usr/bin/env python3
"""Keep one lane of native PyBNF Borghans campaigns busy.

Two modes, sharing one config-substitution, one status TSV and one results reader.

**Deadline mode** (the original, unchanged): full batteries are generation-limited
BIPOP-CMA-ES followed by GNTR, run back to back until a deadline.  When there is no
longer enough time for another full battery, run one explicitly CMA-ES-only tail with a
hard wall-time limit.  The tail is not labelled as a hybrid because PyBNF issue #564
means a budget-expired fit cannot refine.

**Arm mode** (``--template`` + ``--seeds``): run exactly the given seeds against one
arm's template, at one wall-clock budget.  This is what the lanl/PyBNF#563 acceptance
benchmark needs -- four arms across the *same* fixed seeds, so the arms are paired
rather than merely comparable -- and it reuses everything above it.  #564 is fixed
(ADR-0107's ``wall_time_refine_frac`` reserves the refine's share up front), so a
budgeted hybrid arm does now refine; ``Results/method_chain.json`` records whether it
actually did, and that is what the status TSV carries rather than an assumption.

Running
-------
Plain python3 is enough: this driver does not import pybnf, it execs the pybnf entry
point, taken from PYBNF_BIN and falling back to PATH. `.envrc.local` exports PYBNF_BIN;
without it, `pybnf` must be on PATH (activating PyBNF's venv does that).

    python3 -u run_overnight_campaign.py --slot 1 --workers 3 \
        --arm pilot_a1 --template Elowitz_bench_a1_cmaes.conf --seeds 101,102,103 --budget 300
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import subprocess
import time
import shutil
import os
from pathlib import Path


CAMPAIGN = Path(__file__).resolve().parent
# The job directory. This script lives in campaign/, but the model, the .exp data and
# every run happen one level up, and the confs' paths are relative to it.
HERE = CAMPAIGN.parent
BASE = CAMPAIGN / "Elowitz_bench_a1_cmaes.conf"
# Set PYBNF_BIN to the pybnf entry point; otherwise it is resolved on PATH, which is
# what activating the PyBNF environment gives you.
PYBNF = Path(os.environ.get("PYBNF_BIN") or shutil.which("pybnf") or "pybnf")


def replace_key(text: str, key: str, value: object) -> str:
    pattern = rf"(?m)^{re.escape(key)}\s*=.*$"
    replacement = f"{key} = {value}"
    if not re.search(pattern, text):
        raise RuntimeError(f"Template has no {key!r} line")
    return re.sub(pattern, replacement, text)


def make_config(slot: int, run_index: int, seed: int, workers: int,
                tail_budget: int | None, *, template: Path = None,
                arm: str = "overnight", budget: int | None = None) -> tuple[Path, Path, str]:
    template = template or BASE
    # A bare --template name means a template in campaign/, not one in the cwd the
    # driver runs pybnf from (the job directory).
    if not template.is_absolute() and not template.exists():
        template = CAMPAIGN / template.name
    kind = "tail-cmaes" if tail_budget is not None else arm
    stem = f"{arm}_s{slot}_n{run_index:02d}_seed{seed}"
    output = HERE / f"output_{stem}"
    conf = HERE / f"Elowitz_{stem}.conf"
    text = template.read_text()
    text = replace_key(text, "parallel_count", workers)
    text = replace_key(text, "random_seed", seed)
    text = replace_key(text, "output_dir", output.name)
    if budget is not None:
        # An arm template already carries a wall_time_fit line; the deadline-mode template
        # does not, which is what the tail branch below inserts one for.
        text = replace_key(text, "wall_time_fit", budget)
        note = (f"# Acceptance-benchmark run: arm {arm}, seed {seed}, "
                f"wall_time_fit = {budget} s.\n")
    elif tail_budget is not None:
        text = replace_key(text, "max_iterations", 100000)
        text = replace_key(text, "cmaes_restarts", 50)
        text = replace_key(text, "refine", 0)
        marker = f"wall_time_fit = {tail_budget}\n"
        text = text.replace("wall_time_sim = 10\n", marker + "wall_time_sim = 10\n")
        note = ("# Deadline tail: native BIPOP-CMA-ES only. GNTR is deliberately disabled\n"
                "# because a wall-time-expired refinement would be skipped (PyBNF#564).\n")
    else:
        note = ("# Full overnight battery: native BIPOP-CMA-ES followed by native GNTR.\n"
                "# No fit wall clock is set, so the requested refinement is guaranteed to start.\n")
    conf.write_text(note + text)
    return conf, output, kind


def best_objective(output: Path) -> str:
    results = output / "Results"
    candidates = (results / "sorted_params_refine_final.txt",
                  results / "sorted_params_final.txt")
    for path in candidates:
        if path.exists():
            for line in path.read_text().splitlines():
                if line.strip() and not line.startswith("#"):
                    return line.split()[1]
    return "missing"


def method_chain(output: Path) -> tuple[str, str]:
    """``(executed methods, total simulations)`` from ``Results/method_chain.json``.

    The executed chain rather than the requested one, because that distinction is the
    whole reason the file exists: the 15-run baseline asked for ``cmaes,gntr`` and got
    ``cmaes`` (#564), and a benchmark table that reported the request would have been
    wrong about what it measured.
    """
    path = output / "Results" / "method_chain.json"
    if not path.exists():
        return "missing", "missing"
    try:
        doc = json.loads(path.read_text())
    except Exception:
        return "unreadable", "unreadable"
    executed = "+".join(doc.get("executed_methods") or []) or "none"
    total = sum(int(p.get("simulations") or 0) for p in doc.get("phases", []))
    return executed, str(total)


def continuity_norm(output: Path) -> str:
    """The scaled continuity defect at the reported fit, or ``-`` for a run with no
    transcription (arms 1 and 2 never write the file)."""
    path = output / "Results" / "continuity_defects.txt"
    if not path.exists():
        return "-"
    for line in path.read_text().splitlines():
        if line.startswith("scaled_defect_norm_inf\t"):
            return line.split("\t", 1)[1].strip()
    return "-"


def record(status, fields) -> None:
    status.write("\t".join(str(x) for x in fields) + "\n")
    status.flush()


HEADER = ("run", "seed", "kind", "workers", "started", "ended", "seconds", "returncode",
          "best_reduced", "executed", "simulations", "continuity_norm", "output")


def run_one(status, slot, run_index, seed, workers, tail_budget, *, template=None,
            arm="overnight", budget=None) -> str:
    """Run one fit and append its row. Returns the best reduced objective, or ``missing``."""
    conf, output, kind = make_config(slot, run_index, seed, workers, tail_budget,
                                     template=template, arm=arm, budget=budget)
    prefix = conf.with_suffix("").name.replace("Elowitz_", "bnf_")
    stdout_path = HERE / f"{prefix}.out"
    started = time.time()
    started_iso = dt.datetime.now().astimezone().isoformat(timespec="seconds")
    with stdout_path.open("w") as stdout:
        proc = subprocess.run(
            [str(PYBNF), "-c", conf.name, "-o", "-l", prefix, "-L", "critical"],
            cwd=HERE, stdout=stdout, stderr=subprocess.STDOUT, check=False)
    ended = time.time()
    ended_iso = dt.datetime.now().astimezone().isoformat(timespec="seconds")
    best = best_objective(output)
    executed, simulations = method_chain(output)
    record(status, (run_index, seed, kind, workers, started_iso, ended_iso,
                    int(ended - started), proc.returncode, best, executed, simulations,
                    continuity_norm(output), output.name))
    return best


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--slot", type=int, required=True)
    parser.add_argument("--workers", type=int, required=True)
    parser.add_argument("--first-seed", type=int)
    parser.add_argument("--deadline",
                        help="Timezone-aware ISO deadline, e.g. 2026-08-13T08:00:00-06:00")
    parser.add_argument("--template", type=Path,
                        help="Arm mode: the conf template to run (Elowitz_bench_*.conf)")
    parser.add_argument("--arm", default="overnight",
                        help="Arm mode: label used in conf/output/status names")
    parser.add_argument("--seeds",
                        help="Arm mode: comma-separated seeds to run, in order")
    parser.add_argument("--budget", type=int,
                        help="Arm mode: wall_time_fit in seconds for every run")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    status_path = HERE / f"{args.arm}_slot{args.slot}_status.tsv"

    # --- arm mode: exactly these seeds, one budget, paired across arms ---------
    if args.seeds:
        if args.template is None or args.budget is None:
            parser.error("--seeds requires --template and --budget")
        seeds = [int(s) for s in args.seeds.split(",") if s.strip()]
        if args.dry_run:
            for i, seed in enumerate(seeds, 1):
                print(make_config(args.slot, i, seed, args.workers, None,
                                  template=args.template, arm=args.arm,
                                  budget=args.budget)[0])
            return 0
        with status_path.open("a", buffering=1) as status:
            if status_path.stat().st_size == 0:
                record(status, HEADER)
            for i, seed in enumerate(seeds, 1):
                run_one(status, args.slot, i, seed, args.workers, None,
                        template=args.template, arm=args.arm, budget=args.budget)
        return 0

    # --- deadline mode: unchanged ---------------------------------------------
    if args.deadline is None or args.first_seed is None:
        parser.error("deadline mode requires --deadline and --first-seed")
    deadline = dt.datetime.fromisoformat(args.deadline).timestamp()
    # Prior three-worker batteries took 69-78 minutes under simultaneous load.
    # Require two hours before launching another unbounded full battery.
    full_run_guard = 2 * 60 * 60
    # Stop the final wall-time tail ten minutes before 08:00 for PyBNF finalization.
    finish_reserve = 10 * 60

    if args.dry_run:
        conf, output, kind = make_config(args.slot, 1, args.first_seed, args.workers, None)
        print(conf)
        print(output)
        print(kind)
        return 0

    with status_path.open("a", buffering=1) as status:
        if status_path.stat().st_size == 0:
            record(status, HEADER)
        run_index = 1
        seed = args.first_seed
        while True:
            remaining = deadline - time.time()
            if remaining <= finish_reserve:
                break
            tail_budget = None
            if remaining <= full_run_guard:
                tail_budget = int(remaining - finish_reserve)
                if tail_budget < 60:
                    break
            best = run_one(status, args.slot, run_index, seed, args.workers, tail_budget)
            # A tail consumes the remaining search allocation by definition.
            if tail_budget is not None:
                break
            # Do not spin through seeds if a campaign failed before producing results.
            if best == "missing":
                break
            run_index += 1
            seed += 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
