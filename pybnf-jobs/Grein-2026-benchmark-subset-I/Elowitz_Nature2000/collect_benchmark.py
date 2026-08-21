#!/usr/bin/env python3
"""Assemble the lanl/PyBNF#563 acceptance-benchmark table from finished runs.

Reads the per-arm status TSVs ``run_overnight_campaign.py`` wrote in arm mode, puts every
objective through the job's own ``score.py`` (so ``OG`` means exactly what it means for
every other run in this directory -- never a hand-rolled objective), and prints the four
quantities the issue asks for:

    success rate, wall time / model evaluations, continuity norm, final single-shoot
    objective.

Two things it reports that a naive table would not:

* the **executed** method chain per run, not the requested one. The 15-run baseline asked
  for ``cmaes,gntr`` and got ``cmaes`` (#564); a benchmark that reported the request would
  have been wrong about what it measured.
* the **tail**, not just the median. The #563 prototype's paired sweeps moved the tail and
  did not move the median (24-24 over 48 starts, medians tied at every radius), so a table
  summarised by its median would report a null result on the axis the method was never
  claimed to improve.

Usage:  python collect_benchmark.py [--tier A]
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
JSTAR = float((HERE / "jstar.txt").read_text().split()[0])
SOLVED = 1.92                      # chi^2, alpha = 0.05, 1 dof -- the Grein threshold

ARMS = [
    ("bench_a1", "1. BIPOP-CMA-ES"),
    ("bench_a2", "2. BIPOP-CMA-ES + GNTR"),
    ("bench_a3", "3. multistart MS-GNTR"),
    ("bench_a4", "4. BIPOP-CMA-ES + MS-GNTR"),
]


def og(output: Path) -> float | None:
    """The optimality gap, via the job's own score.py -- the same scale every other number
    in this directory is on."""
    proc = subprocess.run([sys.executable, str(HERE / "score.py"), str(output)],
                          cwd=HERE, capture_output=True, text=True)
    if proc.returncode != 0:
        return None
    for line in proc.stdout.splitlines():
        if line.startswith("OPTIMALITY GAP"):
            return float(line.split("=")[-1])
    return None


def rows_for(arm: str) -> list[dict]:
    out = []
    for path in sorted(HERE.glob(f"{arm}_slot*_status.tsv")):
        with path.open() as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row.get("best_reduced") in (None, "", "missing"):
                    continue
                out.append(row)
    return out


def summarize(arm: str, label: str) -> dict | None:
    rows = rows_for(arm)
    if not rows:
        return None
    gaps, seconds, sims, norms, chains = [], [], [], [], set()
    for row in rows:
        gap = og(HERE / row["output"])
        if gap is None:
            continue
        gaps.append(gap)
        seconds.append(int(row["seconds"]))
        try:
            sims.append(int(row["simulations"]))
        except (TypeError, ValueError):
            pass
        chains.add(row.get("executed", "?"))
        if row.get("continuity_norm", "-") not in ("-", "", "missing"):
            norms.append(float(row["continuity_norm"]))
    if not gaps:
        return None
    return {
        "label": label,
        "n": len(gaps),
        "solved": sum(1 for g in gaps if g < SOLVED),
        "best": min(gaps),
        "median": statistics.median(gaps),
        "worst": max(gaps),
        "seconds": statistics.median(seconds) if seconds else float("nan"),
        "sims": statistics.median(sims) if sims else float("nan"),
        "norm": max(norms) if norms else None,
        "chains": ", ".join(sorted(chains)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tier", default="", help="suffix on the arm names, e.g. '_b'")
    args = parser.parse_args()

    print(f"Elowitz_Nature2000 -- lanl/PyBNF#563 acceptance benchmark")
    print(f"J* = {JSTAR:.6f}; solved iff OG < {SOLVED}\n")
    head = ("arm", "n", "solved", "OG best", "OG med", "OG worst", "median s",
            "median sims", "cont. norm", "executed")
    print("| %-26s | %2s | %6s | %8s | %8s | %8s | %8s | %11s | %9s | %-12s |" % head)
    print("|" + "|".join(["-" * w for w in (28, 4, 8, 10, 10, 10, 10, 13, 11, 14)]) + "|")
    any_row = False
    for arm, label in ARMS:
        s = summarize(arm + args.tier, label)
        if s is None:
            print("| %-26s | %2s | %6s | %8s | %8s | %8s | %8s | %11s | %9s | %-12s |"
                  % (label, "-", "not run", "-", "-", "-", "-", "-", "-", "-"))
            continue
        any_row = True
        norm = "-" if s["norm"] is None else ("%.2g" % s["norm"])
        print("| %-26s | %2d | %4d/%-2d | %8.3f | %8.3f | %8.3f | %8.0f | %11.0f | %9s | %-12s |"
              % (s["label"], s["n"], s["solved"], s["n"], s["best"], s["median"],
                 s["worst"], s["seconds"], s["sims"], norm, s["chains"]))
    if not any_row:
        print("\n(no completed runs yet)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
