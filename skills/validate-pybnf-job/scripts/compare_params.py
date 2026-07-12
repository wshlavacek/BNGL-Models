#!/usr/bin/env python
"""Gate 3b (parameter recovery): compare a PyBNF fit's estimates against the paper's reported values.

Published values are read from a small JSON you fill from the paper's table, e.g.:

    { "K1": 9.2e-9, "K2": 4.83e-7, "K1prime": 0.1 }

The fit result is a PyBNF ``sorted_params_final.txt`` (header ``# Simulation Obj <p1> <p2> ...``,
rows sorted best-first) OR any whitespace table with a header naming the parameters; the BEST row
(lowest Obj / first data row) is used. Parameter names are matched after stripping a ``__FREE``
suffix (case-insensitive).

For each shared parameter it prints published, recovered, the log10 ratio, and whether it is within
tolerance (a FACTOR, default 3x -> |log10 ratio| <= log10(3)). Reports how many landed within tol.
Rate-constant fits are often sloppy/non-identifiable -- a few params far off at a good objective is
expected; the reproduced curve (Gate 3a) and the identifiable combinations are the real check.

    PY=~/Code/PyBNF/.venv/bin/python   # any python with numpy
    $PY compare_params.py <published.json> <sorted_params_final.txt> [--factor 3]

Exit 0 on a successful comparison; 1 on a read/parse error.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


def _strip_free(name: str) -> str:
    n = name.strip()
    return n[:-6] if n.lower().endswith("__free") else n


def read_best_fit(path: Path) -> dict[str, float]:
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    header = None
    best = None
    for ln in lines:
        if ln.lstrip().startswith("#"):
            header = ln.lstrip("#").split()
            continue
        if header is None:
            header = ln.split()  # first line was the (uncommented) header
            continue
        best = ln.split()  # first data row = best (sorted best-first)
        break
    if header is None or best is None:
        raise ValueError(f"{path}: could not find a header + a best-fit data row")
    names = [_strip_free(h) for h in header]
    vals = {}
    for h, v in zip(names, best):
        try:
            vals[h] = float(v)
        except ValueError:
            continue  # non-numeric column (e.g. 'Simulation' label)
    return vals


def compare(pub_path: Path, fit_path: Path, factor: float) -> int:
    published = {(_strip_free(k)): float(v) for k, v in json.loads(pub_path.read_text()).items()}
    fit = read_best_fit(fit_path)
    lo = {k.lower(): k for k in fit}
    tol_log = math.log10(factor)

    print(f"ok    published params: {len(published)} from {pub_path.name}")
    print(f"ok    fit best-row params: {len(fit)} from {fit_path.name}  (tol = {factor}x)\n")
    print(f"  {'param':14s} {'published':>13s} {'recovered':>13s} {'log10 ratio':>12s}  within {factor}x?")
    within = shared = 0
    for p, pv in published.items():
        key = lo.get(p.lower())
        if key is None:
            print(f"  {p:14s} {pv:13.4g} {'(not fit)':>13s}")
            continue
        rv = fit[key]
        shared += 1
        if pv > 0 and rv > 0:
            lr = math.log10(rv / pv)
            ok = abs(lr) <= tol_log
            within += ok
            print(f"  {p:14s} {pv:13.4g} {rv:13.4g} {lr:12.2f}  {'yes' if ok else 'NO'}")
        else:
            print(f"  {p:14s} {pv:13.4g} {rv:13.4g} {'n/a (<=0)':>12s}")
    print(f"\n{within}/{shared} shared parameters within {factor}x of the published value.")
    if shared and within < shared:
        print("Note: sloppy/non-identifiable fits legitimately leave some params far off at a good "
              "objective. Judge with Gate 3a (curve reproduction) + an identifiability note, not "
              "this table alone.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("published", type=Path, help="JSON of the paper's reported parameter values")
    ap.add_argument("fit", type=Path, help="PyBNF sorted_params_final.txt (or a headered table)")
    ap.add_argument("--factor", type=float, default=3.0, help="within-tolerance factor (default 3x)")
    args = ap.parse_args()
    for p in (args.published, args.fit):
        if not p.is_file():
            print(f"FAIL  no such file: {p}")
            return 1
    try:
        return compare(args.published, args.fit, args.factor)
    except Exception as exc:  # noqa: BLE001
        print(f"FAIL  {type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
