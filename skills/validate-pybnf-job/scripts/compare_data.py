#!/usr/bin/env python
"""Gate 1 (data provenance): diff a re-digitized/transcribed dataset against a shipped `.exp`.

Both files are read as tables whose FIRST column is the independent variable and whose remaining
columns are observables (named in a header line; a leading '#' and parentheses like ``FL()`` are
tolerated). Delimiter is auto-detected (comma or whitespace). ``_SD`` columns are ignored.

For every observable column present in BOTH files, the re-digitized series is interpolated onto the
shipped file's independent-variable grid (both sorted on x; comparison restricted to the overlap),
and the absolute and relative differences are reported (median + max). A large *systematic*
difference means the wrong panel/series, a missed normalization, or a units error -- not
digitization noise.

    PY=~/Code/PyBNF/.venv/bin/python   # any python with numpy
    $PY compare_data.py <redigitized.csv> <shipped.exp> [--rel-tol 0.15]

Exit 0 always on a successful comparison (PASS/FAIL is a human judgement); 1 on a read/parse error.
With --rel-tol, prints PASS/REVIEW per column as a convenience.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def _norm(name: str) -> str:
    return name.strip().lstrip("#").strip().replace("()", "").lower()


def read_table(path: Path):
    """Return (indep_name, x, {obs_name: y}) from a headered table (csv or whitespace)."""
    lines = [ln.rstrip("\n") for ln in path.read_text().splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"{path}: empty")
    header = lines[0].lstrip("#").strip()
    delim = "," if header.count(",") >= 1 else None  # None => any whitespace
    names = [h for h in (header.split(delim) if delim else header.split()) if h != ""]
    rows = []
    for ln in lines[1:]:
        if ln.lstrip().startswith("#"):
            continue
        parts = ln.split(delim) if delim else ln.split()
        rows.append([float(p) for p in parts])
    data = np.array(rows, dtype=float)
    if data.ndim == 1:
        data = data[None, :]
    cols = {names[i]: data[:, i] for i in range(min(len(names), data.shape[1]))}
    indep_name = names[0]
    x = data[:, 0]
    obs = {n: cols[n] for n in list(cols)[1:] if not _norm(n).endswith("_sd")}
    return indep_name, x, obs


def compare(redig: Path, shipped: Path, rel_tol: float | None) -> int:
    rn, rx, robs = read_table(redig)
    sn, sx, sobs = read_table(shipped)
    print(f"ok    read redigitized {redig.name}: indep='{rn}', {len(rx)} rows, obs={list(robs)}")
    print(f"ok    read shipped     {shipped.name}: indep='{sn}', {len(sx)} rows, obs={list(sobs)}")

    smap = {_norm(k): k for k in sobs}
    shared = [(rk, smap[_norm(rk)]) for rk in robs if _norm(rk) in smap]
    if not shared:
        print(f"FAIL  no observable columns in common (redigitized {list(robs)} vs shipped {list(sobs)})")
        return 0
    ro = np.argsort(rx)
    rx_s = rx[ro]

    allrel = []
    for rk, sk in shared:
        ry = robs[rk][ro]
        so = np.argsort(sx)
        xs, ys = sx[so], sobs[sk][so]
        lo, hi = rx_s.min(), rx_s.max()
        m = (xs >= lo) & (xs <= hi)
        if m.sum() == 0:
            print(f"warn  {sk}: no overlap in independent variable; skipped")
            continue
        interp = np.interp(xs[m], rx_s, ry)
        absd = np.abs(interp - ys[m])
        denom = np.abs(ys[m])
        good = denom > (0.02 * np.nanmax(denom) if np.nanmax(denom) > 0 else 1e-12)
        rel = absd[good] / denom[good]
        allrel.extend(rel.tolist())
        tag = ""
        if rel_tol is not None:
            tag = "  [PASS]" if (len(rel) and np.median(rel) <= rel_tol) else "  [REVIEW]"
        med = np.median(rel) if len(rel) else float("nan")
        mx = rel.max() if len(rel) else float("nan")
        print(f"  {sk:16s} n={m.sum():3d}  abs: med={np.median(absd):.4g} max={absd.max():.4g}  "
              f"rel: med={med:.3f} max={mx:.3f}{tag}")

    if allrel:
        print(f"\nOVERALL rel diff: median={np.median(allrel):.3f}  max={np.max(allrel):.3f}  "
              f"(over {len(allrel)} points with non-trivial magnitude)")
    print("\nInterpret: a few % on a clean linear panel = PASS; ~10-20% off a rasterized log-axis "
          "panel can still PASS if stated; a systematic offset / wrong shape = FAIL (wrong "
          "series/normalization/units). Record the verdict in VALIDATION.md.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("redigitized", type=Path, help="re-digitized / transcribed data (csv or whitespace)")
    ap.add_argument("shipped", type=Path, help="the job's shipped .exp")
    ap.add_argument("--rel-tol", type=float, default=None, help="flag columns whose median rel diff exceeds this")
    args = ap.parse_args()
    for p in (args.redigitized, args.shipped):
        if not p.is_file():
            print(f"FAIL  no such file: {p}")
            return 1
    try:
        return compare(args.redigitized, args.shipped, args.rel_tol)
    except Exception as exc:  # noqa: BLE001
        print(f"FAIL  {type(exc).__name__}: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
