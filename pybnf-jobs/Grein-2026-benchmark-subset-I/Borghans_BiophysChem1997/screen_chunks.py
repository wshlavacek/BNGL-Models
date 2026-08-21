#!/usr/bin/env python3
"""Drive find_oscillating_starts.py in killable chunks, and merge what they find.

One Borghans draw in roughly a hundred sends CVODE into an integration that never returns.
The model class is ``BngsimSbmlModelNoTimeout`` -- it has no wall-clock guard -- and an
in-process ``SIGALRM`` cannot rescue it either: bngsim releases the GIL inside the solve, so
the signal is recorded but the Python handler does not run until the C call returns, which is
the thing that never happens. The first screen burned 41 minutes on one such point.

Process isolation is the only reliable answer. Each chunk is a separate interpreter with its
own seed and a hard subprocess timeout; because the screener writes its keeps incrementally,
a chunk that is killed mid-integration still contributes everything it found before hanging.

Usage:  python screen_chunks.py --want 8 --chunk 150 --timeout 45
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCREEN = HERE / "find_oscillating_starts.py"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--want", type=int, default=8)
    parser.add_argument("--chunk", type=int, default=150, help="draws per chunk")
    parser.add_argument("--timeout", type=float, default=45.0, help="seconds per chunk")
    parser.add_argument("--chunks", type=int, default=40, help="max chunks to try")
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--out", default="oscillating_starts.json")
    args = parser.parse_args()

    merged, screened, killed = [], 0, 0
    for index in range(args.chunks):
        if len(merged) >= args.want:
            break
        part = HERE / f"_chunk_{index:02d}.json"
        part.unlink(missing_ok=True)
        cmd = [sys.executable, "-u", str(SCREEN),
               "--n", str(args.chunk), "--keep", str(args.want),
               "--seed", str(args.seed + 1000 * index), "--out", part.name,
               "--report-every", str(max(args.chunk // 2, 1))]
        try:
            subprocess.run(cmd, cwd=HERE, timeout=args.timeout,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                           check=False)
            outcome = "ok"
        except subprocess.TimeoutExpired:
            killed += 1
            outcome = "killed (hung integration)"
        screened += args.chunk
        found = []
        if part.exists():
            try:
                found = json.loads(part.read_text())
            except Exception:
                found = []
        merged.extend(found)
        print(f"chunk {index:2d}  seed {args.seed + 1000 * index}  {outcome:26s}  "
              f"+{len(found)}  total {len(merged)}", flush=True)

    merged = merged[: args.want]
    (HERE / args.out).write_text(json.dumps(merged, indent=2) + "\n")
    print(f"\n{len(merged)} oscillating starts, from ~{screened} draws across "
          f"{index + 1} chunk(s), {killed} of which hung and were killed")
    for i, item in enumerate(merged, 1):
        p = item["params"]
        print(f"  start {i}: {item['peaks']:2d} peaks, amplitude {item['amplitude']:.2f}, "
              f"relative residual {item['relative_residual']:.2f}, "
              f"scale {p['scale']:.3g}, offset {p['offset']:.3g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
