#!/usr/bin/env python3
"""Rewrite an Elowitz arm's prior box to a given width, keeping everything else identical.

The prior box is the difficulty knob for an optimizer benchmark. Elowitz ships an 8-decade
box (`loguniform_var ... 1e-05 1000`) and this repo has already solved it there to
`OG 0.000175` with gntr -- so at the native width a four-arm comparison risks every arm
solving and discriminating nothing, which is as uninformative as Borghans' every-arm-fails.

Widening is the safe direction and that is why it is the knob used here: it only ADDS volume,
so the solution stays inside the box. Narrowing could exclude the solution, which would not be
calibration but rigging. Every arm is then run at the same width, and the width is reported
with the table -- "solved 6/10" means nothing without saying over what box.

The widening is symmetric in log10 about the native box's own centre, so the native box is
nested inside every wider one.

Running
-------
Plain python3. Needs no PyBNF environment -- it only reads files this job already
produced.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

CAMPAIGN = Path(__file__).resolve().parent
# The job directory. This script lives in campaign/, but the model, the .exp data and
# every run happen one level up, and the confs' paths are relative to it.
HERE = CAMPAIGN.parent
VAR = re.compile(r"^(loguniform_var\s*=\s*(\S+)\s+)(\S+)\s+(\S+)\s*$")

NATIVE_LO, NATIVE_HI = 1e-5, 1e3          # the shipped Elowitz box


def rewrite(text: str, decades: float) -> tuple[str, int]:
    import math
    centre = (math.log10(NATIVE_LO) + math.log10(NATIVE_HI)) / 2.0
    lo = 10.0 ** (centre - decades / 2.0)
    hi = 10.0 ** (centre + decades / 2.0)
    out, n = [], 0
    for line in text.splitlines():
        m = VAR.match(line.strip())
        if m:
            out.append(f"loguniform_var = {m.group(2)} {lo:.6g} {hi:.6g}")
            n += 1
        else:
            out.append(line)
    return "\n".join(out) + "\n", n


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", required=True)
    parser.add_argument("--decades", type=float, required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--budget", type=int, default=300)
    args = parser.parse_args()

    # Templates ship in campaign/; the cwd is the job directory.
    tpl = Path(args.template)
    if not tpl.is_absolute() and not tpl.exists():
        tpl = CAMPAIGN / tpl.name
    text = tpl.read_text()
    text, n = rewrite(text, args.decades)
    for key, value in (("random_seed", args.seed), ("output_dir", args.output_dir),
                       ("wall_time_fit", args.budget)):
        text = re.sub(rf"(?m)^{key}\s*=.*$", f"{key} = {value}", text)
    header = (f"# Calibration: prior box widened to {args.decades:g} decades per axis\n"
              f"# (native is 8). Widening only adds volume, so the solution stays inside.\n")
    (HERE / args.out).write_text(header + text)
    print(f"{args.out}: {n} parameters at {args.decades:g} decades, seed {args.seed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
