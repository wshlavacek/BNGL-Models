#!/bin/bash
# Reproduce the prototype's ONE solve through PyBNF's own job_type = ms.
#
# Same start (seed 3, radius 0.4, the prototype's own start_point), same ladder (8-4-2-1),
# same budget (600 s), sigma profiled in both. If PyBNF lands near OG = -1.2827 the
# implementation is faithful and every earlier experiment of mine was simply pointed at a
# regime of my own invention. If it does not, the implementation differs from the prototype
# and that is the thing to chase.
set -u
cd "$(dirname "$0")"
while pgrep -f "run_paired_oscillating" > /dev/null; do sleep 20; done
echo "=== $(date -Iseconds) prototype seed-3 reproduction, ladder 8-4-2-1, budget 600 ==="
exec python3 -u run_paired_oscillating.py --starts prototype_seed3_start.json \
    --tag proto8 --ms-segments 8 --methods ms,gntr --budget 600 2>> rad_stderr.log
