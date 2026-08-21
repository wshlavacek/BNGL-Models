#!/bin/bash
# Queue the ms_segments = 8 arm behind the running radius-0.4 experiment.
#
# Only `ms` is re-run: gntr has no segment count, so its radius-0.4 results from the
# ms_segments = 4 pass apply unchanged and the pairing is preserved across both arms.
#
# Budget 600 s rather than 300: an 8-4-2-1 ladder has one more rung than 4-2-1, and the
# question is whether m=8 SOLVES, not whether it solves inside 300 s. The ms/4 runs
# converged in 102-126 s, well inside their own budget, so neither arm is budget-limited
# and the comparison stays honest.
set -u
cd "$(dirname "$0")"
while pgrep -f "run_paired_oscillating" > /dev/null; do sleep 30; done
echo "=== $(date -Iseconds) starting ms_segments = 8 arm ==="
exec python3 -u run_paired_oscillating.py \
    --starts radius_starts.json --tag rad8 \
    --ms-segments 8 --methods ms --budget 600
