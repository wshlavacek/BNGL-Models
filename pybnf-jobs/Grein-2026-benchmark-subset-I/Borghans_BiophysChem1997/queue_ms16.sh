#!/bin/bash
# The m=16 rung of the ladder sweep, queued behind the m=8 arm. Same ten radius-0.4 starts,
# so 4 / 8 / 16 are paired and gntr's results carry across all three.
#
# Budget 900 s rather than 600: a 16-8-4-2-1 ladder is five rungs against three, and the
# question is where the trend TURNS OVER, not which config fits in ten minutes. No m=8 run
# hit its 600 s cap (longest was 400 s), so that arm ran to convergence and this one must be
# given the same chance -- otherwise a budget-truncated m=16 would look like a bad ladder
# rather than a slow one. Whether any run hits 900 s is checked afterwards.
set -u
cd "$(dirname "$0")"
while pgrep -f "run_paired_oscillating" > /dev/null; do sleep 30; done
echo "=== $(date -Iseconds) ms_segments = 16 arm ==="
exec python3 -u run_paired_oscillating.py --starts radius_starts.json --tag rad16 \
    --ms-segments 16 --methods ms --budget 900 2>> rad_stderr.log
