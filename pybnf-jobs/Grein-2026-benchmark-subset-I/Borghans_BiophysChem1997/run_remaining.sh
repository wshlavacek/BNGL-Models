#!/bin/bash
# Finish radius-0.4 at the default ladder (starts 6-10, which the displacement guard
# aborted on), then run the ms_segments = 8 arm over all ten.
#
# stderr goes to a FILE, not /dev/null: discarding it to hide CVODE chatter is how the
# guard's abort message was lost and a half-finished run looked like a completed one.
set -u
cd "$(dirname "$0")"
echo "=== $(date -Iseconds) radius-4, starts 6-10 ==="
python3 -u run_paired_oscillating.py --starts radius_starts_6_10.json --tag rad_b \
    --budget 300 2>> rad_stderr.log
echo "=== $(date -Iseconds) ms_segments = 8 arm, all 10 starts ==="
python3 -u run_paired_oscillating.py --starts radius_starts.json --tag rad8 \
    --ms-segments 8 --methods ms --budget 600 2>> rad_stderr.log
echo "=== $(date -Iseconds) done ==="
