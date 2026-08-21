#!/bin/bash
# Reproduce the prototype's ONE solve through PyBNF's own job_type = ms.
#
# Same start (seed 3, radius 0.4, the prototype's own start_point), same ladder (8-4-2-1),
# same budget (600 s), sigma profiled in both. If PyBNF lands near OG = -1.2827 the
# implementation is faithful and every earlier experiment of mine was simply pointed at a
# regime of my own invention. If it does not, the implementation differs from the prototype
# and that is the thing to chase.
#
# Running
#   Run from anywhere; it cd's to the job directory itself. It drives a python script
#   that imports pybnf, so it needs PyBNF's interpreter: PYBNF_PY, exported by
#   .envrc.local, else plain python3 (which fails once the run reaches verify_start).
#
#       campaign/run_prototype_repro.sh
set -u
CAMPAIGN="$(cd "$(dirname "$0")" && pwd)"
# Run from the job directory: the confs name their model and .exp data relative to it.
cd "$CAMPAIGN/.."
while pgrep -f "run_paired_oscillating" > /dev/null; do sleep 20; done
echo "=== $(date -Iseconds) prototype seed-3 reproduction, ladder 8-4-2-1, budget 600 ==="
exec "${PYBNF_PY:-python3}" -u "$CAMPAIGN/run_paired_oscillating.py" --starts prototype_seed3_start.json \
    --tag proto8 --ms-segments 8 --methods ms,gntr --budget 600 2>> rad_stderr.log
