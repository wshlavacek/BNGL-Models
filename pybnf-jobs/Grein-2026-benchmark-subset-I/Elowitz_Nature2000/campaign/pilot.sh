#!/bin/bash
# Headroom pilot: does this benchmark discriminate, or does everything solve / nothing solve?
# Three seeds, arms 1 and 3, short budget. If arm 1 solves 3/3 there is no room for arm 3 to
# show anything; if neither solves at all we are back in Borghans territory and should stop.
#
# Running
#   Run from anywhere; it cd's to the job directory itself. It execs the pybnf entry
#   point, taken from PYBNF_BIN (exported by .envrc.local) and falling back to PATH --
#   bare `pybnf` is not on PATH while this repo's venv is the active one.
#
#       campaign/pilot.sh
set -u
CAMPAIGN="$(cd "$(dirname "$0")" && pwd)"
# Run from the job directory: the confs name their model and .exp data relative to it.
cd "$CAMPAIGN/.."
for arm in a1_cmaes a3_ms; do
  python3 -u "$CAMPAIGN/run_overnight_campaign.py" --slot 1 --workers 3 \
      --arm "pilot_${arm%%_*}" --template "Elowitz_bench_${arm}.conf" \
      --seeds 101,102,103 --budget 300 2>> pilot_stderr.log
done
