#!/bin/bash
# Headroom pilot: does this benchmark discriminate, or does everything solve / nothing solve?
# Three seeds, arms 1 and 3, short budget. If arm 1 solves 3/3 there is no room for arm 3 to
# show anything; if neither solves at all we are back in Borghans territory and should stop.
set -u
cd "$(dirname "$0")"
for arm in a1_cmaes a3_ms; do
  python3 -u run_overnight_campaign.py --slot 1 --workers 3 \
      --arm "pilot_${arm%%_*}" --template "Elowitz_bench_${arm}.conf" \
      --seeds 101,102,103 --budget 300 2>> pilot_stderr.log
done
