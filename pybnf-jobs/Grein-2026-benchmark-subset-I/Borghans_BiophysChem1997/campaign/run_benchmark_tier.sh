#!/bin/bash
# The lanl/PyBNF#563 acceptance benchmark, one tier.
#
# Four arms across the SAME fixed seeds, run one at a time. Sequential on purpose:
# the benchmark reports wall time, and two fits sharing ten cores would measure the
# scheduler rather than the methods -- and unequally, since arms 1/2/4 spend most of
# their budget on three dask workers while arm 3 drives its own search on one.
#
#   ./run_benchmark_tier.sh <tier-suffix> <budget-seconds> <seed,seed,...>
#
# e.g.  ./run_benchmark_tier.sh ''   1500  41,42,43,44,45,46,47,48   # baseline scale
#       ./run_benchmark_tier.sh '_b' 15000 41,42                     # the 10x scaling row
#
# Running
#   Run from anywhere; it cd's to the job directory itself. It execs the pybnf entry
#   point, taken from PYBNF_BIN (exported by .envrc.local) and falling back to PATH --
#   bare `pybnf` is not on PATH while this repo's venv is the active one.
#
#       campaign/run_benchmark_tier.sh <tier>
set -u
CAMPAIGN="$(cd "$(dirname "$0")" && pwd)"
# Run from the job directory: the confs name their model and .exp data relative to it.
cd "$CAMPAIGN/.."

TIER="${1:-}"
BUDGET="${2:-1500}"
SEEDS="${3:-41,42,43,44,45,46,47,48}"

for arm in a1_cmaes a2_cmaes_gntr a3_ms a4_cmaes_ms; do
  label="bench_${arm%%_*}${TIER}"
  echo "=== $(date -Iseconds) arm ${arm} (${label}), budget ${BUDGET}s, seeds ${SEEDS}"
  python3 "$CAMPAIGN/run_overnight_campaign.py" \
      --slot 1 --workers 3 \
      --arm "${label}" \
      --template "Borghans_bench_${arm}.conf" \
      --seeds "${SEEDS}" \
      --budget "${BUDGET}"
done
echo "=== $(date -Iseconds) tier '${TIER}' complete"
python3 "$CAMPAIGN/collect_benchmark.py" --tier "${TIER}"
