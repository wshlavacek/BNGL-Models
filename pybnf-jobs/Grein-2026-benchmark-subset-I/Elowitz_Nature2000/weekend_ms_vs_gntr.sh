#!/bin/bash
# Decide ms vs gntr on Elowitz, paired, over a weekend.
#
# WHY THIS PAIR AND NOTHING ELSE
#   ms drives gntr's OWN trust-region runner (pybnf/shooting/solver.py, "the step math is
#   gntr's, unchanged and not separately tuned" -- ADR-0110), and both resolve their starts
#   through the same GradientOptimizer code, so at a given seed they receive BIT-IDENTICAL
#   start points (verified). The only difference between the two arms is the transcription.
#   Every other comparison available here changes the method as well and cannot attribute a
#   difference to multiple shooting.
#
# WHY THE METRIC IS BINARY
#   Elowitz's landscape is discrete attractors, not a smooth basin (VALIDATION.md):
#       OG ~ 0       the target, reached by ~1 start in 1,000
#       OG 2.4324    the PEtab nominal point, a proven box-constrained local optimum
#       OG 5.8147    the DOMINANT attractor of uniform multistart, 6 of 10 batches
#   So OG among failures indexes WHICH trap a run fell into, not how near it came. A method
#   that reliably reaches 5.81 is reliably trapped, not nearly successful. The criterion is
#   therefore `OG < 1.92`, and OG is recorded but not ranked on.
#
# WHY 100 STARTS AND max_iterations = 500
#   VALIDATION.md's solving recipe, so a start here converges the way one did there. Earlier
#   300 s runs truncated 200 starts into five minutes and told us nothing.
#
# FAIRNESS
#   gntr uses parallel_count = 3; ms drives its own search on one core (a segment is not a
#   PSet evaluation). Equal STARTS is the fair axis and is what is paired; equal wall time is
#   not, and both costs are recorded so the asymmetry is visible rather than hidden.
#
# Append-only logging, one row per run, so the job is safe to stop at any moment and partial
# results are usable. Output dirs are pruned to Results/ after scoring to bound disk.
set -u
cd "$(dirname "$0")"

DEADLINE="${DEADLINE:-}"                 # ISO8601; empty = run until stopped
STARTS="${STARTS:-50}"                   # smaller batches: finer progress, less cap risk
MAXIT="${MAXIT:-500}"
CAP="${CAP:-14400}"                      # per-run safety cap, seconds.
# Generous on purpose: a truncated run silently breaks the comparison (starts stop
# converging), which is exactly what made the earlier 300 s runs uninformative. This is
# a runaway guard, not a budget -- max_iterations is what should end a run.
FIRST_SEED="${FIRST_SEED:-1000}"
LOG="${LOG:-weekend_ms_vs_gntr.tsv}"

[ -s "$LOG" ] || printf "seed\tmethod\tsolved\tOG\tbest_reduced\tsims\tseconds\tstage_trace\n" > "$LOG"

deadline_epoch=0
if [ -n "$DEADLINE" ]; then
  deadline_epoch=$(python3 -c "import datetime,sys;print(int(datetime.datetime.fromisoformat(sys.argv[1]).timestamp()))" "$DEADLINE")
fi

make_conf () {   # $1 method  $2 seed  $3 conf  $4 outdir
  local method="$1" seed="$2" conf="$3" out="$4"
  {
    echo "# Weekend ms-vs-gntr, paired on seed $seed. Identical starts across both methods."
    if [ "$method" = "ms" ]; then
      echo "job_type = ms"
      echo "ms_segments = 4"
      echo "ms_coarsening = 2"
      echo "ms_max_iterations = 25"
      echo "ms_inner_iterations = 50"
      echo "max_iterations = 25"
    else
      echo "job_type = gntr"
      echo "gntr_max_iterations = $MAXIT"
      echo "max_iterations = $MAXIT"
    fi
    echo "population_size = $STARTS"
    echo "refine = 0"
    echo "noise_profiling = 1"
    echo "parallel_count = 3"
    echo "wall_time_sim = 10"
    echo "wall_time_fit = $CAP"
    echo "verbosity = 1"
    echo "sbml_backend = bngsim"
    echo "random_seed = $seed"
    echo "output_dir = $out"
    echo
    sed -n '/^edition = 2/,$p' Elowitz_Nature2000.conf
  } > "$conf"
}

seed=$FIRST_SEED
while true; do
  if [ "$deadline_epoch" -ne 0 ] && [ "$(date +%s)" -ge "$deadline_epoch" ]; then
    echo "=== deadline reached, stopping at seed $seed ==="; break
  fi
  for METHOD in gntr ms; do
    OUT="output_wk_${METHOD}_s${seed}"; CONF="Elowitz_wk_${METHOD}_s${seed}.conf"
    rm -rf "$OUT"
    make_conf "$METHOD" "$seed" "$CONF" "$OUT"
    T0=$(date +%s)
    pybnf -c "$CONF" -o -l "bnf_wk_${METHOD}_s${seed}" -L critical >/dev/null 2>&1
    T1=$(date +%s)
    LINE=$(python3 score.py "$OUT" 2>/dev/null | grep -E "OPTIMALITY GAP|reduced objective|=> ")
    OG=$(echo "$LINE" | grep "OPTIMALITY GAP" | sed 's/.*= *//')
    RED=$(echo "$LINE" | grep "reduced objective" | sed 's/.*= *//' | awk '{print $1}')
    SOLVED=$(echo "$LINE" | grep -q "=> SOLVED" && echo YES || echo no)
    TRACE=$(grep -h "stage trace" "bnf_wk_${METHOD}_s${seed}.out" 2>/dev/null | tail -1 | sed 's/.*stage trace: //')
    SIMS=$(python3 - "$OUT" <<'PY' 2>/dev/null
import json,sys
try:
    d=json.load(open(sys.argv[1]+"/Results/method_chain.json"))
    print(sum(int(p.get("simulations") or 0) for p in d.get("phases",[])))
except Exception: print("?")
PY
)
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$seed" "$METHOD" "$SOLVED" \
      "${OG:-missing}" "${RED:-missing}" "${SIMS:-?}" "$((T1-T0))" "${TRACE:--}" >> "$LOG"
    # Keep the evidence, bound the disk: Results/ is what score.py and any audit need.
    rm -rf "$OUT/Simulations" "$OUT/FailedSimLogs" "$OUT"/*.bp 2>/dev/null
  done
  seed=$((seed+1))
done
echo "=== finished at $(date -Iseconds); $(($(wc -l < "$LOG")-1)) runs logged ==="
