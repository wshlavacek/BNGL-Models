#!/bin/bash
# Five arms x N seeds on Elowitz, native box.
#
# Arms 3 and 5 are the headline pair: ms drives gntr's OWN trust-region runner
# (pybnf/shooting/solver.py, unchanged and not separately tuned), so ms-vs-gntr at the same
# starts and budget isolates what the TRANSCRIPTION adds and nothing else. Every other
# contrast is confounded by a change of method.
set -u
cd "$(dirname "$0")"
SEEDS="${SEEDS:-101 102 103 104 105}"
echo "arm             seed   OG           solved  sims     seconds"
for ARM in a1_cmaes a5_gntr a3_ms a2_cmaes_gntr a4_cmaes_ms; do
  for S in $SEEDS; do
    OUT="output_bench_${ARM}_s${S}"; CONF="Elowitz_run_${ARM}_s${S}.conf"
    rm -rf "$OUT"
    python3 make_box_conf.py --template "Elowitz_bench_${ARM}.conf" --decades 8 \
        --out "$CONF" --output-dir "$OUT" --seed "$S" --budget 0 > /dev/null
    # --budget 0 would blank wall_time_fit; restore whatever the arm's own template set.
    grep -q '^wall_time_fit = 0$' "$CONF" && sed -i '' "s/^wall_time_fit = 0$/$(grep '^wall_time_fit' "Elowitz_bench_${ARM}.conf")/" "$CONF"
    T0=$(date +%s); pybnf -c "$CONF" -o -l "bnf_bench_${ARM}_s${S}" -L critical >/dev/null 2>&1; T1=$(date +%s)
    LINE=$(python3 score.py "$OUT" 2>/dev/null | grep -E "OPTIMALITY GAP|SOLVED|NOT solved")
    OG=$(echo "$LINE" | grep OPTIMALITY | sed 's/.*= *//')
    SOLVED=$(echo "$LINE" | grep -q "=> SOLVED" && echo YES || echo no)
    SIMS=$(python3 -c "
import json,sys
try:
    d=json.load(open('$OUT/Results/method_chain.json'))
    print(sum(int(p.get('simulations') or 0) for p in d.get('phases',[])))
except Exception: print('?')" 2>/dev/null)
    printf "%-15s %4s   %-12s %-6s  %-8s %ss\n" "$ARM" "$S" "${OG:-missing}" "$SOLVED" "$SIMS" "$((T1-T0))"
  done
done
echo done
