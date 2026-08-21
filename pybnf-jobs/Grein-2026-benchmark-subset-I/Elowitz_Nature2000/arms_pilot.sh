#!/bin/bash
# Do the four arms differ at the NATIVE box and a 300 s budget?
#
# Arm 1 is already measured there: 0/3, OG 53.0 / 54.2 / 72.5. So there is headroom -- the
# question is whether any other arm uses it. If one does, that setting IS the discriminating
# regime and no box-widening or budget escalation is needed. Same box, same budget, same
# three seeds as arm 1, so the comparison is paired.
set -u
cd "$(dirname "$0")"
echo "arm      seed   OG          solved   seconds"
for ARM in a2_cmaes_gntr a3_ms a4_cmaes_ms; do
  for S in 101 102 103; do
    OUT="output_pilot_${ARM}_s${S}"; CONF="Elowitz_pilot_${ARM}_s${S}.conf"
    rm -rf "$OUT"
    python3 make_box_conf.py --template "Elowitz_bench_${ARM}.conf" --decades 8 \
        --out "$CONF" --output-dir "$OUT" --seed "$S" --budget 300 > /dev/null
    T0=$(date +%s)
    pybnf -c "$CONF" -o -l "bnf_pilot_${ARM}_s${S}" -L critical > /dev/null 2>&1
    T1=$(date +%s)
    LINE=$(python3 score.py "$OUT" 2>/dev/null | grep -E "OPTIMALITY GAP|SOLVED|NOT solved")
    OG=$(echo "$LINE" | grep OPTIMALITY | sed 's/.*= *//')
    SOLVED=$(echo "$LINE" | grep -q "=> SOLVED" && echo yes || echo NO)
    printf "%-14s %4s   %-11s %-8s %ss\n" "$ARM" "$S" "${OG:-missing}" "$SOLVED" "$((T1-T0))"
  done
done
echo done
