#!/bin/bash
# Calibrate the benchmark's difficulty: how does BIPOP-CMA-ES's success rate decay as the
# prior box widens?
#
# The point is to find the width where arm 1 lands at an INTERMEDIATE success rate. At the
# native 8 decades this repo already solves Elowitz (OG 0.000175, gntr), so every arm would
# solve and the table would discriminate nothing; too wide and we are back in Borghans
# territory with 0/N everywhere. Somewhere between is the only regime where four methods can
# be told apart.
#
# Arm 1 only, because it is the reference the other three are judged against.
#
# Running
#   Run from anywhere; it cd's to the job directory itself. It execs the pybnf entry
#   point, taken from PYBNF_BIN (exported by .envrc.local) and falling back to PATH --
#   bare `pybnf` is not on PATH while this repo's venv is the active one.
#
#       campaign/calibrate_box.sh
set -u
CAMPAIGN="$(cd "$(dirname "$0")" && pwd)"
# Run from the job directory: the confs name their model and .exp data relative to it.
cd "$CAMPAIGN/.."
SEEDS="101 102 103"
BUDGET=300

echo "width  seed   OG          solved   seconds"
for D in 8 12 16 20; do
  for S in $SEEDS; do
    OUT="output_cal_d${D}_s${S}"
    CONF="Elowitz_cal_d${D}_s${S}.conf"
    rm -rf "$OUT"
    python3 "$CAMPAIGN/make_box_conf.py" --template Elowitz_bench_a1_cmaes.conf --decades "$D" \
        --out "$CONF" --output-dir "$OUT" --seed "$S" --budget "$BUDGET" > /dev/null
    T0=$(date +%s)
    "${PYBNF_BIN:-pybnf}" -c "$CONF" -o -l "bnf_cal_d${D}_s${S}" -L critical > /dev/null 2>&1
    T1=$(date +%s)
    LINE=$(python3 score.py "$OUT" 2>/dev/null | grep -E "OPTIMALITY GAP|SOLVED|NOT solved")
    OG=$(echo "$LINE" | grep OPTIMALITY | sed 's/.*= *//')
    SOLVED=$(echo "$LINE" | grep -q "=> SOLVED" && echo yes || echo no)
    printf "%5s  %4s   %-11s %-8s %ss\n" "$D" "$S" "${OG:-missing}" "$SOLVED" "$((T1-T0))"
  done
done
echo "done"
