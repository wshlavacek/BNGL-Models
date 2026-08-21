# Multiple-shooting prototype — lanl/PyBNF#563

A deliberately standalone prototype, outside PyBNF's fit machinery, answering one
question before any of that machinery is built: **does multiple shooting enlarge the
useful convergence region on Borghans?**

It also produced this problem's first solve: `OG = -1.2827`, verified through PyBNF's
own objective and this directory's `score.py`.

## Contents

| file | what it is |
|---|---|
| `msproto.py` | the transcription: segment simulation, data + continuity residuals with analytic Jacobians, the augmented-Lagrangian outer loop, the segment-count homotopy, and the certified single-shoot reconstruction |
| `validate.py` | calibration + correctness gates (see below) — run this first |
| `experiment.py` | the paired convergence-region sweep (perturbations of the PEtab nominal point) |
| `box_experiment.py` | the acceptance-style multistart from the job's own fit box |
| `schedules.py` | the segment-count comparison (`2-1`, `3-1`, `4-2-1`, `8-4-2-1`) |
| `reproduce.py` | re-runs one start and writes its parameter vector |
| `analyze.py`, `analyze_sched.py` | summary tables |
| `solved_seed3_r0.4.json` | the solved parameter vector, with its profiled sigma |
| `verify_gntr.conf` | PyBNF config that re-scores that vector through PyBNF's own objective |
| `sweep.jsonl`, `sched.jsonl`, `box.jsonl` | raw results |

## Running

Everything runs against the parent directory's model and data, on the PyBNF venv:

```bash
BNGSIM_CODEGEN_CACHE_DIR=$PWD/cg ~/Code/PyBNF/.venv/bin/python validate.py
```

`validate.py` is the gate. It checks that the objective reproduces
`nominal_check.json`'s reduced objective, `J_paper` and `OG`; that an *m*-segment
transcription seeded from a continuous trajectory has zero continuity defect and
byte-comparable data residuals to single shooting; and that the analytic Jacobians of
both residual blocks match central finite differences.

The sweeps shard across processes:

```bash
for i in 0 1 2 3; do BNGSIM_CODEGEN_CACHE_DIR=$PWD/cg OMP_NUM_THREADS=1 ~/Code/PyBNF/.venv/bin/python experiment.py sweep.jsonl $i 4 0.1,0.2,0.4,0.8 12 & done; wait
```

To re-run the solved start (deterministic; ~100 s). It writes
`rerun_seed3_r0.4.json` and diffs it against the committed `solved_seed3_r0.4.json`,
which it never overwrites — and note that the solve is **not expected to reproduce**
on a different build of the ODE library:

```bash
BNGSIM_CODEGEN_CACHE_DIR=$PWD/cg ~/Code/PyBNF/.venv/bin/python reproduce.py 3 0.4 600
```

## What it measures, and what it does not

The convergence-region sweep draws starts by perturbing the **PEtab nominal point**.
That is privileged information, so those runs measure *basin size*, not benchmark
performance — the solve reported above is a basin measurement and is not an
acceptance-benchmark result. `box_experiment.py` is the uninformed framing: starts
drawn from the job's own `loguniform_var 0.001 100000` box with nothing seeded.

Every score in every table comes from discarding the auxiliary states and
re-simulating theta with ordinary single shooting, so a run that leaves continuity
unconverged scores as what it actually is.

## Objective

The PEtab problem's own — lognormal (log10-additive Gaussian) noise on
`Ca = Z_state*scale + offset` — with the noise scale **profiled out analytically from
the data term only**, per #563's formulation point 4 and PyBNF #562 / ADR-0108, so
continuity violation can never be absorbed into the reported sigma.
