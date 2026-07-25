# Schwen_PONE2014

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**Setup only** — the job imports, simulates, and scores. Its observables are **log10**
(`lognormal`), and for a log10 observable `-log_likelihood` carries the
change-of-variables Jacobian `Σ log(y_obs·ln10)`, which the paper's Eq. 6 `J*` need
not. The nominal optimality gap is therefore **recorded but not asserted** as a
validation — the same caveat that keeps the other log10 slugs (Elowitz, Borghans) out
of the "objective validated" tier. No optimizer claim is made.

At the PEtab nominal point PyBNF evaluates `−log_likelihood = 943.999`, which sits
just below the tightly-converged reference `J* = 952.422` (the best three Marvin runs
agree to 6+ digits), i.e. the nominal point lands in the reference basin.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `952.4217251306842` |
| paper-scale NLL at the PEtab nominal point | `943.9992974818219` |
| optimality gap at nominal | `-8.422427648862254` |
| scored data points `n` | 286 |
| free parameters `k` | 30 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (ADR-0068),
which handles this problem's estimated noise scales and the log10 (`lognormal`)
objective (EFIM Fisher block for noise scales, ADR-0079/0080/0081). This is a
**default recipe, not a tuned one.**

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`.
This problem exercises **lanl/PyBNF#510**: one experimental condition is measured only
at `t = 0`, which `TimeCourse` previously rejected — failing the entire problem at
config load. It now loads, simulates, and scores. The observables also carry
offset+scale `observableParameters` and per-observable estimated noise scales
(`IR_obs_std`, `std`), imported via the `measurement_params:` sidecars. The run recipe
(`job_type`, `sbml_backend = bngsim`, `wall_time_sim`) is supplied, not recovered.

## Contents

- `Schwen_PONE2014.conf` — the PyBNF job
- `model_Schwen_PONE2014.xml` — SBML model (emitted by the importer)
- `experiment____model1_data*.exp` (+ `_rep2`) — experimental data, 19 conditions
- `experiment____model1_data*_measparams.tsv` — per-measurement observable/noise parameter tables
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `score.py` — scores a run against `J*`

## Running

```bash
pybnf -c Schwen_PONE2014.conf -o
python score.py output
```
