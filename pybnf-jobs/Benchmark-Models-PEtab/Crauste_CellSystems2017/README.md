# Crauste_CellSystems2017

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**Objective validated at the PEtab nominal point** (OG = 0.509 < 1.92). No optimization run has been performed here.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `190.45706550043687` |
| paper-scale NLL at the PEtab nominal point | `190.96591588401623` |
| optimality gap at nominal | `0.5088503835793574` |
| scored data points `n` | 21 |
| free parameters `k` | 12 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through trf's Coleman–Li core, ADR-0068) — handles this problem's estimated noise scale, which plain `trf` refuses. The shipped recipe was
verified to start and run on this problem.

## Contents

- `Crauste_CellSystems2017.conf` — the PyBNF job
- `model_Crauste_CellSystems2017.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment1.exp` — experimental data
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `score.py` — scores a run against `J*`

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`. The
converter preserves both `parameterScale` (lanl/PyBNF#491) and `observableTransformation`
(lanl/PyBNF#499), which plain `petab.v2.petab1to2` drops. The run recipe (`job_type`,
`sbml_backend = bngsim`, `wall_time_sim`) is supplied, not recovered — PEtab specifies a
problem, not a method. `wall_time_sim = 10` caps pathological parameter points; raise it
if valid simulations on your machine are being marked as failures.

## Running

```bash
pybnf -c Crauste_CellSystems2017.conf -o
python score.py output
```
