# Brannmark_JBC2010

**Run cost: `minutes`** — 10,000 evaluations (20 × 500 `gntr`), 22 free parameters.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**Setup only — not fitted.** The job runs and scores correctly, and the PEtab nominal point reproduces the published optimum (`OG = 0.064`, well inside the solved threshold), so the nominal check validates PyBNF's objective against the paper's Eq. 6 NLL.

> **Corrected 2026-08-07** by [lanl/PyBNF#547](https://github.com/lanl/PyBNF/issues/547) (ADR-0104). This problem pre-equilibrates, and the bngsim SBML backend silently dropped `preequilibrate:` — all eight doses simulated identically at the model's authored `insulin_dose_1 = 0.3`, and this README recorded `OG = 1531.44` as "the nominal point is not the published optimum". It is the optimum; the forward model was wrong. A dose-response whose doses all coincide is the signature to check for if this ever recurs.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `141.82485427243665` |
| paper-scale NLL at the PEtab nominal point | `141.88922297826417` |
| optimality gap at nominal | `0.06436870582751908` |
| scored data points `n` | 43 |
| free parameters `k` | 22 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through trf's Coleman–Li core, ADR-0068) — handles this problem's estimated noise scale, which plain `trf` refuses. The shipped recipe was
verified to start and run on this problem.

## Contents

- `Brannmark_JBC2010.conf` — the PyBNF job
- `model_Brannmark_JBC2010.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment__Dose_0___Dose_0.exp`, `experiment__Dose_0___Dose_001.exp`, `experiment__Dose_0___Dose_01.exp`, `experiment__Dose_0___Dose_1.exp`, `experiment__Dose_0___Dose_10.exp`, `experiment__Dose_0___Dose_100.exp` … — experimental data
- `experiment__Dose_0___Dose_100_measparams.tsv`, `experiment__Dose_0___TwoSteps_measparams.tsv` — per-measurement observable/noise parameter tables (ADR-0075)
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
pybnf -c Brannmark_JBC2010.conf -o
python score.py output
```
