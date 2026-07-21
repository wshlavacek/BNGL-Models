# Bertozzi_PNAS2020

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**Setup only — not fitted.** The job runs and scores correctly; the PEtab nominal point is not this problem's published optimum, so no optimality claim is made.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `158.86426270904192` |
| paper-scale NLL at the PEtab nominal point | `178639066015.4067` |
| optimality gap at nominal | `178639065856.54245` |
| scored data points `n` | 22 |
| free parameters `k` | 8 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = cmaes` — CMA-ES with IPOP restarts (ADR-0070/0082) — a global search, chosen because this problem is multimodal or refuses the gradient path. Note: fell back to cmaes: gntr: Error: Condition sets 'I0_' to the value of free parameter 'I0_CA' (a per-condition estimated initial condition, ADR-0076). The gradient path cannot yet route a free parameter that reaches the model through a condition t The shipped recipe was
verified to start and run on this problem.

## Contents

- `Bertozzi_PNAS2020.conf` — the PyBNF job
- `model.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment____u_CA.exp`, `experiment____u_NY.exp` — experimental data
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
pybnf -c Bertozzi_PNAS2020.conf -o
python score.py output
```
