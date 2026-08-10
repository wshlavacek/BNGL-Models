# Weber_BMC2015

**Run cost: `minutes`** — 10,000 evaluations (20 × 500 `gntr`), 36 free parameters.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**Objective validated at the PEtab nominal point** (`OG = -0.0002`, well inside the solved threshold
1.92). **No optimization run has been performed here.** The problem's `nominalValue` point *is* its
published optimum — reproduced to four decimals — so the nominal check validates the import and
PyBNF's objective against the paper's Eq. 6 NLL; it makes no claim about PyBNF's optimizer. This is a
ready-to-run job.

> **Corrected 2026-08-07** by [lanl/PyBNF#547](https://github.com/lanl/PyBNF/issues/547) (ADR-0104). This problem pre-equilibrates, and the bngsim SBML backend silently dropped `preequilibrate:`. Its pre-equilibration condition is all-zeros and matches the model's authored defaults, so both experiments simulated a *completely flat* trajectory — identical to 8 significant figures at every timepoint, including across `t = 24` where `PdBu_time = 24, PdBu_dose = 1` should fire — and this README recorded `OG = 13739.87` as "the nominal point is not the published optimum". It is the optimum; the forward model was wrong. **This slug was previously queued as a tuning candidate on that reading; it is not one.** A trajectory that never moves through an event that should fire is the signature to check for if this ever recurs.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `296.2020024574204` |
| paper-scale NLL at the PEtab nominal point | `296.20179656464325` |
| optimality gap at nominal | `-0.00020589277715998833` |
| scored data points `n` | 135 |
| free parameters `k` | 36 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through trf's
Coleman–Li core, ADR-0068) — handles this problem's estimated noise scale, which plain `trf` refuses.
The shipped recipe was verified to start and run on this problem; it has not been run to completion.
Note that the conf still carries the 20 × 500 placeholder, not the collection's 100 × 1000 working
default — raise it before spending a fitting budget here.

## Contents

- `Weber_BMC2015.conf` — the PyBNF job
- `model_Weber_BMC2015.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment__model1_data1___model1_data2.exp`, `experiment__model1_data1___model1_data2_rep2.exp`, `experiment__model1_data1___model1_data2_rep3.exp`, `experiment__model1_data1___model1_data2_rep4.exp`, `experiment__model1_data1___model1_data2_rep5.exp`, `experiment__model1_data1___model1_data3.exp` … — experimental data
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
pybnf -c Weber_BMC2015.conf -o
python score.py output
```
