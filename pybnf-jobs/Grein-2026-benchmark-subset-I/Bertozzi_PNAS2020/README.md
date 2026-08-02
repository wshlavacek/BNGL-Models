# Bertozzi_PNAS2020

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**SOLVED** — `OG = 5.4e-06` from a from-scratch 20-start `gntr` fit, well inside the
threshold `OG < 1.92`. The PEtab nominal point is *also* this problem's published optimum
(`OG = 5.1e-06`), so the objective is validated independently of the optimizer.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `158.86426270904192` |
| paper-scale NLL at the fitted optimum | `158.8642678` |
| **optimality gap from the fit** | **`5.4e-06`** |
| paper-scale NLL at the PEtab nominal point | `158.86426780179107` |
| optimality gap at nominal | `5.09e-06` |
| scored data points `n` | 22 |
| free parameters `k` | 8 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — the general-objective Fisher/Gauss-Newton trust region (ADR-0068),
PyBNF's fides-analogue and the default for this collection. 20 box-sampled starts,
converged in about a minute.

### This slug used to ship `cmaes`, for two reasons that were both wrong

Until 2026-08-02 this job carried `job_type = cmaes` with a note that the gradient path had
refused it, and a nominal check recording `OG = 1.79e+11`. Two separate defects in PyBNF
produced that picture, and both are now fixed:

- **lanl/PyBNF#531 — the forward model was wrong.** This model derives its only infection
  rate through an SBML `initialAssignment` on a *parameter*, `beta_N = R0_*gamma_/N_`, and
  its conditions set all three of `R0_`, `gamma_` and `N_`. PyBNF's fast simulation path
  never recomputed such a derived parameter, so `beta_N` kept its load-time value: `R0_`
  was **completely inert** and the simulated trajectory was off by six orders of magnitude.
  That is a *scalar-path* defect — the `cmaes` recipe shipped here was fitting a broken
  model, which is where `OG = 1.79e+11` came from.
- **lanl/PyBNF#530 — the gradient path refused it at a boundary that has now moved.**
  `route_experiment` composed the chain rule for a condition-routed free parameter only
  where `d(IC)/d(target)` was a plain 1 (#511). Here `I0_` seeds two species with opposite
  signs (`I_ = I0_`, `S_ = N_ - I0_`), and `R0_`/`gamma_` reach the dynamics only through
  the derived `beta_N` — a case the refusal did not even cover, so those two columns were
  silently wrong rather than refused. Both now carry their real derivatives (ADR-0095), the
  second re-evaluated at each fit point.

`nominal_check.json` here was recomputed after both fixes; the value it carried before that
date is not comparable with anything.

## Contents

- `Bertozzi_PNAS2020.conf` — the PyBNF job
- `model.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment____u_CA.exp`, `experiment____u_NY.exp` — experimental data
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `best_fit_params.txt`, `information_criteria.txt` — the shipped fit's provenance
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
python score.py output          # or: python score.py   (scores the shipped provenance)
```
