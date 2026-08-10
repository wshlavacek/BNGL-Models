# Rahman_MBS2016

**Run cost: `minutes`** — 100,000 evaluations (100 × 1,000 `gntr`) on a 16-reaction ODE model.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**SOLVED** — `OG = 0.000000` from a from-scratch 100-start `gntr` fit in 16 min 01 s. The fit lands
**on** the reference optimum, not merely inside the threshold: `J_paper` and `J*` agree to all seven
reported decimals. The PEtab nominal point is *also* this problem's published optimum
(`OG_nominal = 3.9e-06`), so the objective is validated independently of the optimizer.

This is the collection's **unit-σ** case, and its cleanest fidelity check: every measurement carries
`_SD = 1`, so `Σ log σᵢ` is exactly zero — 23 × `log(1)` — and the restored constant is `(N/2)log(2π)`
alone. The identity `J_paper == −lnL` holds with no σ bookkeeping to get wrong, so a discrepancy here
could only come from the likelihood itself. See `VALIDATION.md`.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `21.153486102003267` |
| paper-scale NLL at the PEtab nominal point | `21.153490000479476` |
| optimality gap at nominal | `3.8984762085192415e-06` |
| scored data points `n` | 23 |
| free parameters `k` | 9 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through trf's
Coleman–Li core, ADR-0068). `population_size = 100`, `max_iterations = 1000` — the collection's
documented working default, not the shipped 20 × 500 placeholder. Rahman was not separately tried at
20 × 500, so it is not itself evidence about the smaller budget.

`tools/fd_check.py` verifies the assembled gradient against central differences across all 9 free
parameters at `2.3e−05`. That matters because a wrong gradient column is silent — the objective stays
correct and the fit merely stops short — which is the failure mode that cost
`Laske_PLOSComputBiol2019` a solved verdict until lanl/PyBNF#534.

## Contents

- `Rahman_MBS2016.conf` — the PyBNF job
- `model_Rahman_MBS2016.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment1.exp` — experimental data
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `score.py` — scores a run against `J*`
- `best_fit_params.txt`, `information_criteria.txt` — the shipped fit's provenance
- `VALIDATION.md` — the full validation against `J*`

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`. The
converter preserves both `parameterScale` (lanl/PyBNF#491) and `observableTransformation`
(lanl/PyBNF#499), which plain `petab.v2.petab1to2` drops. The run recipe (`job_type`,
`sbml_backend = bngsim`, `wall_time_sim`) is supplied, not recovered — PEtab specifies a
problem, not a method. `wall_time_sim = 10` caps pathological parameter points; raise it
if valid simulations on your machine are being marked as failures.

## Running

```bash
pybnf -c Rahman_MBS2016.conf -o
python score.py output
```
