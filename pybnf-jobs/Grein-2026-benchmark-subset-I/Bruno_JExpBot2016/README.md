# Bruno_JExpBot2016

**Run cost: `minutes`** — 10,000 evaluations (20 × 500 `gntr`) on a 6-reaction ODE model.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**SOLVED** — `OG = 1.1e-05` from a from-scratch 20-start `gntr` fit in 41 s, well inside
the threshold `OG < 1.92`. The PEtab nominal point is *also* this problem's published
optimum (`OG = 3.23e-06`), so the objective is validated independently of the optimizer.
See `VALIDATION.md`.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `-46.688197918626756` |
| paper-scale NLL at the fitted optimum | `-46.68818673` |
| **optimality gap from the fit** | **`1.1e-05`** |
| paper-scale NLL at the PEtab nominal point | `-46.688194686350265` |
| optimality gap at nominal | `3.232276490905406e-06` |
| scored data points `n` | 77 |
| free parameters `k` | 13 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — the general-objective Fisher/Gauss-Newton trust region (ADR-0068),
PyBNF's fides-analogue and the default for this collection. 20 box-sampled starts,
converged in 41 s.

This slug shipped `job_type = cmaes` until **lanl/PyBNF#511** (merged #513, 2026-07-23),
with a note that the gradient path had refused it: all 13 free parameters reach the model
only through `condition:` parameter references (`init_b10 = init_b10_1`, a per-condition
estimated initial condition, ADR-0076), and `route_experiment` aborted on any such
perturbation rather than emit a silently-zero column. #511 taught it to compose the chain
rule instead, and the problem moved onto the gradient path and solved. **This README kept
describing the old recipe until 2026-08-02** — the conf itself has said `gntr` since #511.

## Contents

- `Bruno_JExpBot2016.conf` — the PyBNF job
- `model_Bruno_JExpBot2016.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment____model1_data1.exp`, `experiment____model1_data2.exp`, `experiment____model1_data3.exp`, `experiment____model1_data4.exp`, `experiment____model1_data5.exp`, `experiment____model1_data6.exp` — experimental data
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `best_fit_params.txt`, `information_criteria.txt` — the shipped fit's provenance
- `VALIDATION.md` — the full validation against `J*`
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
pybnf -c Bruno_JExpBot2016.conf -o
python score.py output          # or: python score.py   (scores the shipped provenance)
```
