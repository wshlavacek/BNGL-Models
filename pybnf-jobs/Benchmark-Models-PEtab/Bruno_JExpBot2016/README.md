# Bruno_JExpBot2016

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**Objective validated at the PEtab nominal point** (OG = 3.23e-06 < 1.92). No optimization run has been performed here.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `-46.688197918626756` |
| paper-scale NLL at the PEtab nominal point | `-46.688194686350265` |
| optimality gap at nominal | `3.232276490905406e-06` |
| scored data points `n` | 77 |
| free parameters `k` | 13 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = cmaes` — CMA-ES with IPOP restarts (ADR-0070/0082) — a global search, chosen because this problem is multimodal or refuses the gradient path. Note: fell back to cmaes: gntr: Error: Condition sets 'init_b10' to the value of free parameter 'init_b10_1' (a per-condition estimated initial condition, ADR-0076). The gradient path cannot yet route a free parameter that reaches the model through a c The shipped recipe was
verified to start and run on this problem.

## Contents

- `Bruno_JExpBot2016.conf` — the PyBNF job
- `model_Bruno_JExpBot2016.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment____model1_data1.exp`, `experiment____model1_data2.exp`, `experiment____model1_data3.exp`, `experiment____model1_data4.exp`, `experiment____model1_data5.exp`, `experiment____model1_data6.exp` — experimental data
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
pybnf -c Bruno_JExpBot2016.conf -o
python score.py output
```
