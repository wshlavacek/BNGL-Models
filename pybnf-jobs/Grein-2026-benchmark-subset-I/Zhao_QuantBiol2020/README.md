# Zhao_QuantBiol2020

**Run cost: `minutes`** — 100,000 evaluations (100 × 1,000 `gntr`) on a 3-reaction ODE model.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**SOLVED** — `OG = 0.000005` from a from-scratch 100-start `gntr` fit in ≈54 min: the fit lands on
`J*` to five decimal places. The PEtab nominal point is **not** this problem's optimum
(`OG_nominal = 276.12`), so this result comes entirely from the optimizer. The objective is
independently corroborated against upstream's own `simulatedData` (`obj ✓`). See `VALIDATION.md`.

Three of the 28 parameters — `R_Stage_I_Wuhan`, `R_Stage_I_Hubei`, `R_Stage_I_China` — land on their
upper bound of 100. That is the problem's shape rather than a fit pathology: Stage I is the
uncontrolled early-epidemic window, where growth rate saturates in `R` and the likelihood goes nearly
flat. `J*` is matched to 5e−06 with them at the wall, so the reference optimum is the same
constrained one.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `501.22705376318333` |
| paper-scale NLL at the PEtab nominal point | `777.3442620361733` |
| optimality gap at nominal | `276.11720827298996` |
| scored data points `n` | 82 |
| free parameters `k` | 28 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through trf's
Coleman–Li core, ADR-0068) — handles this problem's per-measurement estimated σ, which plain `trf`
refuses. `population_size = 100`, `max_iterations = 1000` — the collection's documented working
default, not the shipped 20 × 500 placeholder.

### The search scale was the whole difficulty

All 28 estimated parameters are `parameterScale = log10` upstream, and every one of them imported as
a *linear* `uniform_var` until **lanl/PyBNF#548** was fixed. Holding everything else constant — same
recipe, same seed, same budget, same objective — and changing only whether they are sampled in log
space:

| configuration | best reduced objective | `OG` |
|---|---|---|
| linear `uniform_var` (before #548), 100 × 1000 | 717.5 and decelerating, run abandoned | ≥ 291 |
| **log `loguniform_var` (after #548), 100 × 1000** | **425.874099** | **0.000005** |

`gamma_*` live on `[1e-08, 1]` with an optimum near 0.05–0.39, and the `sd_*` on `[0.001, 1e5]` with
MLEs of 186–5013; sampled linear-uniform, effectively no box-sampled start begins near the basin. The
defect was silent — the objective, the `obj ✓` oracle check and the finite-difference gradient check
all still passed — and its only symptom was a fit that stalled, which is indistinguishable from a
problem needing a larger budget.

## Contents

- `Zhao_QuantBiol2020.conf` — the PyBNF job
- `model_Zhao_QuantBiol2020.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment____model1_China_Stage_I.exp`, `experiment____model1_China_Stage_II.exp`, `experiment____model1_Hubei_Stage_I.exp`, `experiment____model1_Hubei_Stage_II.exp`, `experiment____model1_Wuhan_Stage_I.exp`, `experiment____model1_Wuhan_Stage_II.exp` … — experimental data
- `experiment____model1_China_Stage_II_measparams.tsv`, `experiment____model1_China_Stage_I_measparams.tsv`, `experiment____model1_Hubei_Stage_II_measparams.tsv`, `experiment____model1_Hubei_Stage_I_measparams.tsv`, `experiment____model1_Wuhan_Stage_III_measparams.tsv`, `experiment____model1_Wuhan_Stage_II_measparams.tsv`, `experiment____model1_Wuhan_Stage_I_measparams.tsv` — per-measurement observable/noise parameter tables (ADR-0075)
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
pybnf -c Zhao_QuantBiol2020.conf -o
python score.py output
```
