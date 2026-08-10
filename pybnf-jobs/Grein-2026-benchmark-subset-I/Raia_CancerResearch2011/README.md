# Raia_CancerResearch2011

**Run cost: `hours`** — 100,000 evaluations (100 × 1,000 `gntr`), 39 free parameters.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**SOLVED** — `OG = 0.000009` from a from-scratch 100-start `gntr` fit, landing on `J*` to five decimal
places; 83 of 100 starts retired on `step is negligible`. At k = 39, n = 205 this is the collection's
largest 🟢 → ✅ conversion.

**The fit beats its own nominal point.** `OG_nominal = 0.78` was already inside the solved threshold,
so this problem could have been "converted" by simply holding the published parameter vector. It was
not: unbiased box-sampled starts found a point ~87,000× closer to `J*` than the published one. That is
the difference between showing the optimizer *holds* a known optimum and showing it *finds* one.

One caveat on the numbers: `sd_pJAK2_rel` lands on its lower bound of `1e-05`. A noise scale driven to
its floor means the corresponding observable is fitted to within the resolution the box permits, and
`J*` is matched with it there — so the reference optimum is the same constrained one. Any
profile-likelihood analysis should widen that box first. See `VALIDATION.md`.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `345.3097672874065` |
| paper-scale NLL at the PEtab nominal point | `346.0897848482896` |
| optimality gap at nominal | `0.7800175608830955` |
| scored data points `n` | 205 |
| free parameters `k` | 39 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through trf's
Coleman–Li core, ADR-0068) — handles this problem's prediction-dependent σ
(`sigma = prediction_formula …`, one of two such problems here with `Armistead`), which plain `trf`
refuses. `population_size = 100`, `max_iterations = 1000` — the collection's documented working
default, not the shipped 20 × 500 placeholder.

This slug is where **lanl/PyBNF#537 / ADR-0100** was found. Its gradient disagreed with central
differences by a factor of two on one column, and was recorded as *not reproducible* after five
attempts — because whether it fired depended on the evaluation point. The cause was an IC-seeding
parameter whose own sensitivity axis *is* the whole derivative, with the seeded contribution summed on
top, so the column read exactly double. Post-fix the worst relative error is `4.85e−05`.

**It also retires the corpus's last `k`-based cost projection.** Issue #38 once projected 47–76 hours
for this slug by extrapolating from k = 39; it ran in well under two hours on ten cores, alongside two
other fits. Cost tracks stiffness and model size, not parameter count.

## Contents

- `Raia_CancerResearch2011.conf` — the PyBNF job
- `model_Raia_CancerResearch2011.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment____model1_data1.exp`, `experiment____model1_data2.exp`, `experiment____model1_data3.exp`, `experiment____model1_data4.exp` — experimental data
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
pybnf -c Raia_CancerResearch2011.conf -o
python score.py output
```
