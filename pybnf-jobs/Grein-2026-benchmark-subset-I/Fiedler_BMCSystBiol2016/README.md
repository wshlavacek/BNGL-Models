# Fiedler_BMCSystBiol2016

**Run cost: `hours`** — 100,000 evaluations (100 × 1,000 `gntr`), 22 free parameters. At 6 h 13 min
this is the most expensive fit in the subset.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**SOLVED, but not saturated** — `OG = 1.003516` from a from-scratch 100-start `gntr` fit in
6 h 13 min, inside the threshold `OG < 1.92`. **The two halves of that sentence should not be
conflated.** The χ² threshold makes `OG = 1.004` statistically indistinguishable from the reference
optimum, so the ✅ is correct as scored — but the optimizer did not find the reference basin:

| point | OG |
|---|---:|
| PEtab nominal (published optimum) | −0.0022 |
| **this fit, 100 × 1000 unbiased** | **1.0035** |
| threshold | 1.92 |

The reference basin demonstrably exists and is reachable — the nominal point sits in it — and 100
box-sampled starts did not find it. Progress was still positive but nearly exhausted at the budget
ceiling (1.113 → 1.034 → 1.004 over the last ~90 minutes, grinding rather than converging). This is
the first problem in the corpus where the collection's 100 × 1000 working default has *not* sufficed,
which makes it the counterexample to the default that `SalazarCavazos` and `Laske` established.

The observables are linear/Gaussian, so the nominal-point evaluation is separately a clean end-to-end
check on the objective: PyBNF's Eq. 6 NLL at the published point reproduces `J*` to ~2e−3. See
`VALIDATION.md`, whose Gate B and gradient sections should be read alongside the ✅.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `-58.58395532491055` |
| paper-scale NLL at the PEtab nominal point | `-58.58618893016114` |
| optimality gap at nominal | `-0.0022336052505949056` |
| scored data points `n` | 72 |
| free parameters `k` | 22 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian
through trf's Coleman–Li core, ADR-0068) — handles this problem's estimated noise
scale (`sigma_pErk`, `sigma_pMek`), which plain `trf` refuses. `population_size = 100`,
`max_iterations = 1000` — the collection's documented working default, and on this problem it is the
one place that default has fallen short. A larger start count is the obvious next experiment.

**This slug carries the corpus's weakest gradient claim.** It is where lanl/PyBNF#535 found seven of
its 22 columns assembled from their initial-condition seed terms alone, several with reversed sign, so
a gradient fit was being steered uphill on them — fixed in ADR-0097. Post-fix, `tools/fd_check.py`
verifies all 22 columns at `3.5e−04`: that passes, but it is the largest residual among the fifteen
clean slugs, an order of magnitude above the rest, and its worst column `tau2` has never converged
under step-size refinement. If this slug's shortfall is ever traced to something other than budget,
the gradient is where to look first.

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`.
This problem exercises **lanl/PyBNF#508**: the per-gel scale factors (`s_pErk_*`,
`s_pMek_*`) are replicate-specific `observableParameters` supplied through the
`measurement_params:` sidecars (`experiment____model1_data*_measparams.tsv`). Before
#508 the replicate dimension was dropped and those bindings were silently lost; the
import now loads, simulates, and scores. The run recipe (`job_type`, `sbml_backend =
bngsim`, `wall_time_sim`) is supplied, not recovered — PEtab specifies a problem, not
a method.

## Contents

- `Fiedler_BMCSystBiol2016.conf` — the PyBNF job
- `model_Fiedler_BMCSystBiol2016.xml` — SBML model (emitted by the importer)
- `experiment____model1_data*.exp` (+ `_rep2`) — experimental data, 3 conditions × 2 replicates
- `experiment____model1_data*_measparams.tsv` — per-measurement observable/noise parameter tables (#508)
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `score.py` — scores a run against `J*`
- `best_fit_params.txt`, `information_criteria.txt` — the shipped fit's provenance
- `VALIDATION.md` — the full validation against `J*`

## Running

```bash
pybnf -c Fiedler_BMCSystBiol2016.conf -o
python score.py output
```
