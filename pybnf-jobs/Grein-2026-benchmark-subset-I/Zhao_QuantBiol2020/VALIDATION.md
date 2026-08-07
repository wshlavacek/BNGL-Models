# VALIDATION — Zhao_QuantBiol2020

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's **lost-search-scale**
case: all 28 of its estimated parameters are `parameterScale = log10` upstream, and every one of them
was imported as a *linear* `uniform_var` until lanl/PyBNF#548 was fixed. It is the second of the two
**per-measurement σ** problems (with `Fiedler`), binding seven `sd_confirmed_infected_*` noise
parameters through `_measparams.tsv` sidecars (ADR-0083).

> **Confidence: 93 / 100.** SOLVED with `OG = 0.000005` from a from-scratch 100-start multi-start, no
> seeding — the fit lands on `J*` to five decimal places. The objective is independently corroborated
> against upstream's own `simulatedData`, and the assembled gradient is verified against central
> differences **in the log sampling space the corrected conf actually searches**. Deductions: the
> model is imported (not re-derived from Zhao et al. 2020); the run is a single seed; and three of
> the 28 parameters sit on a box bound at the optimum (below), which is a property of the problem
> rather than a defect but does mean the reported optimum is a constrained one.

## Gate A — objective fidelity

Linear observables, per-measurement σ, so the restored constant is the `(N/2)log(2π)` form with no
Jacobian. From the fit's own `information_criteria.txt` against the minimized objective:

| term | value |
|---|---|
| PyBNF reduced objective | 425.874099 |
| restored constant (`−lnL − J_reduced`) | 75.352960 |
| `(N/2)·log(2π)`, N = 82 | 75.352960 |
| `J_paper = −log_likelihood` | 501.227058 |

Independently, §2c's oracle recomputes the Eq. 6 NLL at the nominal point straight from upstream's
`measurementData` joined to `simulatedData`, with **no PyBNF in the loop**, and reproduces PyBNF
exactly — this slug carries `obj ✓`.

**One subtlety worth recording, because it cost an afternoon.** PyBNF's *reduced* objective for this
slug already contains the `Σ nⱼ log σⱼ` terms: the reduced objective drops only the
*parameter-independent* per-point constants, and a fitted σ is not parameter-independent. Checking a
recomputed NLL against `reduced_objective` therefore reports a spurious failure on this slug and on
every other estimated-σ slug. `J_paper == −log_likelihood` is the unambiguous scale, and the one
`score.py` uses. See `tools/sigma_profile.py` and its entry in `tools/README.md`.

**Verdict: PASS.**

## Gate B — the gradient is the one the objective implies

The #535 sweep checked this slug and found it clean — but **in linear sampling space**, which is not
the space the corrected conf searches. The gradient is assembled in sampling space `u`, so switching
28 parameters from `uniform_var` to `loguniform_var` changes what is being differentiated, and the
earlier clean result does not carry over. Re-checked at the PEtab nominal point, where
`OG_nominal = 276` makes the gradient large enough for the comparison to be a real test (§2b):

| | worst relative error |
|---|---|
| all 28 columns, log sampling space, `h = 3×10⁻⁴` | **1.95×10⁻⁴** |

Most columns agree to 1e−06 or better; no column is structurally zero. The worst reading,
`R_Stage_II_Hubei` at 1.95e−04, is ordinary finite-difference noise on a stiff epidemic model.

**Verdict: PASS.**

## Gate C — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (**100 starts × 1000 iterations**, `random_seed = 1`,
`sbml_backend = bngsim`, no seeding) converges to `J_pybnf = 425.874099` ⇒
`J_paper = 501.227058` ⇒ **`OG = 0.000005` < 1.92**, against `J* = 501.2270538`. All 100 starts
retired: 75 on `step is negligible`, 25 on `reached max_iterations`. Wall clock ≈ **54 minutes**.

**This slug is the collection's clearest demonstration that a search scale is not cosmetic.** Holding
everything else fixed — same recipe, same seed, same budget, same objective — and changing only
whether the 28 log10 parameters are sampled in log space:

| configuration | best reduced objective | `OG` |
|---|---|---|
| linear `uniform_var` (before lanl/PyBNF#548), 100 × 1000 | 717.5 and decelerating, run abandoned | ≥ 291 |
| **log `loguniform_var` (after #548), 100 × 1000** | **425.874099** | **0.000005** |

The corrected run beat twenty minutes of the linear-scale run's progress in under ninety seconds.
The reason is not subtle once stated: `gamma_*` live on `[1e-08, 1]` with an optimum near 0.05–0.39,
and the `sd_*` on `[0.001, 1e5]` with MLEs of 186–5013. Sampled linear-uniform, a single draw lands
below 5000 with probability 0.05; across 28 such parameters, effectively no box-sampled start begins
anywhere near the basin.

**Verdict: PASS (SOLVED).**

### Three parameters sit on a box bound, and that is the problem's shape

`R_Stage_I_Wuhan`, `R_Stage_I_Hubei` and `R_Stage_I_China` all land on their upper bound of exactly
100. That is not a fit pathology: Stage I is the uncontrolled early-epidemic window, where the
observable is dominated by exponential growth whose rate saturates in `R` — once `R` is large the
likelihood is nearly flat in it, so the MLE runs to whatever ceiling the box provides. The reference
`J*` is matched to 5e−06 with those parameters at the wall, which is the evidence that the reference
optimum is the same constrained one. Any profile-likelihood reading of these three should treat them
as bounded below rather than as point estimates.

### The nominal point is not the optimum, and its distance was overstated

`OG_nominal = 276.12`, so unlike `Brannmark` or `Raia` this problem's PEtab `nominalValue` vector is
not its optimum and the fit was not expected to reproduce it. But that 276 overstates the distance:
the nominal σ are placeholders at 1000, and with the seven estimated σ profiled to their MLEs the
nominal *dynamics* are worth `OG = 135.75`, not 276. Roughly half the apparent gap was the
placeholder, not the model. See the `§` note in the collection README.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`; **no hand corrections**. The conf's parameter
  block was regenerated from a clean re-import after lanl/PyBNF#548; every bound is byte-identical to
  the previous one and only the `*_var` keyword changed.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 100`,
  `max_iterations = 1000`, `wall_time_sim = 10`, `random_seed = 1`.
- k = 28 free parameters (21 model + 7 noise), n = 82 scored points, one observable
  (`observable_confirmed_infected`) across seven conditions — three Wuhan stages, two Hubei, two
  China. Upstream ships a held-out `measurementData_test_*` split which is **not** part of this
  objective; only the training table is scored.

## History — the defect this slug found

**lanl/PyBNF#548 — the log estimation scale.** `petab1to2_preserve_scale` re-injects the
`parameterScale` that `petab.v2.petab1to2` drops, skipping any row that already carries a prior so a
scale petab1to2 already folded into one is not clobbered. That guard cannot be implemented in v2
alone: petab1to2 **materializes** v2's implicit `uniform` default into the converted table whenever
the v1 table has a prior column at all — and Zhao's four prior columns are present and 100% empty.
After conversion the materialized default and a declared `uniform` are the same cell, so the
injection declined and all 28 log10 scales were lost.

**It was silent by construction, and that is the transferable lesson.** The re-injected prior sets
only the search scale and initial sampling; PyBNF's optimiser objective *excludes* the prior. So the
objective was unchanged (`J_paper` at nominal agrees to 13 significant digits before and after the
fix), the `obj ✓` oracle check still reproduced exactly, and the finite-difference gradient check
still passed. Every instrument this corpus uses to catch a bad setup reported green. The only visible
symptom was a fit that descended slowly and stalled — which is indistinguishable from a problem that
needs a larger budget, and was very nearly written up as one.

Fixed in `e008d345`; `Schwen_PONE2014` was the other affected slug, at 24 of 25 parameters.

## Provenance

Run against **bngsim 0.12.1** with PyBNF at `e008d345` (ADR-0103 + ADR-0104 + the #548 fix).

## Bottom line

The lost-search-scale case: 28 log10 parameters against 82 points, solved on `J*` to five decimal
places once they were actually searched in log space. It is also the slug that paid for itself as a
test case — one silent conversion defect affecting any PEtab v1 problem whose parameter table merely
*has* a prior column, which no fixture in the PyBNF suite had caught and which none of this corpus's
own verification gates could see.
