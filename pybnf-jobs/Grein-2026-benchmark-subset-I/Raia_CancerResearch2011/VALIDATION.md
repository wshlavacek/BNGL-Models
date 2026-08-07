# VALIDATION — Raia_CancerResearch2011

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's largest **🟢 → ✅
conversion** (k = 39, n = 205) and one of two problems using **prediction-dependent σ**
(`sigma = prediction_formula …`, with `Armistead`), where the noise scale is a function of the model
prediction rather than a free scalar or a data column.

> **Confidence: 94 / 100.** SOLVED with `OG = 0.000009` from a from-scratch 100-start multi-start,
> no seeding — the fit lands on `J*` to five decimal places, and beats its own nominal point
> (`OG_nominal = 0.78`) by a factor of ~87,000. 83 of 100 starts converged on `step is negligible`.
> Deductions: the model is imported (not re-derived from Raia et al. 2011); the run is a single
> seed; and one of the 39 parameters rests on a box bound at the optimum (below).

## Gate A — objective fidelity

Linear observables with prediction-dependent σ, so the restored constant is the `(N/2)log(2π)` form
with no Jacobian. From the fit's own `information_criteria.txt`:

| term | value |
|---|---|
| PyBNF reduced objective | 156.927377 |
| restored constant (`−lnL − J_reduced`) | 188.382399 |
| `(N/2)·log(2π)`, N = 205 | 188.382399 |
| `J_paper = −log_likelihood` | 345.309776 |

The identity holds to every digit printed. This slug is **not** covered by §2c's independent oracle —
its upstream `measurementData` and `simulatedData` do not join one-to-one — so Gate A here rests on
the internal identity plus Gate C's agreement with `J*`, not on a third-party recomputation. That is
weaker evidence than `Zhao` or `Blasi` carry, and the collection README's `obj ✓` column correctly
leaves this row blank.

**Verdict: PASS (by internal identity; no independent oracle available).**

## Gate B — the gradient is the one the objective implies

This slug is where lanl/PyBNF#537 was found and fixed, and its history is a caution about how a
gradient defect can hide.

| | worst relative error |
|---|---|
| lanl/PyBNF#535 sweep, first pass | a factor-of-two disagreement, recorded as **not reproducible** after five attempts |
| re-run against HEAD after #537 | **4.85×10⁻⁵** |

The factor of two was real. It was an IC-seeding parameter whose own sensitivity axis *is* the whole
derivative, with the seeded contribution being summed on top of it — so the column read exactly
double (ADR-0100, "the parameter and initial-condition sensitivity axes are not independent"). It
presented as irreproducible because whether it fired depended on the evaluation point.

**Verdict: PASS (after #537).**

## Gate C — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (**100 starts × 1000 iterations**, `random_seed = 1`,
`sbml_backend = bngsim`, no seeding) converges to `J_pybnf = 156.927377` ⇒
`J_paper = 345.309776` ⇒ **`OG = 0.000009` < 1.92**, against `J* = 345.3097673`. All 100 starts
retired: 83 on `step is negligible`, 17 on `reached max_iterations`.

**The fit beats its own nominal point, which is the claim §1 asks for.** `OG_nominal = 0.78` — the
PEtab `nominalValue` vector was already inside the solved threshold, so this problem could have been
"converted" by simply holding the published point. It was not: unbiased box-sampled starts found a
point 87,000× closer to `J*` than the published one. That is the difference between showing the
optimizer *holds* a known optimum and showing it *finds* one.

The cost estimate this retires is worth recording. Issue #38 once projected **47–76 hours** for this
slug by extrapolating from `k = 39`. It ran in well under two hours on ten cores, alongside two other
fits. The issue had already retracted that projection in principle — `Elowitz` (k=21) at 8 m 47 s
against `Fiedler` (k=22) at 6 h 13 m is a 40× spread at equal `k` — and this is the confirmation.
**Cost tracks stiffness and model size, not parameter count.**

**Verdict: PASS (SOLVED).**

### One parameter rests on a bound

`sd_pJAK2_rel` lands on its lower bound of `1e-05`. It is a noise scale, and a noise scale driven to
its floor means the corresponding observable is fitted to within the resolution the box permits —
the likelihood keeps rewarding a smaller σ and the bound stops it. `J*` is matched to 9e−06 with it
there, so the reference optimum is the same constrained one. It should be read as bounded-above in
precision rather than as a point estimate, and any profile-likelihood analysis of this problem should
widen that box first.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`; **no hand corrections**. All 39 estimated
  parameters are `parameterScale = log10` upstream and import correctly as `loguniform_var` — this
  slug is **not** affected by lanl/PyBNF#548, whose trigger is a v1 prior column that Raia's
  parameter table does not have.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 100`,
  `max_iterations = 1000`, `wall_time_sim = 10`, `random_seed = 1`.
- k = 39 free parameters, n = 205 scored points.

## Provenance

Run against **bngsim 0.12.1** with PyBNF at `e008d345` (ADR-0100 + ADR-0103 + ADR-0104 + the #548 fix).

## Bottom line

The largest 🟢 conversion in the collection and the last one on §1's checklist: 39 parameters against
205 points, solved on `J*` to five decimal places from unbiased starts, on a problem whose gradient
was silently doubled on one column until lanl/PyBNF#537. It also retires the corpus's last
`k`-based cost projection — 47–76 hours estimated, under two hours actual.
