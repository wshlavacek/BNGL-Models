# VALIDATION — Crauste_CellSystems2017

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's **sparse-data** case:
n = 21 scored points against k = 12 free parameters.

> **Confidence: 88 / 100.** SOLVED with OG = 0.583 from a from-scratch 20-start multi-start, no
> seeding; the fixed-σ objective identity is verified exactly (1.5×10⁻⁷). Deductions: the model is
> imported (not re-derived from Crauste et al. 2017); the run is a single seed; and unlike the other
> solved slugs this one lands *near* rather than *at* the reference optimum — see below, which is a
> statement about the problem, not a defect.

## Gate A — objective fidelity

**Crauste's σ is known and fixed per data point**, supplied by the PEtab measurement table and
imported as `_SD` columns in `experiment____model1_data1.exp`. There is no `noise_model` line and no
free noise parameter; all 12 free parameters are model parameters. (Bruno, Rahman and SalazarCavazos
share this regime — see the noise-regime table in the directory README.)

With σ fixed, the dropped constants include `log σ` as well as `½log(2π)`, so the restored constant is
not the bare `(N/2)log(2π)`:

    −lnL  =  J_pybnf_reduced  +  Σᵢ log σᵢ  +  (N/2)·log(2π)

Checked against the `_SD` columns over all N = 21 scored points:

| term | value |
|---|---|
| `Σ log σᵢ` (from the `.exp` `_SD` columns) | 161.833304 |
| `(N/2)·log(2π)`, N = 21 | 19.297709 |
| predicted restored constant | **181.131013** |
| observed (`−lnL − J_reduced`) | **181.131013** |
| difference | 1.5×10⁻⁷ |

The σᵢ here are large (cell counts), so `Σ log σᵢ` dominates and the restored constant is large and
positive — the mirror image of Bruno, where σᵢ < 1 made it negative. Both are the same identity.

**Verdict: PASS.**

## Gate B — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (20 starts × 500 iterations, box-center + Latin-hypercube seeded by
`random_seed = 1`, `sbml_backend = bngsim`) converges to `J_pybnf = 9.9088609` ⇒
`J_paper = 191.0398740` ⇒ **OG = 0.583 < 1.92**, against `J* = 190.4570655`. Wall time 196 s on
6 cores.

**This one is solved but not saturated, and that is worth stating plainly.** The other solved slugs
reach J\* to 5–6 significant figures; here the gap is 0.58 — comfortably inside the χ² threshold
(1.92, α = 0.05, 1 dof), so the fit is statistically indistinguishable from the reference optimum, but
the optimizer did not land exactly on it. It also did not quite match the PEtab `nominalValue` point,
whose gap is `OG_nominal = 0.509`.

That ordering (fit 0.583 > nominal 0.509, both < 1.92) is the expected signature of a **sparse,
weakly-identified** problem: 12 parameters against 21 points leaves a flat basin floor where many
parameter vectors are near-equivalent in likelihood, and a trust-region method stops when its step
test is satisfied, not when it has crawled to the exact minimum. The verdict is SOLVED by the
benchmark's own criterion; the slug is simply not a precision demonstration the way Armistead or
Bruno are.

**Verdict: PASS (SOLVED).**

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`. Noise imported as fixed per-point σ from the
  measurement table (`_SD` columns); **no hand corrections**.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 20`,
  `max_iterations = 500`, `random_seed = 1`.
- k = 12 free parameters, n = 21 scored points, 4 observables
  (EarlyEffector, LateEffector, Memory, Naive), one experiment.

## Provenance

Run against **bngsim rebuilt from `~/Code/bngsim` at `faf9e6c` (2026-08-01)**. Both that build and the
previous one report version `0.11.35`, so `pip` cannot distinguish them; the twelve intervening
commits touch the Jacobian/codegen/sensitivity paths this fit exercises. An earlier run on the
2026-07-26 build gave the same verdict at the same objective (OG = 0.5828, 226 s), so the result is
stable across both.

## Bottom line

The sparse-data, fixed-σ case: 12 parameters against 21 points, solved at OG = 0.583 in 196 s. It
demonstrates the fixed per-point σ identity a second time (with the opposite sign on `Σ log σᵢ` from
Bruno), and it is the collection's clearest example that "solved" is a statistical criterion — a fit
can be inside the χ² threshold without sitting on the reference optimum.
