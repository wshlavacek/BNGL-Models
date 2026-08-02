# VALIDATION — Armistead_CellDeathDis2024

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's **prediction-dependent
noise** exemplar: σ is not a constant but a formula in the model's own prediction.

> **Confidence: 93 / 100.** SOLVED with OG = 5.8×10⁻⁶ from a from-scratch 20-start multi-start, no
> seeding; the restored constant is the bare `(N/2)log(2π)` to 6 decimal places, confirming the
> reduced-objective bookkeeping for a prediction-dependent σ. Deductions: the model is imported (not
> re-derived from Armistead et al. 2024), and the run is a single seed rather than the paper's ten.

## Gate A — objective fidelity

All four observables carry a **prediction-dependent** noise scale — σ is proportional to the
prediction itself, with an estimated relative coefficient:

```
noise_model Sphinga  = gaussian, sigma = prediction_formula Sphinga*sd_Sphinga_obs
noise_model Cer      = gaussian, sigma = prediction_formula Cer*sd_Cer_obs
noise_model Sphingo  = gaussian, sigma = prediction_formula Sphingo*sd_Sphingo_obs
noise_model S1P      = gaussian, sigma = prediction_formula S1P*sd_S1P_obs
```

This is a strictly harder case than a free scalar σ: the noise scale depends on the parameters *and*
on the simulated trajectory, so its gradient and EFIM Fisher block need the σ-formula chain rule
threaded through the same forward sensitivity as the residual (ADR-0079/0080). That machinery is what
makes `gntr` viable here at all.

The bookkeeping is nonetheless the ordinary one, because every remaining dropped constant is
per-point and parameter-independent:

| term | value |
|---|---|
| `(N/2)·log(2π)`, N = 58 | 53.298435 |
| observed restored constant (`−lnL − J_reduced`) | 53.298435 |

`J_pybnf = −355.214617` ⇒ `−lnL = −301.916182`, and `score.py` reads `−lnL` as the paper's Eq. 6 NLL.

**Verdict: PASS.**

## Gate B — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (20 starts × 500 iterations, box-center + Latin-hypercube seeded by
`random_seed = 1`, `sbml_backend = bngsim`) converges to `J_pybnf = −355.2146166` ⇒
`J_paper = −301.9161821` ⇒ **OG = 5.8×10⁻⁶ < 1.92**, against `J* = −301.9161878`. Wall time 256 s
on 6 cores.

The optimizer found this basin from box-sampled starts; it was **not** seeded at the PEtab
`nominalValue` point (whose own gap is `OG_nominal = 5.8×10⁻⁶`, i.e. the same optimum). This is a
genuine search result, which matters for the benchmark's purpose — the Grein protocol scores what an
optimizer *finds*, not what it is handed.

**Verdict: PASS (SOLVED).**

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`. Four `prediction_formula` noise models emitted
  directly from the PEtab `noiseFormula` column; **no hand corrections**.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 20`,
  `max_iterations = 500`, `wall_time_sim = 10`, `random_seed = 1`.
- k = 14 free parameters (10 model + 4 relative noise coefficients `sd_*_obs`), n = 58 scored points.
- 10 experiments (`experiment____{mutant,wildtype}_rep*.exp`).

## Provenance

Run against **bngsim rebuilt from `~/Code/bngsim` at `faf9e6c` (2026-08-01)**. This matters: both that
build and the previous one report version `0.11.35`, so `pip` cannot distinguish them, and the twelve
intervening commits touch the Jacobian/codegen/sensitivity paths this fit exercises. An earlier run of
this slug on the 2026-07-26 build gave the same verdict (OG = 6.0×10⁻⁶, 265 s), so the result is
stable across both.

## Bottom line

The prediction-dependent-noise exemplar: four observables whose σ is a formula in the prediction
itself, fit by the gradient path to OG = 5.8×10⁻⁶ in 256 s from Latin-hypercube starts. It exercises
the σ-formula chain rule in both the gradient and the EFIM Fisher block (ADR-0079/0080), and confirms
that a prediction-dependent σ leaves the reduced-objective bookkeeping unchanged — the restored
constant is still the bare `(N/2)log(2π)`.
