# VALIDATION — Boehm_JProteomeRes2014

Validation of the PyBNF benchmark job against the **Grein et al. 2026** reference objective. The oracle
here is the benchmark's reference **J\*** (the best Eq. 6 NLL over all optimizer runs on Marvin), not a
primary-source data re-trace: the question is whether **PyBNF's optimizer reaches the benchmark optimum
on PyBNF's own faithful reproduction of the paper's objective**.

> **Confidence: 96 / 100.** Objective fidelity is established three independent ways (below); the fit is
> SOLVED with OG = 0.0012 from a from-scratch multi-start, and the identity is exact to ~1e-7. The
> deduction is the model itself is imported, not re-derived from the Boehm 2014 paper (the benchmark
> collection is the trusted upstream).

## Gate A — objective fidelity (three independent confirmations)

The Boehm objective is the full Gaussian NLL with **estimated** σ (the three `sd_*` are free
parameters). PyBNF's minimized per-point term is `(y−m)²/(2σ²) + log(σ)`; the paper's Eq. 6 per-point
term is that **plus `½log(2π)`**. Over N = 48 scored points the two differ by the constant
`C = 48·½log(2π) = 44.1090496`. Confirmed:

1. **Analytic** — from PyBNF's Gaussian kernel (`pybnf/noise/gaussian.py`): the objective sums
   `data_fit + log_normalizer`; the `½log(2π)` density constant is explicitly dropped from the
   objective and restored only in `log_density`. ⇒ `J_paper = J_pybnf + N·½log(2π)`.
2. **Point evaluation at the published best-fit** — evaluating PyBNF's objective at the PEtab
   `nominalValue` (the published Boehm 2014 fit) gives `J_pybnf(nominal) = 94.1129499`, so
   `J_pybnf(nominal) + C = 138.2219995` vs J\* = 138.2219682 — a residual of **3.1×10⁻⁵**, exactly the
   gap between the published nominal and the benchmark's own optimizer minimum (not an objective
   mismatch).
3. **PyBNF's own reported log-likelihood** — `information_criteria.txt` reports the full normalized
   `lnL = −138.2231909` at the best fit, i.e. `−lnL = J_pybnf + C` to 3.8×10⁻⁸. PyBNF's `lnL` *is* the
   paper's Eq. 6 NLL.

**Verdict: PASS.** `OG = −log_likelihood − J*` is the exact, self-consistent scoring; `score.py`
reproduces it.

## Gate B — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (20 starts × 500 iterations, box-center + Latin-hypercube seeded by
`random_seed = 1`, `sbml_backend = bngsim`) converges to `J_pybnf = 94.1141413` ⇒ `J_paper = 138.2231909`
⇒ **OG = 0.0012 < 1.92**. Every one of the paper's 380 Marvin runs also solved Boehm (min = median),
so the optimum is unambiguous; PyBNF reaches it to 0.0012.

**Verdict: PASS (SOLVED).**

## Configuration

- Import: `petab1to2_preserve_scale` (log10 param scales preserved, PR #491) → `import_job`. Emitted
  `noise_model <obs> = gaussian, sigma = fit sd_<obs>` per observable — correct, since Boehm's
  `observableTransformation = lin`.
- `edition = 2`, `sbml_backend = bngsim` (required for the gradient/sensitivity path), `job_type = gntr`,
  `population_size = 20`, `max_iterations = 500`, `random_seed = 1`.
- No hand corrections needed (linear observable; no import gap).

## Bottom line

A clean pilot: the objective-fidelity relationship (`OG = −lnL − J*`) is nailed three ways, and the
gradient fit reaches the benchmark optimum from scratch (OG = 0.0012). This job validates the whole
import → fit → score pipeline for the linear-observable case.
