# VALIDATION — Perelson_Science1996

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's proof that the
**log10-observable** objective is reproduced exactly, including the change-of-variables Jacobian.

> **Confidence: 95 / 100.** Fidelity for the log10 case is established both analytically and by
> near-exact recovery of J\* (OG = 5×10⁻⁷); the `lognormal` correction is derived and independently
> cross-checked with roadrunner. Deductions: the model is imported (not re-derived from Perelson 1996),
> and the fit uses `cmaes` (the gradient methods fail on this objective — documented).

## Gate A — objective fidelity (log10 + Jacobian)

Perelson's `observableTransformation = log10` means the paper's Eq. 6 NLL is the Gaussian likelihood of
`log10(V)`, transformed back to the observed `V` — so it carries the data Jacobian
`Σ log(V_obs·ln10) = 233.13`. Established two ways:

1. **Independent roadrunner reconstruction** (no PyBNF): scanning `(c, delta)` and profiling σ, the
   full log10 NLL `= N/2·(log(2πσ̂²)+1) + Σ log(V·ln10)` bottoms out at **222.46** (grid-limited) at
   c≈1.79 — matching J\* = 222.281. The *linear* NLL bottoms out at 232.3 and the log10 NLL *without*
   the Jacobian at −9.6; **only log10 + Jacobian matches J\***. This is why the linear import is
   unsolvable (OG ≈ 10) and the `gaussian → lognormal` correction is required.
2. **PyBNF's own log-likelihood** — PyBNF's `lognormal` family reports the full normalized `lnL`
   (matching `scipy.stats.lognorm.logpdf`, Jacobian included) in `information_criteria.txt`:
   `lnL = −222.2807694`. So `−lnL` is exactly the paper's Eq. 6 NLL, and `OG = −lnL − J*`.

**Verdict: PASS.** The `lognormal` correction is necessary and sufficient; `score.py` reads `−lnL`.

## Gate B — the fit reaches the benchmark optimum

Multi-start `cmaes` (24 × ≤500 gen, `random_seed = 1`, `sbml_backend = bngsim`) converges to
`J_paper = −lnL = 222.2807694` ⇒ **OG = 5×10⁻⁷ < 1.92**. Best fit c = 1.861, delta = 0.547, and the
fitted log10 σ = 0.123 (contrast the PEtab linear-scale nominal `sd = 1e5`). All 380 Marvin runs solved
Perelson; PyBNF matches to 7 significant figures.

**Verdict: PASS (SOLVED).**

## Configuration & corrections

- Import: `petab1to2_preserve_scale` → `import_job`.
- **Correction (documented import gap):** `noise_model gaussian → lognormal` for the `log10` observable.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = cmaes` (gradient methods fail on the lognormal
  objective: `gntr` → EFIM `SVD did not converge`; `lbfgs` → `sigma = nan` from the wide `[1e-10,1e10]`
  bounds). `population_size = 24`, `max_iterations = 500`, `random_seed = 1`.

## Bottom line

The log10-observable exemplar: the paper's objective is the log10 Gaussian NLL *with* the data
Jacobian, which PyBNF's `lognormal` `lnL` reproduces exactly (OG = 5×10⁻⁷). The one hand-correction
(`gaussian → lognormal`) is forced by an importer gap and independently justified with roadrunner.
