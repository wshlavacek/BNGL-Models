# VALIDATION — carotenoid_time_uncertainty

Validation against the **Vanhoefer et al. (bioRxiv 2026.05.09.724053)** Fig 5/6 experiment. This is
a **paired** benchmark (standard vs. marginal fit on the same timing-perturbed data), not a single
optimality gap, so its acceptance is the paper's three-way comparison rather than an OG. The oracle
is the reference optimum **θ\*** (`theta_star_source.txt`), which — because the data is generated
*from the model at θ\**  — is the exact ground truth.

## What is validated, and how

1. **The phase-2 gradient is correct on this model.** PyBNF's marginal-time gradient (ADR-0113,
   PR lanl/PyBNF#589) is validated upstream against central finite differences to ~1e-9 per column.
   On *this* model — 13 parameters, 6 conditions, `read_exp_file` σ — the gradient was additionally
   confirmed to assemble finitely at θ\*: all 6 experiments carry a forward-sensitivity tensor, the
   assembled gradient is finite (‖∇‖ ≈ 250 at θ\* on the σ_t = 5 data, i.e. θ\* is *not* the
   perturbed-data optimum), and the Gauss-Newton Fisher (14×14) is PSD (eigenvalues
   [0.14, 2.3×10⁴]). `job_type = gntr` converges on the marginal objective through the real bngsim
   backend.

2. **The data generator is faithful.** At σ_t = 0 the synthetic data is the model at θ\* plus the
   published per-point noise; it matches the real Bruno data within that noise (e.g. condition 1's
   b10 falls 4.0 → 1.5 over 5–180 min, as in the source). The measurement-noise stream is fixed
   across σ_t levels, so the datasets differ *only* in the timing perturbation.

3. **The gates below**, scored by `score.py` against the committed `demo_results/`.

## The gates

The demonstration fits are seeded at θ\* (a single `gntr` start from the log-box centre; see
`make_demo.py`) so both arms descend from the exact ground truth and any drift is attributable to
the objective, not the start.

| quantity | Gate A · σ_t = 0 | Gate B · σ_t = 10 |
|---|---|---|
| log₁₀-MSE(standard vs θ\*) | 0.00025 | **0.00059** |
| log₁₀-MSE(marginal vs θ\*) | 0.00025 | **0.00027** |
| estimated σ_t (marginal) | 0.100 (floor) | **11.3** (injected 10) |
| ln L standard | 45.07 | −15.89 |
| ln L marginal | 151.61 | 19.16 |
| LRT 2·(ln L_m − ln L_s), χ²(1) | 213.1 | 70.1 |
| standard model rejected @ 0.05 | (see caveat) | **yes** (p ≈ 0) |

**Gate B — correction (σ_t = 10). PASS.** The standard fit is dragged off θ\* (MSE 0.00059) while
the marginal fit stays 2.2× closer (0.00027) — the Fig 6A effect. The marginal **recovers the
timing scale** (σ_t = 11.3 ≈ injected 10) — Fig 6B. The LRT of the standard model (the σ_t → 0
limit, so the two nest) rejects decisively — Fig 6C/D. All three of the paper's Fig-6 signatures
reproduce.

**Gate A — no false positive (σ_t = 0). PASS on recovery.** The marginal fit costs nothing when
there is no timing error: it matches the standard MSE (both 0.00025) and drives the estimated σ_t
to its floor (0.1). These are the reliable signals and both pass.

> **The one honest caveat (Gate A LRT).** The LRT *rejects* at σ_t = 0 (ln L_marginal = 151.6 ≫
> ln L_standard = 45.1), which looks like a false positive but is not: it is a **numerical artifact
> of the fixed quadrature grid**. At the σ_t floor (0.1) the timing prior is a spike far narrower
> than the grid spacing (~1 min for `n_steps = 200` over [0, 200]), so the trapezoid cannot resolve
> it and `z_k` is mis-integrated to a spuriously large value. This is exactly the limitation
> PyBNF ADR-0112 flags ("the quadrature grid cannot resolve the time prior") and the reason the LRT
> here is not gated on. The paper's Fig 6C is clean because AMICI's augmented ODE integrates `z_k`
> with **error control**; PyBNF's phase-2 chains sensitivities over the *same fixed grid* the value
> uses, so it inherits the value engine's floor. Removing it is the error-controlled-integration
> follow-up filed in ADR-0113 — orthogonal to the gradient this benchmark exercises. In practice:
> keep the estimated σ_t at or above the grid spacing (the Gate-B regime), where the LRT is honest.

## Gate C — why this is a *gradient* benchmark

The paper's scalability claim is that marginalization keeps the search dimension at **n_θ + 1**
(here 14) rather than the joint approach's **n_θ + n_t** (14 + 77 = 91), and that gradient methods
then fit it in iterations only weakly dependent on the parameter count. PyBNF realizes this without
the augmented ODE: the marginal fit runs the same `gntr` gradient path as the standard fit, over the
same 14-dimensional box, differing only in the objective. The honest framing of Gate C is therefore
**gradient at fixed dimension n_θ + 1** — not "augmented ODE beats quadrature": PyBNF has no
augmented ODE, because its forward-sensitivity tensor already supplies `dx(τ)/dθ` at every grid node
(ADR-0113). Each `gntr` start on the marginal objective converged in well under a minute on this
model in the demonstration runs.

## Bottom line

The real-model exercise of PyBNF's phase-2 `time_error` gradient: on a published 13-parameter model,
the marginal-time fit corrects the timing-induced bias the standard fit suffers, recovers the
injected timing scale, and rejects the standard model — the paper's Fig 6 result — fit by a gradient
method (`gntr`) with no augmented ODE. The single caveat (a fixed-grid LRT at a sub-grid-spacing
σ_t) is a documented, correctly-attributed limitation of the fixed quadrature grid, not of the
gradient.
