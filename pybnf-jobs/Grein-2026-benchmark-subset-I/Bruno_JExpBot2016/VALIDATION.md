# VALIDATION — Bruno_JExpBot2016

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's **fixed-σ** exemplar and
its **per-condition estimated initial condition** exemplar — the first subset-I problem whose free
parameters reach the model only through `condition:` parameter references.

> **Confidence: 94 / 100.** SOLVED with OG = 1.1×10⁻⁵ from a from-scratch 20-start multi-start; the
> objective identity is verified exactly (4.9×10⁻⁷) against the fixed per-point σ carried in the data,
> and the assembled gradient was validated upstream against central finite differences to 2.2×10⁻⁷
> across all 13 free parameters. Deduction: the model is imported (not re-derived from Bruno et al.
> 2016).

## Gate A — objective fidelity

**Bruno does not estimate its noise.** This is the exception to the collection's general rule: every
*other* subset-I problem carries a free `sigma`/`sd_*` parameter, but Bruno's σ is **known and fixed
per data point**, supplied by the PEtab measurement table and imported as `_SD` columns in the six
`experiment____model1_data*.exp` files. All 13 free parameters are model parameters; none is a noise
scale.

With σ fixed, the parameter-independent constants PyBNF drops from the minimized objective include the
`log σ` term as well as `½log(2π)`, so the restored constant is **not** the bare `N/2·log(2π)` that the
estimated-σ slugs show. The identity is

    −lnL  =  J_pybnf_reduced  +  Σᵢ log σᵢ  +  (N/2)·log(2π)

Checked directly against the `_SD` columns over all N = 77 scored points:

| term | value |
|---|---|
| `Σ log σᵢ` (from the `.exp` `_SD` columns) | −150.002162 |
| `(N/2)·log(2π)`, N = 77 | +70.758267 |
| predicted restored constant | **−79.243894** |
| observed (`−lnL − J_reduced`) | **−79.243894** |
| difference | 4.9×10⁻⁷ |

So PyBNF's reported `−lnL` is the paper's Eq. 6 NLL exactly, and `score.py` reads `−lnL`. The negative
restored constant is a consequence of σᵢ < 1 throughout (Gaussian densities > 1, hence a positive
`lnL`), not an anomaly.

**Verdict: PASS.**

## Gate B — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (20 starts × 500 iterations, box-center + Latin-hypercube seeded by
`random_seed = 1`, `sbml_backend = bngsim`) converges to `J_pybnf = 32.55570776` ⇒
`J_paper = −46.68818673` ⇒ **OG = 1.1×10⁻⁵ < 1.92**, against `J* = −46.68819792`. Wall time 41 s
(132 s CPU) on 6 cores.

The recovered vector matches the PEtab `nominalValue` point, whose own optimality gap is
`OG_nominal = 3.2×10⁻⁶` — i.e. the optimizer independently found the published optimum rather than
being seeded at it.

**Verdict: PASS (SOLVED).**

## Why this slug shipped with `cmaes`, and why it no longer does

Bruno was imported on 2026-07-20 with `job_type = cmaes`. That was **not** a judgement about the
landscape — it was a workaround. Every one of its 13 free parameters reaches the model only through a
`condition:` parameter-reference perturbation (a per-condition estimated initial condition, ADR-0076),
and `route_experiment` refused any `is_param_ref` mutation, so the gradient path aborted and the fit
silently fell back to a gradient-free optimizer.

**lanl/PyBNF #511** (merged #513, 2026-07-23 — three days *after* this slug was imported) composes the
chain rule instead of refusing: a condition assignment `target = free_param` gives the referenced free
parameter a `RouteContribution` on the target's own sensitivity column. Because Bruno's `szea` drives
six rate multipliers at once, a `ParamRoute` became a *sum* over contributions; the gradient, the
residual Jacobian and the EFIM/Fisher block all follow. Upstream validated the assembled gradient on
this exact job against central finite differences to 2.2×10⁻⁷.

The `cmaes` recipe was therefore a fossil of a limitation that no longer exists, and the conf now
carries `job_type = gntr` — PyBNF's fides-analogue, which is the method class the Grein study finds
best (`MS+fides` and `MS+CMA-ES` share the podium).

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`. Noise imported as fixed per-point σ from the
  measurement table (`_SD` columns); **no hand corrections**.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 20`,
  `max_iterations = 500`, `wall_time_sim = 10`, `random_seed = 1`.
- 6 conditions imported as 6 `experiment____model1_data*.exp` files with per-experiment `condition:`
  overrides; each condition sets initial amounts and rate multipliers from free parameters.

## Bottom line

The fixed-σ, condition-routed exemplar: 13 free parameters that reach the model *only* through
condition parameter references, fit by the gradient path to OG = 1.1×10⁻⁵ in 41 seconds. It is the
first slug in this collection whose ✅ is owed directly to an upstream gradient-routing fix rather than
to a tuning change, and it demonstrates the objective identity for the one noise regime — known,
fixed, per-point σ — that the rest of subset I does not exercise.
