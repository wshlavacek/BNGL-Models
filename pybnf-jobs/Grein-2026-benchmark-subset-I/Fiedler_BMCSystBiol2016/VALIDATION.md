# VALIDATION — Fiedler_BMCSystBiol2016

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's **replicate
`observableParameters`** case — per-measurement σ bound through `_measparams.tsv` sidecars
(ADR-0083) — and, after this run, its clearest **solved-but-not-saturated** result.

> **Confidence: 76 / 100.** SOLVED with OG = 1.004 from a from-scratch 100-start multi-start, no
> seeding. Inside the χ² threshold, but **the optimizer did not find the reference basin**: the
> problem's own nominal point scores −0.0022, so a fit that reached the published optimum would score
> ≈ 0 and this one is a full 1.0 above it. Deductions: the model is imported (not re-derived from
> Fiedler et al. 2016); the run is a single seed; the result is not saturated (see Gate B); and this is
> the slug whose gradient was wrong until ADR-0097, whose FD residual is still the largest of the
> corpus's clean slugs, and whose `tau2` column has never converged under step-size refinement.

## Gate A — objective fidelity

Noise is **per-measurement σ**, bound via the `_measparams.tsv` sidecars, so every σ is a fitted or
table-supplied quantity and no `Σ log σᵢ` constant survives — the restored constant is `(N/2)log(2π)`
alone:

| quantity | value |
|---|---:|
| PyBNF reduced objective | −123.744014 |
| restored constants | 66.163574 |
| `(N/2)log(2π)`, n = 72 | 66.163574 |
| `J_paper = −log_likelihood` | −57.580439 |
| reference `J*` | −58.583955 |

The import exercises lanl/PyBNF#508 (replicate-specific `observableParameters` / `noiseParameters`
binding): the six `.exp` files come in three condition pairs, each with a `_rep2` sibling, and each
replicate binds its own gel-specific scale (`s_pErk_20140430_gel1` and the seven others).

**Verdict: PASS.**

## Gate B — the fit is inside the threshold but not on the optimum

From-scratch multi-start `gntr` (100 starts × 1000 iterations, box-sampled starts seeded by
`random_seed = 1`, `sbml_backend = bngsim`, no seeding from the nominal point) converges to
`J_pybnf = −123.7440139` ⇒ `J_paper = −57.5804395` against `J* = −58.5839553`:

    OPTIMALITY GAP  OG = 1.003516  <  1.92

Wall time **6 h 13 min 19 s** on 10 cores — by far the most expensive fit in the corpus, and roughly
23× `Rahman_MBS2016` at k = 9 against this slug's k = 22.

**This is solved by the benchmark's criterion and short of the optimum, and the two should not be
conflated.** The χ² threshold (1.92, α = 0.05, 1 dof) makes `OG = 1.004` statistically
indistinguishable from the reference optimum, so the ✅ is correct as scored. But unlike
`Rahman_MBS2016` (0.000000) and `SalazarCavazos_MBoC2020` (2.9×10⁻⁵), which land *on* `J*`, this fit
stopped in a nearby basin. The ordering is the giveaway:

| point | OG |
|---|---:|
| PEtab nominal (published optimum) | −0.0022 |
| **this fit, 100 × 1000 unbiased** | **1.0035** |
| threshold | 1.92 |

The reference basin demonstrably exists and is reachable — the nominal point sits in it — and 100
box-sampled starts did not find it. Progress was still positive but nearly exhausted when the budget
ran out: the best objective moved 1.113 → 1.034 → 1.004 over the run's last ~90 minutes, i.e. grinding
rather than converging. A larger start count is the obvious next experiment; this is the first problem
in the corpus where 100 × 1000 has *not* sufficed, which makes it the counterexample to the working
default that `SalazarCavazos` and `Laske` established.

**Verdict: PASS (SOLVED, not saturated).**

## Gradient — the caveat this slug carries

Fiedler is the slug lanl/PyBNF#535 found **broken**: seven of its 22 columns assembled from their
initial-condition seed terms alone and never entered `sensitivity_params`, several with reversed sign,
so a gradient fit was being steered uphill on them. Its six species are initialised from
`initialAssignment`s that are the closed-form steady state over `K_1, k2 … k11` — the same symbols that
are the rate constants — and ADR-0096's rule dropped their axes on the premise that an IC-seeding id is
absent from the ODE right-hand side. Fixed in **ADR-0097**, which drops an axis only when the model
states the RHS never reads the id.

Post-fix, `tools/fd_check.py` verifies all 22 columns at an interior, bounds-clear point: worst
`|Δ|/‖∇‖∞ = 3.5e−04`. That passes, but it is **the largest residual among the fifteen clean slugs**, an
order of magnitude above the rest, and the worst column is `tau2` — whose own central difference spans
−335 to −348 across step sizes without converging, with the assembled −338 sitting inside that band.
So "the gradient is correct here" is the weakest such claim in the corpus, and if this slug's shortfall
is ever traced to something other than budget, the gradient is where to look first.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`. Noise imported as per-measurement σ via
  `_measparams.tsv` sidecars (ADR-0083); **no hand corrections**.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 100`,
  `max_iterations = 1000`, `wall_time_sim = 30`, `random_seed = 1`.
- k = 22 free parameters, n = 72 scored points, 3 conditions × 2 replicates = 6 experiments.

## Provenance

Run 2026-08-06, 15:52:52 – 22:06:11, starting against **PyBNF `d014272e`** and **bngsim `59c6f38`**.

**PyBNF advanced to `47de23cf` (ADR-0103 — `sbml_atol` derived from the model's own scale rather than
inheriting the backend default `1e-8`) at 16:05, twelve minutes into this run.** Python binds modules
at process start, so the fit executed the code as of 15:52:52, but the overlap is recorded here rather
than assumed harmless. Re-evaluating the recorded best-fit point under `47de23cf` gives `−123.7440183`
against the recorded `−123.7440139` — a shift of 4.5×10⁻⁶, ODE integration noise some six orders below
the threshold, moving OG from 1.003516 to 1.003512. The recorded result stands. Fiedler is also not
among the four slugs whose nominal check the derived tolerance moved (Armistead, Bertozzi, Brannmark,
Giordano).

## Bottom line

The expensive case, and the corpus's second "solved but not saturated" row after
`Crauste_CellSystems2017`: 22 parameters against 72 points, `OG = 1.004` in 6 h 13 m from unbiased
starts, against a reference basin the nominal point proves is reachable. It is the first problem where
the collection's 100 × 1000 working default has fallen short, and — carrying both the corpus's largest
FD residual and a column that has never converged under step refinement — the one whose ✅ should be
read with its Gate B and gradient sections attached.
