# VALIDATION — Laske_PLOSComputBiol2019

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's **mixed noise-regime**
case — 33 of its 42 scored points are natural-log (`lnnormal`), 9 are linear Gaussian — and its
**aliased-parameter** case: a COPASI export in which every rate law reads a `ModelValue_*` alias that
an `initialAssignment` derives, and no source name appears in a rate law at all.

> **Confidence: 92 / 100.** SOLVED with OG = −1×10⁻⁶ from a from-scratch 100-start multi-start, no
> seeding — the fit reaches `J*` itself. The mixed restored-constant identity is verified **exactly**
> against the data, and the assembled gradient is verified against central differences across all 13
> free parameters. Deductions: the model is imported (not re-derived from Laske et al. 2019); the run
> is a single seed; and it needs a larger budget than the collection default, which is a statement
> about the problem's basin structure that one run cannot fully characterise.

## Gate A — objective fidelity

Laske is the collection's only **mixed** regime, so its restored constant is neither of the two
standard forms. Four observables (`R_M5`, `RVSegTot`, `RCSegTot`, `V_rel`) are `lnnormal` — natural-log
normal, imported by lanl/PyBNF#509 — and nine (`IntNucOffset`, `FracNucInt_1..8`) are linear
`gaussian`. Every σ is estimated, so no `Σ log σᵢ` term survives; but the log-transformed points each
contribute a change-of-variables Jacobian `log(y_obs)`:

    −lnL  =  J_pybnf_reduced  +  (N/2)·log(2π)  +  Σ_{lnnormal points} log(y_obs)

Computed straight from the `.exp` columns, with no reference to the run:

| term | value |
|---|---|
| `(N/2)·log(2π)`, N = 42 (all points) | 38.595418 |
| `Σ log(y_obs)` over the 33 `lnnormal` points | 261.089654 |
| **predicted restored constant** | **299.685073** |
| **observed** (`−lnL − J_reduced`) | **299.685073** |

The Jacobian term is nearly seven times the `(N/2)log(2π)` one, which is why the collection README's
noise-regime table — whose "restored constant" column is the σ-source part only — cannot be read as
the whole constant for this slug.

**Verdict: PASS.**

## Gate B — the gradient is the one the objective implies

This slug is where the collection's gradient path was checked against finite differences on a *real*
problem rather than a fixture, and the check found a bug (lanl/PyBNF#534, below). Central differences
of PyBNF's own `loss(u)` against the assembled `gradient(u)`, at the PEtab nominal point, over all 13
free parameters:

| | worst relative error |
|---|---|
| before #534 | **1.0** — `k_syn_R_M` assembled to exactly 0 against ≈ −10.4, at every step size |
| after #534 | **2.5×10⁻³** (at `h = 3×10⁻⁴`; the small-gradient components converge onto the assembled values as `h` grows, which is ordinary FD roundoff on a stiff model) |

**Verdict: PASS (after #534).**

## Gate C — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (**100 starts × 1000 iterations**, `random_seed = 1`,
`sbml_backend = bngsim`) converges to `J_pybnf = -23.631012` ⇒ `J_paper = 276.0540604` ⇒
**OG = −1×10⁻⁶ < 1.92**, against `J* = 276.05406127180015`. The gap is very slightly *negative*, as
`Blasi_CellSystems2016`'s is: PyBNF found a point marginally better than the benchmark's best
recorded run.

**This is the collection's clearest case of the recipe mattering.** Holding the method fixed and
changing only what was broken or budgeted:

| configuration | `OG` |
|---|---|
| 20 × 500, before #531 — derived parameters stale, `k_syn_R_M` and `k_syn_P` inert | 96.7 at nominal, no usable fit |
| 20 × 500, after #531 | 6.76 |
| 100 × 1000, after #531, gradient still missing `k_syn_R_M`'s column (#534) | 0.104 |
| 100 × 1000, after #531 **and** #534 | **−1×10⁻⁶** |

The last row is the interesting one: the same 100 × 1000 budget goes from 0.104 to the reference
optimum purely because one of thirteen gradient columns stopped being identically zero.

**Verdict: PASS (SOLVED).**

### The fitted point is not the nominal point, and should not be

`OG_nominal = 39.85`, so unlike most solved slugs here the PEtab `nominalValue` vector is *not* this
problem's optimum and the fit is not expected to reproduce it. The noise scales make that plain —
their nominals are placeholders:

| parameter | fitted | PEtab nominal |
|---|---|---|
| `sigma_RM5` | 0.29945257 | 1 |
| `sigma_RV` | 0.078306767 | 1 |
| `sigma_RC` | 0.36824068 | 1 |
| `sigma_Vrel` | 0.13080602 | 1 |
| `scale_sigma_FracNucInt` | 0.66512969 | 1 |

The model parameters move correspondingly less (`k_syn_R_M` 30189.6 vs 30600, 1.3%; `k_imp` 0.290 vs
0.296, 2.0%) with two exceptions (`k_syn_R_V` 57%, `k_bind_M1` 15%). That pattern — noise scales
relaxing from placeholder 1 down to their MLE values while most rate constants stay near nominal — is
what a genuine likelihood optimum looks like when the shipped point was never the optimum.
`Blasi_CellSystems2016` documents the same placeholder-σ situation.

## History — two PyBNF fixes, and this slug found the second

**lanl/PyBNF#531 — the forward model.** This is a COPASI export: 27 `ModelValue_*` parameters are
fixed by `initialAssignment`s from the real names (`ModelValue_79 = k_syn_R_M`,
`ModelValue_80 = k_syn_P`, …), every rate law reads the alias, and **no source name appears in a rate
law directly**. PyBNF's fast simulation path recomputed a *species* initial fixed that way but never a
*parameter*, so every alias kept its load-time value. The fitted `k_syn_R_M` and the condition target
`k_syn_P` were both **inert** — and since `experiment____virus_infection_3` is defined by
`k_syn_P = 0`, the no-protein-synthesis experiment simulated identically to
`experiment____virus_infection`. Its nominal check moved from `OG = 96.7` to `39.9` when this was
fixed; the earlier number is not comparable with anything.

**lanl/PyBNF#534 — the gradient column.** Fixing the forward model exposed a routing gap: #530 taught
the router that a *condition target* reaches everything it seeds, but a free parameter bound **by id**
still got only its own axis. For `k_syn_R_M` that axis is identically zero, so its gradient column was
silently zero — the failure mode ADR-0095 exists to prevent, reached through a scope line in that same
ADR. The seed map already held the right term; only the by-id branch had to consult it.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`; **no hand corrections**. Exercises
  lanl/PyBNF#509 (`lnnormal`, natural-log noise) alongside linear Gaussian observables in one problem.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 100`,
  `max_iterations = 1000`, `wall_time_sim = 10`, `random_seed = 1`.
- k = 13 free parameters (8 model + 5 noise scales), n = 42 scored points, 13 observables,
  three conditions/experiments — one of which (`k_syn_P = 0`) only became distinct after #531.

## Provenance

Run against **bngsim 0.12.1** (the released PyPI wheel) with PyBNF at ADR-0094 + ADR-0095 + ADR-0096.
An earlier run of the same recipe on a locally built 0.12.0 reached `OG = 0.104`; the difference is
ADR-0096, not the simulator — the two builds agree to the last digit on this problem's nominal check
and on its 13-column FD comparison.

## Bottom line

The mixed-regime, aliased-parameter case: 13 parameters against 42 points, natural-log and linear
noise in one objective, solved from scratch at the reference optimum. It is also the slug that paid
for itself twice over as a test case — one scalar-path simulation bug and one silently-zero gradient
column, neither of which any fixture in the PyBNF suite had caught.
