# VALIDATION — Weber_BMC2015

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs), corroborated for this slug by §2c's independent
recomputation from upstream's own `simulatedData` tables (`obj ✓`). This is the slug
lanl/bngsim#196 named as **blocked outright**, and the block turned out to be one config line rather
than a missing capability.

> **Confidence: 88 / 100.** SOLVED with `OG = 0.781167` from a from-scratch 100-start multi-start, no
> seeding, in **21 m 25 s**. **Solved but not saturated**, in the `Fiedler_BMCSystBiol2016` sense: the
> published optimum scores `OG_nominal = -0.0002`, and this fit lands 0.78 NLL units above it — inside
> the threshold without reproducing the reference basin. That gap is structural rather than a budget
> shortfall (see Gate C).
> Deductions: the model is imported, not re-derived from Weber et al. 2015; the run is a single seed;
> 24 of 100 start points died on a *second* integration defect (lanl/bngsim#305, below); two fitted
> parameters rest on box bounds; and the recipe needs a hand-set `sbml_atol` that no derivation in
> PyBNF will reach for on its own.

## What was actually blocking this slug

Issue #38 and lanl/bngsim#196 both record Weber as a model that "will not integrate anywhere near"
its nominal point, at every displacement tried including zero. That is true of the runs they
measured, and it is the wrong attribution. Separating the two solves, at the shipped
`sbml_atol = 1e-8`, over the nominal point plus 10 box-sampled points:

| | integrated | wall |
|---|---|---|
| plain forward solve | **7 / 11** | 0.6 s total |
| with the `gntr` forward-sensitivity request applied | **2 / 11** | 80 s total |

The state solve is largely fine and fast. What fails is the **forward-sensitivity** solve — and it
fails *because of the state tolerance*, since CVODES derives its sensitivity tolerances from the
state ones. Every probe in #38 that "failed even at zero displacement" was running through
`tools/fd_check.py`, which applies the sensitivity request; that is why the failure looked
unconditional.

This also explains bngsim#196's own puzzled footnote that the failure "does not reproduce on a bare
model": a bare `Model.from_sbml` + `Simulator.run` has no sensitivity axis to fail on.

## Why the tolerance was unavailable

Weber's seven species sit at `1.24e+02 .. 4.21e+07`, median `4.67e+05`:

| species | nominal y₀ | `rtol * y₀` |
|---|---:|---:|
| `PKDDAGa` | 1.239e+02 | 1.2e-06 |
| `CERT` | 1.608e+05 | 1.6e-03 |
| `PI4K3Ba` | 3.321e+05 | 3.3e-03 |
| `PKD` | 4.665e+05 | 4.7e-03 |
| `PI4K3B` | 1.578e+06 | 1.6e-02 |
| `CERTERa` | 3.195e+07 | 3.2e-01 |
| `CERTTGNa` | 4.208e+07 | 4.2e-01 |

ADR-0103 derives `atol = rtol × median`, which here is **`4.665e-03`** — and then clamps it into
`[1e-16, backend default]`, where the upper end is `1e-8`. The clamp "only ever tightens", so the
derivation computes a usable answer for this model and discards it in favour of one **5.7 decades
tighter**. ADR-0105's per-species vector cannot help either: its entries clamp into
`[scalar_atol, default_atol]`, which for Weber is `[1e-8, 1e-8]`, so the vector is elementwise the
scalar and correctly declines to engage.

That is exactly the tension lanl/bngsim#196 was filed about, and bngsim has since grown both the
per-species vector (bngsim#211) and a trajectory-following `CVodeWFtolerances` atol
(bngsim#213/#258). **Neither is reachable from PyBNF today** — `sbml_atol` is `Optional[float]`, and
PyBNF contains no reference to `TrackingAtol` — so this slug takes the documented scalar off-switch
instead; the clamp is now filed as **lanl/PyBNF#557**, which records that it binds on 10 of the 22
slugs whose nominal state is readable, not just this one. **No upstream change was needed to unblock
it**, which is the headline: #38 listed Weber
under an open upstream capability issue, and it was a config line.

## Gate A — objective fidelity

Mixed estimated + `fix_at` σ (5 fitted, 3 fixed), all-linear observables, so the restored constant is
`(N/2)·log(2π)` plus the fixed-σ terms. Measured at **two very different parameter vectors**:

| point | reduced objective | `J_paper` | difference |
|---|---:|---:|---:|
| PEtab nominal | −69.315679 | 296.201795 | 365.517474 |
| this fit's best | −68.534305 | 296.983169 | 365.517474 |

Identical to all printed digits, which is the check that the dropped constant really is
parameter-independent. So for this slug **`OG = reduced_objective + 69.31548`**, and "solved" means
`reduced < −67.396`.

The tolerance change does not move the objective. `J_paper` at the PEtab nominal point, swept:

| `sbml_atol` | `J_paper` | `OG_nominal` |
|---|---:|---:|
| `1e-8` (old derived value) | 296.2017966 | −0.000206 |
| `1e-6` | 296.2017988 | −0.000204 |
| `1e-5` | 296.2017992 | −0.000203 |
| **`1e-4`** (shipped) | **296.2017950** | **−0.000207** |
| `1e-3` | 296.2018850 | −0.000117 |
| `1e-2` | 296.2018044 | −0.000198 |

Four decades of tolerance move the objective in the 6th decimal. `nominal_check.json` was
regenerated under the shipped conf and still reports the nominal point as the optimum
(`OG_nominal = −0.00020743`, against `−0.00020589` before). This is not an accuracy trade.

**Verdict: PASS.**

## Gate B — the gradient is the one the objective implies

#38 records Weber's gradient as **unverified**, and was right to: the pre-#547 reading of 5.24e−06 was
computed against trajectories where every experiment simulated identically (so it proved nothing), and
the only check that completed afterwards read 5.74e−02 at a tolerance loose enough that FD noise
plausibly dominated. **This is the first check on this slug that both completes and means something.**

At `sbml_atol = 1e-4`, evaluated 3% of each box width off the optimum — the nominal point *is* the
optimum, where both sides of an FD check are noise, which is `tools/README.md`'s own caution:

| `h` | worst relative error | worst column |
|---|---:|---|
| 3e−4 | **2.78e−03** | `s31` |
| 1e−3 | 1.21e−03 | `a22` |
| 3e−3 | 3.60e−03 | `pu5` |

All 36 columns carry magnitude; none is structurally zero. Two things make this FD noise rather than
a defect:

1. **The worst-column identity moves with `h`** while every *assembled* value holds steady
   (`a11`: −2492.85 / −2492.85 / −2492.75; `a31`: −807.346 / −807.342 / −807.332). Per
   `tools/README.md`, an assembled gradient that moved with `h` would mean the evaluation point was
   being clamped; it is not.
2. **The offending columns sit ~6 decades below the gradient norm.** `a22` is −0.0032, `a32` −0.024,
   `s31` −18.5, against a norm dominated by columns near 1e+04 (`pu3` = −10202.8, agreeing to
   8.3e−08). `a22`'s central difference even flips sign across the ladder
   (−0.0091 → +0.0091 → +0.00077), which is what noise on a 3e−03 column looks like.

The assembled gradient is also **invariant to the tolerance** (`a11` reads −2492.85 at `1e-4` and
−2492.75 at `4.665e-3`). What degrades at a looser tolerance is the *objective*, i.e. the FD
reference — which is exactly the quantity a trust-region line search consumes, and the reason `1e-4`
was chosen over looser values that integrate more box points:

| `sbml_atol` | FD worst | FD on `a11` (largest column) |
|---|---:|---:|
| `1e-5` | 1.39e−03 | — |
| **`1e-4`** | 2.78e−03 | **6.85e−06** |
| `1e-3` | 1.61e−02 | — |
| `4.665e-3` | 8.60e−02 | 1.92e−04 |

**Verdict: PASS.**

## Gate C — the fit reaches the benchmark threshold

From-scratch multi-start `gntr` (**100 starts × 1000 iterations**, `random_seed = 1`,
`sbml_backend = bngsim`, `sbml_atol = 1e-4`, no seeding) reaches
`J_pybnf = −68.534305` ⇒ `J_paper = 296.983169` ⇒ **`OG = 0.781167` < 1.92**, against
`J* = 296.2020025`. Wall clock **21 m 25 s** on ten cores. Of the 100 starts, **75 converged on
`step is negligible`**, 24 died with `start point failed to simulate`, and 1 stopped because the
trust-region subproblem could not be factorized (LAPACK failed on the augmented Jacobian).
**None** hit `max_iterations`, so this run was not budget-limited.

### The tolerance is what bought the budget

Measured on 30 box points (seed 11) with the sensitivity request applied — the same thing a start does:

| `sbml_atol` | integrated / 30 | wall |
|---|---:|---:|
| `1e-8` (old derived value) | **6** (20%) | 239 s |
| `1e-5` | 19 (63%) | 64 s |
| **`1e-4`** | **22** (73%) | 25 s |
| `1e-3` | 25 (83%) | 15 s |
| `4.665e-3` | 28 (93%) | 5 s |

The run bears the prediction out: **24 of 100 start points died**, against the 27% the probe forecast. At
the shipped `1e-8` it would have been ~80%, which is why #38's standing instruction was "do not spend
a fitting budget here yet."

### The published optimum is a corner solution

**Five of Weber's 36 parameters sit exactly ON their box bounds at the published optimum** —
`a22`/`a32`/`a33` at their lower bound `1e-4`, `m11` at its upper bound `1e+10`, `pu3` at its upper
bound `1e+08`. #38 lists those five (as the ones `fd_check` had to nudge clear) but does not draw the
consequence: the reference basin is a **corner** of the box, not an interior point, and a box-uniform
multistart under a trust region reaches an interior optimum far more readily than one requiring
several active constraints at once.

This fit recovers **one** of those five coordinates — `m11` at `9.9995e+09`, within 0.005% of its
bound — and leaves the other four in the interior. That is the structural reason for the 0.78 gap.

### It is not a leverage trade — the §2f check

`Schwen_PONE2015` is the caution that a good `OG` can be bought by sacrificing the observable the
source paper publishes. Weber does **not** do that. Per-observable, nominal vs this fit:

| observable | exp | n | −lnL nominal | −lnL fit | gap |
|---|---|---:|---:|---:|---:|
| `yCERTt` | data2 | 3 | 60.8398 | 61.4623 | **+0.6225** |
| `yCERTt` | data3 | 3 | 60.7091 | 60.8712 | +0.1621 |
| `yPI4K3BpRN24` | data2 | 29 | 15.5901 | 15.7470 | +0.1569 |
| `yPKDpN0` | data2 | 9 | 8.2511 | 8.2635 | +0.0124 |
| `yPKDpN25` | data2 | 23 | −16.6212 | −16.6188 | +0.0024 |
| `yPI4K3Bt` | data2 | 3 | 43.8459 | 43.8459 | 0.0000 |
| `yPI4K3Bt` | data3 | 3 | 43.8459 | 43.8459 | 0.0000 |
| `yPKDt` | data3 | 3 | 46.1455 | 46.1455 | 0.0000 |
| `yPKDpN24` | data2 | 9 | 13.4657 | 13.4553 | −0.0104 |
| `yPKDpN0` | data3 | 30 | 46.0641 | 46.0400 | −0.0241 |
| `yCERTpRN24` | data3 | 20 | −25.9340 | −26.0747 | **−0.1407** |
| **total** | | **135** | **296.2018** | **296.9832** | **+0.7814** |

No observable is abandoned: the largest single gap is 0.62 NLL units on a 3-point observable carried
by a fixed σ of 2.4e+08, three observables are reproduced to four decimals, and the fit is *better*
than the published point on three others. **The five fitted σ are recovered to within 0.7%**
(`std_yCERTpRN24` 0.066162 → 0.065685; `std_yPI4K3BpRN24` 0.41422 → 0.41651; `std_yPKDpN0` 0.9741 →
0.97388; `std_yPKDpN24` 1.0803 → 1.079; `std_yPKDpN25` 0.11747 → 0.11747), so the noise model is
recovered and it is the kinetics that sit elsewhere. With 36 parameters against 135 points that is
the expected signature of a flat, structurally non-identifiable manifold, not of a wrong fit.

**Verdict: PASS (SOLVED, not saturated).**

### Two fitted parameters rest on bounds

`m11` lands at `9.9995e+09` against its upper bound `1e+10` (99.995%) — the same bound the published
point occupies, so this is agreement with the reference, not an artefact. `scale_yPI4K3BpRN24` lands
at `999.77` against its upper bound `1000` (99.977%), a 295× increase over its nominal `3.3937`; it is
an observable scale factor compensating a smaller `PI4K3Ba/(PI4K3B + PI4K3Ba)` ratio, and its
observable's likelihood is only 0.157 worse than nominal. Both should be read as bounded rather than
as point estimates, and any profile-likelihood work on this problem should widen those boxes first.

## A second, still-open defect this run exposed

Loosening `atol` removed one failure mode and made a **different** one visible. Counting every failed
integration across the run — not just the 24 that happened to land on a start point and kill it:

| mode | count in this run | moves with `rtol`/`atol`? |
|---|---:|---|
| CVODES `mxstep steps taken before reaching tout` | 4164 stderr lines, 27 `wall_time_sim` trips | **yes** — this is bngsim#196, and `1e-4` fixes it |
| `CVODE made no progress … the step size has collapsed` | 168 occurrences → 84 `SimulationError` | **no** |

Every single step-collapse failure is at `t = 23.999999999999996` (or within 2 ulp of it) — which is
exactly `PdBu_time = 24` and `kb_NB142_70_time = 24`, the discontinuity in this model's
`(time() - PdBu_time) < 0` rate law. bngsim's own error text names the shape: *"typically at a
discontinuity such as an `if(t >= sigma)` rate jump, where `t + h == t`… move the discontinuity onto
an event or a sample time so the integrator can restart across it."*

**This is not bngsim#194, and the root is not missing** — two things I asserted before checking, and
both are wrong. #194 is **closed**, and its subject is a *state*-threshold piecewise ("the state twin of
the GH #72 time roots"); Weber's is a **time** threshold. More to the point, asking the loader directly,
`_collect_time_discontinuity_conditions` over Weber's assignment rules returns
`((time()-PdBu_time)<0)` and `((time()-kb_NB142_70_time)<0)` — **both crossings are registered as
discontinuity triggers.** `(time - p) < 0` survives despite `time` not being bare on either side,
because root registration uses `_make_time_relational_filter`, which admits a relational with *either*
side time-dependent (bngsim GH #259); the bare-symbol requirement in `_clock_threshold_split` governs
the *sensitivity-compensation* path, not root registration.

So what remains is narrower, and is filed as **lanl/bngsim#305**: **a registered time root at which the
integrator still cannot advance.** Two mechanisms are not separated here — a root re-detected at its own
restart point, refining forever just below the crossing; or a post-jump RHS stiff enough at these
box-sampled points that the step genuinely collapses (`u5` goes `0 → PdBu_dose` discontinuously). Only
the first would be a solver bug. The model declares **0 events**, so bngsim's suggested remedy is
available but untaken. Either way the reproducer is this slug's shipped conf: 84 occurrences in 21
minutes across 100 starts.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`; **no hand corrections**. All 36 estimated
  parameters are `parameterScale = log10` upstream and import as `loguniform_var`.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 100`,
  `max_iterations = 1000`, `wall_time_sim = 10`, `random_seed = 1`, **`sbml_atol = 1e-4`**.
- k = 36 free parameters, n = 135 scored points.
- `wall_time_sim` stays at 10. #38 established that 300 s fails identically to 10 s, and that remains
  true: the surviving failures collapse the step size rather than run slowly.
- **One knock-on worth recording:** an explicit `sbml_atol` also becomes the steady-state convergence
  cutoff for this problem's `preequilibrate:` phase — ADR-0105 pairs the two and PyBNF exposes no
  separate key. `1e-4` against states of `1e+02..1e+07` is still a stringent equilibrium, and Gate A's
  nominal `J_paper` — reachable *only* through a correct pre-equilibration, per ADR-0104 — is the
  evidence it converges.

## Provenance

Run against **bngsim 0.12.2** built from source at `114d3b3`, with PyBNF at `095a5a14`.
The installed bngsim core was **stale** on first inspection (compiled binary older than
`src/_bngsim_core.cpp`, which bngsim's own import-time check reports as
`STALE … any correctness verdict drawn from it is a statement about OLD code`); it was rebuilt
before any measurement recorded here was taken. As a cross-check that the rebuild moved nothing,
`Bruno_JExpBot2016` — the corpus's designated regression guard — reproduces its recorded
`J_paper = −46.688194686350265` as `−46.688194686350730`, agreeing to 13 significant digits.

## Bottom line

The slug #38 listed as blocked by an open upstream capability issue, solved by one config line. Its
integration failure was in the **sensitivity** solve rather than the state solve, and the tolerance
that fixes it is the one ADR-0103's own rule computes for this model and then clamps away. `OG` goes
from "no fit has ever been run" to **0.781167**, with the gradient verified for the first time at a
tolerance where the check is a real test. It stops short of the published corner optimum, and the
remaining 24% start mortality is a *second* defect (lanl/bngsim#305) that this run isolated
to `t = 24⁻` — the model's own dose discontinuity, at a crossing bngsim already registers a root for.
