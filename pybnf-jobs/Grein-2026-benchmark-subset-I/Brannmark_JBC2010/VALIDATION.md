# VALIDATION — Brannmark_JBC2010

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs), corroborated for this slug by §2c's independent
recomputation from upstream's own `simulatedData` tables (`obj ✓`). This is the last 🟢 row in the
collection, and the slug issue #38 recorded as the only one with an open upstream issue in front of
it. It did not need one.

> **Confidence: 90 / 100.** SOLVED with `OG = 0.111206` from a from-scratch 100-start multi-start, no
> seeding, in **1 h 41 m 10 s**. **87 of the 100 starts converged on `step is negligible`**, 12 died
> with `start point failed to simulate`, and 1 hit `max_iterations`. The fit lands 0.047 NLL units
> above the published optimum and **recovers all four of the box-bound coordinates that optimum sits
> on**, so this is a reproduction of the reference basin rather than a different basin that happens to
> score.
> Deductions: the model is imported, not re-derived from Brännmark et al. 2010; the run is a single
> seed; 12 of 100 starts died; the recipe needs a hand-set `sbml_atol` that no derivation in PyBNF
> will reach for on its own; one start was still descending when its budget ran out; and the
> *necessity* of that tolerance is **not** established — see "What is not established".

## What was actually blocking this slug — and it is not what was recorded

Issue #38 attributes this slug's failed run to **lanl/bngsim#196**: *"778 of its starts died on CVODES
`mxstep` … a scalar `atol` cannot serve a model that seeds a transient at `1.8e-9` against principal
species at `0.1..10`."* **A scalar serves the state solve perfectly.** Separating the two solves over
30 box points (seed 11), which is what `tools/box_probe.py` exists for:

| | integrated | wall |
|---|---:|---:|
| plain forward solve, derived tolerance | **30 / 30** | 0.9 s |
| with the `gntr` forward-sensitivity request applied | **19 / 30** | 210 s |

Every one of the 30 points integrates, in under a second, at the very tolerance the failure was
attributed to. What dies is the **forward-sensitivity** solve, whose tolerances CVODES scales from
the state ones. Confirmed point-by-point at a single box start, where the split is unambiguous:

| `sbml_atol` | state solve | sensitivity solve |
|---|---|---|
| `3.302e-10` (derived) | **ok, 0.06 s** | **FAIL** — `mxstep`, killed at `wall_time_sim` 10 s |
| `1e-08` | ok, 0.07 s | ok, 0.18 s |
| `1e-07` | ok, 0.07 s | ok, 0.14 s |

This is the same shape as `Weber_BMC2015`, and cleaner: Weber's state solve was merely *mostly* fine
(7 of 11), Brannmark's is perfect.

## Why the tolerance was wrong — and why it is the *mirror* of Weber, not another instance

ADR-0103 derives `atol = rtol × median(y₀)`. Brannmark's nine species span `1.76e-09 .. 1.00e+01`
with median `3.30e-02` (`IRi`), so the derivation returns **`3.302e-10`** — and the backend default
is `1e-08`. The clamp is `[1e-16, default]`, i.e. *only ever tighten*, so unlike Weber the derivation
is **allowed** to apply its answer here. It does, and that is the defect: it tightens this model **1.5
decades below the backend default**, and that tightening is what kills its starts.

ADR-0105's per-species vector clamps each entry into `[scalar_atol, default_atol]` = `[3.302e-10,
1e-08]`, which gives the three principal species their default back and pins the five below the
median at the derived scalar:

| species | nominal y₀ | `rtol × y₀` | vector `atol` |
|---|---:|---:|---:|
| `IRp` | 1.76e-09 | 1.76e-17 | 3.302e-10 |
| `IRiP` | 1.12e-05 | 1.12e-13 | 3.302e-10 |
| `Xp` | 1.58e-04 | 1.58e-12 | 3.302e-10 |
| `IRins` | 1.74e-02 | 1.74e-10 | 3.302e-10 |
| `IRi` | 3.30e-02 | 3.30e-10 | 3.302e-10 |
| `IRSiP` | 1.33e-01 | 1.33e-09 | 1.330e-09 |
| `IRS` | 9.87e+00 | 9.87e-08 | 1.000e-08 |
| `IR` | 9.95e+00 | 9.95e-08 | 1.000e-08 |
| `X` | 1.00e+01 | 1.00e-07 | 1.000e-08 |

The vector is **not buying survival** — 19/30 with it, 19/30 at the flat scalar `3.302e-10` — only
speed (210 s against 303 s). The mortality comes from the five species still pinned at the tight
value, i.e. from ADR-0103's deliberate tightening, which ADR-0105 by construction cannot undo.

**This matters for lanl/PyBNF#557.** That issue predicts Brannmark will be the worked example for its
**ask (b)** — a per-species vector, because "a single number cannot serve `IRp` at `1.8e-09`
alongside `X` at `10`." The measurement says otherwise: a single number serves it fine, and the
vector that already exists does not help. What Brannmark demonstrates is the **mirror of ask (a)** —
Weber shows the derivation refusing to loosen when the model asks it to; Brannmark shows it
tightening below the backend default when it should not. Same rule, two directions, and the second
one has no off-switch short of writing the default back in by hand.

## Gate A — objective fidelity

Four estimated σ, all-linear observables, so the restored constant is `(N/2)·log(2π)` and nothing
else. `noise_model IRS1_P = gaussian, sigma = formula noiseParameter1_IRS1_P` looks like a fifth
regime but is an indirection: the two `_measparams.tsv` sidecars bind it to `sigmaY2Step` (TwoSteps)
and `sigmaY2TimR` (Dose_100), both estimated.

The constant is parameter-independent, measured at **four** parameter vectors spanning six orders of
magnitude in objective — this is the check the shortcut is worth nothing without:

| point | `OG` there | reduced objective | `J_paper` | difference |
|---|---:|---:|---:|---:|
| PEtab nominal | 0.064 | 102.374846217 | 141.889203145 | 39.514356927801 |
| nominal + 3% of each box width | 2384.8 | 2487.155887338 | 2526.670244266 | 39.514356927801 |
| nominal + 10% of each box width | 37373.4 | 37475.737871887 | 37515.252228815 | 39.514356927801 |
| **this fit's best** | 0.111 | 102.421703 | 141.936060 | 39.514357 |

Identical to every printed digit, and `43/2 · log(2π) = 39.5144` independently. So for this slug
**`OG = reduced_objective − 102.310497`**, and "solved" means `reduced < 104.230497`. (Use `score.py`
for anything reported; the shortcut is for watching a run.)

The tolerance change does not move the objective. `J_paper` at the PEtab nominal point:

| `sbml_atol` | `J_paper` | `OG_nominal` |
|---|---:|---:|
| derived vector | 141.889243220 | 0.064389 |
| **`1e-08`** (shipped) | **141.889203145** | **0.064349** |
| `1e-07` | 141.889308498 | 0.064454 |

Stable to 7 significant figures. `nominal_check.json` was regenerated under the shipped conf and
still reports the nominal point as the optimum (`OG_nominal = 0.06434887`, against `0.06436871`
before). This is not an accuracy trade.

**The pre-equilibration risk was the live one and it is addressed.** An explicit `sbml_atol` also
becomes the steady-state convergence cutoff for this problem's `preequilibrate:` phase — ADR-0105
pairs the two and PyBNF exposes no separate key. With `IRp` at `1.76e-09`, a cutoff that declared
equilibrium early would be the ADR-0086/0104 failure mode, and it would be silent. The evidence it
converges is that `J_paper` above is reachable *only* through a correct pre-equilibration (ADR-0104
is what fixed this slug's objective from 1531 to 0.064), and it is unmoved across the sweep.

**Verdict: PASS.**

## Gate B — the gradient is the one the objective implies

This gate is the reason the fitting budget was not spent immediately. ADR-0103 tightened this model
*on purpose*, to give `IRp` usable resolution, and lanl/bngsim#196 records **worst gradient column
`7.38e-05` at `atol = 3.3e-10`**. If loosening to `1e-08` degraded that materially, integrability and
gradient fidelity would be in genuine tension and no scalar would be the right answer.

### First, a methodological correction — the `h` ladder was measuring the wrong thing

`fd_check.py` clips its evaluation point to `[lo + 8h, hi − 8h]` so that no parameter sits on a bound.
**That inset moves with `h`.** Four of Brannmark's 22 parameters sit *exactly* on a box bound at the
PEtab nominal point — `k1d` and `k1f` at the upper `5e+05`, `k1e` at the lower `1e-06`,
`k_IRP_1Step` at the upper `2e+05` — so `--disp 0.03` pushed all four out of the box and the clip
returned them to `bound ∓ 8h`, a *different point on every rung*. The symptom is an **assembled**
gradient that moves with `h`:

| `h` | 1e-4 | 3e-4 | 1e-3 | 3e-3 |
|---|---:|---:|---:|---:|
| assembled `k_IRSiP_2Step` | 904.2 | 924.2 | 994.6 | 1202.1 |

Read naively that is a 33% tolerance-induced gradient error. It is nothing of the sort — it is the
evaluation point sliding. Every number below therefore comes from a point **pinned once**: displaced
in sampling space, inset a fixed 5% of each box width (comfortably clear of `8h` for any `h ≤ 3e-3`),
written to a `param-values.json`, and run with `--disp 0`. The assembled gradient is then **bit-identical
across the whole `h` ladder**, which is what makes "the assembled side holds while the central
difference wanders" a statement about FD noise rather than about clipping. Recorded in
`tools/README.md`.

### The assembled gradient is invariant to the tolerance

Against `3.302e-10` as reference, at `h = 3e-4`, at three points:

| point | ‖g‖ | largest column drift at `1e-08` | at `1e-07` | `1 − cos` at `1e-08` / `1e-07` |
|---|---:|---:|---:|---|
| nominal + 3% | 2.76e+04 | 3.6e−07 ‖g‖ | 2.5e−06 ‖g‖ | 2.1e−13 / 7.0e−12 |
| nominal + 10% | 3.83e+05 | 5.2e−07 ‖g‖ | 2.6e−06 ‖g‖ | 1.6e−13 / 2.3e−12 |
| a box-sampled start | 6.82e+14 | — (reference) | 4.3e−10 ‖g‖ | — / 0.0 |

Three decades of tolerance move the gradient direction by one part in 10¹². **At the box start the
tight tolerance yields no gradient at all** — the sensitivity solve fails, so the reference column is
`1e-08`. That is the whole slug in one row: at a real start point, the tolerance ADR-0103 chose does
not produce a gradient, and the one it clamps away does.

### What degrades is the objective, and by how much

The FD *reference* is a difference of objectives, so it carries the tolerance's noise. Worst relative
error over all 22 columns:

| point | `sbml_atol` | `h`=1e-4 | 3e-4 | 1e-3 | 3e-3 |
|---|---|---:|---:|---:|---:|
| nominal + 10% | `3.302e-10` | 3.91e−04 | **6.92e−05** | 2.22e−04 | 1.28e−04 |
| | **`1e-08`** | 1.62e−03 | **6.55e−04** | 2.74e−04 | 5.37e−05 |
| | `1e-07` | 6.79e−03 | 4.67e−03 | 2.94e−03 | 7.43e−03 |
| a box start | `3.302e-10` | — | — | — | — (no gradient) |
| | **`1e-08`** | 3.77e−02 | 7.18e−03 | 4.71e−03 | 2.47e−03 |
| | `1e-07` | 2.99e−01 | 5.03e−02 | 2.34e−02 | 9.99e−03 |

The `6.92e−05` at the derived tolerance reproduces bngsim#196's recorded `7.38e-05` to the same order,
which is the check that this ladder is measuring the same thing that issue did. **`1e-08` costs about
an order of magnitude of FD precision and `1e-07` costs two**, which is why `1e-08` is shipped:
it recovers 28 of 30 box points (against 19) for a degradation that stays below `1e-03`, and the last
two points are not worth the further 7×. This is a trade against *objective* noise, which is the
quantity a trust-region line search consumes — the same reading as Weber, in the same direction.

The residual flags are FD noise by both of the corpus's tests. **The worst column's identity changes
with `h`** at every tolerance and every point (`k1d` → `k1e` → `k1a` → `k1b` at the derived value;
`k1g` → `k1e` → `k1r` → `k_IRSiP_2Step` at `1e-08`), and **the offenders sit 3 to 12 decades below the
gradient's infinity norm** — `k1e` reads `0.041` against a norm of `2.26e+05`. Critically, the
*derived* tolerance flags the same columns: these are artefacts of differencing a near-zero column,
not of the tolerance.

**Verdict: PASS.**

## Gate C — the fit reaches the benchmark threshold

From-scratch multi-start `gntr` (**100 starts × 1000 iterations**, `random_seed = 1`,
`sbml_backend = bngsim`, `sbml_atol = 1e-8`, **no seeding**) reaches `J_pybnf = 102.421703` ⇒
`J_paper = 141.936060` ⇒ **`OG = 0.111206` < 1.92**, against `J* = 141.82485427243665`. Wall clock
**1 h 41 m 10 s** on ten cores.

| start outcome | count |
|---|---:|
| converged, `step is negligible` | **87** |
| died, `start point failed to simulate` | 12 |
| `reached max_iterations` | 1 |

One start was still descending when its budget ran out, so this run is *marginally* budget-limited —
but 87 of the 100 retired on a negligible step and the best point came from a converged start, so
raising `max_iterations` is not what stands between this result and a better one.

### The tolerance is what bought the budget

Measured on 30 box points (seed 11) with the sensitivity request applied — the same thing a start
does — re-measured at this session's HEAD rather than inherited:

| `sbml_atol` | integrated / 30 | wall |
|---|---:|---:|
| derived vector `[3.3e-10 … 1e-8]` | 19 (63%) | 210 s |
| `3.302e-10` (vector off) | 19 (63%) | 303 s |
| **`1e-08`** (backend default) | **28 (93%)** | 99 s |
| `1e-07` | 30 (100%) | 20 s |
| `1e-06` | 30 (100%) | 5 s |

The run bears the prediction out in the right direction and slightly worse: **12 of 100 starts died**
against the 7% the 30-point sample forecast.

### The error-line count is not the budget loss, and #38 reads it as if it were

This run emitted **331,402** CVODES `mxstep` lines on stderr — *twice* the 164,439 that #38 cites as
evidence of how badly the previous attempt was hurt — while losing only 12 starts and 74 simulations.
The two numbers measure different things: almost every `mxstep` line is a trial point inside a line
search that the trust region rejects and recovers from. A stderr line count is not a budget
measurement, and neither is `ls output/FailedSimLogs | wc -l` (13 files against 12 dead starts here).
The log's own `GNTR start N/100 stopping:` lines are.

Of the 74 `SimulationError`s, 37 are `CVODE integration failed at t=1000000` — the pre-equilibration
horizon — so the surviving failures are concentrated in the equilibration phase rather than the
measured window. **No** failure in this run was a step collapse and **none** tripped `wall_time_sim`,
so `Weber_BMC2015`'s second defect (lanl/bngsim#305) does not appear here.

### The published optimum is a corner, and unlike Weber this fit recovers it

**Four of Brannmark's 22 parameters sit exactly ON their box bounds at the published optimum** —
`k1d` and `k1f` at the upper `5e+05`, `k_IRP_1Step` at the upper `2e+05`, `k1e` at the lower `1e-06`.
That is the same structure that held `Weber_BMC2015` 0.78 NLL units short of its reference basin, and
it is the first thing to check when a multistart lands near but not on a published point.

This fit recovers **all four**: `k1d` at `500000` exactly, `k1f` at `499950` (99.99% of the bound),
`k_IRP_1Step` at `199992` (99.996%), and `k1e` at `1.088e-06` — 0.3% of the box width off its lower
bound in log space. Weber recovered one of five and paid 0.78; Brannmark recovers four of four and
pays 0.047. That is the causal story behind the difference in the two OGs, and it is why this row
carries no `(not saturated)` caveat.

### It is not a leverage trade — the §2f check

`Schwen_PONE2015` is the caution that a good `OG` can be bought by sacrificing the observable the
source paper publishes. Brannmark does not do that. Per experiment/observable, nominal vs this fit:

| experiment | observable | n | −lnL nominal | −lnL fit | gap |
|---|---|---:|---:|---:|---:|
| Dose_100 | `IR1_P` | 12 | 30.0812 | 31.5287 | **+1.4476** |
| Dose_10 | `IRS1_P_DosR` | 1 | 4.5309 | 4.6462 | +0.1153 |
| Dose_1 | `IRS1_P_DosR` | 1 | 4.4184 | 4.5304 | +0.1120 |
| Dose_100 | `IRS1_P_DosR` | 1 | 3.9460 | 4.0570 | +0.1110 |
| Dose_01 | `IRS1_P_DosR` | 1 | 3.8168 | 3.8849 | +0.0681 |
| Dose_3 | `IRS1_P_DosR` | 1 | 4.9217 | 4.9671 | +0.0453 |
| Dose_0 | `IRS1_P_DosR` | 1 | 4.0738 | 4.0022 | −0.0716 |
| Dose_001 | `IRS1_P_DosR` | 1 | 4.2960 | 4.2133 | −0.0827 |
| Dose_100 | `IRS1_P` | 12 | 45.0959 | 44.9308 | −0.1652 |
| TwoSteps | `IRS1_P` | 12 | 36.7084 | 35.1756 | **−1.5328** |
| **total** | | **43** | **141.8892** | **141.9362** | **+0.0470** |

The two largest movements nearly cancel — the fit is 1.45 worse on the `IR1_P` dose response and 1.53
*better* on the `IRS1_P` two-step time course — and the whole net gap is 0.047 NLL units. No
observable is abandoned; the fit beats the published point on 3 of the 10 groups. The four estimated
σ come back within 14% (`sigmaY1TimR` 2.967 → 3.358, `sigmaY2Step` 5.154 → 4.523, `sigmaY2TimR`
10.368 → 10.228, `sigmaYDosR` 17.578 → 18.373), so the noise model is recovered and the residual gap
is kinetic.

Most rate constants land within a factor of 2 of the published values. The one outlier is **`k1g`,
1729 → 33.4 (52×)** — and it is the parameter the gradient analysis already identified as barely
identified: its assembled column reads `−0.89` against a gradient norm of `2.8e+04`, i.e. six decades
down. A 52× move along a direction the objective cannot see is a flat manifold, not a disagreement.

**Verdict: PASS (SOLVED).**

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`; **no hand corrections**. All 22 estimated
  parameters are `parameterScale = log10` upstream and import as `loguniform_var`.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 100`,
  `max_iterations = 1000`, `wall_time_sim = 10`, `random_seed = 1`, **`sbml_atol = 1e-8`**.
- k = 22 free parameters, n = 43 scored points.
- `wall_time_sim` stays at 10; nothing in this run tripped it.
- The conf carries the measurement table above in a comment. `sbml_atol = 1e-8` is the *backend
  default*, so the line reads "do not let ADR-0103 tighten this model below it" — the minimal
  intervention, and the reason no looser value was chosen even though two of them integrate 30/30.

## Provenance

Run against **bngsim 0.12.2** built from source at `ffbf015`, with PyBNF at `095a5a14`. bngsim's
import-time build check reported `fresh` (`built=ffbf0156eea4 HEAD=ffbf0156eea4`) before any
measurement recorded here — a version check is not sufficient for this, and a stale core has bitten
this corpus twice (see `Weber_BMC2015/VALIDATION.md` § Provenance and lanl/PyBNF#558). As a
cross-check, `Bruno_JExpBot2016` — the corpus's designated regression guard — reproduces its recorded
`J_paper = −46.688194686350265` as `−46.688194686350730`, agreeing to 14 significant digits.

## What is **not** established

- **That the tolerance was *necessary*.** At the derived value 19 of 30 box points still integrate
  (63%), and the corpus has solved slugs on less (`Raia` converged on 83/100, `Weber` on 76/100). The
  control — a full 100 × 1000 run at the derived tolerance — was **not** run, so the honest claim is
  that `1e-08` recovers 63% → 93% of the multistart budget and that the fit at `1e-08` solves. Whether
  the derived tolerance would also have solved is open. The `OG = 3.21` figure #38 records is *not*
  evidence either way: it predates ADR-0105 and lanl/PyBNF#553, and it was never scored.
- **That the pre-equilibration cutoff is safe across the whole box.** It is verified at four parameter
  vectors (Gate A), not at every start. A box point with a slower transient could still terminate its
  equilibration early, and the symptom would be a suspiciously *good* objective. This is the first
  thing to re-check if a future run on this slug reports one.
- **Single seed.** `random_seed = 1`, one run. No claim about run-to-run variance.
- **The model is imported, not re-derived** from Brännmark et al. 2010.

## Bottom line

The last 🟢 row, and the one #38 listed as blocked behind an open upstream capability issue, solved by
one config line — `sbml_atol = 1e-8`, the backend default. Its integration failure was in the
**sensitivity** solve rather than the state solve (the state solve integrates 30 of 30 box points in
0.9 s at the tolerance blamed for the failure), and the tolerance that fixes it is the one ADR-0103's
derivation *overrode* — the mirror image of `Weber_BMC2015`, where the same rule refused to loosen.
The gradient is verified for the first time at a point where the check is a real test, and is
invariant across three decades of tolerance to one part in 10¹². `OG` goes from "no fit has ever been
scored" to **0.111206**, from unbiased box-sampled starts, recovering all four coordinates of a
published optimum that sits in a corner of the fit box.
