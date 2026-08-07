# VALIDATION — Schwen_PONE2015

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's **weak-reference** case
and its **log10-dominated** case: `J*` is not converged for this problem, and 993.05 of its 1255.87
restored constant is the change-of-variables Jacobian rather than `(N/2)log(2π)`.

> **Confidence: 71 / 100 — the lowest in the collection, and deliberately so.** The fit is real and
> the arithmetic is sound: `OG = −12.545` from a from-scratch 100-start multi-start, no seeding, and
> the objective reproduces an independent oracle exactly. But **two of the three things an `OG` is
> normally taken to mean do not hold here.** `J*` is an unconverged reference, so clearing the
> threshold is not a demanding test; and the fit does **not** reproduce the kinetics the source paper
> publishes as its figure. Read §Gate C and the callout at the top of `README.md` before quoting this
> number. The 71 is not doubt about the computation — it is the honest distance between "solved" and
> "recovered the published fit".

## Gate A — objective fidelity

Every observable is `log10` (`lognormal`), so `−log_likelihood` carries the change-of-variables
Jacobian `Σ log(y_obs·ln10)` on top of `(N/2)log(2π)`, and here that term **dominates**:

| term | value |
|---|---|
| PyBNF reduced objective | −315.989727 |
| restored constant (`−lnL − J_reduced`) | 1255.866073 |
| of which `(N/2)·log(2π)`, N = 286 | 262.816400 |
| of which the log10 Jacobian | **993.049673** |
| `J_paper = −log_likelihood` | 939.876347 |

That convention is not assumed; it is verified by three solved slugs whose Jacobian is large —
`Perelson` (log10, +233.1287) at `OG = 5e−07`, `Blasi` (ln, −1102.0028) at `−4.3e−07`, `Laske`
(mixed, +261.0897) at `−1e−06`. If Eq. 6 omitted the term, Perelson would miss `J*` by 233 rather
than by 5e−07; a natural-log Jacobian in place of log10 would offset it by `n·log(ln10) = 13.34`,
which the same result excludes.

**Independently confirmed.** §2c's oracle recomputes the NLL at the nominal point straight from
upstream's `simulatedData` + `measurementData`, no PyBNF in the loop, and gets **943.9993** —
matching PyBNF exactly. This slug carries `obj ✓`.

**Verdict: PASS.**

## Gate B — the gradient is the one the objective implies

Covered by the lanl/PyBNF#535 sweep, which found this slug clean. **One caveat, unresolved:** that
check predates lanl/PyBNF#548, so it was performed in *linear* sampling space for the 19 parameters
that were then wrongly linear. The gradient is assembled in sampling space `u`, so it has **not**
been re-verified in the log space this conf now searches. `Zhao_QuantBiol2020`, the other #548 slug,
was re-checked and came back at 1.95e−04; there is no reason to expect worse here, but that is an
expectation rather than a measurement.

**Verdict: PASS (inherited; not re-measured post-#548).**

## Gate C — the fit reaches the benchmark optimum, and what that is worth

From-scratch multi-start `gntr` (**100 starts × 1000 iterations**, `random_seed = 1`,
`sbml_backend = bngsim`, no seeding) converges to `J_pybnf = −315.989727` ⇒
`J_paper = 939.876347` ⇒ **`OG = −12.545379` < 1.92**, against `J* = 952.4217251`. All 100 starts
retired: 64 on `step is negligible`, 36 on `reached max_iterations`. Wall clock ≈ **2 h**. Zero
`mxstep` lines — this model is not wide-dynamic-range and is untouched by lanl/bngsim#196.

**`OG` is negative by 12.5, and that is a statement about the reference, not about the fit.** `J*` is
the best objective over all Marvin runs, so a point 12.5 NLL units below it means the reference is
not converged for this problem (k=30, n=286, log10). The PEtab *nominal* point already scored −8.42.
Gate A's independent oracle is what settles which side is wrong: it reproduces PyBNF exactly, so our
objective is right and `J*` is the outlier. **Clearing a threshold against an unconverged reference
is not a demanding test**, and this row should never be cited as evidence of optimizer strength.

### The fit does not reproduce the paper's published kinetics

`observable_IR2` (high-binding hepatocytes) at the final best point, against the data and against the
PEtab nominal (published) point:

| dose | t | data | **ours** | published |
|---:|---:|---:|---:|---:|
| 100 | 1 → 30 | 4.24 → **9.25** | 5.18 → **4.83** ↓ | 3.50 → **7.74** ↑ |
| 1000 | 1 → 30 | 13.57 → **47.84** | 31.59 → **29.06** ↓ | 16.21 → **48.72** ↑ |
| 10000 | 1 → 30 | 61.28 → **136.19** | 92.23 → **92.04** flat | 73.99 → **112.79** ↑ |

The published point reproduces Fig 11's rising, saturating curves closely. Ours is flat or declining
where the data rise 2–3.5×. `observable_IR1` tracks acceptably in both. Verified at the **final**
best point, not a snapshot; `evaluate_multiple` there returns −315.989705 against the fit's
−315.989727.

**Why, and it is visible in the fitted parameters.** The measurement table binds all 286 points to
just two estimated σ — `IR_obs_std` over the 34 FACS points, `std` over the 252 ELISA points — so
with σ profiled `∂NLL/∂log(RMSⱼ) = nⱼ` and the ELISA carries ≈7.4× the leverage of the panel the
paper publishes. The optimizer took that trade, and the σ show it:

| | nominal | fitted | box |
|---|---:|---:|---|
| `IR_obs_std` (FACS, 34 pts) | 0.047186 | **0.056234** | **pinned at its upper bound** |
| `std` (ELISA, 252 pts) | 0.248324 | 0.186594 | far from either bound (1e−05 … 1000) |

The FACS noise scale ran **into its ceiling** — the fit wanted to call the FACS misfit measurement
noise and the box stopped it — while the ELISA scale tightened. That is the σ-encoding trade made
parametric.

**This is a legitimate likelihood optimum, not a defect.** The objective reproduces the oracle, the
model is properly dose-responsive, and the two-σ grouping is the Benchmarking-Initiative's encoding,
not ours — the paper itself used per-point error estimates ("the shaded error bands correspond to the
estimated error in the data points"). `J*` is defined on the *same* PEtab objective, so the benchmark
comparison is apples-to-apples.

**Verdict: PASS on the benchmark objective (SOLVED). NOT a reproduction of the paper's fit.**

### Five of thirty parameters rest on a box bound

| parameter | value | bound |
|---|---:|---|
| `IR_obs_std` | 0.0562341 | upper (see above — the load-bearing one) |
| `fragments` | 1 | upper |
| `km_nExpID3` | 1e+08 | upper |
| `scaleElisa_nExpID3` | 0.1 | lower |
| `scaleElisa_nExpID4` | 1 | upper |

Five of thirty is high, and it is a further reason to treat this optimum as constrained rather than
interior. Any profile-likelihood reading of this problem should widen those boxes first.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`; **no hand corrections**. Exercises
  lanl/PyBNF#510 (an experiment measured only at `t = 0`, which `TimeCourse` previously rejected).
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 100`,
  `max_iterations = 1000`, `wall_time_sim = 60`, `random_seed = 1`.
- k = 30 free parameters, n = 286 scored points, 4 observables, 19 conditions.
- **19 of its 24 box-sampled parameters were searched linearly until lanl/PyBNF#548** (`e008d345`).
  No fit before 2026-08-07 searched this problem correctly. Its six declared `parameterScaleNormal`
  priors correctly survive as `log-normal`; five genuinely `lin` parameters stay `uniform`.

## Provenance

Run against **bngsim 0.12.2** (released wheel, predating lanl/bngsim#196) with PyBNF at `e008d345`.
Slug renamed from `Schwen_PONE2014` on 2026-08-07 — the paper is 2015; see `README.md`.

## Bottom line

Solved on the benchmark objective at `OG = −12.5`, and the most heavily caveated row in the
collection. Its `J*` is unconverged, its optimum rests on five bounds including a noise scale at its
ceiling, and it does not reproduce the binding kinetics its source paper publishes as Fig 11. All
three are properties of the problem as encoded rather than defects in the fit — and all three are the
reason this file exists.
