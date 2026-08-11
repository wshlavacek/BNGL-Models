# VALIDATION — Elowitz_Nature2000

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). The model is the repressilator of

> Elowitz MB, Leibler S. **"A synthetic oscillatory network of transcriptional regulators."**
> *Nature* **403**, 335–338 (2000). <https://doi.org/10.1038/35002125>

> **Confidence: 88 / 100.** SOLVED with `OG = 0.000175` from a from-scratch multi-start over the
> full PEtab box — **no seeding** — in **54 minutes** and 1,100 unbiased starts. It reproduces `J*`
> to `1.75e−04` and beats its own nominal point (`OG_nominal = 2.4324`). Deductions: the model is
> imported, not re-derived; **3 of 21 parameters rest on a box bound**, one of them (`tau_mRNA` at
> `1e-05`) implausible as biology; and — the reason this is not higher — the entire 2.43 NLL
> improvement over the nominal point is a **4 % reduction in the fitted noise scale**, not a
> qualitatively better trajectory (see Gate C). The two curves are near-indistinguishable by eye.

## Gate A — objective fidelity

A **log10**-transformed observable with **normal** noise and an **estimated** scale
(`sigma = fit sigma`), so the constant restoring PyBNF's reduced objective to Eq. 6 carries a
change-of-variables Jacobian on top of `(N/2)log(2π)`. From the fit's own
`information_criteria.txt`:

| term | value |
|---|---|
| PyBNF reduced objective | −126.163175 |
| restored constant (`−lnL − J_reduced`) | 60.528230 |
| `(N/2)·log(2π)`, N = 58 | 53.298434925871 |
| log10 Jacobian sum, `Σ log(y_obs·ln10)` | 7.229794846622 |
| sum of the two | **60.528229772493** |
| `J_paper = −log_likelihood` | −65.634945 |

The decomposition is reproduced **from upstream data alone** — see the oracle below — and agrees
with the recorded constant `60.52822977249341` to **12 decimals**. Note this is why the collection
README's restored-constant table, which is the σ-source part only, lists `53.298` for this slug and
not `60.528`.

### This slug DOES have a §2c independent oracle, and it reproduces

`tools/sigma_profile.py` reports `joined 0 of 58` for this slug, and an earlier draft of the kickoff
read that as "its upstream rows do not join." **That is wrong**, and this is the first time the
oracle has actually been run here.

The rows join **58 of 58, one-to-one**, on the PEtab identity key
`(observableId, simulationConditionId, preequilibrationConditionId, time)`. The tool additionally
keys on `observableParameters` and `noiseParameters`, which hold different *kinds* of value in the
two tables — parameter **names** in `measurementData` (`background;scale`, `sigma`) against
**resolved numeric values on `parameterScale`** in `simulatedData`
(`-4.98107438218408;-0.279017032524776 `, `-1.14362512938604`, note the trailing space). Those can
never match. This is the same class of over-specification as the `datasetId` bug fixed in `7902faf`,
two columns over.

Recomputing the Eq. 6 NLL at the PEtab nominal point straight from the upstream tables — **no PyBNF
in the loop** — with each point contributing
`0.5log(2π) + log σ + 0.5((log10 y_obs − log10 y_sim)/σ)² + log(y_obs·ln10)`:

| | `J_paper` at the PEtab nominal point |
|---|---|
| oracle, from upstream `simulatedData` | −63.2027999142 |
| PyBNF, same point, same build | −63.2027400751 |
| difference | 5.98e−05 absolute, **9.47e−07 relative** |

That residual is upstream's `simulatedData` being stored to finite precision, not a disagreement.
**Elowitz_Nature2000 earns `obj ✓`.**

> The shared tool is deliberately **not** patched here. Dropping those two columns from its join key
> is *not* the fix — it breaks σ resolution on nine slugs (pandas renames the column
> `noiseParameters_meas`/`_sim` and the resolver stops finding it) and makes `Armistead`, `Fiedler`,
> `Weber`, `Blasi` and `Schwen` over-match. The real fix keys on identity while reading σ from the
> measurement side, and needs revalidating across all 23 slugs. That remains separate work.

## Gate B — the gradient is the one the objective implies

`tools/fd_check.py` at the PEtab nominal point on the current build. §2b recorded this slug FD-clean
at `6.8e−03`. The default step reads red, and **the documented step sweep is what settles it**:

| `h` | worst relative error |
|---|---|
| 1e−4 | 1.36e+00 |
| 3e−4 (default) | 1.00e+00 |
| 1e−3 | 2.66e−01 |
| **1e−2** | **1.40e−02** |

The error falls monotonically as `h` **grows**. That is the roundoff signature `tools/README`
describes — "a column whose FD drifts monotonically *away* from the assembled value as `h` shrinks
is converging at the other end" — the same shape as its `Fiedler` `k2` example. At `h = 1e−2` the
worst column is `1.40e−02`, the same order as the recorded `6.8e−03`. **Clean.**

Two artefacts checked rather than assumed:

- **The assembled gradient itself changes with `h`** (`KM`: 2.99 → 9.83 → 30.39 → 296.9). Expected:
  the tool clips the evaluation point to `[lo + 8h, hi − 8h]`, and with 6 of 21 nominal parameters
  at or essentially on a bound a different `h` genuinely means a different evaluation point. It
  reported moving `KM`, `init_Y_mRNA`, `tps_repr`.
- **`eff_GFP` and `init_GFP_mRNA` assemble to identical values at every `h`** (25.6243 at `h = 1e−2`).
  That is real model symmetry, not a defect — at `h = 1e−2` both match their own central differences
  to `2.4e−04`.

## Gate C — the fit reaches the benchmark optimum

| quantity | value |
|---|---|
| PyBNF reduced objective | −126.163175 |
| `J_paper = −log_likelihood` | −65.634945 |
| reference `J*` | −65.635120 |
| **`OG = J_paper − J*`** | **0.000175** |
| threshold | 1.92 |
| `OG` at the PEtab nominal point | 2.432380 |

`score.py` reads `SOLVED`. The fit reproduces `J*` to `1.75e−04`, i.e. it finds the same optimum the
benchmark's own best run found, from unbiased starts.

### What the 2.43 NLL improvement actually is

Simulating the fit and the nominal point and comparing to all 58 points on the log10 scale the
objective uses:

| | fit | PEtab nominal |
|---|---|---|
| RMS log10 residual | 0.068892 | 0.071842 |
| fitted σ | 0.068891 | 0.071841 |
| max abs log10 residual | 0.231020 | — |

The fitted σ equals the RMS residual to **six decimals** in both cases — exactly the MLE condition
for a Gaussian with an estimated scale, so the objective is internally consistent at both points.

And it accounts for the whole result. With σ estimated, the NLL carries `n log σ`, so the expected
gap is `58 × ln(0.071842 / 0.068892) = 2.4316` against an observed `2.4324 − 0.000175 = 2.4322`.
**The entire improvement is a 4 % reduction in the fitted noise scale.** The two trajectories are
near-indistinguishable by eye; both capture the repressilator's damped oscillations across
`t = 10 … 600`. This is a real, correctly-scored optimum on the benchmark's own objective — and it
is *not* evidence of a qualitatively better description of the biology. That distinction is why the
confidence above is 88 and not higher.

### Three of twenty-one parameters rest on a box bound

`background = 1e−05`, `tau_mRNA = 1e−05`, `tps_repr = 1e−05`, all on the lower bound of
`[1e-05, 1000]`. Every one of the 21 is **inside** the box — checked explicitly, because this
session produced a vivid counter-example (below). `tau_mRNA = 1e−05` is an mRNA lifetime five orders
below anything biological; it should be read as the fit exploiting a flat direction, not as an
estimate.

The fit also sits in a very different region from nominal — `eff_GFP` 224.9 vs `0.00195`, `scale`
`4.75e−05` vs `0.526` — with `scale` and `eff_GFP` trading off against each other, a degeneracy the
observable `GFP*scale + background` makes structural.

### The landscape, and why unbiased search needed ~1,000 starts

Three feasible local optima were characterised this session:

| `OG` | what it is |
|---|---|
| ≈ 0 | `J*`; found by 1 batch in 11 |
| 2.4324 | the PEtab nominal point — **proven** a box-constrained local optimum |
| 5.8147 | the dominant attractor of uniform multistart; 6 of 10 batches |

That the nominal point is itself a constrained optimum is a measurement, not an inference: **both**
`gntr` (trust region) and `lbfgs` (L-BFGS-B) seeded there converge in 30 iterations and stall at
`OG = 2.4304` / `2.4319` with four parameters pinned to bounds. The improving direction points *out*
of the box.

Distribution over 11 independent 100-start batches: `[0, 1.92)` × 1, `[5, 6)` × 6, `[10, 20)` × 2,
`[20, 40)` × 2; median `5.8155`. So the target basin is reached by roughly **1 start in ~1,000**,
which is why the two earlier attempts on this slug (a 20 h `cmaes` run at `OG 5.1051`, a
100 × 1000 `gntr` run at `5.8200`) missed it. The 5.82 they found is the same attractor this run
finds in 6 of 10 batches.

> **A seeded shortcut was tried and is void — do not resurrect it.** A `gntr` polish seeded at the
> nominal point via `logvar` reached `OG = −0.0140`, apparently *beating* `J*`. Three parameters had
> left the box: `init_Y_mRNA = 110709` against an upper bound of 1000, `tps_repr = 4.11e−09` and
> `background = 4.87e−06` against a lower bound of `1e-05`. This is exactly the trap §7 of the
> kickoff warns about — `var`/`logvar` point specs carry no bounds. The box **projection** of that
> point scores `OG = 274.33`, so the result depended entirely on infeasibility. Nothing seeded is
> shipped here; the row's `OG` is from the unbiased fit.

## Configuration

The shipped conf is the importer's `cmaes` recipe and is unchanged. The fit that produced this row
used `gntr` over the same PEtab box:

```
job_type = gntr
population_size = 100        # gradient optimizers use this as the START COUNT
max_iterations = 500
parallel_count = 3
random_seed = 4220003
loguniform_var = <all 21>  1e-05 1000
```

run as eleven independent seeded batches across three concurrent slots, best kept. The solve came
in slot 2's third batch, ~34 minutes in; the run was carried to 54 minutes to characterise the
distribution above.

**Why batches, measured — and it is not what a first pass suggested.** An earlier draft of this file
claimed throughput collapses with `population_size` (119 start-iterations/s at 50 against 9.4/s at
1500). **That was wrong and is retracted.** It measured the first 25 minutes of a 1500-start run,
which is dominated by iteration 0 — where every pathological random start burns the full
`wall_time_sim` before being retired — not by steady-state throughput. Measured properly, over a
300 s window after a 300 s warmup, counting completed start-iterations:

| `population_size` | start-iterations/s |
|---|---:|
| 100 | 28.6 |
| 1500 | **34.9** |

Large populations are, if anything, slightly *faster* per unit work. Two things actually made the
batched run win, and neither is a PyBNF defect:

* **Concurrency.** One `pybnf` process plateaus at ~310% CPU on this problem regardless of
  `parallel_count` (10 and 20 both measured ~310–330%, with 20 marginally slower). Three concurrent
  processes reached ~660–680%, roughly tripling throughput on a 10-core box.
* **Fewer iterations per start.** 500 rather than 1000, which costs nothing here — the historical
  100 × 1000 run reached `5.82` and these batches at 100 × 500 reach `5.81` — so the same wall time
  buys twice the independent draws, which is the only currency that matters on a basin-selection
  problem.

Batching also yields a scored checkpoint every ~12 minutes instead of one result at the end, which
is what made the distribution above measurable at all.

## Provenance

| component | version |
|---|---|
| bngsim | **`ffbf015`**, local **editable** install (not the PyPI wheel) |
| PyBNF | `095a5a14` |
| numpy | 2.5.2 |
| upstream PEtab | `4d20850` (the commit `upstream.json` pins) |

Capability-probed, not version-probed: `_species_value_factors` (#553) present, `_psvs_row_divisor`
(bngsim #161) present, `bngsim.AUTO` present, `BNGSIM_HAS_EVENT_SENS` true, extension
`__build_commit__ = ffbf0156eea4` matching the checkout. The `dff901e → ffbf015` bngsim update this
session carried **nothing in the ODE or sensitivity path** — vendored NFsim, CI, and build tooling
only — and both the Bruno canary (`J_paper = -46.68819468635073`) and this slug's nominal check
(`OG_nominal = 2.432380054142861`) reproduced **bit-identically** across it from cleared codegen
caches.

> `nominal_check.json`'s recorded values predate both, having been produced under bngsim **0.11.35**,
> which is below PyBNF's own `>=0.12.2` pin. Issue #38 flags "whether the collection's nominal checks
> get regenerated under 0.12.2" as an open decision; this row does **not** resolve it. The current
> build reads `OG_nominal = 2.432380054142861` against the recorded `2.4323694626309234`, a
> `1.1e−05` shift that is immaterial against a 1.92 threshold.

### ADR-0105 / lanl/PyBNF#549 is what made this run possible

The previous `gntr` attempt could not evaluate a point near the `cmaes` optimum at all —
`CVODE ... flag=-4` at `t=10`, the sensitivity-augmented system being far stiffer than the plain
trajectory. Per-species `atol` landed after that was measured, and is active by default here
(`sbml_atol` unset). **Zero `flag=-4` occurred in this session** across the pre-flight (10 starts ×
50 iters), the calibration (50 × 2000) and 1,100 production starts. Remaining failures are all at
`t = 0` — pathological random starts — with zero `wall_time_sim` timeouts.

## Bottom line

**SOLVED, unbiased, `OG = 0.000175`.** The fit reproduces the benchmark's `J*` to `1.75e−04` from
1,100 from-scratch starts over the full PEtab box in 54 minutes, and beats its own nominal point.
Gate A is confirmed by an independent oracle that reproduces PyBNF to `9.5e−07` relative — the first
time that oracle has been run for this slug, which also earns it `obj ✓`. Gate B is clean once the
documented step sweep is applied. The honest qualification is Gate C's: the 2.43 NLL margin over the
nominal point is a 4 % reduction in fitted noise, not a visibly better trajectory, and three
parameters rest on a bound.
