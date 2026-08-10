# VALIDATION — Giordano_Nature2020

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This is the collection's **largest** problem on both
axes — k = 50, n = 313 — the SIDARTHE COVID model of

> Giordano G, Blanchini F, Bruno R, Colaneri P, Di Filippo A, Di Matteo A, Colaneri M.
> **"Modelling the COVID-19 epidemic and implementation of population-wide interventions in Italy."**
> *Nature Medicine* **26**, 855–860 (2020). <https://doi.org/10.1038/s41591-020-0883-7>

It is also the slug that produced **lanl/PyBNF#546 / ADR-0103**, and the one whose *apparent*
difficulty was most badly overstated by a placeholder σ.

> **Confidence: 90 / 100.** SOLVED with `OG = 0.135` from a from-scratch 100-start multi-start, no
> seeding, ~2 h 26 m. The result is corroborated in an unusually strong way: a σ-profiling analysis
> performed **before the fit ran** predicted it would need to shrink the residual norm ~10.8× below
> the collection's own reference trajectory, and the completed fit achieved **10.6×**. Deductions:
> the model is imported, the run is a single seed, and 17 of its starts died on CVODES `mxstep`
> (lanl/PyBNF#549) so the effective budget was ~83 starts rather than 100.

## Gate A — objective fidelity

Linear observables with seven estimated scalar σ, so the restored constant is `(N/2)log(2π)` with no
Jacobian:

| term | value |
|---|---|
| PyBNF reduced objective | −3775.834252 |
| restored constant (`−lnL − J_reduced`) | 287.627761 |
| `(N/2)·log(2π)`, N = 313 | 287.627761 |
| `J_paper = −log_likelihood` | −3488.206491 |

**Independently confirmed, and exactly.** §2c's oracle recomputes the Eq. 6 NLL at the nominal point
from upstream's `measurementData` joined to `simulatedData`, with no PyBNF in the loop. The joined
residuals give `Σr² = 9.28e−07`, hence a reduced objective of `4.6397e−07` — matching
`nominal_check.json`'s `4.639740019658988e-07` to every digit, on all 313 rows. That also establishes
something used throughout this file: **upstream's `simulatedData` *is* PyBNF's simulation at the
nominal point.** This slug carries `obj ✓`.

**Verdict: PASS.**

## Gate B — the gradient is the one the objective implies

This slug is where the corpus's worst gradient defect was found, and the history matters because the
first diagnosis was wrong.

Its assembled gradient disagreed with central differences on **41 of 50 parameters**, by up to 26%,
at every step size — and the disagreement partitioned exactly along whether a parameter sat behind a
time gate (41 gated, 9 ungated). The obvious reading is unhandled switching. It is not: bngsim
already registers every `time` inequality as a CVODE root, and this model gets 13 of them.

The cause was the **absolute tolerance**. CVODE weights each state by `rtol·|y| + atol`, so a
constant `atol` declares values beneath it to be noise. Giordano is a population-*fraction* model
whose species sit at `1.7e−8 … 1`, median `3.7e−7`; at the `1e−8` default the absolute term buried
the relative one across the whole early trajectory. The gate correlation is real but incidental — a
gated parameter acts only inside its own stage window, and the earliest windows are where the states
are smallest. Four decades of `rtol` bought nothing; `atol` fixed it.

| | worst relative error |
|---|---|
| before (bngsim default `atol = 1e-8`) | **7.7e−02** |
| after ADR-0103 (`atol` derived from the model's own scale) | **4.5e−04** |

**Verdict: PASS (after lanl/PyBNF#546, ADR-0103).**

## Gate C — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (**100 starts × 1000 iterations**, `random_seed = 1`,
`sbml_backend = bngsim`, no seeding) converges to `J_pybnf = −3775.834252` ⇒
`J_paper = −3488.206491` ⇒ **`OG = 0.135007` < 1.92**, against `J* = −3488.3414981`. Wall clock
**2 h 26 m**. Starts retired: 59 on `reached max_iterations`, 35 on `step is negligible`, 6 on
`start point failed to simulate`.

**The acceptance window here is the tightest in the collection, and it cleared it.** All 7 σ are
estimated over n = 313, so with σ profiled `∂OG/∂log(RMS) = 313` and clearing 1.92 requires matching
the reference residual norm to **0.61%**. The final `OG = 0.135` corresponds to matching it to
**0.043%**.

### The distance was 743.7, not 3776 — and the prediction held

`OG_nominal = 3776` is the number this problem was ranked by, and four-fifths of it was an artefact:
the PEtab nominal σ are placeholders at exactly 1, so `Σ nⱼ log σⱼ ≈ 0` instead of its MLE value.
Profiling the σ at the *nominal dynamics* gives `OG = 743.72` (`tools/sigma_profile.py`, self-check
5.7e−14). That is the honest starting distance, and the descent tracked it:

| stage | reduced | `OG` |
|---|---:|---:|
| PEtab nominal | 0.00 | 3775.97 |
| σ profiled, dynamics unmoved | −3032.25 | 743.72 |
| **final fit** | **−3775.83** | **0.135** |

**The stronger corroboration is a prediction made before the fit ran.** The σ-profiling analysis said
closing 743.7 required shrinking the residual norm by `e^(743.7/313) ≈ 10.8×` relative to upstream's
own reference trajectory. Measured at the final fit, per observable:

| observable | σ at reference | σ at our fit | ratio |
|---|---:|---:|---:|
| `observable_Deaths` | 1.038e−04 | 2.364e−06 | 43.9× |
| `observable_CurrentCases` | 6.370e−05 | 4.353e−06 | 14.6× |
| `observable_ICU` | 1.515e−05 | 1.051e−06 | 14.4× |
| `observable_Recovered` | 2.270e−05 | 2.214e−06 | 10.3× |
| `observable_TotalCases` | 3.114e−05 | 4.598e−06 | 6.8× |
| `observable_DiagHome` | 5.727e−05 | 8.715e−06 | 6.6× |
| `observable_Hospit` | 2.484e−05 | 6.842e−06 | 3.6× |
| **geometric mean** | | | **10.6×** |

Predicted 10.8×, achieved 10.6×. The arithmetic that said this problem was reachable, and by how
much, was right.

**Verdict: PASS (SOLVED).**

### The fit tracks the data, and there is no observable being traded away

Unlike `Schwen_PONE2015`, this fit improves **every** observable relative to the reference trajectory
— the smallest ratio above is 3.6× and the largest 43.9×, all in the same direction. There is no
group being sacrificed for another, which is what a genuine basin looks like rather than a
leverage trade. The seven fitted σ land at `1.05e−06 … 8.72e−06` against data values reaching
`~4e−03`, i.e. the model reproduces the Italy time series to roughly 0.1–0.2% relative.

### Three of fifty parameters rest on a box bound

`eta_22` and `rho_22` at their lower bound of `1e−03`, `kappa_38` at its upper bound of `10`. All
three are stage-gated rates; three of fifty is low for this collection and does not undermine the
optimum, but they should be read as bounded rather than as point estimates.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`; **no hand corrections**. All 50 estimated
  parameters are `parameterScale = log10` upstream and import correctly as `loguniform_var` — this
  slug is the *control* case for lanl/PyBNF#548, since its v1 parameter table has no prior column at
  all and therefore never lost its search scale.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 100`,
  `max_iterations = 1000`, `wall_time_sim = 10`, `random_seed = 1`.
- k = 50 (43 model + 7 noise), n = 313 (one scored experiment `pred1`, 46 timepoints, 7 observables).
  **41 of the 50 sit behind a time gate** — SIDARTHE's rates replicated per NPI stage, suffixed by
  the day the stage begins (`_0`, `_4`, `_12`, `_22`, `_28`, `_38`).

## Known cost: 17 starts lost to integration failure

The run logged **1,976 `mxstep` lines and 17 failed starts**, so the effective budget was ~83 starts,
not 100. This is the same family as `Weber` and `Brannmark`: a population-fraction model whose states
reach `1.7e−08` against others of order 1, served by a single scalar `atol`. ADR-0103's derivation is
what made this problem tractable at all, but it is still one number for a wide range.
**lanl/PyBNF#549** — pass bngsim's per-species `atol` (`CVodeSVtolerances`, bngsim#196) instead of a
scalar — should recover those 17 starts. Worth re-running under it to see whether `OG` improves from
0.135 toward `J*`.

## Provenance

Run against **bngsim 0.12.2** (released wheel, predating bngsim#196) with PyBNF at `e008d345`.

## Bottom line

The collection's largest problem, k = 50 against n = 313, solved at `OG = 0.135` on 100 unbiased
starts — through the tightest acceptance window here (0.61%, cleared at 0.043%). It is also the
clearest case of the corpus's own analysis paying off: `OG_nominal = 3776` had it ranked as the
hardest remaining slug, σ-profiling said the real distance was 743.7 and named the 10.8× residual
shrink required, and the fit delivered 10.6×.
