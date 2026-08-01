# VALIDATION — Blasi_CellSystems2016

Validation of the PyBNF benchmark job against the **Grein et al. 2026** reference objective. The oracle
here is the benchmark's reference **J\*** (the best Eq. 6 NLL over all optimizer runs on Marvin): the
question is whether **PyBNF's optimizer reaches the benchmark optimum on PyBNF's own faithful
reproduction of the paper's objective**.

> **Confidence: 97 / 100.** Objective fidelity is established three independent ways (below); the fit is
> SOLVED with OG = −4.3×10⁻⁷ from a from-scratch multi-start, and it recovers the paper's published
> maximum-likelihood rate constants to four digits without ever being shown them. The deduction is that
> the SBML model is imported rather than re-derived here — though this repository's
> `models/combinatorial_histone_h4_acetylation_blasi2016/` does re-derive the same network from the
> paper independently, and its BNGL steady state agrees with this SBML model to 2×10⁻⁸.

## What makes this problem different

Every one of the 252 measurements is at `time = inf`. There is no time axis and no dose axis — the data
are the **stationary** abundances of the 16 acetylation motifs under one condition, `control`, and they
sum to 1 because the model conserves total H4. That is why Blasi was the last subset-I problem to land:
lanl/PyBNF#509 (`lnnormal` natural-log noise, ADR-0084) was needed for it to *import*, and then config
load still crashed in `TimeCourse` (`OverflowError: cannot convert float infinity to integer`) until
lanl/PyBNF#521 (ADR-0086) taught the importer to emit a steady-state measurement instead of
materializing `t = inf` as a time course.

## Gate A — objective fidelity (three independent confirmations)

The Blasi objective is the Gaussian NLL on **natural-log-transformed** observables with an **estimated**
shared `sigma`. PyBNF minimizes the reduced form and reports the full normalized log-likelihood
separately, so the paper-scale identity `J_paper = −log_likelihood` has to carry both the `½log(2π)`
constants and the log-transform Jacobian `Σ log y_obs`. Over n = 252 points those together are
`−lnL − J_pybnf = −870.430293` (= `252·½log(2π) = 231.5725` minus a Jacobian of `1102.00`). Confirmed:

1. **Independent NumPy re-evaluation.** The whole objective was rebuilt outside PyBNF — steady state
   solved as a null-space problem on the 16×16 generator of the 4-cube, then
   `Σ[½log(2π) + log σ + log y_obs + ½((log y_obs − log y_sim)/σ)²]` — and evaluated at the same
   nominal point. It agrees with PyBNF on **both** scales to **8×10⁻¹³**: reduced
   `227.60339632354192` vs `227.6033963235412`, and `J_paper` `−642.8268965824425` vs
   `−642.8268965824433`. Recorded in `nominal_check.json`.
2. **Backend agreement.** The same nominal point evaluates to `227.6033963235412` on
   `sbml_backend = bngsim` and `227.6033963235455` on `sbml_backend = roadrunner` — a difference of
   4.3×10⁻¹². The steady-state relaxation is not a backend artifact. (The roadrunner check has to be
   run under a metaheuristic `job_type`; `gntr` refuses any backend without forward sensitivities.)
3. **The reference J\* is the paper's own published fit.** This is the strongest check, and it needed
   the observation that the PEtab `nominalValue` point is *not* the optimum. Its eight acetylation rate
   constants are exactly the published MLEs of Table S1/S2 in Blasi et al. (2016), but `nominalValue`
   for the noise scale is a placeholder `sigma = 0.1`, where the published fit profiles `sigma` to its
   maximum-likelihood value `0.2532`. Holding the eight published rate constants and profiling `sigma`
   gives `J_paper = −1090.5610757` against `J* = −1090.5618246`, i.e. **OG = 7.5×10⁻⁴** — inside the
   solved threshold, with the residual explained by the four-digit rounding of the published constants.
   So PyBNF's objective, the benchmark's `J*`, and the paper's published maximum-likelihood fit are all
   the same number.

**Verdict: PASS.** `OG = −log_likelihood − J*` is the exact, self-consistent scoring; `score.py`
reproduces it.

## Gate B — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (20 starts × 500 iterations, box-center + Latin-hypercube seeded by
`random_seed = 1`, `sbml_backend = bngsim`) converges in about 25 s to `J_pybnf = −220.1315318` ⇒
`J_paper = −1090.561825` ⇒ **OG = −4.3×10⁻⁷ < 1.92**. The gap is negative at the level of the printed
precision of `log_likelihood`: PyBNF lands on the benchmark optimum, marginally past the best of the
paper's 380 Marvin runs. The landscape is benign — 365 of those 380 runs also reached it, and the 15
that did not are every Nelder-Mead run but one and six of the ten pyswarm runs — so the optimum is
unambiguous.

**Verdict: PASS (SOLVED).**

## Gate C — parameter recovery

The fit was started from the box with no knowledge of the published values, and it recovers all eight
published rate constants to their full four published digits (`best_fit_params.txt` vs Table S2 rank 1;
the table is in `README.md`). The recovered `sigma = 0.253210` equals the profile maximum-likelihood
value `sqrt(Σresidual²/n)` of the published fit. This is a stronger statement than the objective match:
the optimum is *identified*, not merely *attained*.

**Verdict: PASS.**

## Configuration

- Import: `petab1to2_preserve_scale` (parameter scales preserved, lanl/PyBNF#491) → `import_job`, from
  upstream `Benchmark-Models-PEtab` commit `4d20850`. Emitted
  `noise_model <obs> = lnnormal, sigma = fit sigma` per observable — correct for
  `observableTransformation = log` with `noiseDistribution = normal`.
- `edition = 2`, `sbml_backend = bngsim` (required for the gradient/sensitivity path),
  `job_type = gntr`, `population_size = 20`, `max_iterations = 500`, `random_seed = 1`.
- One hand correction: the importer's dead `noise_model x_k5k16 = lnnormal, sigma = read_exp_file _SD`
  line is deleted. `observable_k5k16` is declared in the PEtab observables table but has **zero**
  measurement rows — the K5K16 motif is below the LC-MS quantification limit — so no `.exp` file has an
  `x_k5k16` column and the line is never consulted. Dropping it keeps the conf from implying a 16th
  fitted observable. Harmless either way; worth filing upstream as an importer nit.

## Bottom line

The last of the 23 subset-I problems, and the cleanest recovery in the collection: PyBNF imports a
purely steady-state PEtab problem, fits it in 25 s, reaches the benchmark optimum
(OG = −4.3×10⁻⁷), and lands on the published maximum-likelihood parameters of Blasi et al. (2016) to
four digits. It also closes the loop on what `J*` means for this problem — it is the paper's own fit,
which the nominal-point check could not have shown, because that point pins `sigma` at a placeholder.
