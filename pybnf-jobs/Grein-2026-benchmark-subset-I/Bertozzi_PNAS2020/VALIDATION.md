# VALIDATION — Bertozzi_PNAS2020

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's **derived-parameter**
case: its only infection rate is fixed by an SBML `initialAssignment` on a *parameter*,
`beta_N = R0_*gamma_/N_`, from three quantities the conditions set.

> **Confidence: 90 / 100.** SOLVED with OG = 5.4×10⁻⁶ from a from-scratch 20-start multi-start; the
> estimated-σ objective identity is verified to 14 significant figures, and the fit independently
> reproduces the PEtab `nominalValue` point, which is itself the published optimum
> (`OG_nominal = 5.09×10⁻⁶`). Deductions: the model is imported (not re-derived from Bertozzi et al.
> 2020); the run is a single seed; and two parameters land ~1.5% off nominal along a shallow
> direction (below).

## Gate A — objective fidelity

Bertozzi is an **estimated scalar σ** problem — `sd_I_CA` and `sd_I_NY` are free parameters, two of
the eight — so the restored constant is the bare `(N/2)log(2π)`:

    −lnL  =  J_pybnf_reduced  +  (N/2)·log(2π)

| term | value |
|---|---|
| `J_pybnf_reduced` at the nominal point | 138.64762007128826 |
| `−lnL` = `J_paper` | 158.86426780179107 |
| observed restored constant | **20.216647730502814** |
| predicted `(N/2)log(2π)`, N = 22 | **20.2166477305028** |
| difference | 1.4×10⁻¹⁴ |

**Verdict: PASS.**

## Gate B — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (20 starts × 500 iterations, `random_seed = 1`,
`sbml_backend = bngsim`) converges to `J_pybnf = 138.6476200551601` ⇒ `J_paper = 158.8642678` ⇒
**OG = 5.4×10⁻⁶ < 1.92**, against `J* = 158.86426270904192`. Wall time about one minute on 6 cores.

The PEtab `nominalValue` point evaluates to `OG_nominal = 5.09×10⁻⁶`, so for this problem the nominal
point **is** the published optimum and the fit found it from scratch. The two agree on six of the
eight parameters to between 5×10⁻⁶ and 6×10⁻⁵:

| parameter | fitted | PEtab nominal | rel. diff |
|---|---|---|---|
| `I0_CA` | 232.91036 | 232.9026 | 3.3×10⁻⁵ |
| `I0_NY` | 3631.9544 | 3631.9735 | 5.3×10⁻⁶ |
| `R0_CA` | 4.9992874 | 4.9272168 | **1.5×10⁻²** |
| `R0_NY` | 1.5 | 1.5 | 0 |
| `gamma_CA` | 0.038328067 | 0.039032488 | **1.8×10⁻²** |
| `gamma_NY` | 0.15974815 | 0.15974647 | 1.1×10⁻⁵ |
| `sd_I_CA` | 171.54239 | 171.552 | 5.6×10⁻⁵ |
| `sd_I_NY` | 638.74077 | 638.73702 | 5.9×10⁻⁶ |

`R0_CA` and `gamma_CA` are the exception, and the reason is a property of the problem rather than of
the fit — the two are **not independently identified**, only their product is. Both series cover the
early exponential phase, where `S ≈ N` and

    dI/dt ≈ (beta_N·N − gamma_)·I = gamma_·(R0_ − 1)·I

so `R0_` and `gamma_` trade off along that combination. It is conserved to 2.7×10⁻⁵ even where the
individual parameters differ by 1.5%:

| | fitted | PEtab nominal | rel. diff |
|---|---|---|---|
| `gamma_CA·(R0_CA − 1)` | 0.1532849554 | 0.1532890426 | 2.7×10⁻⁵ |
| `gamma_NY·(R0_NY − 1)` | 0.0798740750 | 0.0798732350 | 1.1×10⁻⁵ |

The fit simply slid along that valley until `R0_CA` reached its upper bound (5) while the likelihood
barely moved; `R0_NY = 1.5` sits at its **lower** bound in both. The objective difference is the whole
story: `5.4×10⁻⁶` fitted against `5.09×10⁻⁶` nominal — the fit is 3×10⁻⁷ *worse*, i.e. the two points
are indistinguishable.

**Verdict: PASS (SOLVED).**

## History — this slug was wrong twice, in two different ways

Until 2026-08-02 this job shipped `job_type = cmaes` with `OG_nominal = 1.79×10¹¹` and a `⚪ setup
only` status. Both numbers were artifacts, and neither had anything to do with CMA-ES.

**lanl/PyBNF#531 — the forward model.** PyBNF's fast simulation path (the cached-clone + `set_param`
route) recomputed a species initial fixed by an `initialAssignment`, but never a *parameter* fixed
that way. `beta_N` therefore kept the value it had at load, computed from the SBML file's own
`R0_ = 0.1, gamma_ = 0.1, N_ = 1` — `0.01` instead of `7.6×10⁻⁹`:

    R0_ 3.0 → 3.3, max trajectory delta : 0          # R0_ is completely inert
    engine beta_N (fast path)           : 0.01
    engine beta_N (reload)              : 8.34e-09

Two of the eight free parameters had an identically flat objective and the trajectory was off by six
orders of magnitude. This is a **scalar-path** defect: it corrupted the objective for every
`job_type`, so the `cmaes` recipe shipped here was optimizing a broken model. The reload path
(`_needs_structural_reload`) was always correct, and is the `rtol=0, atol=0` parity oracle the fix is
tested against.

**lanl/PyBNF#530 — the gradient column.** `route_experiment` composed the chain rule for a free
parameter reaching the model only through a `condition:` parameter reference, but only where
`d(IC)/d(target)` was a plain `1` (#511); anything else was an honest refusal. Bertozzi needed the
whole generalization at once:

- `I0_` seeds **two** species with opposite derivatives — `I_ = I0_` (+1, written as a legal
  one-argument `<times/>`) and `S_ = N_ - I0_` (−1);
- `R0_` and `gamma_` seed the **parameter** `beta_N`, with derivatives `gamma_/N_` and `R0_/N_` that
  are not numbers at all until the fit vector is known — and `gamma_` is *also* a rate constant the
  ODE reads directly, so its column is the sum of two terms.

The second of those was never covered by the refusal, because it is not a species initial condition:
those two columns were **silently wrong** rather than refused, which is why #530 could not be landed
without #531. Both now carry their real derivatives (ADR-0095), the point-dependent ones re-evaluated
at each evaluated PSet.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`; **no hand corrections**.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 20`,
  `max_iterations = 500`, `wall_time_sim = 10`, `random_seed = 1`.
- k = 8 free parameters (6 model + 2 noise scales), n = 22 scored points, 2 observables
  (`y_I_CA`, `y_I_NY`, both reading the single model column `I_` — the ADR-0077 shape), two
  conditions, two experiments.

## Provenance

Run against **bngsim 0.12.0** (`~/Code/bngsim` at `f4b24ac`) with PyBNF at ADR-0094 + ADR-0095. The
`nominal_check.json` here was recomputed after both fixes; the value it carried before 2026-08-02 is
not comparable with anything.

## Bottom line

The derived-parameter case, and the collection's sharpest reminder that a shipped recipe can be a
symptom rather than a choice. Both the `cmaes` recipe and the `1.79×10¹¹` nominal gap were downstream
of a simulation bug, not of the problem's difficulty: with the forward model correct and the
condition-routed gradient columns complete, it solves on the gradient path in about a minute, from
scratch, to six significant figures.
