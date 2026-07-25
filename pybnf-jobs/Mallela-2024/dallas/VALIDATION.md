# VALIDATION — Mallela-2024/dallas

Primary-source validation of the PyBNF job `pybnf-jobs/Mallela-2024/dallas/`.
Confidence is **earned from the gate evidence below**.

> **Confidence: 85 / 100** — every input is the authors' own (model, data, and the
> published MAP from `examples/Vax_and_Variants/Dallas/`), the model file differs from the
> independently verified library model only by the removed actions block, and the gradient fit
> **improves on the published MAP**
> (neg_bin NLL **4875.96 → 4762.00**, Δ = -113.96). The deduction is that this is a
> *2-parameter local refinement*, not a reproduction of the authors' full 22-parameter
> posterior: the remaining parameters — every switch time among them — are pinned, and one of the
> two blockers below had to be fixed upstream in PyBNF for a gradient fit of this model class to
> run at all.

Primary source (in the PyBNF repo; the authors' own PyBioNetFit setup):
- Paper: Mallela A, Chen Y, Lin YT, Miller EF, Neumann J, He Z, Nelson KE, Posner RG, Hlavacek WS.
  "Impacts of Vaccination and Severe Acute Respiratory Syndrome Coronavirus 2 Variants Alpha and
  Delta on Coronavirus Disease 2019 Transmission Dynamics in Four Metropolitan Areas of the United
  States." *Bull Math Biol* 2024; 86:31. DOI 10.1007/s11538-024-01258-4.
- Author files used: `~/Code/PyBNF/examples/Vax_and_Variants/Dallas/` — `Dallas.bngl`,
  `Dallas.conf` (the legacy `am` / `neg_bin_dynamic` adaptive-MCMC sampler), `Dallas.exp` (daily
  NYT case data, t = 0..648), and `Output/adaptive_files/MLE_params.txt` (the published MAP,
  == Table 1).

"The paper's result" for this job = **the fit of the Dallas MSA daily detected-case counts**
(panel A of Fig. 2) at the published MAP.

---

## Gate 0 — Materials inventory

| needed | present? | path / note |
|---|---|---|
| paper PDF | ✅ | `dev/papers/Mallela2024/BullMathBiol24.pdf` |
| author model/job/data files | ✅ | `Vax_and_Variants/Dallas/` — `.bngl` + `.conf` + `.exp` |
| authors' fit result (MAP) | ✅ | `Dallas/Output/adaptive_files/MLE_params.txt` (22 params + r) |

**Verdict:** PASS — a complete, self-consistent authors' PyBioNetFit setup.

## Gate 1 — Data provenance

| `.exp` | source | method | units | diff vs. author file | verdict |
|---|---|---|---|---|---|
| `dallas.exp` | authors' `Dallas.exp` | copied verbatim | daily new detected cases (integer counts) | byte-identical | PASS |

t = 0..648, day 0 = 2020-01-21, day 648 = 2021-10-30. The series carries a few **negative** days
(a state revising its cumulative count downward); PyBNF's count-domain guard scores those as 0, and
`make_reproduction.py` reproduces that convention exactly.

**Verdict:** PASS — byte-level provenance to the authors' own data file.

## Gate 2 — Model fidelity

`dallas.bngl` is the curated library model
`models/covid19_vaccination_and_variants_mallela2024_dallas.bngl` with only the
*simulation* actions removed. That library model is itself independently verified (its notebook
re-implements Appendix Eqs. 1–40 in SciPy and matches BioNetGen to ≤2.3e-5 on every observable).

| aspect | authors' model | our `dallas.bngl` | verdict |
|---|---|---|---|
| social distancing | 5 periods (n = 4); Ptau/Lambdatau step at σ, τ₁..τ<sub>4</sub> | identical | match |
| variants | Ytheta ×y₁ at θ₁ (Alpha), ×y₂ at θ₂ (Delta); Utheta1/Utheta2 gate S<sub>V</sub> reinfection | identical | match |
| vaccination | empirical μ(t) table → V₁..V₆ → S<sub>V,1..4</sub> | identical | match |
| `counter()` seed / rate | `counter()=0`, `0->counter() 1` (t = sim time) | identical | match |
| S₀ | 7,573,136 | identical | match |
| free parameters | 22 `X__FREE` + `r__FREE` | ids declared; conf frees only `beta`, `fD`; r fixed | scoped (deliberate) |
| fit observable | `fDCs_Cum` | identical | match |
| network cap | none needed (finite: 42 species / 88 reactions) | none | match |
| actions block | generate + simulate | removed (synthesized from conf) | expected (edition-2) |

**Independent structural check:** the generated `reference/dallas.net` is **identical** (modulo
comment lines) to the library model's `models/.../reference/…_dallas.net`, and the
SBML conversion passes the bngsim L0–L3 faithfulness ladder (below).

**Verdict:** PASS.

## Gate 3 — Reproduce, and improve on, the paper's parameters (the headline)

`make_reproduction.py` simulates through **BNG2.pl** (not bngsim — so the figure reproduces without
the gradient toolchain) and scores the model exactly as PyBNF does: a negative-binomial −logpmf per
day with the prediction taken as the **mean**, `prob` clipped as PyBNF clips it, negative counts
skipped, on the cumulative→incident differenced series.

| metric | model @ published MAP | model @ L-BFGS-B fit |
|---|---|---|
| **neg_bin NLL (the fit objective)** | **4875.96** | **4762.00**  (Δ -113.96, −2.3 %) |
| median \|rel err\| vs. raw daily counts | 33.2 % | 34.4 % |
| median \|rel err\| vs. the 7-day rolling mean | 28.6 % | 26.6 % |
| peak of the 7-day mean, vs. data | +23.6 % | -1.8 % |
| cumulative reported cases, vs. data | +18.0 % | +4.0 % |
| `beta` | 0.304173 | 0.310366 |
| `fD` | 0.481101 | 0.400519 |

Why the two relative-error rows: the model predicts a **smooth mean**, while the NYT daily counts
carry a strong day-of-week reporting cycle, revision spikes, and the occasional negative day — so
the 7-day rolling mean is the fair target for a *shape* comparison, while the NLL (scored on the
raw counts, as the fit itself is) is the metric the job actually optimizes. Both are reported so
neither flatters the result.

**Reading the numbers.** FD drops 17 % and beta rises 2 %, pulling the model's total case count from +18 % to +4 % of the data and its peak 7-day mean from +24 % to within 2 %.

**Verdict:** PASS — the model tracks every surge in sequence across ~4–5 orders of magnitude, and
the 2-parameter gradient refinement lowers the paper's own objective. Because the published MAP is
an MCMC-sampled *mode* of a 22-parameter posterior rather than an optimum, an improvement
here is expected; it is evidence the gradient path works, not a correction to the paper.

## Native-only verification (not PEtab)

`neg_bin` + `location = mean` + `cumulative` are outside the PEtab v2 subset → native-only:

1. **Tier-1** (`scripts/check_conf.py`): edition 2, `lbfgs` resolves, data bound, 2 free params by
   id, no `__FREE`. **PASS.**
2. **Native-only guard:** `pybnf.petab.export_job` raises `NotImplementedError` — *"the whole-fit
   noise model is mean-centered (location = mean); PEtab v2 takes the prediction as the
   distribution median for every noise family"*. **PASS.** (No PEtab round-trip is run, by design.)
3. **Real bngsim gradient fit** (`pybnf -c dallas.conf`): converges to a finite objective
   **4762.00**; a few L-BFGS-B line-search trials at the box edges still return `inf` (CVODES
   stalls when a trial β puts a surge on top of a switch time) and are absorbed by
   `max_failed_simulations`. **PASS.**
4. **Reproduction:** Gate 3.

## Gate 4 — SBML / SED-ML / COMBINE exports

`reference/` carries `dallas.net`, `dallas.xml` (SBML L3V2), `dallas.sedml` (SED-ML L1V3) and
`dallas.omex`, produced by `bngsim.convert.net_to_omex(..., gate="full")`.

| level | check | verdict |
|---|---|---|
| L0 | syntactic validity (libsbml consistency) | pass |
| L1 | structural equivalence (species/reaction counts + topology) | pass |
| L2 | round-trip identity (max scale-relative \|Δ dy/dt\| = 3.5e-18) | pass |
| L3 | numerical equivalence — all **27,258** species·time cells, max rel \|Δ\| = 0 | pass |
| L4 | symbolic equivalence | inconclusive (piecewise/table kinetics — the checker punts) |

The SED-ML carries **this job's** protocol (a uniform time course over t = 0..648, 649 points,
CVODE, atol/rtol 1e-7) supplied explicitly, because the job-form `.bngl` has no actions block to
transcribe.

> **Conversion gotcha.** The 67-branch empirical vaccination table `v_rate()` becomes an
> equally deep MathML `<piecewise>`. Raise `sys.setrecursionlimit` and run the conversion on a
> thread with an enlarged `threading.stack_size` (100k / 512 MB were used for all four slugs), or
> the SBML reload overflows the stack — it bites hardest on New York City (353 branches, a 3.2 MB
> SBML) and Phoenix (341).

## Toolchain requirement

A gradient fit of **this class of model** — `if()`-gated rate laws plus a count likelihood — did
not run on stock PyBNF. Two independent defects, both fixed upstream in PyBNF
(`fix: gradient-fit a piecewise model with a count likelihood`), were required:

1. **Over-requested output-sensitivity selectors.** `bngsim_model/net_model.py` asked bngsim for an
   `expression:` sensitivity for **every** global function whenever `print_functions` was on.
   bngsim differentiates function bodies symbolically and refuses any body carrying an `if()`
   (lanl/bngsim#198), raising `ValueError` if such a selector is requested. **All 14** of this
   model's functions are `if()` chains, so *every* simulation died and the objective was `inf`
   everywhere. The request is now filtered to the differentiable expressions; the scored column
   (`fDCs_Cum`) is an observable and is unaffected.
2. **Divide-by-zero in the negative-binomial derivative.** `d(data_fit)/d(mean) =
   r (mean − obs)/(mean (r + mean))` is `0/0 → nan` at a mean-centered prediction of exactly 0 —
   which is what this model predicts for the first ~15 days, before `t0` gates transmission on.
   The `nan` reached the optimizer as a `nan` parameter proposal. The mean is now floored at the
   image of the clip `data_fit` already applies to `prob`, so the slope is finite and consistent
   with the value being scored.

Both are ordinary bugs with unit tests, not workarounds; the full PyBNF suite (3304 tests) passes
with them.

## Divergence & corrections

- Scope vs. paper: **narrower by design** — a 2-parameter gradient refinement of the published MAP,
  not a reproduction of the authors' 22-parameter adaptive-MCMC posterior. Every other
  parameter is pinned at the published value.
- No switch time is free: their gradient is identically zero without bngsim PR #50 (issue #48).
- The dispersion is **fixed** at the published r = 3.021658179294529 rather than fit as a nuisance, so the
  likelihood surface is exactly the one the paper's r defines.
- `uniform_var` bounds are narrowed to the MAP ±50 % (not the authors' verbatim 0–20 / 0–1),
  because a gradient job starts at the **box center**; the authors' bounds start at β = 10, where
  CVODES stalls and the objective is `inf`. This is the single most important conf choice here.
- No change to the science: S₀, every fixed constant, the data, and the model structure are the
  authors' own.

## Bottom line

A small, fast, fully verified gradient job that does something the paper's own sampler did not: it
finds a strictly better point of the paper's own objective in two parameters, in a few minutes,
from the paper's own starting estimate. It also pins down exactly why gradient fitting a piecewise
epidemiological model was previously impossible in PyBNF. Most valuable next step: with bngsim
PR #50's switch-time sensitivities wired through PyBNF, free `t0`/`t_sigma`/`t_tau*`/`theta*` and
re-run — that is the fit this model class really wants.
