# Grein 2026 benchmark, subset I — PyBNF fitting jobs

PyBNF fitting jobs for the **23 subset-I problems** of the **Grein et al. 2026** optimizer benchmark,
each carrying that benchmark's reference objective. These are the first data points placing PyBNF on
that leaderboard.

> Grein T, Penas DR, Weindl D, Lakrisenko P, Banga JR, Hasenauer J.
> **"A benchmark collection for optimizer evaluation in systems biology"** (working title).
> *bioRxiv* 2026.07.11.737731 — 33 optimizers × 30 PEtab problems, >1.5M core-hours.
> Data: `Benchmarking-Initiative/Benchmark-Models-PEtab` (the problems, PEtab v1) +
> `ICB-DCM/optimizer-benchmark-2026-suppl-code-and-data` (reference objectives `best_fx_marvin.csv`).

Grein et al. cover 30 problems; this directory covers subset I, which is 23 of them.

Unlike the hand-built BNGL jobs elsewhere in `pybnf-jobs/`, these are **PEtab-imported SBML** jobs
(`edition = 2`, `sbml_backend = bngsim`): ODE models converted from the PEtab v1 collection through
PyBNF's importer.

## This directory is not a vendored copy

It was named `Benchmark-Models-PEtab` until 2026-07-31, after the upstream repository, which made it
read as a checkout of that repository that we had no business editing. It never was one. There is no
submodule, and **not one upstream PEtab problem file lives here** — no `*.yaml`, no
`measurementData_*.tsv`, `observables_*.tsv`, `experimentalCondition_*.tsv` or `parameters_*.tsv`. A
PEtab problem *is* those tables. What is here is a corpus of **PyBNF jobs derived from** them, plus a
verbatim copy of each problem's SBML model.

### What is upstream and what is ours

Counting every file inside the 23 slug directories:

| category | files | bytes | |
|---|---:|---:|---|
| SBML model, verbatim from upstream | 23 | 1,681,108 | **copied** — 7% of bytes, 6% of files |
| `.exp` and `_measparams.tsv` — PyBNF-format *translations* of the upstream measurement tables, emitted by the importer | 185 | 141,974 | derived |
| `jstar.txt` — one number, transcribed from the ICB-DCM suppl repository | 23 | 433 | derived |
| `.conf`, `score.py`, `nominal_check.json`, `README.md`, `VALIDATION.md`, `best_fit_params.txt`, `information_criteria.txt` | 152 | 22,128,239 | **written here** |
| total | 383 | 23,951,754 | |

**Not one file here is copied except the 23 SBML models**, which are now 7% of the bytes and 6% of the
files — so the directory is 94% locally authored by file count and 93% by byte count. **Editing this
directory is normal**: the locally authored files are the deliverable, and adding a job necessarily
edits the coverage matrix below. Nothing upstream is ever modified.

> Earlier revisions of this table read `350 files / 2,164,145 bytes`, with the SBML models called out
> as *"80% of bytes"* and the locally written files as `119 / 340,630`. That was true before the ✅ rows
> started shipping their fits: the 20 `best_fit_params.txt` files are 21.7 MB on their own — each is
> the run's full sorted 5001-row parameter table — which is what moved the SBML share from 80% of bytes
> to 7%. Counted over `git ls-files` only; the `output/` trees and `bnf_*.log` files a run leaves behind
> are gitignored and are not part of this.

### Upstream pin

`upstream.json` pins the exact upstream commit the 23 SBML files were taken from —
`4d2085084b289f6215a95475b1ee639fd7d42283` (2026-07-21) — and records, per file, a sha256 of its
LF-normalized content plus the raw sha256 on each side. Verifying is one uniform rule with no
per-file exceptions: re-fetch the pinned commit and assert, for every model,

    sha256(local LF-normalized) == sha256(upstream LF-normalized) == sha256_lf

All 23 pass as of 2026-07-31.

The comparison is on LF-normalized content because `.gitattributes` stores `*.xml` as LF here
whatever upstream ships. One upstream file, `Armistead_CellDeathDis2024`, is CRLF, so its
`sha256_raw_upstream` differs from `sha256_raw_local`; `upstream.json` flags that in
`upstream_line_endings` and `models_with_crlf_upstream`, and it is expected rather than drift. The
`*.xml text eol=lf` rule means a future re-copy of a CRLF upstream file is normalized on `git add`
instead of silently landing with foreign line endings.

`upstream.json` also records where `jstar.txt` comes from: `data/best_fx_marvin.csv` in
`ICB-DCM/optimizer-benchmark-2026-suppl-code-and-data` — **not** `suppl/data/best_fx_marvin.csv`,
which is the path the preprint's own text implies and which 404s.

## The scoring: optimality gap (OG)

    OG = −log_likelihood − J*        "solved" iff  OG < 1.92     (χ², α = 0.05, 1 dof)

where **J\*** = `min` over all optimizer runs on the *Marvin* cluster of the paper's Eq. 6 Gaussian
negative log-likelihood (`data/best_fx_marvin.csv`, column `fx_best`). Each slug ships
`score.py`, which computes exactly this.

## Objective fidelity

Subset I is **not one noise regime**. Earlier revisions of this file asserted that every problem
estimates its σ as a free `sigma`/`sd_*` parameter; that is wrong, and the correction matters because
the restored constant — the bridge between what PyBNF minimizes and the paper's Eq. 6 NLL — differs
by regime. Classified from each conf's `noise_model` lines and `.exp` columns:

| regime | n | slugs | restored constant |
|---|---:|---|---|
| **estimated scalar σ** (`sigma = fit p`) | 12 | Bertozzi, Blasi, Boehm, Borghans, Brannmark, Elowitz, Giordano, Laske, Okuonghae, Perelson, Schwen, Sneyd | `(N/2)log(2π)` |
| **fixed per-point σ** (data `_SD` columns) | 4 | Bruno, Crauste, Rahman, SalazarCavazos | `Σ log σᵢ + (N/2)log(2π)` |
| **estimated, prediction-dependent σ** (`sigma = prediction_formula …`) | 2 | Armistead, Raia | `(N/2)log(2π)` |
| **per-measurement σ** (`sigma = formula noiseParameterN_*`, bound via `_measparams.tsv` sidecars, ADR-0083) | 2 | Fiedler, Zhao | `(N/2)log(2π)` |
| **mixed estimated + `fix_at`** | 1 | Weber | `(N/2)log(2π)` + the fixed-σ terms |
| **unit σ** (`objective = sos`, no noise model) | 2 | Oliveira, Smith | `(N/2)log(2π)` |

The objective is the *full* Gaussian NLL (Eq. 6, with the `log(2πσ²)` normalizer) in every regime, not
a bare sum-of-squares. The fixed-per-point identity is verified numerically against the `.exp` `_SD`
columns in two slugs' `VALIDATION.md` — Bruno to 4.9×10⁻⁷ (where σᵢ < 1 makes `Σ log σᵢ` negative) and
Crauste to 1.5×10⁻⁷ (where σᵢ ≫ 1 makes it large and positive).

**The restored-constant column above is the σ-source part only.** A slug whose observables are
log-transformed carries the change-of-variables Jacobian on top of it, and for six of them that term
*dominates* — so the column is not the whole restored constant for those, and reading it as such
will mislead. Measured at each nominal point (`J_paper − reduced` against `(N/2)log(2π)`):

| slug | scale | restored | of which `(N/2)log(2π)` | of which Jacobian |
|---|---|---:|---:|---:|
| `Blasi_CellSystems2016` | ln | −870.4303 | 231.5725 | −1102.0028 |
| `Borghans_BiophysChem1997` | log10 | 114.7780 | 102.0022 | 12.7758 |
| `Elowitz_Nature2000` | log10 | 60.5282 | 53.2984 | 7.2298 |
| `Laske_PLOSComputBiol2019` | lin/ln | 299.6851 | 38.5954 | 261.0897 |
| `Perelson_Science1996` | log10 | 247.8317 | 14.7030 | 233.1287 |
| `Schwen_PONE2015` | log10 | 1255.8661 | 262.8164 | 993.0497 |

`Laske_PLOSComputBiol2019` is the one to read first, because it is the **mixed** case: 33 of its 42
points are natural-log (`lnnormal`) and 9 are linear Gaussian, so its Jacobian is `Σ log(y_obs)` over
the `lnnormal` points *only*. Computed straight from the `.exp` columns that is 261.089654, which with
`(N/2)log(2π) = 38.595418` gives 299.685073 — the observed restored constant to every digit. Its
`VALIDATION.md` records the check.

### Why `objective = sos` is faithful for Oliveira and Smith

These two ship a plain sum of squares and no `noise_model` line, which looks like a fidelity break. It
is not, and the reason is worth recording because it is not obvious from the conf alone. **Both
problems specify σ ≡ 1 upstream**, checked against the pinned commit:

- `Smith_BMCSystBiol2013` — `observables_*.tsv` carries the literal `noiseFormula = 1.0` for all nine
  observables.
- `Oliveira_NatCommun2021` — `noiseFormula = noiseParameter1_*`, bound by the measurement table to
  `sd_cumulative_cases` / `sd_cumulative_deaths`, which `parameters_*.tsv` gives `nominalValue = 1`
  and **`estimate = 0`** — fixed, not estimated.

With σ ≡ 1 uniformly, Eq. 6 is `½Σr² + (N/2)log(2π)`. Minimizing `sos` (`Σr²`) and minimizing that NLL
differ by a positive affine map, so they have the **same argmin** — the optimizer lands on the same
point either way. (Uniformity is what makes this safe: with per-observable σ ≠ 1 the two objectives
would weight residuals differently and the argmins would diverge.)

Scoring is unaffected regardless, because `score.py` never reads the `sos` value: `−lnL` in
`information_criteria.txt` is computed from the **noise model's** per-point `log_density` (ADR-0056),
independent of the fitting objective. Oliveira confirms the arithmetic exactly — `J_paper − reduced =
110.27262 = 120 × ½log(2π)`. For Smith the same constant (56.97) is invisible only because it falls
below the ulp of its nominal objective, 6.85×10³².

PyBNF **minimizes a reduced objective** that drops the parameter-independent per-point
constants — `½log(2π)`, and (for a log10 observable) the change-of-variables Jacobian
`Σ log(y_obs·ln10)` — because they do not affect the argmin. It then **reports the full normalized
log-likelihood** at the best fit in `information_criteria.txt` (matching `scipy.stats.norm.logpdf` /
`lognorm.logpdf`), restoring every dropped constant. Therefore

    J_paper  ==  −log_likelihood

exactly, for **both** linear and log10 observables.

**This identity is now corroborated across the collection, not just argued.** For twelve problems the
PEtab `nominalValue` point is the published optimum, and evaluating PyBNF's objective there
reproduces the paper's `J*` to within the solved threshold — seven to ~10⁻³ or better (four of those
to ~10⁻⁵), plus `Fiedler_BMCSystBiol2016` at ~2e−3 (see `nominal_check.json` in each slug). That is
an end-to-end check of the whole imported chain: SBML model → simulation → observable formulas →
noise model → objective.

**The log10 Jacobian convention is settled, and it matches.** An earlier revision of this section, and
the `nominal_check.json` of every log10 slug, hedged that `−log_likelihood` carries a change-of-variables
Jacobian the paper's Eq. 6 `J*` "need not", so a log10 nominal check was recorded but not asserted. That
hedge predated any log10 slug being *fitted*, and three solved slugs now settle it — each with a large
Jacobian, each landing on `J*`:

| slug | scale | Jacobian | OG **from a fit** |
|---|---|---:|---:|
| `Perelson_Science1996` | log10 | +233.1287 | 5.0e−07 |
| `Blasi_CellSystems2016` | ln | −1102.0028 | −4.3e−07 |
| `Laske_PLOSComputBiol2019` | lin/ln | +261.0897 | −1.0e−06 |

Were the term absent from Eq. 6, Perelson would miss `J*` by 233 and Blasi by 1102. A *natural-log*
Jacobian in place of log10 would offset Perelson by `n·log(ln10) = 13.34`, which the same result
excludes. Both sides carry the identical term, so a log10 nominal OG is as much a validation as a
linear one.

> **Two nominal checks were recomputed on 2026-08-02** and their old values should not be used.
> `Bertozzi_PNAS2020` moved from `OG = 1.79e+11` to `5.09e−06` — its nominal point *is* the
> published optimum after all — and `Laske_PLOSComputBiol2019` from `96.7` to `39.9`. Both had been
> evaluated against a forward model that **lanl/PyBNF#531** has since fixed: a parameter fixed by an
> SBML `initialAssignment` was never recomputed when its dependencies changed, so Bertozzi's
> `beta_N = R0_*gamma_/N_` and Laske's 27 COPASI-style `ModelValue_*` aliases all kept their
> load-time values. That silently pinned any fitted or condition-set quantity reaching the dynamics
> only through them. It is a **scalar-path** defect, so it affected every `job_type`, not just the
> gradient ones. No other subset-I model derives a parameter this way — checked by scanning all 23
> SBML files for an `initialAssignment` whose symbol is a non-rule-governed parameter.

## Coverage: all 23 subset-I problems

| slug | run cost | J\* | scale | k | n | optimizer | OG | obj ✓ | status |
|---|---|---|---|---|---|---|---|---|---|
| `Armistead_CellDeathDis2024` | `minutes` | −301.9161878 | lin | 14 | 58 | gntr | 5.8e−06 |   | ✅ **solved** |
| `Bertozzi_PNAS2020` | `minutes` | 158.8642627 | lin | 8 | 22 | gntr | 5.4e−06 |   | ✅ **solved** |
| `Blasi_CellSystems2016` | `minutes` | −1090.5618246 | ln | 9 | 252 | gntr | −4.3e−07 | ✓ | ✅ **solved** |
| `Boehm_JProteomeRes2014` | `minutes` | 138.2219682 | lin | 9 | 48 | gntr | 0.0012 | ✓ | ✅ **solved** |
| `Bruno_JExpBot2016` | `minutes` | −46.6881979 | lin | 13 | 77 | gntr | 1.1e−05 | ✓ | ✅ **solved** |
| `Crauste_CellSystems2017` | `minutes` | 190.4570655 | lin | 12 | 21 | gntr | 0.583 | ✓ | ✅ **solved** |
| `Fiedler_BMCSystBiol2016` | `hours` | −58.5839553 | lin | 22 | 72 | gntr | 1.004 |   | ✅ **solved** (not saturated) |
| `Laske_PLOSComputBiol2019` | `hours` | 276.0540613 | lin/ln | 13 | 42 | gntr | −1e−06 | ✓ | ✅ **solved** |
| `Perelson_Science1996` | `minutes` | 222.2807689 | log10 | 3 | 16 | cmaes | 5e−7 |   | ✅ **solved** |
| `Rahman_MBS2016` | `minutes` | 21.1534861 | lin | 9 | 23 | gntr | 0.000000 | ✓ | ✅ **solved** |
| `Raia_CancerResearch2011` | `hours` | 345.3097673 | lin | 39 | 205 | gntr | 9e−06 |   | ✅ **solved** |
| `SalazarCavazos_MBoC2020` | `minutes` | 366.8615730 | lin | 6 | 18 | gntr | 2.9e−05 | ✓ | ✅ **solved** |
| `Sneyd_PNAS2002` | `minutes` | −319.7923458 | lin | 15 | 135 | gntr | 1.4e−5 |   | ✅ **solved** |
| `Schwen_PONE2015` ‖ | `hours` | 952.4217251 | log10 | 30 | 286 | gntr | −12.55 ¶ | ✓ | ✅ **solved** |
| `Elowitz_Nature2000` | `hours` | −65.6351201 | log10 | 21 | 58 | cmaes | 2.43 † |   | ⚪ setup only |
| `Borghans_BiophysChem1997` | `hours` | −132.0084765 | log10 | 23 | 111 | cmaes | 48.7 † | ✓ | ⚪ setup only |
| `Zhao_QuantBiol2020` | `minutes` | 501.2270538 | lin | 28 | 82 | gntr | 5e−06 | ✓ | ✅ **solved** |
| `Brannmark_JBC2010` | `minutes` | 141.8248543 | lin | 22 | 43 | gntr | 0.064 † | ✓ | 🟢 objective validated ‡ |
| `Giordano_Nature2020` | `hours` | −3488.3414981 | lin | 50 | 313 | gntr | 0.135 | ✓ | ✅ **solved** |
| `Weber_BMC2015` | `minutes` | 296.2020025 | lin | 36 | 135 | gntr | 0.781 | ✓ | ✅ **solved** (not saturated) ‡ |
| `Okuonghae_ChaosSolitonsFractals2020` | `minutes` | 373.5476580 | lin | 16 | 92 | gntr | 0.0012 |   | ✅ **solved** |
| `Oliveira_NatCommun2021` | `minutes` | 7904.9343174 | lin | 12 | 120 | gntr | 0.011 |   | ✅ **solved** |
| `Smith_BMCSystBiol2013` ✱ | `hours` | 20922.1642440 | lin | 25 | 62 | gntr | 0.502 | ✓ | ✅ **solved** |

`k` = free parameters, `n` = scored data points.
**† = optimality gap at the PEtab nominal point, not from a fit.** Only the twenty ✅ rows report an
OG from an actual optimization run.

**"(not saturated)" = the fit cleared the threshold without reaching the best point the problem is
known to have.** Two rows say this and they say it for different reasons, which is why the phrase is
not a single caveat. `Fiedler_BMCSystBiol2016` ran out of budget while still descending (1.113 → 1.034
→ 1.004 over its last ~90 min) short of a basin its own nominal point proves is there.
`Weber_BMC2015` converged — 75 of its 100 starts retired on `step is negligible`, none on
`max_iterations` — but to an interior optimum 0.78 NLL units above its published one, because
**five of its 36 parameters sit exactly on their box bounds at that published optimum**
(`a22`/`a32`/`a33` at `1e-4`, `m11` at `1e10`, `pu3` at `1e8`). The reference basin is a *corner* of
the box, which a box-uniform multistart under a trust region reaches far less readily than an interior
point; this fit recovers one of the five (`m11`, to within 0.005% of its bound). Neither row is a
leverage trade in the `Schwen_PONE2015` ¶ sense — Weber's five fitted σ come back within 0.7% of their
published values, and it *beats* the published point on three of its eleven observable/experiment
groups. So a `(not saturated)` OG is a sound benchmark number and a weak claim about recovering
published kinetics; read the slug's `VALIDATION.md` before quoting it as the latter.

**`Weber_BMC2015` is the one row whose recipe needs a hand-set ODE tolerance**, and removing it makes
the problem unfittable rather than merely slower. Its conf carries `sbml_atol = 1e-4` because ADR-0103
derives `atol = rtol × median(y₀)` = `4.665e-3` for this model and then clamps it to `1e-8` — the
derivation "only ever tightens", so it discards its own answer for one **5.7 decades tighter** — and
ADR-0105's per-species vector clamps into `[scalar_atol, default_atol]` = `[1e-8, 1e-8]` here, so it
is elementwise the scalar and correctly declines to engage. At the derived value **only 6 of 30
box-sampled starts integrate** under the sensitivity request the gradient path applies, against 22 of
30 at `1e-4`. Two things worth carrying to any future slug: the failure is in the forward-**sensitivity**
solve rather than the state solve (the plain forward model manages 7 of 11 box points at the derived
tolerance, in 0.6 s), and the tolerance is a pure trade against *objective* noise — the assembled
gradient is invariant across the sweep, while the FD reference degrades from 2.78e−03 at `1e-4` to
8.60e−02 at `4.665e-3`, which is what a line search consumes. Full sweep in that slug's
`VALIDATION.md`. Filed upstream as **lanl/PyBNF#557**, which records that this clamp binds on **10 of
the 22** slugs whose nominal state is readable — `Crauste`, `Laske`, `Okuonghae`, `Oliveira`,
`Perelson`, `Rahman`, `Raia`, `Smith`, `Weber`, `Zhao` — though nine of the ten already solve at the
clamped value, so the blast radius of changing the default is why it is an opt-in request.

**✱ = this row — nominal point *and* fit — was produced on a stack no other row shares.** It is two
differences, not one, and both are load-bearing.

*The objective.* `Smith_BMCSystBiol2013` is the only slug of the 23 whose SBML declares
`hasOnlySubstanceUnits` species *and* gives its compartments a size other than 1 (down to `6.4e-14`).
Until lanl/PyBNF#553 (merged 2026-08-10 as `85b36a96`) the bridge handed observable formulas each
such species as `amount / compartment_size`, inflating all 62 scored points by up to 1.6e13 and
putting this row's OG at `6.9e+32`. That number was an artefact of the defect, not a distance from
the optimum; the honest nominal-point figure is `8.7e+05`. Every other row was produced before that
fix and is unaffected — re-evaluating all 23 at their nominal points moves only this one, the other
22 bit-identical.

*The gradient.* The `OG = 0.502` fit additionally required a **local, editable `bngsim` working copy
at `dff901e`**, not the released PyPI wheel. Reaction 7 (`R5f`) is `per_species_volume_scaling`, and
without bngsim#160/#161 the whole model declines the analytic `df/dp`, all 25 columns fall back to
CVODES difference quotients, and every start times out to `inf`. PyBNF's capability gate cannot
detect this — `BNGSIM_HAS_EVENT_SENS` is a **version floor of exactly 0.12.2**, which the PyPI 0.12.2
wheel passes while lacking the commit. Probe capability, never version.

Reproducing this row therefore needs PyBNF at or past `85b36a96` **and** a bngsim carrying
bngsim#161 (`'_psvs_row_divisor' in inspect.getsource(bngsim._codegen)`). It was run against PyBNF
`095a5a14` and bngsim `dff901e`; see the slug's `VALIDATION.md` for the full provenance, including
why `nominal_check.json` still records the wheel's value (they differ by 1.12e-10 relative).

**¶ = solved on the benchmark objective, but *not* a reproduction of the source paper's fit.** Two
things hold at once for `Schwen_PONE2015` and both must be quoted together. Its `J*` is **not
converged** — the PEtab nominal point already scores 8.42 NLL units better, and §2c's independent
oracle reproduces PyBNF exactly at 943.9993, so the reference is the outlier — which makes clearing
the threshold an undemanding test. And its 286 points are bound to just **two** estimated σ
(`IR_obs_std` over 34 FACS points, `std` over 252 ELISA points), where the source paper used
per-point error estimates; with σ profiled, `∂OG/∂log(RMSⱼ) = nⱼ`, so the ELISA assay carries ≈7.4×
the leverage of the panel the paper publishes as its figure. The fit takes that trade: it beats the
published parameter vector on objective (−315.99 vs −311.87) while `observable_IR2` comes out flat or
declining where the data rise 2–3.5×, and its FACS noise scale runs **into its upper bound**. None of
that is a defect — the encoding is upstream's and `J*` is defined on the same objective, so the
benchmark comparison is sound — but a good `OG` here is evidence about the **PEtab objective**, not
about recovering the published kinetics. See that slug's `README.md` callout and `VALIDATION.md`.

**‖ = renamed locally; upstream and Grein et al. still say `Schwen_PONE2014`.** The paper is
[Schwen et al. 2015, *PLOS ONE* 10(7):e0133653](https://doi.org/10.1371/journal.pone.0133653) —
published 2015-07-29, and the model's own SBML `<notes>` say 2015 — so the local slug was corrected
to `Schwen_PONE2015` on 2026-08-07. **Nothing upstream was renamed**: the SBML file is still
`model_Schwen_PONE2014.xml`, `upstream.json`'s `upstream_path` still points at
`Benchmark-Models/Schwen_PONE2014/`, and `jstar.txt` still comes from the `Schwen_PONE2014` row of
`best_fx_marvin.csv`. Because the local slug is therefore **not** a reliable join key for this one
problem, both `nominal_check.json` and `upstream.json` record an explicit `upstream_slug`, and
`tools/sigma_profile.py` reads that field rather than inferring the upstream directory from the
local name. Any future tool that joins to upstream by slug must do the same.

**§ = this row's `OG_nominal` is inflated by placeholder nominal σ, and is not a difficulty
ranking.** `OG_nominal` evaluates the objective at the PEtab `nominalValue` vector — *every*
parameter, the estimated noise parameters included. Where those nominal σ are placeholders rather
than fitted values, the number is dominated by `Σ nⱼ log σⱼ` sitting far from its MLE, and says
almost nothing about how far the *dynamics* are from the optimum. `tools/sigma_profile.py` computes
the honest version: hold every non-noise parameter at nominal, set each estimated σ to its MLE
`σⱼ = √(Σⱼ r²/nⱼ)`, and re-score. That is the best OG reachable without moving the dynamics at all.

| slug | nominal σ | `OG_nominal` | `OG` σ-profiled | inflation |
|---|---|---:|---:|---:|
| `Giordano_Nature2020` | 1 (placeholder) | 3776 | **743.72** | 3032 |
| `Zhao_QuantBiol2020` | 1000 (placeholder) | 276.12 | **135.75** | 140.4 |
| `Brannmark_JBC2010` | already at MLE | 0.06437 | 0.0641 | — (none) |

This matters because issue #38 orders the remaining ⚪ tuning candidates "roughly by nominal-point
distance, i.e. plausibly by difficulty". For a placeholder-σ slug that ordering is measuring the
placeholder, not the problem. It is the same pattern `Laske`'s and `Blasi`'s `VALIDATION.md` already
describe — σ nominals of exactly 1 relaxing to their MLEs during a fit — read here as a correction to
a *ranking* rather than as a note about one slug.

> **The tool self-checks, and the check is the load-bearing part.** Substituting the nominal σ back
> in must reproduce `nominal_check.json`'s `J_paper`; that residual is reported next to the
> inflation, and a profiled number is only meaningful when it is orders smaller. `Giordano` reads
> `5.7e−14` against an inflation of 3032, and `Zhao` `4.7e−06` against 140.4 — both decisive.
> `Brannmark` reads `2.9e−04` against an inflation of `2.9e−04`, i.e. **the same order**, which is
> the tool correctly reporting that it cannot resolve an effect this small; the honest reading is
> that Brannmark's nominal σ are already at their MLE and there is no inflation to correct.
> `Laske_PLOSComputBiol2019` reads `4.9e+02` against `5.5e+02` and is excluded outright — its
> upstream `simulatedData` is not the nominal-point trajectory, so the join computes residuals
> against the wrong point. **Check against `J_paper`, never against `reduced_objective`:** PyBNF
> drops only the *parameter-independent* per-point constants, so `Σ nⱼ log σⱼ` stays inside the
> reduced objective. Comparing against `reduced_objective` reports a spurious failure on every
> estimated-σ slug.

Coverage is the same as the `obj ✓` column's and for the same reason — it is the same oracle. Of the
23 slugs, 9 are computable; the other 14 skip for a reported reason (no `simulatedData`, or rows that
will not join one-to-one), and 5 of the 9 have no estimated σ at all, which makes profiling a no-op.

**obj ✓ = the objective has been checked against an independent oracle.** The Eq. 6 NLL is
recomputed at the nominal point straight from the upstream PEtab tables — `simulatedData` (the
collection's own reference simulation) joined to `measurementData`, with the declared
`observableTransformation` and nominal σ — with **no PyBNF in the loop**, and compared to what PyBNF
reports at the same point. `✓` = reproduces PyBNF (14 slugs) — twelve of them exactly, plus `Smith`
at 3.1e−09 and `SalazarCavazos` at 2.0e−08 relative, which is integrator tolerance rather than a
digit-for-digit match. `✗` = disagrees, i.e. a defect — **no row carries one now**; the two that did
(`Brannmark`, `Weber`) were fixed by lanl/PyBNF#547 and now reproduce. Blank = not checked: either
upstream ships no `simulatedData` (`Bertozzi`, `Okuonghae`, `Oliveira`), or the rows could not be
joined (`Fiedler`, `Raia`), or the checker cannot key them (`Elowitz`, below), or the checker's own σ
handling is the doubtful half rather than PyBNF's (`Armistead` and `Sneyd`, where PyBNF matches `J*`
to 0.0000 and 0.0006, and `Perelson`).

> **`Elowitz_Nature2000` is a checker limitation, not a data property, and it is the same bug as
> Smith's two columns over.** Its rows join **58 of 58, one-to-one** on the PEtab identity key. The
> checker reports `joined 0 of 58` because it *also* keys on `observableParameters` and
> `noiseParameters`, which hold different kinds of value on the two sides — parameter **names** in
> `measurementData` (`background;scale`, `sigma`) against the **resolved numeric values** in
> `simulatedData` (`-4.98107438218408;-0.279017032524776 `, `-1.14362512938604`, log10 scale). Those
> can never match. Dropping the two columns is **not** the fix: verified across all 23, it breaks σ
> resolution on nine slugs — once a column is not a join key, the merge renames it
> `noiseParameters_meas`/`_sim` and the resolver stops finding it — and makes `Armistead`, `Fiedler`,
> `Weber`, `Blasi` and `Schwen` over-match. The fix is to key on identity while reading σ from the
> measurement side. Until then `Elowitz`'s oracle is **untested, not absent**, and this row's blank
> should not be read as evidence either way. `Fiedler` and `Raia` are different: their rows genuinely
> do not join one-to-one on identity keys at all.

> **`Smith_BMCSystBiol2013` moved out of the "could not be joined" list on 2026-08-11, and the reason
> is worth keeping.** Its rows always joined; the checker's key did not. `datasetId` is a PEtab
> *visualization grouping label*, not part of a measurement's identity, and upstream tags 13 of
> Smith's 62 rows inconsistently across the two tables — all ten of figure 2C and all three of
> figure 2D read `fig2C`/`fig2c`/`fig2D` in `measurementData` and `fig2A` in `simulatedData`.
> Including the column joined 62 − 13 = 49 and the slug was written off. On the identity keys
> `(observableId, simulationConditionId, time)` it joins 62 of 62, one-to-one. `tools/sigma_profile.py`
> no longer keys on `datasetId`; the change was diffed across all 23 slugs and moves only this one.

**This column is orthogonal to the status column, and that distinction is the point.** A ⚪ row with
`✓` is a job whose objective is known good and whose nominal point simply is not the optimum
(`Zhao`, 136 away once its placeholder σ are profiled out — see § — and genuinely unrun). A ⚪ row
with no mark is a job nobody has checked. Conflating
those two is exactly how `Weber_BMC2015` sat as an ordinary ⚪ "ready-to-run" job whose objective was
wrong by 13,740, one command away from being handed a multi-start budget. **That is not a
hypothetical** — it is what happened, and this column exists because of it.

**‡ = these two rows were `⚠️ blocked` until 2026-08-07, and their old OGs should not be used.** Both
use **pre-equilibration**, which the bngsim SBML backend silently dropped: it neither applied the
conditions nor ran the unmeasured equilibration phase, so the scored run simulated the model exactly
as authored. `Brannmark`'s eight dose experiments came out byte-identical and `Weber`'s trajectory
was flat across the timepoint where its dose fires, and the objective was wrong before any optimizer
started. They showed `OG_nominal` of `1.5e+03` and `1.4e+04`, which were recorded as "the nominal
point is not the optimum" — it is. **lanl/PyBNF#547 (ADR-0104)** fixed it, and both now land on the
independent NLL recomputed from upstream's own `simulatedData` tables: `Weber` at `−0.0002` and
`Brannmark` at `+0.064`. Checked one level deeper for `Brannmark`, whose observable formulas are
simple enough to evaluate by hand, its trajectory now agrees with that reference simulation to
`1.4e−05` relative across all 43 measured points, where before 42 of the 43 were off by more than
0.1% (worst: a factor of 8.5). Exactly 2 of the 23 slugs pre-equilibrate and exactly those 2 were
affected; every other slug checked that way reproduced PyBNF exactly throughout.

Four status levels, and the difference matters:

- ✅ **solved** — a real fit was run and reached `OG < 1.92`.
- 🟢 **objective validated** — no fit was run, but the problem's `nominalValue` point *is* its
  published optimum, and PyBNF's objective there lands within the solved threshold of `J*`. This
  validates the import and the objective; it makes no claim about PyBNF's optimizer.
- ⚪ **setup only** — the job imports, simulates, and scores correctly, but its nominal point is not
  the published optimum, so nothing about optimality is claimed. These are ready-to-run jobs.
- ⚠️ **blocked** — the job imports and runs, but PyBNF's objective for it is *wrong*, so neither its
  OG nor any fit against it means anything until the upstream defect lands. Not a ready-to-run job:
  do not spend a fitting budget here. **No slug is in this state**; `Brannmark` and `Weber` were,
  from 2026-08-06 until lanl/PyBNF#547 landed on 2026-08-07. The level is kept because the failure it
  names is silent by nature and will recur.

### Coverage is now complete

Every import blocker has landed: **lanl/PyBNF#492–#496**, then **#508 (Fiedler, replicate
`observableParameters`), #509 (Laske, natural-log `lnnormal` noise), #510 (Schwen, `t=0`-only
experiment)**, and finally **#521 (Blasi, steady-state `time = inf` measurements, ADR-0086)**. That
took importable coverage from 15/23 to **23/23**, all of which also run.

`Blasi_CellSystems2016` was the last to land and needed two fixes rather than one. All 252 of its
measurements are stationary abundances at `t = inf` — no time axis, no dose axis — with observables
log-transformed under normal noise. #509 (`lnnormal`) let it *import*; config load then still crashed
in `TimeCourse` (`OverflowError: cannot convert float infinity to integer`), because the importer
materialized `time = inf` as a time course instead of using PyBNF's steady-state support
(`steady_state=>1` ParamScan, ADR-0046). #521 closed that gap, and the problem now fits in about 25 s
and **solves** — recovering the published maximum-likelihood rate constants of Blasi et al. (2016) to
four digits. It is also the one slug whose nominal point is *not* its optimum despite carrying the
published parameter values, because PEtab's `nominalValue` for its noise scale is a placeholder; see
its `VALIDATION.md`.

## Optimizer choice

`job_type` per slug was chosen by **running each candidate**, not by assumption:

- **`gntr`** — the general-objective Fisher/Gauss-Newton trust region (EFIM Hessian fed through
  `trf`'s Coleman–Li core, ADR-0068). This is the method that handles an **estimated noise scale**;
  plain `trf` refuses it (`trf` needs "an exact least-squares residual … with a fixed noise scale").
  It is PyBNF's fides-analogue and the default here. It works on log10/`lognormal` objectives too —
  the EFIM Fisher block for noise scales landed in ADR-0079/0080/0081.
- **`cmaes`** with IPOP restarts (ADR-0070, restart trigger fixed in ADR-0082) — used where the
  gradient path genuinely refuses the problem, and for the three strongly multimodal problems
  (Borghans, Elowitz, Okuonghae), where a local method from a few starts lands in a local basin.

**No** gradient-path refusal remains. The last one was `Smith_BMCSystBiol2013`, and it is worth
recording how it went away, because the entry stood here after it had stopped being true.

The model contains discrete events (state-dependent jumps), across which forward output
sensitivities went stale, so bngsim refused them and lanl/PyBNF#461 hoisted that refusal into a
blanket pre-flight gate (`_require_differentiable_dynamics`). bngsim applies the event's own jump
now, and PyBNF lifts the gate at `BNGSIM_HAS_EVENT_SENS` — which is a **version floor of exactly
0.12.2**, not a capability probe. So the refusal ended at the 0.12.2 release, and this table was
already stale before anyone tried `gntr` here.

Two consequences worth keeping, since the version string cannot distinguish the builds:

- A build that reports `0.12.2` passes the gate whether or not it carries bngsim#160/#161, which
  emits the analytic sensitivity RHS for a cross-compartment reaction. Smith's reaction 7 is
  `per_species_volume_scaling`; without that commit the whole model declines the analytic `df/dp`
  and all 25 columns fall to CVODES difference quotients, which blows the wall clock and returns
  `inf` for every start. The gate says yes and the fit still cannot run, with nothing pointing at
  why.
- `gntr` is the recipe here now. It was blocked by a second, unrelated defect until 2026-08-10 —
  see the `nominal_check.json` note and lanl/PyBNF#553 below.

The ADR-0076 condition-routing line is now **clear**, in two steps:

- `Bruno_JExpBot2016` left the list at **lanl/PyBNF #511** (merged #513, 2026-07-23), which taught
  `route_experiment` to compose the chain rule for a free parameter that reaches the model only
  through a `condition:` parameter reference. It fits on the gradient path and **solves**
  (`OG = 1.1×10⁻⁵` in 41 s). Its `VALIDATION.md` records the history.
- `Bertozzi_PNAS2020` left it at **lanl/PyBNF #530** (ADR-0095, 2026-08-02), which stopped assuming
  that seed derivative is a plain `1`. Bertozzi needed all of it: `I0_` seeds two species with
  opposite signs (`I_ = I0_`, `S_ = N_ - I0_`), and `R0_`/`gamma_` reach the dynamics only through
  the derived `beta_N = R0_*gamma_/N_`, whose derivatives are expressions re-evaluated at each fit
  point. Note that the last of those was never covered by the refusal at all — those two columns
  were silently wrong rather than refused, which is why #530 could not be landed without #531.
  It now **solves** (`OG = 5.4×10⁻⁶`) and its conf carries `job_type = gntr`.

Every shipped conf was verified to start and complete a tiny run.

**These recipes are reasonable defaults, not tuned ones — and that is a real limitation.** A full
`gntr` run of the shipped `SalazarCavazos_MBoC2020` conf (20 starts × 500) converged to `OG = 10.2`,
i.e. *worse* than that problem's own nominal point (`OG = 0.33`): 20 box-sampled starts were not
enough to find the reference basin. `Laske_PLOSComputBiol2019` is the worked example in the other
direction, and it is now **solved**: the collection-default 20 × 500 reaches only `OG = 6.76` on it,
while 100 × 1000 — the budget its conf now carries — reaches the reference optimum itself. Expect to tune
`population_size` / `max_iterations`, or to switch to `cmaes` with IPOP restarts, before treating any
⚪ or 🟢 row as a statement about PyBNF's optimizers. The twenty ✅ rows are the only ones where a fit
was actually driven to `OG < 1.92`.

`SalazarCavazos_MBoC2020` is the second worked example, and the sharper one: at 20 × 500 it reaches
`OG = 10.2` — *worse* than its own nominal point — and at 100 × 1000 it lands on `J*`. **100 × 1000 is
now the working default for this collection**; a conf still carrying 20 × 500 is an untuned placeholder.
The one problem where that budget has not sufficed is `Fiedler_BMCSystBiol2016`, solved at `OG = 1.004`
without reaching a reference basin its own nominal point proves is there.

## Import + fit pipeline (reproduce)

```bash
# 1. import — parameterScale (#491) AND observableTransformation (#499) are preserved,
#    both of which plain petab.v2.petab1to2 drops. No hand-editing of noise_model is needed.
python -c "from pybnf.petab import petab1to2_preserve_scale, import_job; import tempfile; \
  v2 = petab1to2_preserve_scale('<id>/<id>.yaml', tempfile.mkdtemp()); \
  import_job(v2, 'out', job_type='gntr')"
# 2. set the run recipe: sbml_backend = bngsim (REQUIRED for gradients), wall_time_sim, random_seed
# 3. fit + score
pybnf -c <id>.conf -o
python score.py output          # OG = -log_likelihood - J*
```

`wall_time_sim` caps pathological parameter points, where CVODE would otherwise grind for seconds
before failing; a timed-out simulation is simply a failed one, so the optimum is unchanged. Raise it
if valid simulations on your machine are being marked as failures.

## Per-slug contents

`<id>.conf` (the runnable fit) · the SBML model (verbatim) · `experiment*.exp` (data) ·
`*_measparams.tsv` (per-measurement observable/noise parameter tables, where the problem uses them) ·
`jstar.txt` (the reference J\*) · `nominal_check.json` (the nominal-point evaluation) · `score.py` ·
`README.md`. The twenty solved slugs additionally ship `best_fit_params.txt`,
`information_criteria.txt`, and `VALIDATION.md` from their fits.
