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
| SBML model, verbatim from upstream | 23 | 1,681,108 | **copied** — 80% of bytes, 7% of files |
| `.exp` and `_measparams.tsv` — PyBNF-format *translations* of the upstream measurement tables, emitted by the importer | 185 | 141,974 | derived |
| `jstar.txt` — one number, transcribed from the ICB-DCM suppl repository | 23 | 433 | derived |
| `.conf`, `score.py`, `nominal_check.json`, `README.md`, `VALIDATION.md`, `best_fit_params.txt`, `information_criteria.txt` | 116 | 323,552 | **written here** |
| total | 347 | 2,147,067 | |

The SBML files dominate the byte count and nothing else, which is why the directory *looks* mostly
vendored while being 93% non-copied by file count. **Editing this directory is normal**: the
locally authored files are the deliverable, and adding a job necessarily edits the coverage matrix
below. Nothing upstream is ever modified.

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

### Why `objective = sos` is faithful for Oliveira and Smith

These two ship a plain sum of squares and no `noise_model` line, which looks like a fidelity break. It
is not, and the reason is worth recording because it is not obvious from the conf alone. **Both
problems specify σ ≡ 1 upstream**, checked against the pinned commit:

- `Smith_BMCSystBiol2013` — `observables_*.tsv` carries the literal `noiseFormula = 1.0` for all seven
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

**This identity is now corroborated across the collection, not just argued.** For ten problems the
PEtab `nominalValue` point is the published optimum, and evaluating PyBNF's objective there
reproduces the paper's `J*` to within the solved threshold — six to ~10⁻³ or better (four of those
to ~10⁻⁵), plus `Fiedler_BMCSystBiol2016` at ~2e−3 (see `nominal_check.json` in each slug). That is
an end-to-end check of the whole imported chain: SBML model → simulation → observable formulas →
noise model → objective. (For the log10 slugs the nominal check carries the change-of-variables
Jacobian and is recorded but not asserted as a validation.)

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

| slug | J\* | scale | k | n | optimizer | OG | status |
|---|---|---|---|---|---|---|---|
| `Armistead_CellDeathDis2024` | −301.9161878 | lin | 14 | 58 | gntr | 5.8e−06 | ✅ **solved** |
| `Bertozzi_PNAS2020` | 158.8642627 | lin | 8 | 22 | gntr | 5.4e−06 | ✅ **solved** |
| `Blasi_CellSystems2016` | −1090.5618246 | ln | 9 | 252 | gntr | −4.3e−07 | ✅ **solved** |
| `Boehm_JProteomeRes2014` | 138.2219682 | lin | 9 | 48 | gntr | 0.0012 | ✅ **solved** |
| `Bruno_JExpBot2016` | −46.6881979 | lin | 13 | 77 | gntr | 1.1e−05 | ✅ **solved** |
| `Crauste_CellSystems2017` | 190.4570655 | lin | 12 | 21 | gntr | 0.583 | ✅ **solved** |
| `Fiedler_BMCSystBiol2016` | −58.5839553 | lin | 22 | 72 | gntr | −0.0022 † | 🟢 objective validated |
| `Perelson_Science1996` | 222.2807689 | log10 | 3 | 16 | cmaes | 5e−7 | ✅ **solved** |
| `Rahman_MBS2016` | 21.1534861 | lin | 9 | 23 | gntr | 3.9e−06 † | 🟢 objective validated |
| `Raia_CancerResearch2011` | 345.3097673 | lin | 39 | 205 | gntr | 0.78 † | 🟢 objective validated |
| `SalazarCavazos_MBoC2020` | 366.8615730 | lin | 6 | 18 | gntr | 0.326 † | 🟢 objective validated |
| `Sneyd_PNAS2002` | −319.7923458 | lin | 15 | 135 | gntr | 1.4e−5 | ✅ **solved** |
| `Laske_PLOSComputBiol2019` | 276.0540613 | lin/ln | 13 | 42 | gntr | 39.9 † | ⚪ setup only |
| `Schwen_PONE2014` | 952.4217251 | log10 | 30 | 286 | gntr | −8.42 † | ⚪ setup only |
| `Elowitz_Nature2000` | −65.6351201 | log10 | 21 | 58 | cmaes | 2.43 † | ⚪ setup only |
| `Borghans_BiophysChem1997` | −132.0084765 | log10 | 23 | 111 | cmaes | 48.7 † | ⚪ setup only |
| `Zhao_QuantBiol2020` | 501.2270538 | lin | 28 | 82 | gntr | 276 † | ⚪ setup only |
| `Brannmark_JBC2010` | 141.8248543 | lin | 22 | 43 | gntr | 1.5e+03 † | ⚪ setup only |
| `Giordano_Nature2020` | −3488.3414981 | lin | 50 | 313 | gntr | 3.8e+03 † | ⚪ setup only |
| `Weber_BMC2015` | 296.2020025 | lin | 36 | 135 | gntr | 1.4e+04 † | ⚪ setup only |
| `Okuonghae_ChaosSolitonsFractals2020` | 373.5476580 | lin | 16 | 92 | cmaes | 4.7e+05 † | ⚪ setup only |
| `Oliveira_NatCommun2021` | 7904.9343174 | lin | 12 | 120 | gntr | 9.6e+06 † | ⚪ setup only |
| `Smith_BMCSystBiol2013` | 20922.1642440 | lin | 25 | 62 | cmaes | 6.9e+32 † | ⚪ setup only |

`k` = free parameters, `n` = scored data points.
**† = optimality gap at the PEtab nominal point, not from a fit.** Only the eight ✅ rows report an
OG from an actual optimization run.

Three status levels, and the difference matters:

- ✅ **solved** — a real fit was run and reached `OG < 1.92`.
- 🟢 **objective validated** — no fit was run, but the problem's `nominalValue` point *is* its
  published optimum, and PyBNF's objective there lands within the solved threshold of `J*`. This
  validates the import and the objective; it makes no claim about PyBNF's optimizer.
- ⚪ **setup only** — the job imports, simulates, and scores correctly, but its nominal point is not
  the published optimum, so nothing about optimality is claimed. These are ready-to-run jobs.

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

**One** gradient-path refusal remains, and it has nothing to do with condition routing:

| slug | refusal | fixable by |
|---|---|---|
| `Smith_BMCSystBiol2013` | the model contains **discrete events** (state-dependent jumps); forward output sensitivities go stale across a jump, so bngsim cannot supply a gradient there (`_require_differentiable_dynamics`, lanl/PyBNF #461) | nothing in the ADR-0076 line — this needs event-aware sensitivities |

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
enough to find the reference basin. `Laske_PLOSComputBiol2019`, run once on the corrected forward
model, went the other way — `OG = 6.76`, well in from its nominal `39.9`, but still outside the
threshold. Neither run is shipped, and neither slug's status changes. Expect to tune
`population_size` / `max_iterations`, or to switch to `cmaes` with IPOP restarts, before treating any
⚪ or 🟢 row as a statement about PyBNF's optimizers. The eight ✅ rows are the only ones where a fit
was actually driven to `OG < 1.92`.

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
`README.md`. The four solved slugs additionally ship `best_fit_params.txt`,
`information_criteria.txt`, and `VALIDATION.md` from their fits.
