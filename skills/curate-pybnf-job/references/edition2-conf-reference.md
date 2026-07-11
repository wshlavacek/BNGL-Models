# Edition-2 (new-era) PyBNF `.conf` reference

Authoritative config surface for authoring an edition-2 job setup, grounded in the
PyBNF source (`~/Code/PyBNF`, v1.6.0). Every claim cites `file:line`. Items marked
**[E2]** require `edition = 2`; **[LEGACY]** items are edition-1-only and are
*rejected* under edition 2 — never use them here. When in doubt, open the cited
source or a tutorial `.conf` under `~/Code/PyBNF/examples/tutorial/`.

## Contents
1. Canonical skeleton & syntax
2. `edition` gate
3. Run selector: `job_type` (+ allowed codes)
4. Objective surface: `objective` / `noise_model` / `profile_objective`
5. Free parameters: `*_var` families and the `parameter:` record
6. Noise models (families, verbs)
7. Block directives: `model:` / `condition:` / `experiment:` / `observable:` / `normalization`
8. General / execution options
9. Method-specific keys (per `job_type`)
10. Comment idioms

---

## 1. Canonical skeleton & syntax

General syntax (`parse.py:105-484`): one directive per line; `key = value` for
scalars; `directive: name, field: value, ...` for block directives (colon-delimited,
fields in any order, only the required field mandatory); `#` starts a comment
(`parse.py:108`); blanks/comment-only lines are skipped (`parse.py:506`); keywords are
case-insensitive.

Canonical order (matches every real-world conf and the tutorial):

```conf
output_dir = output       # relative -> lands in <slug>/output/ (run pybnf from the job folder)
edition = 2
bngl_backend = bngsim
model: <name>.bngl

# conditions (only if the experiment perturbs / pre-equilibrates)
condition: <c>, perturbations: <param> = <val>

# experiments (each binds one .exp; this is what makes exp_data non-empty)
experiment: <name>, [preequilibrate: <c0>,] [condition: <c>,] [type: parameter_scan,] [method: ssa|nf,] data: <file>.exp

job_type = de
objective = sos            # or chi_sq when the .exp has _SD columns

max_iterations = 50
population_size = 12
verbosity = 2

# free params bind to model ids by NAME (ADR-0034) — no __FREE alias
loguniform_var = <model_param_id> <lo> <hi>
```

---

## 2. `edition` gate  **[E2]**

`edition = 2` opts into the modern language (`edition.py:1-27,38`). Absence ⇒ implicit
edition 1 (legacy). Edition 2 requires PyBNF ≥ 1.5.0 (`edition.py:46-49`). Every E2
feature below is guarded by `require_edition(ed, 2, feature)` which raises an
explanatory error, not a silent reinterpret, if used without `edition = 2`
(`edition.py:107-126`).

---

## 3. Run selector: `job_type = <code>`  **[E2]**

`job_type` names the algorithm; there is **no default** — you must set it
(`config.py:510-515`). Using legacy `fit_type` under edition 2 is an error
(`config.py:504-509`). Codes (`registry.py:63-90`):

| Code | Family | Meaning |
|---|---|---|
| `de` | optimizer | Differential Evolution (island/synchronous) — the workhorse for paper fits |
| `ade` | optimizer | Asynchronous DE |
| `pso` | optimizer | Particle Swarm |
| `ss` | optimizer | Scatter Search |
| `sa` | optimizer | Simulated Annealing *(deprecated)* |
| `sim` | optimizer | Nelder–Mead Simplex (also the default `refine` polisher) |
| `powell` | optimizer | Powell |
| `cmaes` | optimizer | CMA-ES |
| `trf` | optimizer | Trust-Region Reflective least-squares (**gradient; needs `bngl_backend = bngsim`**) |
| `lbfgs` | optimizer | L-BFGS-B (**gradient; needs bngsim**) |
| `profile_likelihood` | optimizer | Profile-likelihood confidence intervals |
| `mh` / `pt` / `am` / `dream` / `p_dream` / `hmc` | sampler | Bayesian posteriors (MCMC/HMC) |
| `check` | checker | Model checking against BPSL `.prop` constraints (no fitting) |

For a first paper fit, prefer a global metaheuristic (`de`, `ss`, `pso`) — that is
what the 2019-paper corpus uses (`examples/real-world/*/*.conf`).

---

## 4. Objective surface  **[E2]**

Under edition 2, `objfunc` is rejected (`config.py:2131-2138`). Set **exactly one** of
`objective`, `noise_model`, `profile_objective` (`config.py:2139-2153`).

### `objective = <token>` — the common case (`objective.py:1457-1468`)

| token | Meaning | Use when |
|---|---|---|
| `sos` | plain sum of squares (Gaussian, σ=1) | `.exp` has **no** `_SD` columns |
| `chi_sq` | Gaussian, σ read per point from the `_SD` columns | `.exp` **has** `_SD` columns |
| `chi_sq_dynamic` | Gaussian, σ fit as one free param | one global unknown noise level |
| `sod` | sum of absolute differences (Laplace, scale=1) | robust to outliers |
| `norm_sos` | per-point relative error (Gaussian, relative σ) | data spans orders of magnitude |
| `ave_norm_sos` | column-mean-normalized | mixed-scale observables |
| `lognormal` / `laplace` | log-Gaussian / Laplace shortcuts | — |

`neg_bin`, `lognormal`, `expression`, `callable`, and the inline analytical targets
(`banana`, `gaussian`, …) also exist but are **not PEtab-exportable** (except `sos`,
`chi_sq`, `sod`, `norm_sos`, `ave_norm_sos`, and per-observable Gaussian/Laplace) — see
`petab-compliance.md`. Keep real-world job setups inside the exportable subset.

### `noise_model = <family>, …` — assemble a likelihood (see §6). `profile_objective = kl|wasserstein` — column-joint shape objectives (not exportable).

---

## 5. Free parameters

### 5a. `<family>_var = <id> <p1> [<p2>] [b|u]`

Binds by **name** to the model parameter `<id>` (ADR-0034 — no `__FREE` marker;
`config.py:643-644,700-701`). `log`-prefixed variants search on a log10 scale. Bounded
families take an optional trailing `b`/`u` reflect-bounds flag. Families
(`priors/*.py`):

| keyword | prior | args (`p1 p2`) |
|---|---|---|
| `uniform_var` / `loguniform_var` | Uniform box (bounded) | `lower upper` |
| `normal_var` / `lognormal_var` | Normal | `mean sd` |
| `laplace_var`, `cauchy_var`, `gumbel_var`, `logistic_var` (+ `log`*) | 2-param | `location scale` |
| `gamma_var`, `inv_gamma_var`, `weibull_var` (+ `log`*) | 2-param, support >0 | `shape scale` |
| `beta_var` / `logbeta_var` | Beta, support [0,1] | `alpha beta` |
| `exponential_var`, `rayleigh_var`, `half_normal_var`, `half_cauchy_var`, `chisquare_var` (+ `log`*) | 1-param | one value |
| `var` / `logvar` | no prior (Simplex start point) | value [step] |

Most paper fits use `uniform_var` (linear scaling factors) or `loguniform_var` (rate
constants over decades), bracketing the published value — exactly the receptor/igf1r/
egfr_ode confs.

### 5b. `parameter:` record (fully labeled, for priors)  **[E2]**

`parameter: <id>[, prior: <family>][, parameter_scale: lin|log10|ln][, <field>: <n> …][, lower: <n>, upper: <n>][, initial_value: <n>]` (`config.py:2309-2410`). Use when
you need a prior a positional keyword can't spell (e.g. `student_t`) or explicit
truncation bounds. `parameter_scale: log` is rejected as ambiguous — use `log10`/`ln`.

Noise/observable nuisances (e.g. `sigma`, `scale`, `b_C`) are declared like any free
parameter (`config.py:2261-2283`).

---

## 6. `noise_model` families & verbs  **[E2]**

```
noise_model [<observable>] = <family>, <param> = <verb> [<arg>][, …][, location = mean|median][, cumulative]
```

Omit `<observable>` ⇒ whole-fit default; include it ⇒ per-observable override that
layers over the default (`objective.py:1479-1488`; example `10_per_observable_noise`).

Families (`objective.py:583-590`): `normal`/`gaussian` (param `sigma`), `lognormal`
(`sigma`), `laplace` (`scale`), `student_t` (`sigma`, `df`=4), `neg_bin`
(`dispersion`). **PEtab-exportable:** `normal`/`gaussian` and `laplace`. `neg_bin`,
`lognormal`, `student_t` are not (see `petab-compliance.md`).

Sigma-source verbs (`objective.py:620-662`): `read_exp_file _SD` (per-point from data),
`fix_at <n>`, `fit <name>` (estimate a free param), `relative [cv]`, `column_mean`,
`formula <expr>`. Example: `noise_model = normal, sigma = read_exp_file _SD` (identical
to `objective = chi_sq`).

---

## 7. Block directives

### 7a. `model: <file>[, <file>…]`  **[E2]**
Pure declaration, **no data binding**; repeatable (multi-model joint fits). Stem =
modelId. Extensions `.bngl` / `.xml` / `.ant` / `.target` (`parse.py:159,167-179`).
Never use legacy `model = m.bngl : data.exp`.

### 7b. `condition: <name>, perturbations: <var op val>[, …]`  **[E2]**
A named parameter state (PEtab Condition). `op` ∈ `= * / + -` (`=` absolute; others
relative to nominal). `model: <file>` sub-field required only with multiple models
(`config.py:1218-1258`). Example: `condition: withligand, perturbations: Ligand_isPresent = 1`.

### 7c. `experiment: <name>, …, data: <file>.exp`  **[E2]**
Binds a `.exp` and synthesizes the simulation (`parse.py:355-428`, `config.py:1260+`).
Fields (any order, only `data:` required):

| field | meaning |
|---|---|
| `data: <f…>` | comma list of `.exp` (replicates) / `.con` / `.prop` (BPSL constraints) files. Mixing `.exp` + `.prop` = data fusion; `.prop` alone = constraint-only (then `t_end:` is required). BPSL makes the job native-only — see `bpsl-constraints.md` |
| `condition: <c>` | apply a named Condition (omitted ⇒ wildtype) |
| `preequilibrate: <c0>` | equilibrate to steady state under `<c0>`, unmeasured, before measuring; state carries over (ADR-0052) |
| `type: time_course \| parameter_scan` | overrides the data-driven inference |
| `method: ode \| ssa \| pla \| nf` | simulator (default `ode`) |
| `t_end: <t>` | scan measurement time (absent ⇒ steady state) / constraint-only endpoint |
| `t_start`, `n_steps` | constraint-only time-course grid |
| `equil_t_end: <t>` | **required** for `method: nf` pre-equilibration (NFsim has no steady-state solve) |
| `gml: <n>`, `complex: 0\|1` | NFsim options (only `method: nf`) |

**Data-type inference:** a `.exp` whose first column is `time` ⇒ time course; any other
first column ⇒ `parameter_scan` over that column (a model parameter). A scan with no
`t_end` runs to steady state. Examples:
`experiment: timecourse, data: timecourse.exp`;
`experiment: doseresponse, type: parameter_scan, t_end: 1200, data: doseresponse.exp`;
`experiment: receptor, preequilibrate: noligand, condition: withligand, data: receptor.exp`.

### 7d. `observable:` — two forms  **[E2]**
1. Column rename: `observable: <entity>, column: <header>` — map a data column to a
   model observable name.
2. Measurement model: `observable: <id>, formula: <PEtab expr>` — a post-simulation
   `observableFormula` (ADR-0036), e.g. `observable: obs_scaled, formula: scale * Obs_B`.
   Both are PEtab-exportable.

### 7e. `normalization <obs> = init|peak|zero|unit`  **[E2]**
Per-observable prediction transform (ADR-0053). **NOT PEtab-exportable** — avoid it in
real-world job setups (the shipped `igf1r.conf` uses `normalization = init` and is
therefore not exportable). If the data is genuinely relative, prefer encoding the
reference in the model/observable instead.

---

## 8. General / execution options (`config_schema.py`)

| key | default | meaning |
|---|---|---|
| `output_dir` | `pybnf_output` | results dir; **relative to where you invoke pybnf**, not the conf's location — use a plain `output` and run pybnf from inside the job folder, so results land in `<slug>/output/` |
| `bngl_backend` | `auto` | **set `bngsim`** (new-era default; required for gradient fits) |
| `delete_old_files` | `1` | clear stale output |
| `verbosity` | *(required)* | 0–3 |
| `max_iterations`, `population_size` | *(required)* | search budget |
| `refine` | `0` | polish best fit with a local optimizer (`refine = 1`) |
| `wall_time_sim` | auto (3600 s) | per-simulation timeout — raise for big networks |
| `max_failed_simulations` | `100` | abort threshold |
| `smoothing` | `1` | average N stochastic replicates (SSA/NF) |
| `normalization = <type>` | — | whole-fit form (legacy-compatible; still not PEtab-exportable) |
| `bootstrap` | `0` | bootstrap replicates for uncertainty |
| `random_seed` | `None` | reproducibility |

`sbml_backend` (`roadrunner`/`bngsim`), `sbml_integrator`, cluster keys, etc. exist but
are rarely needed for a BNGL paper fit.

---

## 9. Method-specific keys (only the common ones)

Each `job_type` narrows to its own schema; unrecognized keys warn.
- **DE (`de`):** `mutation_rate` 0.5, `mutation_factor` 0.5, `stop_tolerance` 0.002,
  `de_strategy` `rand1`; island keys `islands` 1, `migrate_every` 20, `num_to_migrate` 3.
- **PSO:** `cognitive` 1.5, `social` 1.5, `particle_weight` 0.7.
- **Scatter Search (`ss`):** `local_min_limit` 5.
- **Refinement:** `refine = 1` + `refine_method = sim|powell|cmaes|trf|lbfgs`;
  `simplex_max_iterations` bounds the polish.
- **Samplers (`dream`, `pt`, …):** `burn_in`, `sample_every`, `credible_intervals`, etc.
  (`samplers/base.py:27-60`).

Full per-method tables: `optimizers/*.py`, `samplers/*.py` in the PyBNF source.

---

## 10. Comment idioms (match the corpus)

Every real-world/tutorial conf opens with a banner: a title line naming the biology +
`job_type`, prose on the model/data/protocol, the source paper, and a
`pybnf -c <conf>` run line — then the directives, grouped with `### … ###` section
headers (`General Options`, `Fitting Options`, `Experiments`, `Parameters`). Cite the
paper (PMCID/DOI) and which figure/table each `.exp` came from. See
`examples/real-world/egfr_ode/egfr_ode.conf` for the model to imitate.
