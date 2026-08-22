# tools — regenerating and verifying a slug

Six scripts used to produce and check the numbers this corpus ships. All are run from the
repository with PyBNF available (`uv run --project ~/Code/PyBNF --extra tests python
tools/<script>.py …`).

They exist because each was originally written ad hoc, and each has a non-obvious failure mode
that silently produces a *plausible wrong answer* rather than an error. Those are written down here.

> **Run them from the PyBNF project, not from this one.** `uv run` resolves its project from the
> working directory, so a bare `uv run --extra tests` inside `pybnf-jobs/` finds *this* repo's
> `pyproject.toml`, which defines no such extra, and dies with `Extra 'tests' is not defined`
> rather than with anything about the script. `--project ~/Code/PyBNF` is not optional.

## `nominal_check.py` — regenerate a slug's `nominal_check.json`

```bash
python tools/nominal_check.py <slug-dir> <upstream-petab-dir> [--write]
```

Evaluates PyBNF's objective at the PEtab `nominalValue` point and puts it on the paper's Eq. 6 scale
the same way a real run does — through `likelihood_information_criteria`, which is what writes
`information_criteria.txt` at the end of a fit. Without `--write` it prints the new values and the
old ones for comparison.

Keys the script does not recompute are carried forward from the existing file, so anything hand-added
there survives. `optimizer` is the exception: it is re-derived from the conf's `job_type` on every
run, because as a carried-forward key it went stale twice — `Okuonghae` and `Bruno` both moved from
`cmaes` to `gntr` and kept reporting `cmaes` until 2026-08-10. Any future field that mirrors the conf
should be derived the same way rather than carried.

`<upstream-petab-dir>` is the problem's directory in a checkout of
`Benchmarking-Initiative/Benchmark-Models-PEtab` at the commit `upstream.json` pins. That checkout is
deliberately **not** vendored here (see the collection README) and is gitignored if you keep it under
`pybnf-jobs/`. Both PEtab layouts are handled (`parameters.tsv` and `parameters_<slug>.tsv`).

> **Validate any change to this script against `Bruno_JExpBot2016`.** Bruno's model has no derived
> parameter and no measurement-model layer, so its recorded value is stable across every fix this
> corpus has seen. The script must reproduce `J_paper = -46.688194686350265` **exactly**. If it does
> not, the script is wrong, not the slug.

> **Gotcha: a Result can only be scored once.** A measurement-model layer (ADR-0036) materializes its
> observable columns onto the `Data` in place, so scoring the same `Result` a second time collides
> with the columns the first pass added. The script simulates twice for this reason — a real fit
> re-simulates for its information criteria for exactly the same reason.

## `fd_check.py` — verify a job's assembled gradient against finite differences

```bash
python tools/fd_check.py <slug-dir | job.conf> [param-values.json] [-h STEP]
```

Mirrors `GradientOptimizer.gradient_at` — builds one routing per scored `(model, suffix)`, applies
the union sensitivity request, simulates once, assembles the gradient in sampling space `u` — then
compares against central differences of the same objective in the same space. Every free parameter is
printed, so a **structurally zero** column stands out from a merely small one.

This is what found **lanl/PyBNF#534**, the first time it was pointed at a real model rather than a
fixture: a free parameter that reached the trajectory only by seeding another entity assembled to
exactly `0` against a central difference of `-10.4`.

> **Gotcha: score the Result before assembling.** The script calls `evaluate_multiple` on the
> simulated `Result` first. Skip that and the objective's `_pset_values` is unbound and its
> measurement layer unmaterialized, so **every** free parameter living in an observable formula or a
> noise scale reads back as exactly `0.0` — which looks precisely like a product bug and is not one.

> **Gotcha: sweep the step size before believing a red row.** On a stiff model `h = 1e-7` is pure
> roundoff. The default is `3e-4`, which gave the best agreement on Laske. A real defect does **not**
> move with `h`; FD noise does. `k_syn_R_M` read `0` at every step size — that is what a real one
> looks like.
>
> Sweep **upward** as well. On these models the small-`h` end is the roundoff end: a column whose FD
> drifts monotonically *away* from the assembled value as `h` shrinks is converging at the other end.
> Fiedler's `k2` reads 216.07 / 212.04 / 210.02 at `h` = 1e-4 / 3e-4 / 1e-3 against an assembled
> 210.95, and lands on it at 1e-2.

> **Gotcha: a gradient check is only a test where the gradient is large.** Twelve of the eighteen
> `gntr` slugs have `OG_nominal` ≈ 0 — their PEtab `nominalValue` point *is* the optimum, so the true
> gradient there is zero and "assembled vs central difference" compares two quantities that are both
> noise. Every one of those slugs flags red at the nominal point, and every one of those flags is an
> artefact; the slugs that came back clean were exactly the ones far from the optimum (Blasi at
> OG 448, Zhao at 276, Laske at 40). Read `nominal_check.json`'s `OG_nominal` first, and for a slug
> at its optimum re-evaluate somewhere the gradient has magnitude — displacing a few percent of each
> box width along sampling space is enough. This is what separated `Fiedler_BMCSystBiol2016`'s real
> defect (lanl/PyBNF#535) from the eleven artefacts around it.

> **Gotcha: an experiment must travel with its `data_key`.** `gradient_at` appends
> `(sim, exp, routing, suffix)`; this script used to append the first three only. The visible
> consequence is mild — a column under an analytic per-series `scale` (ADR-0066, lanl/PyBNF#533) is
> *refused* rather than assembled, which is loud. The invisible one is not: `_iter_scored_points`
> repoints `objective._scale_factors` at the current experiment **only when `data_key` is supplied**,
> so without it `_scale_factors` keeps whatever the last `evaluate()` inside `evaluate_multiple` left
> there, and every experiment's residual is scored with the last experiment's profiled `c*`. Nothing
> errors, and the ratio between two experiments' `c*` presents as a clean multiplicative factor on the
> affected columns — indistinguishable by eye from a product bug. Only bites a slug that declares
> `normalization … = scale`. Fixed; found while auditing this script for lanl/PyBNF#537.

> **Gotcha: both sides must include the constraint penalty.** `loss_at` calls
> `evaluate_multiple(..., config.constraints)`, which adds `cset.total_penalty(...)` — so the finite
> difference has always differentiated the penalty. The assembled side called
> `assemble_gaussian_gradient` alone and never its sibling `assemble_constraint_gradient`, which is
> the other half of what `gradient_at` does. On a slug with active constraints every constraint-touched
> column therefore read red, by construction. Only bites a slug with a `.con` / `.prop` file — no PEtab
> import in this corpus has one, which is why it went unnoticed. Fixed; found alongside the above.

> **Gotcha: keep every parameter clear of its box bound.** `FreeParameter.set_value` *clamps* an
> out-of-box value, so a parameter sitting on a bound has `f(u+h) == f(u)`: its central difference is
> a half-step, or exactly zero when both sides clamp. That reads as a dead gradient column and is
> purely an artefact of where the point was chosen. Four of the five remaining red rows in the #535
> sweep were this, and nothing else. Clip the evaluation point to `[lo + 8h, hi - 8h]` and say so
> when a parameter had to be moved.

> **Gotcha: that same `8h` inset makes an `h` ladder run at a *different point* on every rung —
> whenever `--disp` pushes a bound-sitting parameter out of the box.** The clip is
> `[lo + 8h, hi - 8h]`, so its position depends on `h`; a parameter whose displaced value lands
> outside the box comes back to a bound offset that *moves with the step size*. Every other column
> then shifts too, because the evaluation point did. On `Brannmark_JBC2010` four of 22 parameters
> (`k1d`, `k1e`, `k1f`, `k_IRP_1Step`) sit **exactly on** a bound at the PEtab nominal point, so
> `--disp 0.03` sent all four out and the clip returned them to `bound ∓ 8h`. The symptom is an
> *assembled* gradient that moves with `h` — `k_IRSiP_2Step` read 904 / 924 / 995 / 1202 across
> `h` = 1e-4 … 3e-3 — which `Weber_BMC2015`'s `VALIDATION.md` correctly names as the signature of a
> clamped point, and which reads as a 33% tolerance-induced gradient error if you do not check.
> **The fix is to pin the point rather than to re-derive it per rung:** compute the displaced vector
> once, inset it by a fixed fraction of each box width (0.05 is comfortably clear of `8h` for any
> `h ≤ 3e-3`), write it to a `param-values.json`, and run the ladder with `--disp 0`. The assembled
> gradient then reads *identical* across the whole ladder, which is what makes "the assembled side
> holds while the central difference wanders" a statement about FD noise rather than about clipping.

`param-values.json` is a `{parameterId: value}` map; without it the script uses each parameter's
median-quantile value, which may be a meaningless point for the model. To evaluate at the PEtab
nominal point:

```python
import csv, json
rows = csv.DictReader(open('<upstream>/parameters_<slug>.tsv'), delimiter='\t')
json.dump({r['parameterId']: float(r['nominalValue']) for r in rows
           if str(r.get('estimate', '1')).strip() == '1'}, open('nominal.json', 'w'))
```

## `sigma_profile.py` — how far is the nominal point once its estimated σ are free?

```bash
python tools/sigma_profile.py <upstream-root> <slug-dir>...
```

`nominal_check.json`'s `OG_nominal` evaluates the objective at the PEtab `nominalValue` vector —
**every** parameter, the estimated noise parameters included. Many slugs here ship placeholder σ
nominals (`Giordano` 1, `Zhao` 1000), and for those the number is dominated by `Σ nⱼ log σⱼ` sitting
far from its MLE. It is then not a statement about the dynamics, which is how it gets read: issue #38
orders the remaining ⚪ candidates "roughly by nominal-point distance, i.e. plausibly by difficulty".

This holds every non-noise parameter at nominal, sets each estimated σ to `√(Σⱼ r²/nⱼ)`, and
re-scores. No PyBNF and no simulation — the trajectory is upstream's own `simulatedData`, the same
oracle §2c uses, so this inherits exactly that coverage. Measured on the corpus: `Giordano` 3776 →
**743.72**, `Zhao` 276.12 → **135.75**, `Brannmark` unchanged.

> **The self-check is the load-bearing part, so it is printed rather than thresholded.** Substituting
> the nominal σ back in must reproduce `nominal_check.json`'s `J_paper`. That residual is reported in
> the last column, and a profiled number means something only when it is orders below the inflation
> it claims. `Giordano` reads `5.7e-14` against 3032. `Brannmark` reads `2.9e-04` against an inflation
> of `2.9e-04` — the same order, i.e. the tool saying it cannot resolve an effect that small, and the
> right conclusion is that Brannmark has no inflation rather than that it has one of `2.9e-04`.
> `Laske` reads `4.9e+02` against `5.5e+02`: its `simulatedData` is not the nominal-point trajectory
> at all, so every residual is against the wrong point. A threshold would have turned all three into
> the same uninformative red.

> **Gotcha: check against `J_paper`, never against `reduced_objective`.** PyBNF's reduced objective
> drops only the *parameter-independent* per-point constants; `Σ nⱼ log σⱼ` depends on a fitted σ and
> therefore stays inside it. The first version of this script compared to `reduced_objective` and
> reported a confident failure on `Brannmark`, `Laske` **and** `Zhao` — three slugs, one wrong
> baseline. `J_paper == -log_likelihood` is the unambiguous scale and the one `score.py` uses.

> **Gotcha: a one-to-one join is not optional.** `measurementData` and `simulatedData` must join on
> the full key (observable, condition, pre-equilibration condition, time, observable/noise
> parameters, dataset). Where they do not, the script refuses rather than guessing — a duplicated
> join silently multiplies the residual sum, which would read as a real distance. Normalizing the key
> columns matters too: the two tables disagree on dtype wherever a cell is sometimes blank, and on
> formatting between `1` and `1.0`.

## `box_probe.py` — how many box-sampled starts integrate, and which solve fails

```bash
python tools/box_probe.py <slug-dir> [-n N] [--sens] [--atol A] [--seed S]
```

Draws `N` points from the conf's own sampling distribution and simulates each, reporting how many
survive. **Run it twice, with and without `--sens`.** The difference isolates the forward-sensitivity
solve from the state solve, and on two slugs that difference *was* the answer:

| slug | plain forward solve | with the gradient path's sensitivity request |
|---|---|---|
| `Weber_BMC2015` | 7 / 11, 0.6 s | **2 / 11, 80 s** |
| `Brannmark_JBC2010` | 30 / 30, 1.0 s | **19 / 30, 210 s** |

Both had been recorded as lanl/bngsim#196 — "a scalar `atol` cannot serve a model spanning ten decades".
Both are really the *sensitivity* solve, which CVODES scales from the state tolerances. `fd_check.py`
always applies the sensitivity request, so a probe run through it cannot separate the two and reads as
an unconditional failure.

> **Gotcha: `--atol` writes a temporary conf rather than poking the config.** The tolerance keys are read
> when `load_config` *builds* the model, so assigning `config.config['sbml_atol']` afterwards is a silent
> no-op — and the symptom is a sweep whose rows are all identical, which reads as "the tolerance changes
> nothing". The header line prints `_config_atol` so you can see it moved.

> **Gotcha: an explicit `--atol` turns the ADR-0105 per-species vector OFF.** `sbml_atol` is that
> mechanism's documented off-switch, so passing the same scalar the derivation would produce is not a
> null comparison — it is vector-versus-scalar. On `Brannmark` both give 19/30, but the vector is 1.4×
> faster. Compare against the no-`--atol` run to price the vector.

> **Gotcha: `--sens` is the number that predicts start mortality**, because every gradient `job_type`
> applies the union sensitivity request on each evaluation. It has held: 22/30 here forecast the 24 dead
> starts out of 100 in Weber's real run.

> **Gotcha: a dead point is not always a tolerance problem.** Points that fail both with and without
> `--sens`, and fail fast, are bad parameter points. Weber has two such in eleven and they stay dead at
> every tolerance.

## `linear_scope.py` — which observable-layer parameters could a profile actually remove

```bash
python tools/linear_scope.py <slug-dir | corpus-root>... [--json out.json]
```

Static classification, no simulation. For every declared free parameter that reaches an
observable formula — directly, or bound per data row through a measurement-params table — it
reports whether the formula is affine in it, whether it is a *pure* multiplicative factor, what
space its observable's noise family scores its residual in, how many series it is tied across,
and whether any noise source also reads it.

Written for lanl/PyBNF#572, which proposed profiling an observable's `scale`/`offset` out of the
search by variable projection. The distinction the script exists to draw is the one that decided
that issue: variable projection is a **linear least-squares** identity, so it needs the family's
residual to be `d - y` with `y` affine in the coefficients. On a `lognormal` / `lnnormal` family
the residual is `log d - log(a*s + b)`, which is affine in `log a` only when `b == 0`. The output
column is therefore `closed_form ∈ {linear-lsq, log-geomean, none}`, not a yes/no.

Corpus-wide result (all 23 slugs): **47** affine observable-layer parameters in 9 slugs —
19 `linear-lsq`, 1 `log-geomean`, 27 `none`. Every coupled `(scale, offset)` pair in the corpus
is on a log family and has no closed form; of 8 offset-role parameters exactly one
(`Laske`'s `Int_nuc_off`) does.

> **Gotcha: PEtab's parser builds symbols with assumptions, so `sympy.Symbol(name)` is not the
> symbol in the expression.** `sp.Symbol('scale') in expr.free_symbols` is `False` for a `scale`
> that is plainly there, `diff` against the bare symbol is identically `0`, and **every affinity
> test then passes vacuously** — the script reports the entire corpus as linear, Michaelis-Menten
> denominators included, and looks perfectly plausible doing it. Resolve symbols by *name* against
> `expr.free_symbols`. This is `_symbol()`, and it is the single most dangerous line in the file.

> **Gotcha: a row-varying token is a parameter too, and the placeholders are not interchangeable.**
> `Brannmark`'s `k_IRSiP_1Step`/`k_IRSiP_2Step` never appear in an `observable:` line; they are
> bound per row to `observableParameter1_IRS1_P`. And `Schwen`'s `observable_Insulin` is
> `op1 + op2*g(·; op3, op4)` — two placeholders affine, two inside a Michaelis-Menten denominator.
> Testing "are all this observable's placeholders affine" rejects the whole observable and
> silently drops the two parameters that are exactly the partial separability #572 is about.
> Test per placeholder and keep the binding keyed by it.

> **Gotcha: a `sos` slug has no sigma, and `_spec_for` will hand you a Gaussian anyway.**
> `Smith_BMCSystBiol2013` is `objective = sos`; the family it reports is a structural default and
> only its `ln_base` means anything. Smith's residual *is* `sim - exp` on the linear scale, so the
> classification stands — but the sigma-weighting question does not arise there, and reading the
> reported family as a noise model would be wrong.

> **Gotcha: "affine in each" is not "affine in the pair".** `Schwen`'s `scale*(IR1 + IR1in + offset)`
> is affine in `scale` and affine in `offset` and quadratic in the two together. The joint test
> (all second partials, cross terms included) is what separates a genuine two-column span from a
> reparametrization, and the `reparam` flag records the difference.

## `linear_profile.py` — what the landscape looks like once those parameters are profiled

```bash
python tools/linear_profile.py <slug-dir> [--params a,b] [-n N] [--noise-profiling]
                               [--point LABEL=vals.json ...] [--json out.json]
```

Simulates once per point, then minimizes **PyBNF's own `evaluate_multiple`** over the named
parameters with θ held fixed, and prints the searched score, the profiled score, and the
no-dynamics (`flat`) reference side by side. It runs **no fit** — a `Borghans` point is one
simulation plus a few hundred re-scorings of it.

The profile is numeric on purpose. #572's closed form does not exist on four of the six slugs it
names, and a numeric profile is the exact conditional optimum whatever the family is, so the
landscape question can be answered separately from the closed-form question. ADR-0108 pinned its
own sigma closed form against a numeric minimization of the reported objective for the same
reason: re-deriving the algebra and measuring *that* confirms the algebra, not the landscape.

Measured on `Borghans_BiophysChem1997` (76 box draws, `--noise-profiling`): the searched
objective spans `[-160.9, +339.7]`; the profiled objective spans `[-167.15, -165.98]`, with
**46 of 76 draws landing exactly on** the no-dynamics score `-165.982113` and 39 of 76 profiling
the scale below `1e-6` — i.e. discarding the dynamics. See ADR-0123 in PyBNF.

> **Gotcha: a Result can only be scored once.** The measurement layer materializes its observable
> columns onto the `Data` in place (ADR-0036) and `_add_column` raises rather than overwriting.
> The inner minimization scores one Result hundreds of times, so every scoring takes a `deepcopy`
> of the simulated data. This one fails loudly, which is the good case.

> **Gotcha: `FreeParameter.set_value` clamps to the box, and leaving the box is the whole point.**
> Carrying a trial value in a `PSet` silently returns the bound instead, and the measured
> "profiled" landscape is then a bounded-search landscape wearing the wrong name. The inner loop
> hands `evaluate_multiple` a duck-typed `(name, value)` sequence — all it reads is
> `{p.name: p.value for p in pset}` — and never a `PSet`.

> **Gotcha: the flat reference is only a reference line when sigma is profiled.** With a searched
> free sigma it is evaluated at whatever sigma the current draw carries and moves by *hundreds* of
> objective units from row to row — a different number each row under the same name. The tool
> recomputes it at every point and prints the spread for exactly this reason; a spread that is not
> ~0 means the run needed `--noise-profiling` and every comparison against that column is void.

> **Gotcha: pin the non-intercept coefficients to exactly `0`, not to their lower bound.** A
> `loguniform_var` scale bottoms out at `1e-3`, and on a trajectory spanning ten decades
> `1e-3 * s` is not remotely constant. The first version used the bound and produced a "flat line"
> that varied with the draw, which reads as sigma contamination and is not.

> **Gotcha: for one or two coefficients, grid before you polish.** The profiled surface of a
> log-family offset is not convex in `log10(offset)`, and a bare local solver started at the drawn
> value converges to whichever side of the data it started on — tens of objective units, larger
> than the whole effect being measured.

> **Gotcha: `noise_profiling = 1` is all-or-nothing and several slugs refuse it.** `Brannmark`
> refuses (`IRS1_P` takes sigma from a `PerMeasurementFormulaSigma`); `Weber` and `Borghans`
> accept. That is not a limitation of this script, it is ADR-0108's rule, and it is load-bearing
> for the comparison: on `Weber`, switching sigma profiling *on* takes the observable-scale
> profile's benefit from thirteen orders of magnitude to nothing.

> **Gotcha: `$VAR` holding several arguments does not word-split in zsh.** Building a
> `--point a --point b` list in a shell variable and passing it unquoted sends argparse one giant
> argument and gets `unrecognized arguments`. Write the flags out, or use an array.

> **The self-check that makes the numbers trustworthy: the profiled score can never exceed the
> searched score**, because the searched value is in the set being minimized over. It held at
> every one of the ~250 points measured. Check it before believing any row.

### `--closed-form`: checking #572's variable projection against the numeric profile

Builds `Phi` by **evaluating the model's own prediction** at basis coefficient vectors — with every
profiled coefficient at `0` the aligned prediction is `Phi.0`, and setting coefficient `j` to `1`
gives `Phi.e_j` — then solves `c* = (Phi^T W Phi)^-1 Phi^T W d` and scores it. Formula-agnostic, no
new sensitivity machinery, and the `pred(0) != 0` case reports itself rather than being absorbed
into the intercept.

Measured on `../../Synthetic-2026-linear-observable` (80 box draws): the closed form and the numeric
profile agree to **1.5e-14 relative** at the 53 points where both can reach the same optimum. At the
other 27 the projection returns a **negative `scale`** and scores a median 14.5 objective units
better — because it is unconstrained by construction and does not know `scale` was declared
`loguniform`. That is the finding, not a discrepancy.

> **Gotcha: run it with sigma SEARCHED, never with `--noise-profiling`.** `W = diag(1/sigma^2)` has
> to be the same matrix on both sides of the comparison. Under noise profiling sigma moves with
> every evaluation, so the closed form solves a different weighted problem than the one the numeric
> profile converges to, and the two disagree for a reason that is an error in neither.

> **Gotcha: per-point fit weights are not on the `aligned_prediction_data` seam.** A job that
> declares none (weight 1 everywhere) is unaffected; a job with weights would need them folded into
> `W` before this comparison means anything.

> **Gotcha: `aligned_prediction_data` returns `None` unless every scored point is a linear-scale
> Gaussian.** That is not a limitation to work around — it is exactly #572's precondition, so the
> inapplicable case reports itself instead of producing a number.
