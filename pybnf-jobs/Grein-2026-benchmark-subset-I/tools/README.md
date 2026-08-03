# tools — regenerating and verifying a slug

Two scripts used to produce and check the numbers this corpus ships. Both are run from the
repository with PyBNF available (`uv run --project ~/Code/PyBNF python tools/<script>.py …`).

They exist because both were originally written ad hoc, and both have a non-obvious failure mode
that silently produces a *plausible wrong answer* rather than an error. Those are written down here.

## `nominal_check.py` — regenerate a slug's `nominal_check.json`

```bash
python tools/nominal_check.py <slug-dir> <upstream-petab-dir> [--write]
```

Evaluates PyBNF's objective at the PEtab `nominalValue` point and puts it on the paper's Eq. 6 scale
the same way a real run does — through `likelihood_information_criteria`, which is what writes
`information_criteria.txt` at the end of a fit. Without `--write` it prints the new values and the
old ones for comparison.

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

> **Gotcha: keep every parameter clear of its box bound.** `FreeParameter.set_value` *clamps* an
> out-of-box value, so a parameter sitting on a bound has `f(u+h) == f(u)`: its central difference is
> a half-step, or exactly zero when both sides clamp. That reads as a dead gradient column and is
> purely an artefact of where the point was chosen. Four of the five remaining red rows in the #535
> sweep were this, and nothing else. Clip the evaluation point to `[lo + 8h, hi - 8h]` and say so
> when a parameter had to be moved.

`param-values.json` is a `{parameterId: value}` map; without it the script uses each parameter's
median-quantile value, which may be a meaningless point for the model. To evaluate at the PEtab
nominal point:

```python
import csv, json
rows = csv.DictReader(open('<upstream>/parameters_<slug>.tsv'), delimiter='\t')
json.dump({r['parameterId']: float(r['nominalValue']) for r in rows
           if str(r.get('estimate', '1')).strip() == '1'}, open('nominal.json', 'w'))
```
