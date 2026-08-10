# Smith_BMCSystBiol2013

**Run cost: `hours`** — 100,000 evaluations (100 × 1,000 `gntr`) over a 367-reaction network — the subset's most expensive product of budget and model size.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**Setup only — not fitted.** The job runs and scores correctly; the PEtab nominal point is not
this problem's published optimum, so no optimality claim is made.

> **Requires lanl/PyBNF#553** (merged 2026-08-10 as `85b36a96`). All 133 species here are declared
> `hasOnlySubstanceUnits="true"` and the compartments run from `8.3e-12` down to `1e-13`. Before
> that fix the bridge handed observable formulas each species as `amount / compartment_size`, so
> every scored point came out up to 1.6e13 too large and the nominal point evaluated to
> `6.85e+32`. **That figure appeared in this file and in `nominal_check.json` and was an artefact,
> not a distance from the optimum.** The tell was that `PTP_activ__2D` — the one observable that
> is a ratio within a single compartment, where the volume cancels — was accurate while every
> scaled-amount observable was off by 10–14 orders. Multiplying the old predictions by their
> compartment volume returns the measured data (`PI3K_activity` 515 vs 515, `IRSYp` 220 vs 220,
> `Glucose_uptake` 229 vs 230).
>
> The practical consequence was that the problem was **unfittable**: absorbing `1/V ≈ 1e13` needs
> the `sc_*` scale parameters to reach ~1e-15 against a box of `[1e-4, 1e4]`, and in a 100-start
> `gntr` run six of the nine pinned hard at the lower bound while the objective stalled at 1.5e+27.
> Smith is the only slug of the 23 with both preconditions, so no other slug's numbers moved.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `20922.16424399946` |
| paper-scale NLL at the PEtab nominal point | `888198.438390848` |
| optimality gap at nominal | `867276.2741468486` |
| scored data points `n` | 62 |
| free parameters `k` | 25 |

With `objective = sos` and σ fixed at 1, the constant restoring the reduced objective to the
paper's Eq. 6 scale is `(n/2)log(2π) = (62/2)log(2π) = 56.974189`, so **`OG < 1.92` means a
reduced objective at or below `20867.110055`**. (`nominal_check.json` computes that constant as
`0.0` only because it vanished against the old `6.85e+32`; the identity is confirmed on
`Oliveira_NatCommun2021`, the corpus's other unit-σ `sos` slug, whose restored constant is
exactly `(120/2)log(2π) = 110.272624`.)

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — 100 starts × 1000 iterations, the collection default. This slug was
`cmaes` until 2026-08-10 because **the gradient path used to refuse it**, and both halves of
that refusal are now gone.

**The pre-flight gate.** The model contains discrete events (state-dependent jumps), across
which forward output sensitivities went stale, so bngsim refused them and lanl/PyBNF#461
hoisted the refusal into a blanket gate (`_require_differentiable_dynamics`). bngsim applies
the event's own jump now, and PyBNF lifts the gate at `BNGSIM_HAS_EVENT_SENS` — a **version
floor of exactly 0.12.2**, not a capability probe. The refusal therefore ended at the 0.12.2
release; this README asserted it for longer than it was true.

**What the version floor cannot see.** Because it is a floor, any build reporting `0.12.2`
passes it — including one without bngsim#160/#161, which emits the analytic sensitivity RHS
for a cross-compartment reaction. Reaction 7 (`R5f`) here is `per_species_volume_scaling`;
without that commit the whole model declines the analytic `df/dp`, all 25 columns fall back to
CVODES difference quotients, and every start times out and returns `inf`. The gate says yes and
the fit still cannot run. Measured on a build carrying it: 0 decline warnings, 0
difference-quotient fallbacks, and the event condition rooted and jumped in flight
(`'PI345P3>pip3_basal'`).

**And a second, unrelated blocker.** Even with a usable gradient this problem could not be
fitted until lanl/PyBNF#553 — see Status. The gradient itself is sound: a finite-difference
check of the assembled gradient against central differences of the objective, at the nominal
point where the gradient has real magnitude, agrees to a worst relative error of **8.3e-05**
across all 25 parameters, with no structurally zero columns.

## Contents

- `Smith_BMCSystBiol2013.conf` — the PyBNF job
- `model_Smith_BMCSystBiol2013.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment____figure2B__0_0__1_3em10.exp`, `experiment____figure2B__0_0__1_5em08.exp`, `experiment____figure2B__0_0__1_5em09.exp`, `experiment____figure2B__0_0__1_5em10.exp`, `experiment____figure2B__0_0__1_5em11.exp`, `experiment____figure2B__0_0__1em06.exp` … — experimental data
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `score.py` — scores a run against `J*`

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`. The
converter preserves both `parameterScale` (lanl/PyBNF#491) and `observableTransformation`
(lanl/PyBNF#499), which plain `petab.v2.petab1to2` drops. The run recipe (`job_type`,
`sbml_backend = bngsim`, `wall_time_sim`) is supplied, not recovered — PEtab specifies a
problem, not a method. `wall_time_sim = 10` caps pathological parameter points; raise it
if valid simulations on your machine are being marked as failures.

## Running

```bash
pybnf -c Smith_BMCSystBiol2013.conf -o
python score.py output
```
