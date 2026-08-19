# Weber_BMC2015

**Run cost: `minutes`** — 100 starts × 1000 iterations (`gntr`), 36 free parameters; **21 m 25 s
measured** on ten cores.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**✅ SOLVED — `OG = 0.781167`** (threshold 1.92), from a from-scratch multi-start `gntr` run
(100 starts × 1000 iterations, `random_seed = 1`, no seeding) in **21 m 25 s**. **75 of the 100 starts
converged on `step is negligible`**, 24 died with `start point failed to simulate`, and 1 stopped on a
LAPACK factorization failure in the trust-region subproblem; **none hit `max_iterations`**, so the run
was not budget-limited. See `VALIDATION.md`.

**Solved but not saturated**, in the `Fiedler_BMCSystBiol2016` sense. The PEtab `nominalValue` point
*is* the published optimum (`OG_nominal = -0.0002`, reproduced to four decimals), and this fit lands
0.78 NLL units above it. That is structural rather than a budget shortfall: **five of the 36
parameters sit exactly on their box bounds at the published optimum** (`a22`/`a32`/`a33` at `1e-4`,
`m11` at `1e10`, `pu3` at `1e8`), so the reference basin is a *corner* of the box. The fit recovers
one of those five (`m11`, to within 0.005% of its bound) and leaves the rest interior. It is **not** a
leverage trade — the five fitted σ come back within 0.7% of their published values, three observables
reproduce to four decimals, and the fit is *better* than the published point on three others.

> **Corrected 2026-08-19 — the tolerance was never what made this slug fittable.** The account
> below stood for nine days and is wrong; it is kept because the wrong version is the memorable one.
>
> Until 2026-08-11 this slug would not fit at PyBNF's derived tolerance: 6 of 30 sensitivity-applied
> box points integrated against 22 at `sbml_atol = 1e-4`, ~80% of starts died, and the conclusion
> drawn — here, in `VALIDATION.md`, and in [lanl/PyBNF#557](https://github.com/lanl/PyBNF/issues/557)
> — was that ADR-0103's only-ever-tighten clamp had made the model unfittable.
>
> **It was [lanl/bngsim#305](https://github.com/lanl/bngsim/issues/305)**: a registered
> time-discontinuity root the integrator could never *reach*, so the run wedged one ulp below
> `t = 24` — this model's `PdBu_time` / `kb_NB142_70_time` crossing. Closed 2026-08-11; the fix's own
> changelog measures this slug and reports its step count roughly halving. The bngsim build every
> measurement above was taken on (`114d3b3`, Aug 10 16:55) is a **provable git ancestor** of that fix
> (`7b9140b`, Aug 11 12:27).
>
> Re-measured 2026-08-19 on bngsim 0.13.0, same probe, same 30 box points, same sensitivity request:
> **every tolerance arm integrates 30/30** within 7.2 s of each other — unset, `1e-08`, `1e-04`, and
> PyBNF's new `auto` / `tracking`.
>
> **What survives.** The clamp #557 describes is real and unchanged: this model's own scale asks for
> `4.665e-03` and ADR-0103 hands it `1e-08`, with ADR-0105's vector clamping into `[1e-8, 1e-8]` and
> correctly declining. #557 shipped the opt-in for it ([ADR-0114](https://github.com/lanl/PyBNF/blob/main/docs/adr/0114-an-only-ever-tighten-clamp-is-a-no-regression-rule-rather-than-a-property-of-the-model-so-sbml-atol-gains-an-opt-in-that-lets-the-derivation-loosen-and-one-that-follows-the-trajectory.md)),
> so `sbml_atol = auto` now reaches that answer. **This slug deliberately does not use it** — see the
> conf comment: `auto` costs ~7x the finite-difference disagreement at the default step size, because
> what a looser tolerance costs is objective *smoothness*, which is what a trust-region line search
> consumes. `1e-4` earns its place on integrator cost (78641 CVODE steps → 36154 over 20 box points)
> and gradient agreement, not on integrability.

> **Corrected 2026-08-07** by [lanl/PyBNF#547](https://github.com/lanl/PyBNF/issues/547) (ADR-0104). This problem pre-equilibrates, and the bngsim SBML backend silently dropped `preequilibrate:`. Its pre-equilibration condition is all-zeros and matches the model's authored defaults, so both experiments simulated a *completely flat* trajectory — identical to 8 significant figures at every timepoint, including across `t = 24` where `PdBu_time = 24, PdBu_dose = 1` should fire — and this README recorded `OG = 13739.87` as "the nominal point is not the published optimum". It is the optimum; the forward model was wrong. **This slug was previously queued as a tuning candidate on that reading; it is not one.** A trajectory that never moves through an event that should fire is the signature to check for if this ever recurs.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `296.2020024574204` |
| **paper-scale NLL at this fit's best point** | **`296.983169`** |
| **optimality gap from the fit** | **`0.781167`** |
| paper-scale NLL at the PEtab nominal point | `296.20179502327426` |
| optimality gap at nominal | `-0.00020743414614798894` |
| scored data points `n` | 135 |
| free parameters `k` | 36 |

The two nominal-point rows were regenerated 2026-08-10 under the conf's explicit `sbml_atol = 1e-4`;
they previously read `296.20179656464325` / `-0.00020589277715998833`. Four decades of tolerance move
the objective in the 6th decimal — see `VALIDATION.md` Gate A for the full sweep.

For this slug `J_paper - reduced_objective = 365.517474` exactly, at both the nominal point and the
fit's best point, so **`OG = reduced_objective + 69.31548`** and "solved" means `reduced < -67.396`.

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through trf's
Coleman–Li core, ADR-0068) — handles this problem's estimated noise scale, which plain `trf` refuses.
The shipped recipe is the collection's 100 × 1000 working default (raised from the 20 × 500 import
placeholder on 2026-08-10) and has been run to completion: 75 of its 100 starts converged on
`step is negligible` and none hit `max_iterations`.

The conf also carries an explicit **`sbml_atol = 1e-4`**. It is no longer what makes this problem
fittable — lanl/bngsim#305 was, and it is fixed (see the correction above) — but it is still the right
value, and the conf comment carries the measurements: half the integrator steps of the derived
tolerance, and the best finite-difference gradient agreement of any arm at the default step size.
Do not replace it with `sbml_atol = auto` without re-reading that comment; `auto` is measurably worse
here on exactly the quantity a trust-region line search consumes.

## Contents

- `Weber_BMC2015.conf` — the PyBNF job
- `model_Weber_BMC2015.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment__model1_data1___model1_data2.exp`, `experiment__model1_data1___model1_data2_rep2.exp`, `experiment__model1_data1___model1_data2_rep3.exp`, `experiment__model1_data1___model1_data2_rep4.exp`, `experiment__model1_data1___model1_data2_rep5.exp`, `experiment__model1_data1___model1_data3.exp` … — experimental data
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `best_fit_params.txt` — the fit's sorted final parameter sets (`output/Results/sorted_params_final.txt`)
- `information_criteria.txt` — `log_likelihood` / AIC / BIC / AICc at the best fit
- `VALIDATION.md` — the four-gate validation of this result
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
pybnf -c Weber_BMC2015.conf -o
python score.py output
```
