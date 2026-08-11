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

> **Unblocked 2026-08-10, and not the way this slug's history predicted.** Issue #38 and
> [lanl/bngsim#196](https://github.com/lanl/bngsim/issues/196) both recorded Weber as unable to
> integrate anywhere near its nominal point, at every displacement tried including zero, and #38's
> standing instruction was "do not spend a fitting budget here yet." **No upstream change was needed.**
> Two corrections came out of actually measuring it:
>
> 1. **The failure was in the forward-*sensitivity* solve, not the state solve.** At the shipped
>    tolerance the plain forward model integrates 7 of 11 box points in 0.6 s; with the `gntr`
>    sensitivity request applied it manages 2 of 11 in 80 s. CVODES derives its sensitivity tolerances
>    from the state ones, which is why loosening the state `atol` rescues it — and why bngsim#196's own
>    footnote that this "does not reproduce on a bare model" was right without knowing the reason: a
>    bare `Simulator.run` has no sensitivity axis to fail on.
> 2. **The fix is `sbml_atol = 1e-4`, and the reason it was unavailable is a clamp.** ADR-0103 derives
>    `atol = rtol × median(y₀)`, which for this model is `4.665e-3`; the derivation "only ever
>    tightens", so it discards its own answer for `1e-8`, **5.7 decades tighter**. ADR-0105's
>    per-species vector clamps into `[scalar_atol, default_atol]` = `[1e-8, 1e-8]` here, so it is
>    elementwise the scalar and correctly declines to engage. The conf carries the reasoning and the
>    measured sweep.
>
> **A second defect is still open and cost this run 24 of its 100 starts.** Loosening `atol` removed the
> `mxstep` failures and made a distinct one visible: 84 simulations died with
> `CVODE made no progress … the step size has collapsed`, every one at `t = 23.999999999999996` —
> exactly `PdBu_time = 24` / `kb_NB142_70_time = 24`, the discontinuity in this model's
> `assignmentRule u5 = piecewise(0, (time - PdBu_time) < 0, PdBu_dose)`, with `h ≈ eps·t`. **This is
> not `lanl/bngsim#194`** (closed, and about *state* thresholds) **and the root is not missing**: the
> loader registers `((time()-PdBu_time)<0)` and `((time()-kb_NB142_70_time)<0)` as discontinuity
> triggers. It is a registered time root the integrator cannot advance past — filed as
> [lanl/bngsim#305](https://github.com/lanl/bngsim/issues/305). See `VALIDATION.md` for what is and is
> not established. The tolerance clamp behind item 2 above is
> [lanl/PyBNF#557](https://github.com/lanl/PyBNF/issues/557).

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

The conf also carries an explicit **`sbml_atol = 1e-4`**, which is what makes this problem fittable at
all. Do not remove it without reading the comment above it: at the value PyBNF derives on its own,
80% of box-sampled starts fail to integrate under the sensitivity request the gradient path needs.

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
