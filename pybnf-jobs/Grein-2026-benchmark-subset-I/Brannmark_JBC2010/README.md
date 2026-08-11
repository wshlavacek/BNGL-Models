# Brannmark_JBC2010

**Run cost: `hours`** — 100 starts × 1000 iterations (`gntr`), 22 free parameters; **1 h 41 m 10 s
measured** on ten cores.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**✅ SOLVED — `OG = 0.111206`** (threshold 1.92), from a from-scratch multi-start `gntr` run
(100 starts × 1000 iterations, `random_seed = 1`, no seeding) in **1 h 41 m 10 s**. **87 of the 100
starts converged on `step is negligible`**, 12 died with `start point failed to simulate`, and 1 hit
`max_iterations`. See `VALIDATION.md`.

**It recovers the published basin, not merely the threshold.** The PEtab `nominalValue` point *is*
the published optimum (`OG_nominal = 0.0643`), and **four of the 22 parameters sit exactly on their
box bounds there** — `k1d`/`k1f` at the upper `5e+05`, `k_IRP_1Step` at the upper `2e+05`, `k1e` at
the lower `1e-06`. That corner structure is what held `Weber_BMC2015` 0.78 NLL units short of its
reference basin; this fit lands **all four** (`k1d` exactly, `k1f` at 99.99% of its bound,
`k_IRP_1Step` at 99.996%, `k1e` 0.3% of a box width off) and pays only 0.047. It is **not** a leverage
trade: the four estimated σ come back within 14%, and the fit *beats* the published point on 3 of its
10 experiment/observable groups.

> **Unblocked 2026-08-11, and #38's diagnosis of it was wrong in a specific way.** This issue recorded
> the slug's failed run under [lanl/bngsim#196](https://github.com/lanl/bngsim/issues/196) — "a scalar
> `atol` cannot serve a model that seeds a transient at `1.8e-9` against principal species at
> `0.1..10`". **No upstream change was needed**, and a scalar serves it fine:
>
> 1. **The failure is in the forward-*sensitivity* solve, not the state solve.** At the derived
>    tolerance the plain forward model integrates **30 of 30** box points in 0.9 s; with the `gntr`
>    sensitivity request applied it manages 19 of 30 in 210 s. CVODES derives its sensitivity
>    tolerances from the state ones. At a single box start the split is unambiguous: state solve ok in
>    0.06 s, sensitivity solve dead on `mxstep`.
> 2. **The fix is `sbml_atol = 1e-8` — the backend default — and this slug is the *mirror* of Weber,
>    not another instance of it.** ADR-0103 derives `atol = rtol × median(y₀)`, which here is
>    `3.302e-10`; the clamp "only ever tightens", so unlike Weber the derivation is *allowed* to apply
>    its answer, and does — **1.5 decades below the backend default**. That tightening is what kills
>    the starts. ADR-0105's per-species vector clamps into `[3.302e-10, 1e-08]` and so cannot undo it
>    (19/30 either way; it buys speed, not survival). So Brannmark is **not** the worked example for
>    [lanl/PyBNF#557](https://github.com/lanl/PyBNF/issues/557)'s ask (b) that issue predicted — it is
>    the mirror of ask (a). The conf carries the reasoning and the measured sweep.
>
> **A stderr line count is not a budget measurement.** #38 cites "164,439 `mxstep` lines" as evidence
> of how badly the previous attempt was hurt. This *successful* run emitted **331,402** of them while
> losing 12 starts and 74 simulations — almost every one is a trial point a trust region rejects and
> recovers from. Use the log's `GNTR start N/100 stopping:` lines; `ls output/FailedSimLogs | wc -l`
> is also a proxy (13 files against 12 dead starts here).

> **Corrected 2026-08-07** by [lanl/PyBNF#547](https://github.com/lanl/PyBNF/issues/547) (ADR-0104). This problem pre-equilibrates, and the bngsim SBML backend silently dropped `preequilibrate:` — all eight doses simulated identically at the model's authored `insulin_dose_1 = 0.3`, and this README recorded `OG = 1531.44` as "the nominal point is not the published optimum". It is the optimum; the forward model was wrong. A dose-response whose doses all coincide is the signature to check for if this ever recurs.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `141.82485427243665` |
| **paper-scale NLL at this fit's best point** | **`141.936060`** |
| **optimality gap from the fit** | **`0.111206`** |
| paper-scale NLL at the PEtab nominal point | `141.8892031447126` |
| optimality gap at nominal | `0.06434887227595709` |
| scored data points `n` | 43 |
| free parameters `k` | 22 |

The two nominal-point rows were regenerated 2026-08-11 under the conf's explicit `sbml_atol = 1e-8`;
they previously read `141.88922297826417` / `0.06436870582751908`. Three decades of tolerance move the
objective in the 7th significant figure — see `VALIDATION.md` Gate A for the sweep.

For this slug `J_paper - reduced_objective = 39.514356927801` exactly, verified at four parameter
vectors spanning `OG` 0.06 to 37,373, so **`OG = reduced_objective - 102.310497`** and "solved" means
`reduced < 104.230497`. That shortcut is for watching a run; report `score.py`.

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through trf's
Coleman–Li core, ADR-0068) — handles this problem's estimated noise scale, which plain `trf` refuses.
`population_size = 100`, `max_iterations = 1000` — the collection's documented working default, run to
completion: 87 of its 100 starts converged on `step is negligible` and 1 hit `max_iterations`.

The conf also carries an explicit **`sbml_atol = 1e-8`**, which is what makes this problem fittable.
Do not remove it without reading the comment above it: it is the *backend default*, and at the value
PyBNF derives on its own — 1.5 decades tighter — only 63% of box-sampled starts integrate under the
sensitivity request the gradient path needs, against 93% here.

## Contents

- `Brannmark_JBC2010.conf` — the PyBNF job
- `model_Brannmark_JBC2010.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment__Dose_0___Dose_0.exp`, `experiment__Dose_0___Dose_001.exp`, `experiment__Dose_0___Dose_01.exp`, `experiment__Dose_0___Dose_1.exp`, `experiment__Dose_0___Dose_10.exp`, `experiment__Dose_0___Dose_100.exp` … — experimental data
- `experiment__Dose_0___Dose_100_measparams.tsv`, `experiment__Dose_0___TwoSteps_measparams.tsv` — per-measurement observable/noise parameter tables (ADR-0075)
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
pybnf -c Brannmark_JBC2010.conf -o
python score.py output
```
