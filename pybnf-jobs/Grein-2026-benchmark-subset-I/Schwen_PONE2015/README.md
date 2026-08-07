# Schwen_PONE2015

**Run cost: `hours`** — 100,000 evaluations (100 × 1000 `gntr`), 30 free parameters.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

> Schwen LO, Schenk A, Kreutz C, Timmer J, Bartolomé Rodríguez MM, Kuepfer L, Preusser T.
> **"Representative Sinusoids for Hepatic Four-Scale Pharmacokinetics Simulations."**
> *PLOS ONE* **10**(7): e0133653 (2015). <https://doi.org/10.1371/journal.pone.0133653>

The PEtab problem is the paper's **third** proof-of-concept application — cellular insulin
binding, internalization and extraction in two entities of mouse hepatocytes (*low-binding* and
*high-binding*), fitted to flow-cytometry (FACS) and ELISA data. It is the model of **Fig 11** and
Eq 21/22, not the full four-scale framework the paper's title describes.

---

> ## ⚠️ READ FIRST — this problem's σ encoding is not the paper's, and it changes what a good `OG` means
>
> **The paper estimated an error for every data point. The PEtab problem collapses that to two free
> σ.** Fig 11's caption says the shaded bands "correspond to the **estimated error in the data
> points**"; the PEtab measurement table binds every row to one of exactly two noise parameters:
>
> | free σ | covers | points |
> |---|---|---:|
> | `IR_obs_std` | `observable_IR1` + `observable_IR2` + `observable_IRsum` (all FACS) | **34** |
> | `std` | `observable_Insulin` (all ELISA) | **252** |
>
> Both are estimated (`estimate = 1`, `parameterScale = log10`), so with σ profiled out each group
> contributes `nⱼ·log σⱼ + nⱼ/2` and **`∂NLL/∂log(RMSⱼ) = nⱼ`**. The ELISA assay therefore carries
> **≈7.4× the leverage** (252 vs 34) of the FACS panel that Fig 11 publishes as *the* fit.
>
> **The consequence, measured at the final best point.** The published parameter vector is **not**
> the MLE of this objective. The completed fit reaches `−315.9897` against the published point's
> `−311.8668` — a better objective — while `observable_IR2` comes out flat or declining at all
> three doses, where the data rise 2–3.5× and the published point tracks them closely:
>
> | dose | `IR2` at t=1 → t=30 | data | published | ours |
> |---:|---|---:|---:|---:|
> | 100 | | 4.24 → 9.25 | 3.50 → 7.74 | 5.18 → 4.83 |
> | 1000 | | 13.57 → 47.84 | 16.21 → 48.72 | 31.59 → 29.06 |
> | 10000 | | 61.28 → 136.19 | 73.99 → 112.79 | 92.23 → 92.04 |
>
> **This is not a defect and not ours to fix.** The objective reproduces the independent oracle
> exactly (943.9993), the model is properly dose-responsive, and the encoding is the
> Benchmarking-Initiative's. `J*` is defined on the *same* PEtab objective, so the benchmark
> comparison is apples-to-apples and a `✅` here is sound.
>
> **The fitted σ show the trade directly.** `IR_obs_std` (FACS, 34 points) ran **into its upper
> bound** at 0.056234 — the fit wanted to call the FACS misfit measurement noise and the box
> stopped it — while `std` (ELISA, 252 points) *tightened* from 0.248324 to 0.186594, nowhere near
> either of its own bounds.
>
> **What it does mean:** a good `OG` on this slug is evidence about the **PEtab objective**, not
> evidence that the paper's Fig 11 kinetics were reproduced. Any write-up must say which. Combined
> with `J*` being unconverged here (below), Schwen is the collection's weakest link between "solved"
> and "recovers the published fit".
>
> **What is *not* a discrepancy:** that 88% of the points are ELISA. The paper fitted both assays —
> *"A selection of the FACS data is shown Fig 11. In addition, the dose- and time dependency of
> insulin depletion in the medium was evaluated by an ELISA."* Fig 11 plots a **selection** of one
> assay, not the objective. Our 19 conditions match that sentence exactly: FACS dose-response at
> 15 min (`data1`–`data4`), FACS time courses at three doses (`data5`–`data7`), ELISA dose and time
> dependency (`data8`–`data19`).

---

## The slug is named 2015; upstream says 2014

This directory was `Schwen_PONE2014` until 2026-08-07. **The paper is 2015** — received 2015-04-17,
published 2015-07-29, and `journal.pone.0133653` is a 2015 DOI. The model's own SBML `<notes>` say
so too: *"PEtab implementation of the model from Schwen et al. (2015), PLoS One, 10, e0133653"*.

**Only the local slug was renamed.** Everything that is upstream's to name still says 2014, because
changing it would break the join to the collection this corpus exists to be compared against:

| | name | why |
|---|---|---|
| this directory, `*.conf` | `Schwen_PONE2015` | ours to name; the year should be right |
| `model_Schwen_PONE2014.xml` | **unchanged** | verbatim upstream file, pinned by sha256 in `upstream.json` |
| `upstream.json` → `upstream_path` | **unchanged** | `Benchmark-Models/Schwen_PONE2014/...` is upstream's path |
| Grein et al.'s own tables, `best_fx_marvin.csv` | **2014** | the leaderboard key; `jstar.txt` comes from it |

`nominal_check.json` and `upstream.json` both record `upstream_slug: Schwen_PONE2014`, and
`tools/sigma_profile.py` reads that field rather than guessing from the directory name — the local
slug is deliberately **not** a reliable join key for this one problem.

## Status

**✅ SOLVED — `OG = −12.545379`**, 100 × 1000 unbiased starts, ~2 h, 64 of 100 starts
converging on a negligible step. See `VALIDATION.md`, and read the callout above before quoting
that number: the reference it beats is unconverged, and the fit does not reproduce Fig 11.

Earlier revisions of this file carried two claims that are now known to be wrong; both are corrected
here rather than quietly dropped.

**Correction 1 — the log10 objective is not a caveat.** This file used to say the nominal gap was
"recorded but not asserted" because a log10 observable carries the change-of-variables Jacobian
`Σ log(y_obs·ln10)` that Eq. 6's `J*` need not, and that this kept the log10 slugs out of the
"objective validated" tier. That convention is now **settled and verified** by three solved slugs
whose Jacobian is large: `Perelson` (log10, +233.1287) at `OG = 5e−07`, `Blasi` (ln, −1102.0028) at
`−4.3e−07`, and `Laske` (mixed, +261.0897) at `−1e−06`. If Eq. 6 omitted the term, Perelson would
miss `J*` by 233 rather than by 5e−07. This slug carries **`obj ✓`**.

**Correction 2 — `J*` is not "tightly converged"; it is the outlier.** This file used to describe
`J* = 952.422` as tightly converged and the nominal point as landing "just below" it. In fact the
nominal point scores `−log_likelihood = 943.999`, i.e. **8.42 NLL units better than the reference**,
which is 4.4× the solved threshold *in the wrong direction*. Recomputing that NLL straight from
upstream's own `simulatedData` + `measurementData`, with no PyBNF in the loop, gives **943.9993 —
matching PyBNF exactly**. So our setup is right and the reference is not converged for this problem
(k=30, n=286, log10). **`J*` here is a weak reference, not a target.**

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `952.4217251306842` |
| paper-scale NLL at the PEtab nominal point | `943.9992974818219` |
| optimality gap at nominal | `−8.422427648862254` |
| scored data points `n` | 286 |
| free parameters `k` | 30 |

## The search scale was wrong until 2026-08-07

**19 of this problem's 25 log10 parameters were imported as linear `uniform_var`** and searched on a
linear box. `petab1to2` materializes PEtab v2's implicit `uniform` default into the converted
parameter table whenever the v1 table has a prior column at all — Schwen's is present and populated
for exactly 6 of 30 rows — and the re-injection then declined to overwrite what looked like a
declared prior. Fixed as [lanl/PyBNF#548](https://github.com/lanl/PyBNF/issues/548); the conf's
parameter block was regenerated from a clean re-import. Its six real `parameterScaleNormal` priors
correctly survive as `log-normal`, and its five genuinely `lin` parameters stay `uniform`.

**No fit before that date searched this problem correctly.** `Zhao_QuantBiol2020` was the other
affected slug, at 28 of 28.

## How the 286 points divide

| observable | points | share of points | governing σ |
|---|---:|---:|---|
| `observable_Insulin` (ELISA depletion) | 252 | 88.1% | `std` |
| `observable_IR1` (low-binding) | 15 | 5.2% | `IR_obs_std` |
| `observable_IR2` (high-binding) | 15 | 5.2% | `IR_obs_std` |
| `observable_IRsum` (weighted, 0.605/0.395) | 4 | 1.4% | `IR_obs_std` |

These are shares of **scored points**, not of the likelihood — the two coincide only because each
group carries its own free σ, which is the point made in the callout at the top of this file. Fig 11's
bottom panel is the `IR1`/`IR2` time courses at three doses: **30 of 286 points**, against 252 for the
assay the figure does not show.

The weights `0.605`/`0.395` in `observable_IRsum` are the measured low-binding cell fraction; the
paper reports `η_l,obs = 0.606`.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (ADR-0068), which handles
this problem's estimated noise scales and the log10 (`lognormal`) objective (EFIM Fisher block for
noise scales, ADR-0079/0080/0081). `population_size = 100`, `max_iterations = 1000` — the
collection's documented working default, not the shipped 20 × 500 placeholder.

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`.
This problem exercises **lanl/PyBNF#510**: one experimental condition is measured only
at `t = 0`, which `TimeCourse` previously rejected — failing the entire problem at
config load. It now loads, simulates, and scores. The observables also carry
offset+scale `observableParameters` and per-observable estimated noise scales
(`IR_obs_std`, `std`), imported via the `measurement_params:` sidecars. The run recipe
(`job_type`, `sbml_backend = bngsim`, `wall_time_sim`) is supplied, not recovered.

## Contents

- `Schwen_PONE2015.conf` — the PyBNF job
- `model_Schwen_PONE2014.xml` — SBML model, **verbatim upstream, deliberately not renamed**
- `experiment____model1_data*.exp` (+ `_rep2`) — experimental data, 19 conditions
- `experiment____model1_data*_measparams.tsv` — per-measurement observable/noise parameter tables
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `score.py` — scores a run against `J*`

## Running

```bash
pybnf -c Schwen_PONE2015.conf -o
python score.py output
```
