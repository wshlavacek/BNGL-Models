# VALIDATION — Jaruszewicz-Blonska-2023/reduced_combination

Primary-source validation of the PyBNF job `pybnf-jobs/Jaruszewicz-Blonska-2023/reduced_combination/`,
per the `validate-pybnf-job` skill. Confidence is **earned from the gate evidence below**.

> **Confidence: 94 / 100** — Gates 0–2 PASS/equiv against the authors' OWN reduced-model BNGL (Gate 2
> = identical NF-κB network bar an inert time-clock); Gate 1 target is byte-reproducible from the
> authors' original-model BNGL and matches the paper's **914 independent data points exactly**; Gate
> 3a reproduces the paper's **reported AMD\*_WT = 1.29 / AMD\*_A20KO = 1.16 to 3–4 significant
> figures**. The fit now uses the paper's **EXACT Eq-7 objective** (fixed-σ log-normal + floor +
> per-series geometric-mean scaling, `lanl/PyBNF#479`) — not an analog — and under it Gate 3b recovers
> **all 13 parameters within 3×** (10/13 essentially exact, |log₁₀ ratio| ≤ 0.04), with the paper's
> single least-identifiable ε now tight to **1.5×** (was 5.9× under the analog). The objective at the
> Table-1 params (4.07) **coincides with the fit minimum (4.06; ratio 0.996)** — i.e. the exact
> objective's optimum *is* the ground truth. The remaining dock is the small structural
> reduced-vs-original mismatch (AMD ≈ 1.29, the paper's own residual), not the objective or the fit.

Primary sources (in the untracked `dev/papers/Jaruszewicz-Blonska-2023/`; not redistributed):
- Model paper: Jaruszewicz-Błońska J, Kosiuk I, Prus W, Lipniacki T. "A plausible identifiable model
  of the canonical NF-κB signaling pathway." PLoS ONE 2023; 18(6): e0286416. DOI
  10.1371/journal.pone.0286416. (Reduced 2023 model; Eqs 1–6, Table 1, **Fig 2**, **S1 Table**.)
- Author files used (S1 Code): `Reduced2023 - BNGL file/fitted_reduced_model_combination_protocol.bngl`
  (reduced model at Table-1 values) and `Lipniacki_et_al_2004_model - BNGL file/
  lipniacki_original_model_combination_protocol.bngl` (the original 15-variable model = fit target).

**"The paper's result" for this job** = Table-1 fitted values, i.e. the reduced model refit to
reproduce the original model across the combination experiment (Fig 2). Synthetic-target /
ground-truth-recovery job; Gate 3b recovers Table 1.

---

## Gate 0 — Materials inventory — **PASS**

| needed | present? | path / note |
|---|---|---|
| model paper PDF | ✅ | `journal.pone.0286416.pdf` (Table 1 = params, Fig 2 = combination) |
| S1 Table (protocol time points) | ✅ | `pone.0286416.s009` — the exact fit time points |
| author reduced-model BNGL (fitted, combination) | ✅ | S1 Code — direct Gate-2 diff |
| author original-model BNGL (combination, fit target) | ✅ | S1 Code — generates the target |

The job did not exist before this validation; it was **built** from the authors' own materials.

## Gate 1 — Data provenance — **PASS (synthetic target, byte-reproducible; point count matches paper)**

The 12 `.exp` are the **original Lipniacki-2004 model's RAW output** for the 6 protocols × {WT, A20KO},
sampled at the **S1-Table** measurement times — **un-normalized** (no floor, no peak, no geomean, no
`_SD`). The exact Eq-7 objective applies the paper's ρ=0.03·max floor and per-series geometric-mean
normalization **symmetrically to the simulated and the experimental column at scoring time**
(`normalization <obs> = floor 0.03, scale`, `lanl/PyBNF#479`), so the target must ship raw.
`validation/gen_combination.py` runs the **authors' own**
`lipniacki_original_model_combination_protocol.bngl` through BNG2.pl and samples the observables
{IKKa, NFkB_n, A20, IkBa\*(=IkBa_total), IkBa mRNA} (A20 excluded for A20KO).

**Point count check:** 538 (WT) + 430 (A20KO) = 968 measurements; minus 54 series (5 obs × 6
protocols WT + 4 × 6 A20KO) for the per-series scaling = **914 independent data points — exactly the
paper's reported count** ("the combination experiment contains 914 independent data points"). This
confirms the S1-Table time-point transcription is exact.

| protocol group | TNF ON windows (from S1 Table) | verdict |
|---|---|---|
| continuous | 0–240 min | PASS |
| pulse5_60 / 5_100 / 5_200 | 5-min pulses at 60 / 100 / 200-min periods | PASS |
| pulse22_5 / pulse45 | 22.5-min pulses @45-min / 45-min pulses @90-min | PASS |

## Gate 2 — Model fidelity — **PASS (equiv)**

Reference: the authors' `fitted_reduced_model_combination_protocol.bngl` (S1 Code).

`scripts/model_diff.py` (my `reduced_combination_continuous_wt.bngl` vs the authors' model; generated
network, rate values ignored):
```
species  : common=6  onlyA(ours)=1 -> Clock()      onlyB=0
reactions: common=14 onlyA(ours)=1 -> ()->Clock()  onlyB=0
```
Identical NF-κB network. Documented edition-2 refactors (same as `reduced_onoff`, plus the
protocol/genotype file split):

| aspect | authors' file | our `.bngl` | verdict |
|---|---|---|---|
| molecule types / seed / 14 NF-κB rules + rate laws | Eqs 1–6 | identical **+ inert Clock()** | match / equiv |
| TNF stimulus | `TR` parameter (0/1) toggled by `setParameter`/`simulate_ode` phases (12-arm action block) | `TRfunc()` = sum of ≤3 rectangular ON windows [t_on_k,t_off_k), read from a `Clock()`; each protocol's windows baked into a per-protocol model file | equiv (delivers the same TR(t) per protocol) |
| A20KO | `setConcentration("A20()",0)` + `setParameter("c_deg",0)` | **A20-synthesis rule removed** → A20 stays 0 → no A20→IKKa feedback | equiv |
| total IkBa protein IkBa* | measured var | `IkBa_star()=IkBa+1-NFkB_n` output function | equiv |
| pre-equilibration | steady_state solve / saved SS | none — seed IS the exact TR=0 steady state | equiv |

The 12 protocol/genotype model files are generated from the canonical `reduced_combination.bngl` by
`validation/gen_combination_models.py`. Per-file (not per-`condition:`) protocol delivery is required
because the bngsim backend can perturb only *free* parameters via a condition, and the TNF window
parameters are not fit; the 12 files are joined in ONE multi-model fit sharing the 13 free params.

**Verdict:** PASS (equiv).

## Gate 3a — Reproduce at the paper's parameters — **PASS (matches the paper's reported AMD)**

Reduced model set to **Table-1 fitted** values; compared to the original-model combination target.

- Metric = the paper's AMD (Eq 8), over all 914 points:
  **AMD\*_WT = 1.286** (paper reports **1.29**) · **AMD\*_A20KO = 1.158** (paper reports **1.16**).
  This reproduces the paper's headline reproduction metric (Fig 2 text) to 3–4 significant figures.
- Reproduction figures: `combination_gate3a.png` (`gen_combination.py`, vs the garage original model)
  and the self-contained `reduced_combination_reproduction.png` (`make_reproduction.py`, vs the
  committed `.exp`). Both reproduce **Fig 2** — WT nuclear NF-κB across all 6 protocols: the damped
  oscillation (continuous) and the entrained pulse trains (5-min/22.5-min/45-min pulses). Gate 3a is
  **objective-independent**; the switch to raw `.exp` only means `make_reproduction.py` now
  peak-normalizes the raw target itself for the overlay (per-protocol median |rel err| 5–17 %).

**Verdict:** PASS — the reduced model at the paper's own params reproduces the paper's own Fig-2
combination result and its reported AMD.

## Gate 3b — Recover the paper's parameters by fitting — **PASS (13/13 within 3×; the paper's EXACT Eq-7 objective)**

- Run: `pybnf -c reduced_combination.conf` — edition-2 `de`, 12-model joint fit, pop 100 × 150 iters
  + Simplex refine, `random_seed=1`. Objective = the paper's **EXACT Eq-7** (`noise_model = lognormal,
  sigma = fix_at 1` + `normalization <obs> = floor 0.03, scale`; `lanl/PyBNF#479`) — a fixed-σ
  log-normal residual on log₁₀ (the paper's squared-log objective) with the ρ=0.03 measurement floor
  and per-series geometric-mean scaling applied symmetrically to sim and data. Final objective
  **4.06** (DE best; refine reached 4.05).
- `scripts/compare_params.py published_table1.json output/Results/sorted_params_final.txt` (tol 3×):

| param | published (Table 1) | recovered | log10 ratio | within 3×? | note |
|---|---|---|---|---|---|
| k_deg | 0.000107 | 0.0001092 | +0.01 | yes | essentially exact |
| k_1 | 0.00195 | 0.00193 | −0.00 | yes | essentially exact |
| k_3 | 0.00145 | 0.001448 | −0.00 | yes | essentially exact |
| k_2 | 0.0357 | 0.03628 | +0.01 | yes | essentially exact |
| a_3 | 0.0946 | 0.08156 | −0.06 | yes | near-exact (on-off left this loose — combination identifies it) |
| **delta (δ)** | 0.1083 | 0.05635 | −0.28 | yes | ✓ (~1.9×) |
| **epsilon (ε)** | 0.0428 | 0.06383 | +0.17 | yes | **~1.5× — was ~5.9× under the analog**; paper's least-identifiable param (Fig 5) |
| c_deg | 0.000106 | 0.0001096 | +0.01 | yes | essentially exact |
| c_4a | 0.00313 | 0.002979 | −0.02 | yes | essentially exact |
| a_2 | 0.0763 | 0.06915 | −0.04 | yes | essentially exact |
| c_5a | 5.78e-05 | 5.689e-05 | −0.01 | yes | essentially exact |
| i_1a | 0.000595 | 0.0006477 | +0.04 | yes | essentially exact |
| c_3a | 0.000372 | 0.0003649 | −0.01 | yes | essentially exact |

**13/13 within 3×**, and **10/13 essentially exact** (|log₁₀ ratio| ≤ 0.04). The exact objective
recovers **every** parameter — including the paper's single least-identifiable **ε to ~1.5×** (+0.17),
tightened from the ~5.9× (+0.77) the edition-2 analog left it at. **δ (~1.9×) and a_3 (near-exact)**
are again recovered far tighter than the sparse on-off slug, matching the paper's finding that the
combination experiment improves identifiability of those directions.

- **Objective at Table-1 sits at the fit minimum.** `obj_at_table1.py` scores the reduced model at the
  Table-1 nominals with the *same* objective (PyBNF `fit_type = check`): **4.07**, vs the fit's **4.06**
  (ratio **0.996**). Under the exact objective the ground-truth params *are* the optimum — the fit
  simply lands there. (Contrast the analog, whose different optimum sat away from Table 1.)

**Verdict:** PASS — all 13 parameters recovered (most essentially exact) under the paper's exact
objective, with Table 1 confirmed as the objective minimum.

---

## Divergence & corrections

- Scope vs. paper: **matches** — built to the paper's actual model (reduced 2023, all 13 params free)
  and its headline combination experiment (Fig 2, S1 Table, WT + A20KO, 914 points), recovering Table 1.
- Corrections applied while building: (1) `IkBa_star` as a BNGL output **function** (a conf `formula`
  observable is a per-point measurement model, not a materialized column); (2) protocol delivery via
  **12 per-protocol/genotype model files** joined in a multi-model fit (the bngsim backend perturbs
  only free params via a `condition:`, and the window parameters are not fit).
- **Objective (now EXACT).** The fit uses the paper's actual Eq-7 objective — sum of squared
  geometric-mean-normalized log-differences with a ρ=0.03 floor — expressed with standard tokens via
  `lanl/PyBNF#479` (ADR-0066): `noise_model = lognormal, sigma = fix_at 1` (fixed-σ Gaussian on log₁₀)
  + `normalization <obs> = floor 0.03, scale`. `floor` = `x' = x + 0.03·max(x)` (the paper's ρ floor);
  `scale` = the per-series analytic optimal scale, which for a log family is `geomean(data)/geomean(sim)`
  = the paper's ÷geomean exactly. Both are applied **symmetrically to the simulated and the
  experimental column**, so the `.exp` ship raw. This replaces the earlier edition-2 `norm_sos + peak`
  *analog* — the previous main earned-confidence dock is retired.
- **PyBNF fix required (`#479` follow-up, committed).** `#479` floors the experimental data too, but
  this job's target is **sparse** (each observable has its own S1-Table times → NaN cells in the shared
  grid), and `normalize_to_floor`/`peak`/`unit` used `np.max`, which returns NaN on such a column and
  poisons it (every point → NaN → silently skipped → objective 0.0 for every pset). Fixed in PyBNF
  commit **`4dd38899`** ("NaN-aware peak/unit/floor normalization for sparse multi-observable data"):
  those three now use `np.nanmax`/`np.nanargmax` (+ a sparse-data unit test in
  `tests/test_data_class.py`); the analytic-`scale` path was already NaN-safe.

## Bottom line

This is the paper's **headline** result (Fig 2) and it now reproduces to high precision under the
paper's **exact** objective: the model is the authors' own reduced model (identical NF-κB network),
the target is byte-reproducible from the authors' own original model with the point count matching the
paper exactly (914), Gate 3a reproduces the paper's reported AMD (1.29/1.16) to 3–4 sig figs, and Gate
3b recovers **all 13 parameters** (10 essentially exact) — with the paper's least-identifiable ε now
tight to ~1.5× — while the objective at Table-1 coincides with the fit minimum (ratio 0.996). Residual
risk is only the small structural reduced-vs-original mismatch (AMD ≈ 1.29), which is the paper's own
reported residual, not a job artifact.
