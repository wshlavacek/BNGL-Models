# VALIDATION — Jaruszewicz-Blonska-2023/reduced_combination

Primary-source validation of the PyBNF job `pybnf-jobs/Jaruszewicz-Blonska-2023/reduced_combination/`,
per the `validate-pybnf-job` skill. Confidence is **earned from the gate evidence below**.

> **Confidence: 90 / 100** — Gates 0–2 PASS/equiv against the authors' OWN reduced-model BNGL (Gate 2
> = identical NF-κB network bar an inert time-clock); Gate 1 target is byte-reproducible from the
> authors' original-model BNGL and matches the paper's **914 independent data points exactly**; Gate
> 3a reproduces the paper's **reported AMD\*_WT = 1.29 / AMD\*_A20KO = 1.16 to 3–4 significant
> figures**; Gate 3b recovers 12/13 parameters within 3× (most near-exact), with only ε (the paper's
> single least-identifiable parameter) loose. Docked mainly for the edition-2-native objective being
> an *analog* of the paper's exact geometric-mean/log objective.

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

The 12 `.exp` are the **original Lipniacki-2004 model's** output for the 6 protocols × {WT, A20KO},
sampled at the **S1-Table** measurement times and peak-normalized (with the paper's ρ=0.03·max
floor). `validation/gen_combination.py` runs the **authors' own**
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
  oscillation (continuous) and the entrained pulse trains (5-min/22.5-min/45-min pulses).

**Verdict:** PASS — the reduced model at the paper's own params reproduces the paper's own Fig-2
combination result and its reported AMD.

## Gate 3b — Recover the paper's parameters by fitting — **PASS (12/13 within 3×; only the paper's least-identifiable ε loose)**

- Run: `pybnf -c reduced_combination.conf` — edition-2 `de`, 12-model joint fit, pop 100 × 150 iters
  + refine, `random_seed=1` (converged in ~18 min, final objective **109.3**). Objective **norm_sos +
  normalization=peak** (edition-2-native analog of the paper's Eq-7 log objective; native-only).
- `scripts/compare_params.py published_table1.json output/Results/sorted_params_final.txt` (tol 3×):

| param | published (Table 1) | recovered | log10 ratio | within 3×? | note |
|---|---|---|---|---|---|
| k_deg | 0.000107 | 0.000113 | +0.02 | yes | near-exact |
| k_1 | 0.00195 | 0.00185 | −0.02 | yes | near-exact |
| k_3 | 0.00145 | 0.00119 | −0.08 | yes | ✓ |
| k_2 | 0.0357 | 0.0264 | −0.13 | yes | ✓ |
| a_3 | 0.0946 | 0.0562 | −0.23 | yes | ✓ (on-off left this loose — combination identifies it) |
| **delta (δ)** | 0.1083 | 0.187 | +0.24 | yes | ✓ (~1.7×; on-off left this at ~3×) |
| **epsilon (ε)** | 0.0428 | 0.250 | +0.77 | NO (~5.9×) | paper's least-identifiable param (Fig 5) |
| c_deg | 0.000106 | 0.000110 | +0.02 | yes | near-exact |
| c_4a | 0.00313 | 0.00359 | +0.06 | yes | near-exact |
| a_2 | 0.0763 | 0.0544 | −0.15 | yes | ✓ |
| c_5a | 5.78e-05 | 6.86e-05 | +0.07 | yes | ✓ |
| i_1a | 0.000595 | 0.000911 | +0.18 | yes | ✓ |
| c_3a | 0.000372 | 0.000331 | −0.05 | yes | near-exact |

**12/13 within 3×**, most within 1.3× (6 near-exact within ~1.15×). Only **ε** drifts (~5.9×) — the
single parameter the paper's identifiability analysis flags as least identifiable across all protocols
(Fig 5B/5C/5D). Critically, the richer combination experiment recovers **δ (~1.7×) and a_3 (~1.7×)**
far tighter than the sparse on-off slug (which left both near/beyond 3×), exactly matching the paper's
finding that the combination experiment improves identifiability of those directions.

**Verdict:** PASS — 12/13 parameters recovered (most near-exact); only the paper's own
least-identifiable parameter ε remains loose (~5.9×).

---

## Divergence & corrections

- Scope vs. paper: **matches** — built to the paper's actual model (reduced 2023, all 13 params free)
  and its headline combination experiment (Fig 2, S1 Table, WT + A20KO, 914 points), recovering Table 1.
- Corrections applied while building: (1) `IkBa_star` as a BNGL output **function** (a conf `formula`
  observable cannot be peak-normalized → `inf`); (2) the paper's ρ=0.03·max floor on the target
  (original IKK_a is 0 at rest → `norm_sos` ÷0); (3) protocol delivery via **12 per-protocol/genotype
  model files** joined in a multi-model fit (the bngsim backend perturbs only free params via a
  `condition:`, and the window parameters are not fit).
- Objective/normalization is an **edition-2-native analog** of the paper's exact geometric-mean/log
  objective (native-only) — the main earned-confidence dock.

## Bottom line

This is the paper's **headline** result (Fig 2) and it reproduces to high precision: the model is the
authors' own reduced model (identical NF-κB network), the target is byte-reproducible from the
authors' own original model with the point count matching the paper exactly (914), Gate 3a reproduces
the paper's reported AMD (1.29/1.16) to 3–4 sig figs, and Gate 3b recovers 12/13 parameters (most
near-exact) with only the paper's own least-identifiable ε loose. Residual risk: the objective is an
analog, not the paper's exact objective. Highest-value next step: implement the paper's exact
objective (geometric-mean-normalized log sum-of-squares) as a native/callable PyBNF objective and
confirm ε tightens.
