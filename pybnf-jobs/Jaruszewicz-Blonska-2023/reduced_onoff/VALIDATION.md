# VALIDATION — Jaruszewicz-Blonska-2023/reduced_onoff

Primary-source validation of the PyBNF job `pybnf-jobs/Jaruszewicz-Blonska-2023/reduced_onoff/`, per
the `validate-pybnf-job` skill. Confidence is **earned from the gate evidence below**.

> **Confidence: 89 / 100** — Gates 0–2 PASS/equiv against the authors' OWN reduced-model BNGL (Gate 2
> = identical NF-κB network bar an inert time-clock); Gate 1 target is byte-reproducible from the
> authors' original-model BNGL; Gate 3a reproduces the original on-off dynamics (AMD\*≈1.4). The fit
> now uses the paper's **EXACT Eq-7 objective** (fixed-σ log-normal + floor + per-series
> geometric-mean scaling, `lanl/PyBNF#479`) — not an analog — recovering **11/13 parameters within 3×**
> (up from 10/13), with the identifiable directions essentially exact and **ε and a_3 both tightened**
> vs the analog. The 2 that still drift (δ, ε) are exactly the paper's documented least-identifiable
> parameters (Fig 5); on this deliberately sparse 41-point protocol the objective minimum sits ~16 %
> below the value at Table 1 (0.38 vs 0.45) because the sloppy a_3/δ/ε manifold admits a
> marginally-better-fitting point — the information-poor protocol, exactly as the paper describes.
> Docked for that residual non-identifiability of δ (the richer combination slug recovers all 13).

Primary sources (in the untracked `dev/papers/Jaruszewicz-Blonska-2023/`; not redistributed):
- Model paper: Jaruszewicz-Błońska J, Kosiuk I, Prus W, Lipniacki T. "A plausible identifiable model
  of the canonical NF-κB signaling pathway." PLoS ONE 2023; 18(6): e0286416. DOI
  10.1371/journal.pone.0286416. (Reduced 2023 model; Eqs 1–6, Table 1, Table 2, Fig 8, S6 Fig.)
- Data paper: same (the fit target is the paper's ORIGINAL Lipniacki-2004 model, S1 Code).
- Author files used (S1 Code, `pone.0286416.s015.zip`): `Reduced2023 - BNGL file/
  fitted_reduced_model_on_off_protocol.bngl` (the reduced model at Table-1 fitted values) and
  `Lipniacki_et_al_2004_model - BNGL file/lipniacki_original_model_on_off_protocol.bngl` (the
  original 15-variable model = the fit target).

**"The paper's result" for this job** = Table-1 (fitted) parameter values, i.e. the reduced model
refit to reproduce the original model. The on-off protocol (Table 2) is the paper's proposed
identifiability-optimal experiment; S6 Fig shows the reduced fitted on-off WT dynamics. This is a
**synthetic-target / ground-truth-recovery** job (like Gupta fceri): the "data" is the original
model's output; Gate 3b must recover Table 1.

---

## Gate 0 — Materials inventory — **PASS**

| needed | present? | path / note |
|---|---|---|
| model paper PDF | ✅ | `journal.pone.0286416.pdf` |
| SI figures/tables/text | ✅ | `pone.0286416.s001–s014` (Table 1/2 in main PDF; S1/S2 Text = model math) |
| author reduced-model BNGL (fitted) | ✅ | S1 Code `fitted_reduced_model_on_off_protocol.bngl` — enables direct Gate-2 diff |
| author original-model BNGL (fit target) | ✅ | S1 Code `lipniacki_original_model_on_off_protocol.bngl` |

The job did not exist before this validation; it was **built** from the authors' own materials
(divergence policy: "if missing, build"). No `NEEDED.md` existed for this paper; nothing is missing.

## Gate 1 — Data provenance — **PASS (synthetic target, byte-reproducible)**

The two `.exp` are the **original Lipniacki-2004 model's RAW** on-off output (the paper's fit target),
sampled at the Table-2 measurement times — **un-normalized** (no floor, no peak, no geomean, no `_SD`).
Provenance is exact and reproducible: `validation/gen_onoff.py` runs the **authors' own**
`lipniacki_original_model_on_off_protocol.bngl` through BNG2.pl, samples observables {IKKa, NFkB_n,
A20, IkBa*(=IkBa_total), IkBa mRNA} at Table-2 times (protocol t = 3600 s pre-equilibration offset),
and writes the raw target.

| `.exp` | source | measured variables × times | normalization | verdict |
|---|---|---|---|---|
| `reduced_onoff_wt.exp` | original model, WT | 5 vars, Table-2 times (27 measurements) | none — raw output | PASS |
| `reduced_onoff_a20ko.exp` | original model, A20KO | 4 vars (no A20), Table-2 times (23 measurements) | none — raw output | PASS |

Total 50 measurements → 41 independent after per-series scaling — **matches the paper's "50 in silico
measurements / 41 independent data points"** for the on-off protocol exactly. The paper's Methods
recipe (`x' = x + 0.03·max(x)`, then divide by the geometric mean) is now applied **by the objective
at scoring time**, symmetrically to sim and data (`normalization <obs> = floor 0.03, scale`;
`lanl/PyBNF#479`) — so the target ships raw. See the Gate 3b objective note.

## Gate 2 — Model fidelity — **PASS (equiv)**

Reference: the authors' `fitted_reduced_model_on_off_protocol.bngl` (S1 Code).

`scripts/model_diff.py` (generated network, rate values ignored):
```
species  : common=6  onlyA(ours)=1 -> Clock()      onlyB=0
reactions: common=14 onlyA(ours)=1 -> ()->Clock()  onlyB=0
```
The 6 NF-κB species and all 14 NF-κB reactions are **identical** to the authors' model. Documented
functional-equivalence refactors (edition-2, same idioms as Lin-2021/nyc_multiphase):

| aspect | authors' file | our `.bngl` | verdict |
|---|---|---|---|
| molecule types / seed species | IKK(st~n~a), IkBa, IkBa_mRNA, A20, NFkB; (1,0,0,0,0,0) | same **+ inert Clock()** (seed 0) | equiv |
| 14 NF-κB reaction rules + rate-law structure | Eqs 1–6 | identical | match |
| TNF stimulus | `TR` parameter (0/1) × k_1, k_2, set by `setParameter` phases in a `begin actions` block | `TRfunc()=if(simtime<7200,1,0)` × k_1, k_2, with `simtime` from the Clock (rule `0->Clock() 1`) | equiv (delivers the same TR(t)) |
| total IkBa protein IkBa* | (measured var) | `IkBa_star()=IkBa+1-NFkB_n` output function (NOT a conf `formula` observable — that cannot be peak-normalized) | equiv |
| pre-equilibration | 1 h TR=0 phase | none — seed IS the exact TR=0 steady state (verified: all six ODEs = 0 at (1,0,0,0,0,0)) | equiv |
| initial conditions / units | non-dimensional | same | match |

**Verdict:** PASS (equiv) — identical published NF-κB model; the only structural addition is the
inert time-clock that delivers the TNF stimulus in an edition-2 (conf-synthesized-time-course) form.

## Gate 3a — Reproduce at the paper's parameters — **PASS**

Model set to the **Table-1 fitted** values (the shipped `.bngl` nominals). `make_reproduction.py`
(and `validation/gen_onoff.py`) simulate the reduced model at those params and compare to the
original-model on-off target.

- Metric = the paper's Average Multiplicative Distance (AMD, Eq 8): **AMD\*_WT = 1.41**, **AMD\*_A20KO
  = 1.32** (the paper reports numeric AMD only for the combination experiment, 1.29/1.16; the on-off
  values are of the same order — the on-off has far fewer, higher-information points).
- Overlay (`onoff_gate3a.png`, `reduced_onoff_reproduction.png`): the reduced model tracks the
  original for all observables in both genotypes — both NF-κB peaks (TNF-on spike + TNF-off rebound),
  the IκBα* degradation/resynthesis, the mRNA pulse. The one visibly looser variable is WT A20 (the
  reduced A20 rises faster: the reduction eliminated A20 mRNA — a documented, paper-accepted artifact).
- Gate 3a is **objective-independent**; the switch to raw `.exp` only means `make_reproduction.py` now
  peak-normalizes the raw target itself for the overlay (per-obs median |rel err|: IKK_a 4 %, tIkBa
  2 %, NF-κB 24 %, IκBα* 45 %, A20 77 % — the A20-mRNA artifact).

**Verdict:** PASS — the reduced model at the paper's own params reproduces the original on-off
dynamics (consistent with S6 Fig).

## Gate 3b — Recover the paper's parameters by fitting — **PASS (11/13 within 3×; the paper's EXACT Eq-7 objective; sloppy directions per the paper)**

- Run: `pybnf -c reduced_onoff.conf` — edition-2 `de`, pop 100 × 200 iters + Simplex refine,
  `random_seed=1`. Objective = the paper's **EXACT Eq-7** (`noise_model = lognormal, sigma = fix_at 1`
  + `normalization <obs> = floor 0.03, scale`; `lanl/PyBNF#479`) — fixed-σ log-normal residual on
  log₁₀ with the ρ=0.03 floor and per-series geometric-mean scaling applied symmetrically to sim and
  data. Final objective **0.380**.
- `scripts/compare_params.py published_table1.json output/Results/sorted_params_final.txt` (tol 3×):

| param | published (Table 1) | recovered | log10 ratio | within 3×? | note |
|---|---|---|---|---|---|
| k_1 | 0.00195 | 0.00220 | +0.05 | yes | identifiable ✓ |
| k_3 | 0.00145 | 0.00139 | −0.02 | yes | ✓ |
| c_3a | 0.000372 | 0.000403 | +0.04 | yes | ✓ |
| i_1a | 0.000595 | 0.000491 | −0.08 | yes | ✓ |
| c_deg | 0.000106 | 0.000149 | +0.15 | yes | ✓ |
| c_5a | 5.78e-05 | 8.49e-05 | +0.17 | yes | ✓ |
| k_deg | 0.000107 | 6.95e-05 | −0.19 | yes | ✓ |
| k_2 | 0.0357 | 0.0198 | −0.26 | yes | ✓ |
| c_4a | 0.00313 | 0.00584 | +0.27 | yes | ✓ |
| a_2 | 0.0763 | 0.1565 | +0.31 | yes | ✓ |
| **a_3** | 0.0946 | 0.1973 | +0.32 | yes | now within 3× (was ~3× / +0.48 under the analog); couples to δ |
| **epsilon (ε)** | 0.0428 | 0.159 | +0.57 | NO (~3.7×) | tightened from ~9.5× under the analog; paper: least identifiable (Fig 5) |
| **delta (δ)** | 0.1083 | 0.0121 | −0.95 | NO (~9×) | paper: least identifiable (Fig 5); trades off with a_3/ε |

**11/13 within 3×** (up from 10/13 under the analog). The identifiable directions are recovered
tightly (|log₁₀ ratio| ≤ 0.19 for 7/13). The 2 that drift, **δ and ε**, are exactly the parameters the
paper's own linear + Monte Carlo identifiability analysis flags as least sensitive: "the parameters ε
and δ appear as two out of three least sensitive parameters for all protocol sets" (Results, Fig
5B/5C). Relative to the edition-2 analog, the exact objective **tightens ε (from ~9.5× to ~3.7×) and
a_3 (into 3×)**; on this sparse protocol δ instead drifts low (a_3, δ, ε co-vary along the sloppy
manifold — a_3 enters only as `a_3·δ` and `a_3·(1−NFkBn)`, Eqs 3, 5), the exact trade-off the paper
documents. The richer combination slug recovers **all three** tightly.

- **Objective at Table-1 vs the fit minimum.** `obj_at_table1.py` scores the reduced model at the
  Table-1 nominals with the *same* objective (PyBNF `fit_type = check`): **0.452**, vs the fit's
  **0.380** (ratio 0.84). On this deliberately sparse 41-point protocol the empirical minimum sits
  ~16 % below the ground-truth value — the sloppy a_3/δ/ε manifold admits a marginally-better-fitting
  point at 41 samples. (On the data-rich combination slug the two coincide, ratio 0.996.)

**Verdict:** PASS — identifiable parameters recovered essentially exactly under the paper's exact
objective; the 2 remaining sloppy directions (δ, ε) are precisely the paper's documented
least-identifiable parameters for this information-poor protocol.

---

## Divergence & corrections

- Scope vs. paper: **matches** — the job was built (it did not exist) to the paper's actual model
  (reduced 2023, Eqs 1–6, all 13 params free) and protocol (Table 2 on-off, WT + A20KO), recovering
  Table 1. No simplification.
- Corrections applied while building: `IkBa_star` is a BNGL output **function** (a conf `formula`
  observable is a per-point measurement model, not a materialized column). Recorded in the
  `.bngl`/`.conf` comments.
- **Objective (now EXACT).** The fit uses the paper's actual Eq-7 objective — sum of squared
  geometric-mean-normalized log-differences with a ρ=0.03 floor — via `lanl/PyBNF#479` (ADR-0066):
  `noise_model = lognormal, sigma = fix_at 1` + `normalization <obs> = floor 0.03, scale`, applied
  symmetrically to sim and data (so the `.exp` ship raw). This replaces the earlier `norm_sos + peak`
  *analog* — the previous main earned-confidence dock is retired.
- **PyBNF fix required (`#479` follow-up, committed).** `#479` floors the experimental data too, but
  this job's target is **sparse** (each observable has its own Table-2 times → NaN cells), and
  `normalize_to_floor`/`peak`/`unit` used `np.max`, which returns NaN on such a column and poisons it
  (objective 0.0 for every pset). Fixed in PyBNF commit **`4dd38899`** (NaN-aware peak/unit/floor
  normalization for sparse multi-observable data): those three now use `np.nanmax`/`np.nanargmax`
  (+ a sparse-data unit test); the analytic-`scale` path was already NaN-safe.

## Bottom line

Solid: the model is the authors' own reduced model (identical NF-κB network), the fit target is
byte-reproducible from the authors' own original model, the point count matches the paper (41
independent), Gate 3a reproduces the original dynamics, and Gate 3b — now under the paper's **exact**
objective — recovers the identifiable parameters essentially exactly (11/13 within 3×, ε and a_3
tightened vs the analog), while the 2 drifting directions (δ, ε) are precisely the paper's documented
sloppy directions for this information-poor protocol. Residual risk: on the sparse 41-point on-off,
δ remains non-identifiable and the objective minimum sits ~16 % below the value at Table 1. The
data-rich [`reduced_combination`](../reduced_combination/) slug recovers all 13 parameters and places
Table 1 exactly at the objective minimum — the on-off is the identifiability-optimal *design* but
still the information-poorer of the two experiments, exactly as the paper frames it.
