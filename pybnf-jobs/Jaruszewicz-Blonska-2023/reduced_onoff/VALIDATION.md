# VALIDATION — Jaruszewicz-Blonska-2023/reduced_onoff

Primary-source validation of the PyBNF job `pybnf-jobs/Jaruszewicz-Blonska-2023/reduced_onoff/`, per
the `validate-pybnf-job` skill. Confidence is **earned from the gate evidence below**.

> **Confidence: 86 / 100** — Gates 0–2 PASS/equiv against the authors' OWN reduced-model BNGL (Gate 2
> = identical NF-κB network bar an inert time-clock); Gate 1 target is byte-reproducible from the
> authors' original-model BNGL; Gate 3a reproduces the original on-off dynamics (AMD\*≈1.4); Gate 3b
> recovers 10/13 parameters within 3× at a good objective, and the 3 that drift (a_3, δ, ε) are
> exactly the paper's documented least-identifiable parameters. Docked for the edition-2-native
> objective being an *analog* of the paper's exact geometric-mean/log objective (not the same), and
> for ε drifting ~9×.

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

The two `.exp` are the **original Lipniacki-2004 model's** on-off output (the paper's fit target),
sampled at the Table-2 measurement times and normalized. Provenance is exact and reproducible:
`validation/gen_onoff.py` runs the **authors' own** `lipniacki_original_model_on_off_protocol.bngl`
through BNG2.pl, samples observables {IKKa, NFkB_n, A20, IkBa*(=IkBa_total), IkBa mRNA} at Table-2
times (protocol t = 3600 s pre-equilibration offset), and writes the target.

| `.exp` | source | measured variables × times | normalization | verdict |
|---|---|---|---|---|
| `reduced_onoff_wt.exp` | original model, WT | 5 vars, Table-2 times (27 measurements) | peak + ρ=0.03·max floor (per series, union grid) | PASS |
| `reduced_onoff_a20ko.exp` | original model, A20KO | 4 vars (no A20), Table-2 times (23 measurements) | same | PASS |

Total 50 measurements → 41 independent after per-series scaling — **matches the paper's "50 in silico
measurements / 41 independent data points"** for the on-off protocol exactly. Normalization = the
edition-2-native analog of the paper's Methods recipe (`x' = x + 0.03·max(x)`, then divide by a
per-series constant); the paper uses the geometric mean, this job uses the trajectory peak (both
remove the reduced-vs-original unit difference; see Gate 3b note).

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

**Verdict:** PASS — the reduced model at the paper's own params reproduces the original on-off
dynamics (consistent with S6 Fig).

## Gate 3b — Recover the paper's parameters by fitting — **PARTIAL/PASS (identifiable params + sloppy directions per the paper)**

- Run: `pybnf -c reduced_onoff.conf` — edition-2 `de`, pop 100 × 150 iters + refine, `random_seed=1`.
  Objective **norm_sos + normalization=peak** (edition-2-native analog of the paper's Eq-7 log
  objective; native-only). Final objective **7.24** (< the objective **17.1** at the exact Table-1
  params — a different-objective optimum sits lower, as for tlbr/igf1r/Rukhlenko).
- `scripts/compare_params.py published_table1.json output/Results/sorted_params_final.txt` (tol 3×):

| param | published (Table 1) | recovered | log10 ratio | within 3×? | note |
|---|---|---|---|---|---|
| k_1 | 0.00195 | 0.00220 | +0.05 | yes | identifiable ✓ |
| k_2 | 0.0357 | 0.0233 | −0.19 | yes | ✓ |
| k_3 | 0.00145 | 0.00239 | +0.22 | yes | ✓ |
| k_deg | 0.000107 | 0.000315 | +0.47 | yes | ✓ |
| c_deg | 0.000106 | 0.000138 | +0.11 | yes | ✓ |
| c_3a | 0.000372 | 0.000289 | −0.11 | yes | ✓ |
| c_4a | 0.00313 | 0.00498 | +0.20 | yes | ✓ |
| c_5a | 5.78e-05 | 9.5e-05 | +0.22 | yes | ✓ |
| a_2 | 0.0763 | 0.0469 | −0.21 | yes | ✓ |
| i_1a | 0.000595 | 0.000689 | +0.06 | yes | ✓ |
| **a_3** | 0.0946 | 0.0312 | −0.48 | NO (~3×) | least-identifiable (couples to δ) |
| **delta (δ)** | 0.1083 | 0.347 | +0.51 | NO (~3.2×) | paper: least identifiable (Fig 5) |
| **epsilon (ε)** | 0.0428 | 0.00444 | −0.98 | NO (~9.5×) | paper: least identifiable (Fig 5) |

**10/13 within 3×.** The 3 that drift are exactly the parameters the paper's own linear + Monte Carlo
identifiability analysis flags as least sensitive / least identifiable: "the parameters ε and δ
appear as two out of three least sensitive parameters for all protocol sets" (Results, Fig 5B/5C),
and a_3 enters the reduced model only as the products `a_3·δ` and `a_3·(1−NFkBn)` (Eqs 3, 5) so a_3
and δ trade off. The identifiable directions are recovered tightly (|log10 ratio| ≤ 0.22 for 8/13).

**Verdict:** PARTIAL/PASS — identifiable parameters recovered; the sloppy directions (a_3, δ, ε)
drift exactly as the paper predicts. Scored as the paper's own identifiability result, not as a
failure (cf. tlbr/igf1r policy).

---

## Divergence & corrections

- Scope vs. paper: **matches** — the job was built (it did not exist) to the paper's actual model
  (reduced 2023, Eqs 1–6, all 13 params free) and protocol (Table 2 on-off, WT + A20KO), recovering
  Table 1. No simplification.
- Corrections applied while building: (1) `IkBa_star` moved from a conf `formula` observable to a
  BNGL output **function** — a `formula` observable is a per-point measurement model, not a
  materialized column, and `normalization=peak` silently yields `inf` on it. (2) Added the paper's
  ρ=0.03·max floor to the target (the original IKK_a is exactly 0 at rest → `norm_sos` divides by 0).
  Both are recorded in the `.bngl`/`.conf` comments and here.
- Objective/normalization is an **edition-2-native analog** of the paper's exact geometric-mean/log
  objective (which is not in the PEtab-exportable subset); this is the main earned-confidence dock.

## Bottom line

Solid: the model is the authors' own reduced model (identical NF-κB network), the fit target is
byte-reproducible from the authors' own original model, the point count matches the paper (41
independent), Gate 3a reproduces the original dynamics, and Gate 3b recovers the identifiable
parameters while the drifting ones are precisely the paper's documented sloppy directions. Residual
risk: the objective is an analog, not the paper's exact objective, and ε is recovered only to ~9×.
Highest-value next step to raise confidence: implement the paper's exact objective (sum of squared
geometric-mean-normalized log-differences) as a native/callable PyBNF objective and confirm ε/δ/a_3
tighten (or a Monte Carlo ensemble like the paper's, rather than a single fit).
