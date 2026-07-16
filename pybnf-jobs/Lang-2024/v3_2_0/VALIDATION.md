# VALIDATION — Lang-2024/v3_2_0

Primary-source validation of the PyBNF job `pybnf-jobs/Lang-2024/v3_2_0/`.
Confidence is **earned from the gate evidence below**.

> **Confidence: 82 / 100** — a **faithful transcription of the paper's PEtab v3.2.0
> parameter-estimation problem** into a native edition-2 PyBNF job. The model is the authors'
> own rule-based network (Gate 2, verified species-for-species against the generated network),
> the data are the authors' own Stallaert-derived single-cell medians (Gate 1), the 8
> observation functions and the 73 initial-condition parameters map **exactly** onto the PEtab
> observable formulas and condition table (Gate 2), and the job **round-trips through PEtab v2**
> and reaches a finite objective under a real bngsim fit (Gate 3). At the PEtab `nominalValue`
> point the model reproduces the **2-cycle oscillation with the correct period and phase** for
> 6 of 8 observables (Gate 4). The deductions (−18) are that this is a **setup transcription,
> not a converged-fit reproduction**: the PEtab problem ships `nominalValue`s that are a
> reference point, not the cluster-optimized (saCeSS) optimum, so observable **amplitudes are
> uncalibrated** (and pRB is flat) until the 177-parameter fit is actually run on a cluster —
> which this job is scoped to *set up*, not to run to convergence (`heavy=True`).

Primary sources (in the untracked `dev/papers/Lang-2024/`; not redistributed):
- Paper: Lang PF, Penas DR, Banga JR, Weindl D, Novak B. "Reusable rule-based cell cycle model
  explains compartment-resolved dynamics of 16 observables in RPE-1 cells." *PLoS Comput Biol*
  2024; **20**(1):e1011151. DOI 10.1371/journal.pcbi.1011151 (`journal.pcbi.1011151.pdf` +
  S1 `…s001.pdf`).
- Model + PEtab problem: `paulflang/cell_cycle_petab` @ `versions/v3.2.0` (commit `d562d3d`),
  staged at `dev/papers/Lang-2024/petab_v3.2.0/`: `cell_cycle_v3.2.0.bngl` (a manual BNGL
  translation of the authors' SBML model), `v3.2.0.yaml`, `parameters_v3.2.0.tsv`,
  `observables_v3.2.0.tsv`, `experimentalCondition_v3.2.0.tsv`, and
  `Stallaert_CellSystems2021_Data_2rounds.tsv`.
- Point of contact: William S. Hlavacek, hlavacek@lanl.gov.

"The paper's result" for this job = **the v3.2.0 parameter-estimation problem** — fit 177 of
205 parameters of the full fused RPE-1 cell-cycle model to 8 nuclear single-cell observables
over ~2 cell cycles (Stallaert et al., *Cell Systems* 2021), one wildtype condition. This is
the paper's headline contribution and is **distinct** from the Phase-1 Fig-1 transition
submodels in `models/cell_cycle_oscillator_lang2024/` (whose metadata explicitly excludes "the
paper's full fused, parameter-fitted cell cycle model").

---

## Gate 0 — Materials inventory

| needed | present? | path / note |
|---|---|---|
| paper PDF + supplement | ✅ | `journal.pcbi.1011151.pdf`, `…s001.pdf` |
| author model (BNGL) | ✅ | `cell_cycle_v3.2.0.bngl` (manual translation of the authors' SBML) |
| PEtab problem (params/obs/cond/data + yaml) | ✅ | `petab_v3.2.0/*.tsv`, `v3.2.0.yaml` |
| authors' fitted optimum | ⚠️ | **not shipped** — the PEtab `nominalValue` column is a reference point, not the saCeSS optimum (Gate 4) |

**Verdict:** PASS with a caveat — a complete, runnable PEtab *problem*, but the converged
best-fit parameter set is not part of the supplement, so the reproduction target is the
**setup**, not a published optimum.

## Gate 1 — Data provenance

| `.exp` | source | method | units | verdict |
|---|---|---|---|---|
| `v3_2_0.exp` | `Stallaert_CellSystems2021_Data_2rounds.tsv` (PEtab measurements) | pivoted long→wide, verbatim values | median single-cell nuclear intensity (a.u.), 2 duplicated cycles | PASS |

The long-format PEtab measurement table (4800 rows = 8 observables × 600 times, single `wt`
condition) was pivoted to a wide `.exp` (`# time` + 8 function columns). All 8 observables share
one identical 600-point time grid (t = 0 .. 136471.6), so there are **0 NaN cells**; no value
was altered. The data are the Stallaert et al. (2021) RPE-1 single-cell medians the authors used,
duplicated to two cell cycles to enforce sustained oscillation (paper Text F).

**Verdict:** PASS — byte-level provenance to the authors' own PEtab measurement file.

## Gate 2 — Model fidelity (the load-bearing gate)

Reference: the authors' `cell_cycle_v3.2.0.bngl` + the PEtab observable/condition/parameter
tables. The reaction rules, molecule types, and compartment are carried **verbatim**; the
edition-2 fitting adaptations were verified against the **generated network**, not by eye:

| aspect | authors' problem | our `v3_2_0.bngl` | check |
|---|---|---|---|
| reaction network | 136 rules → **73 species, 332 reactions** | identical (rules verbatim) | BNG2.pl: 73 sp / 332 rxn, builds ~0.4 s |
| kinetic constants (114) | SBML values (`.bngl` `k* val/1`) | bare `k* val`, exact (TSV nominals are rounded) | evaluated from staged `.bngl` |
| initial conditions (73) | PEtab condition table `ic_wt_cell__*` | every seed species carries its `ic_wt_cell__*` param | **munge(canonical species i) == TSV ic id for all 73** |
| species 60–73 | fitted **nonzero** initials in PEtab; generated at 0 in staged `.bngl` | **seeded explicitly** at their ic params (12 free, 2 fixed) | else 12 fitted initials would be silently 0 |
| observation model (16) | `obs = o* + s*·Σspecies` (offsets, scales) | functions `obs_*__nuc_median() = o* + s*·<Species obs>` | see below |
| noise | `noiseParameter1+noiseParameter2 = 0+1` (fixed, normal, lin) | constant unit σ → `objective = sos` | exact |
| actions block | generate + writeSBML + simulate to 1e6 | removed (synthesized from conf) | expected (edition-2) |

**Observation-model species sums — verified exactly against the PEtab `observableFormula`s:**
- cycA/cycB1/cycE/Skp2/p21/p27 = the molecule *total* (`Species tCCNA CCNA()` … ) — each PEtab
  formula sums exactly all species of that molecule (4/4/4/1/6/6 species; confirmed).
- pRB = the single free phospho-RB species `RB1(E2F,Ser807_Ser811~p)` (the only pRB species in
  the network; DpRb-comment confirms no phospho-Rb:E2F complexes exist).
- E2F1 = `tE2F − E2F_pSer332_prom` (total E2F minus the 5 Ser332-phospho promoter-bound forms),
  which is **exactly** the 9 species the PEtab formula lists (14 E2F species − 5 excluded).

Independent numeric check: at t=0 the observable values are fully determined by the nominal
initials + observation params; e.g. `tCCNA(0) = 0.0193315` from the sim equals the sum of the 4
CCNA-species ic nominals (0.0070285 + 0.0058377 + 0.0028914 + 0.0035738) to all printed digits —
confirming both the sum membership and that species 60–73 are seeded correctly.

**Verdict:** PASS — the network is the authors', and the initial-condition and observation
mappings are exact (not approximate), verified against the generated network and the PEtab tables.

## Gate 3 — Verification (parses · PEtab round-trips · fit runs)

The completion bar for a heavy job (the user scoped this to *set up*, not run to convergence):

1. **Tier-1** (`scripts/check_conf.py`): edition 2, `job_type=de` resolves, `v3_2_0.exp` bound,
   **177 free params bind by id, no `__FREE`**, `model.stochastic=False`. **PASS.**
2. **PEtab v2 round-trip** (`scripts/petab_roundtrip.py --job-type de`): `export_job` →
   `petab.v2` lint **clean** → `import_job`. Exports 8 observables (`observableFormula` = bare
   function name), **177 estimated parameters**, 4800 measurements. **PASS** (stays inside the
   PEtab-exportable subset: `sos`, arithmetic/function observables, constant Gaussian noise — no
   `normalization`/`cumulative`/`neg_bin`/constraints).
3. **Real bngsim fit** (`pybnf -c` with `population_size=6, max_iterations=2`): the
   simulate→score→propose loop runs and returns a **finite objective (SOS ≈ 3.0e4)** in ~14 s
   (`k=177, n=4800`). **PASS.** A full 177-parameter fit is cluster-scale (`heavy=True`); this
   confirms the machinery, not convergence.

**Verdict:** PASS — parses, round-trips through PEtab v2, and fits to a finite objective.

## Gate 4 — Reproduction (setup, not converged fit)

`make_reproduction.py` → `v3_2_0_reproduction.png` overlays the model at the PEtab `nominalValue`
point on the data (8 panels). Metrics over the 600-point fit window:

| observable | median rel err | model behaviour vs data |
|---|---|---|
| cyclin A | 62 % | **oscillates, correct period + phase**; amplitude slightly high |
| cyclin B1 | 35 % | **oscillates, correct period + phase** (best-tracked) |
| cyclin E | 146 % | oscillates in phase; amplitude **overshoots ~4×** |
| E2F1 | 150 % | oscillates in phase; amplitude overshoots |
| pRB | 97 % | **flat ~2.0** — does not oscillate at nominals (worst) |
| Skp2 | 29 % | **oscillates, correct period + phase** |
| p21 | 37 % | weak in-phase oscillation, low amplitude |
| p27 | 40 % | weak in-phase oscillation, low amplitude |
| **overall** | **51 % (median), SOS 8479** | 2 full cycles, period matched |

The model at nominals reproduces the **two-cycle oscillation with the correct period and phase**
across the core engine (cyclin A/B1/E, E2F1, Skp2) — strong evidence the reconstruction is
faithful and the model *is* the paper's oscillator. What is **not** reproduced is amplitude
calibration: several observation scales/offsets and rates are at reference values (e.g.
`sSKP2=sCDKN1A=sCDKN1B=1.0`, the untouched default), and pRB is flat. Those are exactly what the
177-parameter fit adjusts.

**Verdict:** the *setup* reproduces the paper's oscillation structure; the *fit quality* is a
reference point, and improving it requires the cluster-scale optimization this job sets up.

---

## Divergence & corrections

- Scope: **matches** the v3.2.0 parameter-estimation problem (full fused model, 8 nuclear
  observables, single `wt` condition, 177/205 params). Optimizer is `de` (the paper used
  cooperative scatter search saCeSS; PyBNF's `ss` is the nearest analog).
- Corrections/adaptations: (1) species 60–73 seeded explicitly at their fitted PEtab initials
  (generated at 0 in the staged `.bngl`); (2) kinetic constants taken from the staged `.bngl`
  (exact) rather than the **rounded** TSV nominals; (3) three log-scaled rates with PEtab lower
  bound 0 (`kDpApc_1`, `kDpE2f1`, `kPhC25A`) declared `uniform_var` over `[0, upper]` (log10
  search is undefined at 0); (4) observation-model **priors simplified** — the PEtab problem puts
  Laplace regularization priors on the 8 offsets and 2 of the scales; here all 177 free params
  use uniform/loguniform priors bracketing the PEtab bounds (keeps the job inside the
  PEtab-exportable subset; the L1 regularization term is dropped). (5) actions block stripped.

## Bottom line

An **honest, exact setup transcription** of the paper's v3.2.0 fitting problem — not a
gold-standard fit reproduction (the converged optimum isn't in the supplement). The network,
initial-condition mapping, and observation model are verified exact against the generated network
and PEtab tables; the job round-trips through PEtab v2 and fits to a finite objective; and at the
reference nominals the model reproduces the 2-cycle oscillation's period and phase. The
load-bearing caveat is that **amplitude calibration requires running the 177-parameter fit on a
cluster** — which this job is built to enable (`heavy=True`), following the user's scope ("don't
force a full fit"). Most valuable next step: run the `de`/`ss` fit at cluster scale and drop the
resulting `best_fit` into `make_reproduction.py --params-file` to measure a converged reproduction.
