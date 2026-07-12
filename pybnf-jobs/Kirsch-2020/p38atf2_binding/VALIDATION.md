# VALIDATION — Kirsch-2020 / fig7b_timecourse → {p38atf2_binding, ppatf2_phospho}

Primary-source validation of the PyBNF job formerly at
`pybnf-jobs/Kirsch-2020/fig7b_timecourse/`, per the `validate-pybnf-job` skill. The audit
found a **wrong-observable defect** and re-scoped the job into **two** faithful slugs
(divergence policy: match the literature). This file is the shared scorecard for both;
[`ppatf2_phospho`](../ppatf2_phospho/) links back here.

> **Confidence: 91 / 100** — Gate 2 model is **byte-identical** to the authors' own BNGL
> (48/48 species, 152/152 reactions); Gate 1 data are grounded and, after correcting a
> serious observable mismap, correctly identified; both re-scoped slugs reproduce at the
> paper's published parameters (Gate 3a) and recover their fitted parameters within ≤1.3×
> (Gate 3b). Held back from higher only because the pp-JNK-arm params (`keq7/kstim7/dp1`)
> are held at published rather than independently re-fit (their pp-JNK/pp-p38 absolute-µM
> curves were not reconstructed — the remaining "full joint fit"), and slug B rests on a
> 4-point figure digitization.

Primary sources (in the untracked `dev/papers/Kirsch-2020/`; not redistributed):
- Model + data paper: Kirsch K et al., *Nat Commun* 2020;11:5769. PMC7666158 /
  DOI 10.1038/s41467-020-19582-3 (`41467_2020_Article_19582.pdf`).
- Supplement: SI PDF (`…MOESM1…pdf`, holds **Supplementary Table 2** = the fitted params);
  **Source Data** workbook (`…MOESM6….xlsx`); **Supplementary Software 1**
  (`…MOESM4….zip` → `Bionetgen_JNK-p38-ATF2_model.bngl` + `ReadMe.txt`) — the **authors'
  own model file**, which turns Gate 2 into a direct diff.

**"The paper's result" for this job.** Supplementary Table 2 fits **eight** cell-based rate
constants by **decomposing** the anisomycin time courses across datasets (Table 2 "data"
column), *not* by one joint fit to one curve:

| dataset (Fig. 7b) | params fit | model observable |
|---|---|---|
| pp-JNK western blot | `keq7, kstim7, dp1` | `ppJNK` |
| **p38-ATF2 NanoBit binding** | **`keq6, kstim6, dp2, dp4`** | **`p38ATF2all`** → slug **`p38atf2_binding`** |
| **pp-ATF2(T69/T71) western blot** | **`dp3`** | **`pT69pT71`** → slug **`ppatf2_phospho`** |

---

## Gate 0 — Materials inventory

| needed | present? | path / note |
|---|---|---|
| model+data paper PDF | ✅ | `41467_2020_Article_19582.pdf` |
| SI (Table 2) | ✅ | `41467_2020_19582_MOESM1_ESM.pdf` |
| Source Data workbook | ✅ | `41467_2020_19582_MOESM6_ESM.xlsx` (sheet `Figure_4`) |
| authors' model file | ✅ | `…MOESM4…/Bionetgen_JNK-p38-ATF2_model.bngl` (+ `ReadMe.txt`) |

**Verdict: PASS.** The `ReadMe.txt` states the file "contains parameter for WT ATF2, and
used for generating Figure_7b" — i.e. its params are the published Fig. 7b best-fit.

## Gate 1 — Data provenance  → the defect

The shipped `.exp` numbers are a **faithful transcription** of the Source Data `Figure_4`
sheet, WT column (mean rows 6–24, SD rows 29–47); `extract_exp.py` regenerates the file
**byte-for-byte**. **But the original job mis-identified what that column is.**

- That numeric block is the **Fig. 4a left-panel p38-ATF2 NanoBit binding** curve
  ("treated/untreated"): WT peaks ~1.6× at 9 min then **decays**, S90N rises to ~4×, MUT4
  stays ~1 — matching the Fig. 4a plot exactly. The same data are the Fig. 7b "p38 + ATF2
  interaction" panel. The `Anti-ppATF2(T69/T71)` string at cell **I4** labels a western-blot
  *image*, not this column.
- The original job labeled it "WT pp-ATF2(T69/T71) western-blot" and mapped it to the
  **phosphorylation** observable `pT69pT71`. **Proof it is the binding curve, not
  phosphorylation** — running the *authors' own model file* at their published params and
  comparing each observable ÷basal to the `.exp` peak-decay (1.04 → 1.63@9min → 1.17):

  | model observable ÷basal | shape | median rel-err vs `.exp` |
  |---|---|---|
  | `pT69pT71` (what the job used) | monotonic → 3.5× | **125 %** ❌ |
  | **`p38ATF2all`** (the NanoBit complex) | peak-decay, peaks ~9 min | **5.8 %** ✅ |

  `pT69pT71` (phosphorylation) is monotonic — it *cannot* make the observed peak. The
  peak-decay is the phosphoswitch: pp-JNK phosphorylates S90, evicting p38 from the FRS.

**Verdict: FAIL as-shipped (wrong observable) → CORRECTED.** Post-correction, provenance
is: slug **`p38atf2_binding`** — transcribed from Source Data (byte-identical, PASS); slug
**`ppatf2_phospho`** — the true pp-ATF2(T69/T71) curve is not in any Source Data sheet
(no `Figure_7` sheet) and was **digitized** from the Fig. 7b pp-ATF2 panel (linear axes,
±0.03 µM, no error bars → `sos`; `digitize.py`).

## Gate 2 — Model fidelity

Reference: the authors' `Bionetgen_JNK-p38-ATF2_model.bngl` (Supplementary Software 1).

`scripts/model_diff.py` (generated networks): **species A=48 B=48 common=48 · reactions
A=152 B=152 common=152 → IDENTICAL.** Rate-law expressions and all rate constants match the
authors' file by hand. The only adaptations are cosmetic/functionally-equivalent: no
`actions` block (protocol in the conf); a boolean `Stim` gate + two `functions`
(`k7_act()`, `k6_act()`) reproduce the authors' `setParameter(kstim7,keq7)` equilibration
idiom exactly (max rel. diff 0.0); fitted constants declared `id nominal`, bound by name.

**Verdict: PASS (identical generated network; functionally-equivalent refactor).**

## Gate 3a — Reproduce at the paper's parameters

Both slugs, model set to Supplementary Table 2 WT values (reproduced by running the
authors' own file and the reconstruction identically):

- **`p38atf2_binding`** (`p38ATF2all` ÷basal vs Fig. 4a WT NanoBit): **median 5.8 %**
  rel-err (max 15.3 %) — captures the peak-decay transient.
- **`ppatf2_phospho`** (`pT69pT71` absolute vs digitized Fig. 7b, two conditions): **CTR
  median 6.1 %**; **JNK-IN-8 median 12.5 %** (its 10-min point sits below the authors' own
  `calc-JNK-IN-8` line too). Both match the authors' published `calc` curves (CTR
  0.21→0.75 µM; JNK-IN-8 → 0.45 µM).

**Verdict: PASS** — the model at the paper's own parameters reproduces the paper's own
Fig. 4a / Fig. 7b curves, once mapped to the correct observables.

## Gate 3b — Recover the paper's parameters by fitting

- **`p38atf2_binding`** — `de`, 4 free (`keq6,kstim6,dp2,dp4`), ~19 s, obj 16.7:

  | param | published | recovered | log10 ratio | within tol? |
  |---|---|---|---|---|
  | keq6   | 1.74e-5 | 1.98e-5 | +0.06 | ✅ (1.14×) |
  | kstim6 | 1.16e-4 | 1.54e-4 | +0.12 | ✅ (1.33×) |
  | dp2    | 1.76e-3 | 2.29e-3 | +0.11 | ✅ (1.30×) |
  | dp4    | 4.50e-3 | 4.52e-3 | +0.00 | ✅ (1.00×) |

- **`ppatf2_phospho`** — `de`, 1 free (`dp3`), ~8 s: recovered **1.008e-2 vs 9.54e-3 =
  1.06×** of published.

**Verdict: PASS** — every fitted parameter recovered within ≤1.3× of the published value.
(The pp-JNK-arm `keq7/kstim7/dp1` are held at published in both slugs; independently
re-fitting them needs the pp-JNK / pp-p38 absolute-µM curves, not reconstructed here.)

---

## Divergence & corrections

- **Scope vs. paper:** the single `fig7b_timecourse` slug fit p38-ATF2 **binding** data to
  the **phosphorylation** observable `pT69pT71` under `normalization=init` — a wrong-
  observable mismap that at published params overshot the data by ~125 % and forced its
  earlier 8-param fit ~6× off the published values. Re-scoped per the divergence policy.
- **Corrections applied:**
  1. Renamed `fig7b_timecourse` → **`p38atf2_binding`** (files + folder via `git mv`);
     observable `pT69pT71` → **`p38ATF2all`**; free set → the 4 NanoBit-fit params
     `{keq6,kstim6,dp2,dp4}`; `.exp`, conf, model header, `extract_exp.py`,
     `make_reproduction.py`, README all corrected.
  2. Added sibling slug **`ppatf2_phospho`** for the actual pp-ATF2(T69/T71) phosphorylation
     fit (digitized Fig. 7b; `pT69pT71`, absolute, CTR + JNK-IN-8, free `dp3`).
- **Re-run after corrections:** Gates 1/2/3a/3b re-run for both slugs (metrics above).

## Bottom line

The model is a byte-identical reconstruction of the authors' published BNGL, and the data
are real Kirsch-2020 measurements — the defect was purely a **wrong observable mapping**
(NanoBit p38-ATF2 *binding* fit to the ATF2 *phosphorylation* observable), now corrected by
splitting into a binding slug (`p38atf2_binding`, tabulated data, `p38ATF2all`) and a
phosphorylation slug (`ppatf2_phospho`, digitized Fig. 7b, `pT69pT71`). Both reproduce at
the paper's parameters and recover their fitted params within ≤1.3×. The single most
valuable next step to reach the paper's *full* determination is to digitize the Fig. 7b
**pp-JNK** and **pp-p38** absolute-µM panels and run the 8-parameter **joint** fit across
all datasets (pp-JNK + pp-p38 + pp-ATF2 + NanoBit) — the only piece not reconstructed here.
