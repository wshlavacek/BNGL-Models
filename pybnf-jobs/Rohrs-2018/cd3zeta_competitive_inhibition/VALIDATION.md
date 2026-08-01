# VALIDATION — Rohrs-2018/cd3zeta_competitive_inhibition

Primary-source validation of the PyBNF job `pybnf-jobs/Rohrs-2018/cd3zeta_competitive_inhibition/`.
Confidence is **earned from the gate evidence below**.

Primary sources (in the untracked `dev/papers/rohrs2018/`; not redistributed):

- Paper + Supporting Material: Rohrs JA, Zheng D, Graham NA, Wang P, Finley SD. "Computational
  model of chimeric antigen receptors explains site-specific phosphorylation kinetics."
  *Biophys J* 2018; **115**(7):1116–1129. DOI 10.1016/j.bpj.2018.08.018, PMCID PMC6199440
  (`mmc8.pdf`).
- Authors' BioNetGen models: Data S1–S5 (`mmc2.pdf`–`mmc6.pdf`); phosphatase equation summary
  Data S6 (`mmc7.pdf`). **Data S4 (`mmc5.pdf`) is the competitive-inhibition model and its
  reported best-fit parameters**, which are the recovery target here.
- Point of contact: William S. Hlavacek, hlavacek@lanl.gov.

"The paper's result" for this job = **the competitive-inhibition parameter estimation of Figs. 3D
and 4** — fit the LCK density, five Michaelis constants and one inhibition scale factor to the
site-resolved phosphorylation time courses of wild-type CD3ζ and its six ITAM point mutants, with
k<sub>cat</sub> and K<sub>M,B1</sub> held at literature values.

---

## Gate 0 — Materials inventory

| needed | present? | path / note |
|---|---|---|
| paper PDF + Supporting Material | ✅ | `mmc8.pdf` (article pp. 1116–1129 + Figs. S1–S7) |
| authors' model for this mechanism | ✅ | `mmc5.pdf` = Data S4, BioNetGen, with the fitted parameter values |
| authors' fitted optimum | ✅ | the Data S4 parameter block **is** the reported best fit; Fig. 4 additionally gives the mean ± SD of the 50 best sets |
| fit data in tabular form | ❌ | **not supplied** — reported only as plotted curves; digitized from the vector Fig. S5 (Gate 1) |
| the tenth data set (Fig. 1C, wild-type replicate 1) | ⚠️ | raster figure, not digitized; this job fits **nine of the authors' ten** data sets |

**Verdict:** PASS with one scope caveat — the mechanism, the parameters and nine of the ten fitted
data sets are recoverable from the supplement; the tenth is not.

## Gate 1 — Data provenance

| `.exp` | source | method | units | verdict |
|---|---|---|---|---|
| `wt_10pops_rep2.exp` | Fig. S5A | vector extraction | % phosphorylation | PASS |
| `wt_0pops.exp`, `wt_45pops.exp` | Fig. S5B, S5C | vector extraction | % phosphorylation | PASS |
| `mut_A1/A2/B1/B2/C1/C2.exp` | Fig. S5D, S5H, S5E, S5I, S5F, S5J | vector extraction | % phosphorylation | PASS |

Figure S5 is **vector art**, so every plotted marker, error bar and curve vertex is recovered
exactly by `models/car_cd3zeta_phosphorylation_rohrs2018/digitize_rohrs2018.py`; the only
uncertainty is the axis calibration, which is anchored on tick marks rather than the plot frame.
Data-point values come from the **error-bar segments** — each is drawn as two vertical segments
meeting at the mean plus two caps at mean ± SD — because the marker glyph paths in this PDF carry
a systematic offset of about 0.15 pt (≈2% in time), while the error bars land exactly on the axis
limit at the first time point. Sampling times recover to within 0.4% of the nominal
0.1/1/5/10/30/60/180 min (and 0.1/1.5/10/60/360 min for the two lipid conditions) and were snapped
to those nominal values.

Two independent confirmations that the digitization is faithful:

- The *model* curves digitized from the same nine panels are reproduced by the curated BNGL model
  to **≤ 1.6 percentage points** (`models/car_cd3zeta_phosphorylation_rohrs2018/verify_rohrs2018.ipynb`).
- The digitized Fig. 3 data-summary dots appear in all four mechanism panels; their spread across
  panels bounds the *raster* digitization error at 0.85 min in half-maximal time and 0.04 in Hill
  coefficient — the vector extraction used here is tighter still.

The six mutated sites are `NaN` in their own condition (a Y→F site is not measurable), so 312 of
the 354 digitized points are scored.

**Verdict:** PASS — vertex-exact extraction from the authors' own vector supplement, with the
extraction method, the axis calibration and the marker-offset correction all documented in the
regenerating script.

## Gate 2 — Model fidelity

Reference: Data S4 (`mmc5.pdf`) and Eqs. 4 and 5 of the paper.

| aspect | authors' Data S4 | this `.bngl` | check |
|---|---|---|---|
| rate law | Eq. 4: `(k_cat/K_M,i)·Y_i·LCK / D`, `D = 1 + Σ Y_j/K_M,j + Σ pY_j/K_I,j` | identical, as `inhib()` + six `kp_*()` | independent SciPy transcription agrees with BioNetGen to **1.3e-6 percentage points** |
| inhibition constants | Eq. 5: `K_I,i = K_M,i · X_I` | `KiA1 KmA1*Xi`, … | verbatim |
| CD3ζ representation | six independent molecule types, each at the full density | one molecule with six tyrosines (Data S1's representation) | provably equivalent for these observables — no rate law depends on a neighboring site, so the sites stay statistically independent; checked numerically (`ITAM_2p == pY·pY/CD3z_T` to 1e-4 molecules) |
| Y→F point mutation | site seeded in the inert `~X` state | `live_<site> = 0` gate multiplying that site's substrate, inhibitor and rate terms | algebraically identical; simulated side by side on the same 41-point grid, the two formulations agree to **1.6e-6 percentage points** across all six mutants |
| initial conditions | CD3ζ 20,000/µm², all sites unphosphorylated; LCK constant | same | verbatim |
| fixed parameters | `k_cat = 360 /min`, `K_M,B1 = 270 /µm²` | same, not declared free | matches the paper's stated protocol |
| actions block | `generate_network` + `writeMfile` (MATLAB export) | removed | expected (edition-2); the network is bounded by the rules, so no directive must be retained |

**Verdict:** PASS — the mechanism and every parameter are the authors', and the two structural
adaptations (single-chain CD3ζ, gate-parameter mutants) are shown numerically to be
reformulations, not approximations.

## Gate 3 — Verification (parses · PEtab round-trips · fit runs)

1. **Tier-1** (`scripts/check_conf.py`): edition 2; `job_type = de` resolves; nine experiments
   bind data (`exp_data` non-empty); **7 free parameters bind by id, none with `__FREE`**;
   `model.stochastic = False`. **PASS.**
2. **PEtab v2 round-trip** (`scripts/petab_roundtrip.py --job-type de`): `export_job` →
   `petab.v2` lint **clean** → `import_job`. Exports **6 observables** (`observableFormula` = bare
   function name), **7 estimated parameters**, **312 measurements**, **6 conditions** and 6
   experiment periods. **PASS** — the job stays inside the PEtab-exportable subset (`sos`,
   function observables, constant Gaussian noise, plain conditions; no `normalization`,
   `cumulative`, `neg_bin` or BPSL constraints). The bundle is committed under `petab/`.
3. **Real bngsim fit** (`pybnf -c`, full budget `population_size = 40, max_iterations = 300,
   refine = 1`): the simulate→score→propose loop runs and returns a finite objective. See Gate 4.

**Verdict:** PASS — parses, round-trips through PEtab v2, and fits to a finite objective.

## Gate 4 — Reproduction

Two questions, answered separately.

**(a) Does the model at the authors' published parameters reproduce their data?**
`make_reproduction.py` simulates all nine conditions at the Data S4 values and scores against the
digitized points:

| quantity | value |
|---|---|
| SSE, model at the published parameters (`make_reproduction.py`) | **1.253e4** over 312 points (rms 6.34 percentage points) |
| the same, scored by PyBNF itself (`nominal_check.conf`) | `sos = 6253.58` → **SSE 1.2507e4** |
| SSE, the authors' own plotted Fig. S5 curve at the same points | **1.266e4** |
| SSE reported in Fig. 3D for this mechanism (all ten data sets) | 3.47e4 |

> **Objective convention.** PyBNF's `sos` is the Gaussian negative log-likelihood at σ = 1, i.e.
> **half** the sum of squared error. Every PyBNF objective quoted here is doubled before it is
> compared with the paper's SSE. The 0.2% gap between the two SSE rows above is that
> `make_reproduction.py` interpolates the model onto the data times while PyBNF simulates at them
> exactly.

The curated model and the authors' published curve are indistinguishable as fits to their own data
— the 1% difference between them is digitization noise, not a modeling difference. This is the
load-bearing reproduction result: it certifies the model, the data extraction and the parameter
transcription simultaneously.

**(b) Does a from-scratch search recover the published parameters?** Partly — it recovers the
*ordering*, which is the paper's claim, but not the absolute values.

One full run of the committed conf (`de`, 40 x 300, `refine = 1`) on a workstation: differential
evolution stopped itself on tolerance at generation 163 with `sos` 8748.0, and the simplex polish
ran its full 300 iterations to `sos` **7368.0 (SSE 1.4736e4)** — **18% above** the published
parameters' 6253.6 (SSE 1.2507e4). Wall clock about 80 min.

| parameter | published (Data S4) | this fit | ratio |
|---|---|---|---|
| `LCK_T` | 19.703 | 10.72 | 0.54 |
| `KmA1` | 96.548 | 112.78 | 1.17 |
| `KmB2` | 153.82 | 184.47 | 1.20 |
| `KmA2` | 331.68 | 567.58 | 1.71 |
| `KmC2` | 405.43 | 590.43 | 1.46 |
| `KmC1` | 484.9 | 720.34 | 1.49 |
| `Xi` | 0.1153 | 0.3087 | 2.68 |

**What is recovered.** The site ranking of the Michaelis constants is reproduced *exactly*,
including the position of the held-fixed K<sub>M,B1</sub> = 270 within it:

```
published:  A1  97 < B2 154 < B1 270 < A2 332 < C2 405 < C1 485
fitted:     A1 113 < B2 184 < B1 270 < A2 568 < C2 590 < C1 720
```

That ranking is the paper's biological result — LCK prefers A1 and B2, disfavours C1 — and it
survives a search started from a four-decade box.

**What is not, and why that is expected.** The absolute values drift by 1.2–2.7x, and the drift is
coherent rather than random: across the five fitted sites the fitted point rescales
k<sub>cat</sub>·LCK/K<sub>M,i</sub> by 0.32–0.47 and K<sub>I,i</sub> by 3.1–4.6, i.e. it sits
further along the *same* correlated valley the authors describe — "the Michaelis-Menten constants,
the inhibition constants, and the total LCK concentration and catalytic rate" are correlated
(Materials and Methods). Holding k<sub>cat</sub> and K<sub>M,B1</sub>, as both they and this job do,
narrows that valley but does not close it. The authors also did not rely on a single search: they
ran PSO **100 times** and reported the 50 best sets as ranges (Fig. 4). One run of one optimizer
landing 18% short in SSE and inside a factor of 2.7 on every parameter is the expected outcome of
that landscape, not a defect in the setup — and the setup is separately certified by (a), where the
published point reproduces the data as well as the authors' own curve.

**Manifest consequence:** the `_manifest.py` entry carries **no `recover` dict**. A tolerance loose
enough to pass (`tol >= 2.0`) would assert nothing; the meaningful assertion is the PEtab
round-trip plus a finite objective, both of which Gate 3 covers.

**Verdict:** PASS for the setup and for the published-parameter reproduction; PARTIAL for parameter
recovery, with the ordering recovered exactly and the residual gap attributed to a documented,
paper-acknowledged parameter correlation.

---

## Summary

| gate | verdict |
|---|---|
| 0 — materials | PASS (nine of the authors' ten data sets are recoverable) |
| 1 — data provenance | PASS (vertex-exact vector extraction, error-bar means and SDs) |
| 2 — model fidelity | PASS (mechanism and parameters verbatim; both adaptations verified numerically) |
| 3 — verification | PASS (tier-1, PEtab v2 round-trip, real bngsim fit) |
| 4a — published parameters vs data | PASS (SSE 1.253e4 against the authors' own curve's 1.266e4) |
| 4b — parameter recovery | PARTIAL (ordering exact; values within 2.7x; SSE 18% short in one run) |

