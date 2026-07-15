# VALIDATION — Lin-2021/nyc_multiphase

Primary-source validation of the PyBNF job `pybnf-jobs/Lin-2021/nyc_multiphase/`.
Confidence is **earned from the gate evidence below**.

> **Confidence: 90 / 100** — a **gold-standard reproduction**. The model, data, and fitted
> parameters are all the authors' own (`examples/COVID19forecasting_aMCMC/BigApple/`), the model
> ships the **correct NYC census S0**, and the **published MAP reproduces the data** (median 15.9 %
> rel err over t=50..151; peak within 9 %; cumulative within 2 %), capturing the broad NYC plateau
> that the single-phase model cannot. This is the strongest kind of job — same class as
> `Salazar-Cavazos-2019` and `Erickson-2019/igf1r`. The 10-point deduction is the neg-bin/single-
> trajectory residual (a few early-count and peak points) and the 10-parameter fit's inherent
> sloppiness (parameter recovery by re-fitting is not identifiable, though the reproduction is solid).

Primary source (in the PyBNF repo; the authors' own PyBioNetFit forecasting setup):
- Paper: Lin YT, et al. "Daily Forecasting of Regional Epidemics of Coronavirus Disease with
  Bayesian Uncertainty Quantification, United States." *Emerg Infect Dis* 2021; 27(3):767–778.
  DOI 10.3201/eid2703.203364.
- Author files used: `~/Code/PyBNF/examples/COVID19forecasting_aMCMC/BigApple/` — `m1.bngl` (the
  two-phase model, **with the correct S0 = 19,216,182**), `m1.conf` (the legacy `am` /
  `neg_bin_dynamic` adaptive-MCMC sampler), `m1.exp` (daily NYT case data, t = 0..170), and
  `adaptive_files/MLE_params.txt` (the published MAP — 10 values).

"The paper's result" for this job = **the fit of the NYC MSA daily detected-case counts** with the
two-phase (n=1) model at the published MAP.

---

## Gate 0 — Materials inventory

| needed | present? | path / note |
|---|---|---|
| paper PDF | ✅ | Lin 2021 (`dev/papers/.../20-3364-combined.pdf`) |
| author model/job/data files | ✅ | `BigApple/` — `m1.bngl` + `m1.conf` + `m1.exp` |
| authors' fit result (MAP) | ✅ | `BigApple/adaptive_files/MLE_params.txt` (10 values) — **reproduces** (Gate 3) |

**Verdict:** PASS — a **complete, self-consistent** authors' PyBioNetFit setup (contrast the
single-phase DataS1 files, which had a placeholder S0). Gold-standard.

## Gate 1 — Data provenance

| `.exp` | source | method | units | diff vs. author file | verdict |
|---|---|---|---|---|---|
| `nyc_multiphase.exp` | authors' `m1.exp` | copied verbatim | daily new detected cases (integer counts) | byte-identical | PASS |

The data are the authors' MSA-aggregated NYT daily case counts (t = 0..170). The values at t=50..53
(106, 108, 202, 128) are identical to the single-phase `nyc.exp` — the **same NYC daily data**, just
extended to the full t=0..170 window. No value changed.

**Verdict:** PASS — byte-level provenance to the authors' own data file.

## Gate 2 — Model fidelity

Reference: the authors' own `m1.bngl`. Molecule types, seed species, observables, rate-law functions
(the two-phase `Lambda`/`P` switching), reaction rules, and fixed constants carried **verbatim**.
Edition-2 adaptations:

| aspect | authors' model | our `nyc_multiphase.bngl` | verdict |
|---|---|---|---|
| two-phase social distancing | σ=t0+t_delta, τ₁=t0+t_delta+t_delta2; (p0,λ0)/(p1,λ1) | identical | match |
| `counter()` seed / rate | `counter()=1`, `0->counter() 1` | identical (t = 1 + sim_time) | match |
| fixed constants | Lin-2021 region-independent | identical (verbatim) | match |
| free parameters | 10 `X__FREE` (9 model + r) | 9 declared as ids; `r` deleted (→ conf `r_disp`) | equiv |
| fit observable | `fDCs_Cum` (cumulative counter) | identical | match |
| **`S0`** | **`19,216,182` (correct census)** | identical | match |
| network cap | bare `generate_network` (finite) | none needed (~35 species / 61 reactions) | match |
| actions block | generate + simulate to t=186 | removed (synthesized from conf) | expected (edition-2) |

Independent check: `counter()` seeds at 1 and grows at rate 1 (observable `t` = 1 + sim time); the
two-phase switches fire at σ ≈ 32.9 and τ₁ ≈ 143.1 under both BNG2.pl and bngsim (CVode logs the
discontinuities during the fit). BNG2.pl generates 61 reactions in ~0.07 s.

**Verdict:** PASS — the only substantive change is the edition-2 free-parameter rename; the network
is the authors' own.

## Gate 3 — Reproduce at the paper's parameters (the headline)

- Params used: the 10 values in `BigApple/adaptive_files/MLE_params.txt` (= the `.bngl` nominals):
  t0=32.83, t_delta=0.073, t_delta2=110.22, beta=2.033, lambda0=0.098, p0=0.872, lambda1=5.181,
  p1=0.736, fD=0.115, r=10.27.
- Reproduction: `make_reproduction.py` → `nyc_multiphase_reproduction.png` (daily new detected
  cases, model vs. data, linear + log).

| metric | value |
|---|---|
| median \|rel err\|, window t=50..151 (same as single-phase slug) | **15.9 %** |
| median \|rel err\|, full t=0..170 (ey>0) | 18.4 % (inflated by tiny early-count points) |
| model peak | 14,435/day @ day 77 (2020-04-07) — data 15,925 @ day 73 (within 9 %) |
| cumulative t=50..151 | model 498,528 vs. data 487,332 (**+2 %**) |

On both scales the model tracks the data across ~5 orders of magnitude, and — the key win over the
single-phase model — it **captures the broad plateau** (days ~68–85, daily ~10–14k) instead of a
single sharp spike. The two social-distancing phases (a strong phase 1 from ~day 33, a milder phase
2 from ~day 143) are what let it hold the plateau and then track the long decline.

**Verdict:** PASS — the authors' own published MAP reproduces the NYC first wave. Gold-standard.

## Gate 3b — Why this matters for the single-phase `nyc/` slug

This slug is the **control** that diagnosed the single-phase slug's failure. The single-phase DataS1
model at its published MLE overshoots ~4× because (a) its shipped `S0` was a placeholder (Alabama's
4.9M, not NYC's 19.2M) and (b) one social-distancing phase cannot hold the NYC plateau. Here, with
the **correct S0** and the **two-phase** model, the published MAP reproduces at 16 %. So the `nyc/`
discrepancy is fully explained — a placeholder + model-complexity issue, not a data problem.

## Native-only verification (not PEtab)

`neg_bin` + mean-centering are outside the PEtab v2 subset → native-only (like `Erickson-2019/igf1r`):
1. **Tier-1** (`scripts/check_conf.py`): edition 2, `de` resolves, data bound, 10 free params by id,
   no `__FREE`. **PASS.**
2. **Native-only guard** (`scripts/petab_roundtrip.py`): `export_job` raises `NotImplementedError`.
   **PASS.**
3. **Real bngsim fit** (`pybnf -c nyc_multiphase.conf`, pop 60 × 80 + refine): converges to a
   finite objective (neg_bin nll **1018.99**); the `counter()`/phase-change switches fire under
   bngsim. Notably the **published MAP scores nll 1003.5** (full t=0..170) — *better* than this de
   run, i.e. the de search could not beat the published params in a modest budget, confirming the
   MAP is near-optimal. (The 10-parameter landscape is sloppy/multimodal; the reproduction at the
   published MAP, not re-fit parameter recovery, is the validation.)
4. **Reproduction:** Gate 3 (median 15.9 %).

The likelihood mirrors the authors' `neg_bin_dynamic`: a whole-fit
`noise_model = neg_bin, dispersion = fit r_disp, location = mean` + a per-observable override adding
`cumulative` on `fDCs_Cum`.

## Divergence & corrections

- Scope vs. paper: **matches** the two-phase NYC fit; a `de` optimizer replaces the authors'
  adaptive-MCMC sampler (as requested).
- Corrections: `r r__FREE` deleted (NB dispersion → conf `r_disp`); actions block stripped;
  `wall_time_sim = 30` added (a high-β trial at S0 ~19.2M is stiff at the phase discontinuities).
  No change to the science — S0 and every fixed constant are the authors' own.

## Bottom line

The strongest kind of job: the authors published a complete, self-consistent PyBioNetFit setup
(correct census S0, model, data, and a MAP that reproduces). This two-phase slug both adds a
gold-standard NYC COVID reproduction **and** explains the sibling single-phase slug's
non-reproduction (placeholder S0 + insufficient model complexity). Most valuable next step: none for
correctness; optionally add the adaptive-MCMC (`am`) UQ run or other MSAs' multi-phase fits.
