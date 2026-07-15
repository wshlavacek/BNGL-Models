# VALIDATION — Lin-2021/nyc

Primary-source validation of the PyBNF job `pybnf-jobs/Lin-2021/nyc/`.
Confidence is **earned from the gate evidence below**.

> **Confidence: 78 / 100** — the model is the authors' own compartmental model (Gate 2), the
> data are the authors' own NYT-derived daily case counts (Gate 1), and a PyBNF `de` fit
> reproduces the NYC first-wave epidemic across ~5 orders of magnitude (median 16 % rel err,
> Gate 3). This is **not** a gold-standard authors-reproduce-exactly port (unlike
> `Salazar-Cavazos-2019` or `Erickson-2019/igf1r`): the shipped DataS1 model files carry a
> **placeholder `S0` and `sigma`**, so the shipped `MLE_params.txt` does **not** reproduce the
> data (it overshoots ~4×, Gate 3b). `S0` was corrected to NYC's census population and the
> reproduction uses an **independent PyBNF `de` fit** rather than the published MLE. The
> deductions (−22) are that inconsistency plus the single-phase model's structural residual (a
> sharper/earlier peak than the data's broad plateau).

Primary sources (in the untracked `dev/papers/Lin-2021-Mallela-2023/`; not redistributed):
- Model + data paper: Mallela A, Lin YT, Hlavacek WS. "Differential contagiousness of
  respiratory disease across the United States." *Epidemics* 2023; 45:100718.
  DOI 10.1016/j.epidem.2023.100718 (open access).
- Model lineage: Lin YT, et al. "Daily Forecasting of Regional Epidemics of Coronavirus Disease
  with Bayesian Uncertainty Quantification, United States." *Emerg Infect Dis* 2021;
  27(3):767–778. DOI 10.3201/eid2703.203364.
- Author files used: `DataS1/PyBioNetFit/New York-Newark-Jersey City, NY-NJ-PA/` —
  `…​.bngl` (the single-phase model, **with placeholder S0/σ**), `…​.conf` (the legacy `am` /
  `neg_bin_dynamic` adaptive-MCMC sampler), `…​.exp` (daily NYT case data, t = 50..151), and
  `…_output/adaptive_files/MLE_params.txt` (the shipped MLE — see Gate 3b). A byte-identical
  model copy also lives at `~/Code/PyBNF/examples/Mallela2022MSAs/New York-…​/`.

"The paper's result" for this job = **the fit of the NYC MSA daily detected-case counts over
the first wave** (model day t = 50..151, i.e. 2020-03-11..2020-06-20; day 0 = 2020-01-21).

---

## Gate 0 — Materials inventory

| needed | present? | path / note |
|---|---|---|
| model + data paper PDF | ✅ | Mallela 2023 (`1-s2.0-S1755436523000543-main.pdf`) + Lin 2021 (`20-3364-combined.pdf`) |
| author model/job files | ✅ | `DataS1/PyBioNetFit/New York-…/` — `.bngl` + `.conf` + `.exp` |
| authors' fit result | ⚠️ | `adaptive_files/MLE_params.txt` present **but does not reproduce** (Gate 3b) |

**Verdict:** PASS with a caveat — the authors shipped a complete PyBioNetFit *setup*, but it is
a **job-setup template** with un-substituted per-MSA fixed parameters, not a completed fit that
round-trips.

## Gate 1 — Data provenance

| `.exp` | source | method | units | diff vs. author file | verdict |
|---|---|---|---|---|---|
| `nyc.exp` | authors' `…​.exp` | transcribed verbatim | daily new detected cases (integer counts) | values identical; only the function-observable header gained parens (`CumNum_detected_cases_Cum()`) | PASS |

The data are **daily new detected cases** (non-monotonic — rises to 15,925 at t≈73 = 2020-04-03,
falls to ~669 at t=151), the authors' MSA-aggregated NYT county case counts. No value changed.

**Verdict:** PASS — byte-level provenance to the authors' own data file.

## Gate 2 — Model fidelity

Reference: the authors' own `…​.bngl`. The molecule types, seed species, observables, rate-law
functions, reaction rules, and fixed epidemiological constants (Mallela 2023 Table 2) are carried
**verbatim**. Edition-2 fitting adaptations:

| aspect | authors' model | our `nyc.bngl` | verdict |
|---|---|---|---|
| compartments / rules / functions | SEIR-type + `~M`/`~P`/`~Q` + `counter()` time | identical | match |
| fixed epidemiological constants | Table 2 values | identical (verbatim) | match |
| free parameters | 6 `X__FREE` (t0, beta, lambda0, p0, fD, r) | 5 declared as ids; `r` deleted (it is only the NB dispersion → conf `r_disp`) | equiv |
| observable ↔ measured quantity | `CumNum_detected_cases_Cum()` (a function) | identical function | match |
| network cap | bare `generate_network` (finite) | none needed (finite, **31 species / 62 reactions**) | match |
| actions block | generate + simulate to t=152 | removed (synthesized from conf) | expected (edition-2) |
| **`S0`** | **`4903185` (placeholder)** | **`19216182` (corrected — see Gate 3b)** | **corrected** |
| `sigma` | `63` | `63` (kept; ≈ 2020-03-24 NY PAUSE) | kept |

Independent checks (BNG2.pl and bngsim):
- **Counter/time mechanism:** `counter()` == simulation time exactly (`max|t−counter| = 0`); the
  `if(t≥sigma)` social-distancing switch fires at day 63 (protected `SP` = 0 before t=63, grows
  after) under **both** BNG2.pl and the in-process bngsim backend (CVode logs the t=63
  discontinuity during the fit). This resolves the one real edition-2 risk (the `if()`/counter
  behaviour under bngsim).

**Verdict:** PASS on the reaction network; the one substantive change beyond the edition-2 rename
is correcting the placeholder `S0` (next gate).

## Gate 3a — Reproduce the fit (independent `de`, corrected S0)

- Nominals in `nyc.bngl` = an independent PyBNF `de`+refine fit (`nyc.conf`, population 30 ×
  60 iterations): `t0=26.633, beta=0.538, lambda0=9.160, p0=0.592, fD=0.883, r_disp=9.887`
  (obj = neg_bin nll **839.88**).
- Reproduction: `make_reproduction.py` → `nyc_reproduction.png` (daily new detected cases, model
  vs. data, linear + log).

| metric (fit window t=50..151, 102 pts) | value |
|---|---|
| median \|rel err\| | **16.3 %** (max 203 % at the earliest low-count points) |
| model peak | 22,492/day @ day 67 (2020-03-28) |
| data peak | 15,925/day @ day 73 (2020-04-03) |
| cumulative over window | model 554,464 vs. data 487,332 (+14 %) |

On the **log scale** the model tracks the data across ~5 orders of magnitude — the exponential
growth (t=50..63), the peak, and the long declining tail (t=80..151). On linear scale the
**single-phase (n=0)** model fits a sharper, ~6-days-earlier peak than the data's broad
plateau (days 70–90); that is the model's structural limit, not a setup error (the paper selects
n=0/1/2 per MSA; NYC's broad plateau is what a second social-distancing phase would capture).

**Verdict:** PASS — the fit reproduces the NYC first-wave epidemic; residuals are the
single-phase peak/plateau mismatch.

## Gate 3b — The shipped MLE does NOT reproduce (the headline finding)

The QUICKSTART recon assumed `MLE_params.txt` is a published MAP that reproduces the data (as in
`egfr_simpull`). **It is not**, and the cause is a **placeholder in the shipped DataS1 files**:

- The shipped `.bngl` carries `S0 = 4,903,185` — **Alabama's 2019 population**, *identical across
  every MSA* (verified: NYC, Detroit, Appleton all ship the same S0), and `sigma = 63` shared too.
  The paper states S0 was set **per-MSA from census data** and σ as the **per-MSA day cumulative ≥
  200** (Mallela 2023 Sec 2.5) — neither is in the shipped files.
- With the shipped placeholders, the published MLE (t0=27.37, β=0.671, λ₀=0.105, p0=0.752,
  fD=0.262, r=11.21) scores **neg_bin nll ≈ 6.1 million** on the NYC data (median rel err 108 %,
  peak 70,428/day = **4.4× the data**). The same test on other MSAs is likewise bad (Detroit 676 %
  rel err, Appleton 94 %) — **systematic, not NYC-specific.**
- An (S0, σ) sensitivity sweep at the published MLE finds **no** combination that reproduces:
  even at NYC's real census S0 (19.2M) with an earlier σ, the best case is nll ≈ 2,300 / 65 % rel
  err — still far worse than the `de` optimum (~840 / 16 %). The MLE params are simply inconsistent
  with the DataS1 `.exp` under this model.

The scoring is trustworthy: our independent neg_bin nll computation reproduces PyBNF's own
objective exactly (860.64 at the pre-refine `de` best), so these numbers are PyBNF's, not an
artifact of a re-implementation.

**Resolution (per Bill):** correct `S0` to NYC's 2019 census MSA population (**19,216,182**), keep
`sigma = 63` (= 2020-03-24, matching NY's 2020-03-22 "PAUSE" stay-at-home order; the pre-window
cumulative needed to derive the ≥200 day exactly is not in `nyc.exp`), and use the independent
`de` fit (Gate 3a) as the reproduction rather than the non-reproducing published MLE.

**Verdict:** the published MLE is **not reproducible** from the shipped files; documented and
worked around, not hidden.

---

## Native-only verification (not PEtab)

`neg_bin` (and mean-centering) were removed from PEtab v2, so this job is **native-only** (like
`Erickson-2019/igf1r`). Verified accordingly:

1. **Tier-1** (`scripts/check_conf.py`): edition 2, `job_type=de` resolves, `nyc.exp` bound, 6
   free params (5 model ids + `r_disp` nuisance) bind by id, no `__FREE`. **PASS.**
2. **Native-only guard** (`scripts/petab_roundtrip.py --job-type de`): `export_job` raises
   `NotImplementedError` (mean-centering / neg_bin outside the PEtab-exportable subset). **PASS**
   (do NOT expect a round-trip).
3. **Real bngsim fit** (`pybnf -c nyc.conf`): converges to a finite objective (neg_bin nll
   839.88); the `counter()`/day-63 switch fires correctly under bngsim. **PASS.**
4. **Reproduction:** Gate 3a (median 16.3 % rel err).

The conf's likelihood is the edition-2 successor to the authors' `neg_bin_dynamic`: a whole-fit
`noise_model = neg_bin, dispersion = fit r_disp, location = mean` (mean centering is MANDATORY —
the ODE gives the mean count) **plus** a per-observable override adding the `cumulative` flag that
differences the running `CumNum_detected_cases_Cum()` into daily incidence. (`cumulative` is
per-observable ONLY — a whole-fit `cumulative` is rejected by the parser, so the QUICKSTART's
single-line block was corrected to this two-line form.)

## Divergence & corrections

- Scope vs. paper: **matches** the single-phase per-MSA fit, but a `de` optimizer replaces the
  authors' adaptive-MCMC sampler (as requested — cheaper than the Bayesian posterior).
- Corrections applied: (1) `S0` 4,903,185 → 19,216,182 (placeholder → NYC census); (2) nominals =
  independent `de` fit, not the non-reproducing MLE; (3) `r r__FREE` deleted (NB dispersion →
  conf `r_disp`); (4) actions block stripped; (5) `wall_time_sim = 30` added (a high-β trial at the
  census S0 is very stiff at the σ discontinuity and must be failed fast, not left to grind).
- Re-run after setup: tier-1 PASS, native-only guard PASS, `de` fit + reproduction as above.

## Bottom line

An **honest reconstructed fit**, not a gold-standard reproduction. The model and data are the
authors' own, the network and time mechanism are verified, and a `de` fit reproduces the NYC
first-wave epidemic across ~5 orders of magnitude (median 16 %). The load-bearing caveat — worth
surfacing to the authors — is that the **shipped DataS1 PyBioNetFit files use a placeholder `S0`
(and `σ`), so the shipped `MLE_params.txt` does not reproduce any MSA's data**; this job corrects
`S0` and substitutes an independent fit. Most valuable next step: obtain the authors' real per-MSA
`S0`/`σ` (and the multi-phase n=1/2 config) to test whether the published MLE then reproduces.
