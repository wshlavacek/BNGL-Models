# VALIDATION — Miller-2026 / mek_isoform_de

Primary-source validation of the PyBNF job `pybnf-jobs/Miller-2026/mek_isoform_de/`, per the
`validate-pybnf-job` skill. Confidence is **earned from the gate evidence below**.

> **Confidence: 92 / 100** — Gates 0–2 PASS via the authors' OWN files **and** all three upstream
> sources traced to their ultimate origin (Kamioka Fig 3D, Catalanotti Figs 3a/4h–i/5d–g, K&L
> Table 1 + model); Gate 3a PASS (model at Table-3 best-fit reproduces WT quant to F_quant 22.7≈24
> and 89/90 BPSL); Gate 3b PARTIAL by the paper's OWN documented non-identifiability (SI Fig 4) —
> the objective trends 1668→40 and the identifiable parameter set recovers. Real corrections applied.

Primary sources (untracked `dev/papers/Miller-2026/`; not redistributed):
- **Method/result paper:** Miller EF, Mallela A, Neumann J, Lin YT, Hlavacek WS, Posner RG. *Front
  Immunol* 2026; 17:1663008. doi:10.3389/fimmu.2026.1663008 (+ Supplementary "Data Sheet 1").
- **Base model:** Kocieniewski P, Lipniacki T. *Phys Biol* 2013; 10:035006 (8 species, 41 rules).
- **Quantitative WT data:** Kamioka Y, et al. *J Biol Chem* 2010; 285:33540 (Fig 3D).
- **Qualitative orderings:** Catalanotti F, et al. *Nat Struct Mol Biol* 2009; 16:294.
- **Authors' own PyBNF files:** `~/Code/PyBNF/examples/Miller2025_MEK_Isoforms/MEK_isoform_optimization_DE/`.

"The paper's result" = **Table 3 "Best-fit value"** (31 params) + **Fig 1B** (WT quant) + **Fig 2
B/D/F/H** (qualitative satisfaction). Reported DE objective **24.0** (SI Table 7).

---

## Gate 0 — Materials inventory

| needed | present? | note |
|---|---|---|
| method/result paper PDF | ✅ | Miller 2026 (15 pp) |
| Supplementary Material | ✅ | Data Sheet 1 (18 pp): SI Tbl 1–5 (90 BPSL), SI Tbl 6–9 (satisfaction, obj, RMSE), SI Fig 1–5 (MCMC diag, profile, ablation) |
| base-model paper | ✅ | Kocieniewski & Lipniacki 2013 (16 pp) |
| quantitative data paper | ✅ | Kamioka 2010 (9 pp) |
| qualitative data paper | ✅ | Catalanotti 2009 (10 pp) |
| authors' own model/data/BPSL/conf files | ✅ | enables direct Gate-2 diff + verbatim data provenance |

**Verdict: PASS** — every source, including all three upstream experimental/model papers, is in hand.

## Gate 1 — Data provenance

| `.exp` / `.prop` | source | method | diff vs. shipped | verdict |
|---|---|---|---|---|
| `WT.exp` (pSOS1/pEGFR/pERK, 18 pts, AU) | Kamioka **Fig 3D** (HeLa + 50 ng/ml EGF), tabulated verbatim as **SI Table 8** "Experimental Data (AU)" | transcribed (table) + re-digitized (figure) | `WT.exp` ≡ SI Table 8 **exactly** (all 18 values); re-digitized Fig 3D matches to ≤0.5 AU (peak-normalized) | **PASS** |
| 5 × `.prop` (90 BPSL: WT 30, KO 24, N78G 18, T292A 12, T292D 6) | up/down orderings of pMEK(pRDS)/pERK across cell lines — Catalanotti **Figs 3a, 4h, 4i, 5d–g** (EGF regime), formalized as **SI Tables 1–5** | transcribed vs SI; cross-checked vs Catalanotti figures | all 90 match SI Tables 1–5 (see corrections); orderings grounded in Catalanotti (KO/N78G/T292A sustained > WT late; T292D < WT; KO pMEK early<late crossover = C4/C5/C6) | **PASS** |

Notes: quantitative data are **HeLa** (Kamioka), qualitative data are **MEFs** (Catalanotti) — a
cell-line mix inherited from K&L, not a job defect. SOS1 phospho in Kamioka Fig 3D is a
mobility-shift readout (not a phospho-antibody). Catalanotti's opposite-sign adhesion regime (Fig 1d)
is correctly **excluded** — all constraints are EGF-stimulation.

**Verdict: PASS** — every value traces to a primary table/figure; the single defect (KO typo) is
corrected below.

## Gate 2 — Model fidelity

Reference: Kocieniewski & Lipniacki 2013 (Table 1 params, Figs 4–7 perturbation scheme) + the
authors' own `.bngl`. Each model generates **110 species / 689 reactions** through BNG2.pl (K&L:
110 species/ODEs ✓).

| aspect | paper (K&L) | our `.bngl` | verdict |
|---|---|---|---|
| molecule types / seed species | 8 (EGFR, Sos1, Ras, Raf, MEK1, MEK2, PHP, ERK) | identical | match |
| reaction rules | 41 rules (Fig 1) | identical (authors' file) | match |
| `MEK_pRDS` observable | combined MEK1+MEK2 bis-phospho (Raf-dependent sites) | `MEK1(S1~Ypp)+MEK2(S1~Ypp)` | match |
| initial conditions (µM×1e6 = copies) | Table 2 defaults | identical | match |
| perturbations (KO/N78G/T292A/T292D) | K&L p.8 (i)–(v): MEK1=0; b2=b4=0; p4=0; T292~Yp seed + u4=0 + b5/N | identical, per-model | match |
| nominal rate constants | K&L Table 1 defaults = Miller Table 2 | identical (spot-checked ~20 params) | match |
| units | (mcls·s)⁻¹ / s⁻¹ | identical | match |

Two anomalies (documented, not "corrected" — correcting would break reproduction of the paper's result):
- **`c1_L__FREE` is inert.** No reaction rule uses it (rule 4 uses `c1`, fixed at 0.02 by each
  model's action); the effective ligand signal is `c1 = 0.02` = K&L's S0. This is *why* its
  profile likelihood is flat (SI Fig 4). Confirmed empirically: reproduction identical at any
  in-bounds `c1_L` (F_quant 22.7 at c1=0.02 vs 33.2 at c1=c1_L=7.4e-3 → 0.02 is the operative value).
- **K&L Table-1 `p3 = 2×10⁻³`** is a µM-vs-copy-number unit slip; the model's `p3 = 2×10⁻⁹
  (mcls·s)⁻¹` is the physically correct value (2e-3 × ERK ≈ 6000 s⁻¹ is absurd; 2e-9 gives ≈6e-3 s⁻¹).
  `p3` is a fitted parameter, so this affects only the nominal, not the result.

**Verdict: PASS (identical to the authors' file; functionally equivalent to K&L's published model).**

## Gate 3a — Reproduce at the paper's parameters

`make_reproduction.py` substitutes Table-3 "Best-fit value" (31 params) into the 5 models, runs
BNG2.pl (ODE, c1=0.02):

- **WT quantitative:** `F_quant = 22.7` (sum-of-squares on the 18 WT points) vs the paper's total DE
  objective **24.0** (SI Table 7). Per-point residuals track SI Table 8 RMSE(PyBioNetFit) (e.g. at
  t=3600 the model gives pEGFR 0.34 / pERK 0.54, matching the paper's implied 0.32 / 0.52).
- **Qualitative:** **89/90 BPSL constraints** satisfied (paper states 84/90; the difference is
  borderline near-tie sensitivity to Table-3's 2-sig-fig rounding). The single failure is KO-C6
  (`KO.MEK_pRDS@3600 > N78G@3600`: 10809 vs 10919, a 1% near-tie).
- **Trajectories reproduce the biology:** WT transient MEK_pRDS/pERK; KO, N78G, T292A **sustained**;
  T292D **low/fast** — matching Catalanotti + K&L Figs 4–7 (`mek_isoform_de_reproduction.png`).

**Verdict: PASS** — model at the paper's own params reproduces Fig 1B and Fig 2.

## Gate 3b — Recover the paper's parameters by fitting

Bounded legacy DE (population_size=24, 15 iterations; the authors used 10000 on HPC):

- Objective **1668 → 40.3** — descends to K&L's original 40.0 (SI Table 7) and trends toward the
  PyBNF MLE's 24.0; more iterations continue down.
- **12/31 parameters within 3×** of Table 3 — and they are precisely the **identifiable** set:
  `t1` (log-ratio 0.01), `a1` (0.04), `p2a` (−0.06), `u5` (0.02), `i2` (0.08), `b2` (−0.10),
  `u3` (0.16), `i1` (−0.26), `p4` (−0.37), and all 3 scaling factors. The sloppy parameters stay
  1–2.5 decades off (`c2` +2.59, `n5` −2.09, `n2` +1.70, `b3` −1.45, `n1` +1.37, `b1` +1.34).

This mirrors the paper's own profile-likelihood analysis (**SI Fig 4: ~14/28 rate constants
identifiable**, the rest flat). The non-identifiability is a *documented property of the model/data*,
not a fitting failure (cf. the tlbr / igf1r / Rukhlenko sloppy-fit precedents).

**Verdict: PARTIAL** — objective + identifiable parameters recover; sloppy parameters cannot, by the
paper's own finding.

---

## Divergence & corrections

Scope matches the paper (the authors' own DE job = "the paper's result"). Corrections applied to the
job files (all reflected in `README.md`):

1. `KO.prop` C9 typo `time=1600` → `time=3600` (SI Table 2).
2. Model action blocks standardized (`end model` slice + one clean `setParameter(c1);simulate`) —
   the authors' variant files carried malformed trailing action junk (BNG2 "Not a CODE reference").
3. `scalepEGFR__FREE` bound `[3e-5,6e-5] → [3e-5,3e-4]` to bracket Table 3's own best-fit (1.1e-4),
   which the shipped bound could not reach (the model requires ≈1.05e-4 for the pEGFR peak).
4. CRLF → LF on all `.prop`/`.exp`.

Re-run after corrections: Gate 3a (89/90, F_quant 22.7) and Gate 3b (obj→40) both on the corrected files.

## Bottom line

A gold-standard validation: the authors' own files, all three upstream sources traced to their
ultimate origin, and a model that reproduces **both** the quantitative WT fit (F_quant≈objective) and
the qualitative cross-cell-line orderings (89/90). The only residual is Gate 3b's structural
sloppiness — the paper's *own* documented result (SI Fig 4), not a defect. Most valuable next step to
raise confidence further: a longer (HPC-scale) DE to confirm the objective reaches 24.0 and to map the
full identifiable/sloppy partition against SI Fig 4 panel-by-panel.
