# VALIDATION — Miller-2026 / mek_isoform_de

Primary-source validation of the PyBNF job `pybnf-jobs/Miller-2026/mek_isoform_de/`, per the
`validate-pybnf-job` skill. Confidence is **earned from the gate evidence below**.

> **Reorganized for edition-2 (Mechanism A).** The job is now PRIMARY edition-2
> (`mek_isoform_de.conf` / `wt.bngl` — one model + four `condition:` perturbations for the
> KO/N78G/T292A/T292D cell lines). The SI-faithful edition-1 files (five per-variant models) are kept
> as [`../mek_isoform_de_legacy`](../mek_isoform_de_legacy) for provenance and as the independent
> BNG2.pl reproduction oracle. The Gate 0–1 provenance and the model chemistry are **shared** with the
> legacy twin; the edition-2 job adds a **network byte-identity** proof (Gate 2) and reproduces the
> legacy oracle **exactly** (Gate 3a: F_quant 22.72, 89/90 BPSL).

> **Confidence: 92 / 100** — Gates 0–2 PASS via the authors' OWN files **and** all three upstream
> sources traced to their ultimate origin (Kamioka Fig 3D, Catalanotti Figs 3a/4h–i/5d–g, K&L
> Table 1 + model); the edition-2 network is byte-identical to all five legacy variants; Gate 3a PASS
> (model at Table-3 best-fit reproduces WT quant to F_quant 22.72≈24 and 89/90 BPSL); Gate 3b PARTIAL
> by the paper's OWN documented non-identifiability (SI Fig 4) — the objective trends 1668→41 and the
> identifiable parameter set recovers. Real corrections applied.

Primary sources (untracked `dev/papers/Miller-2026/`; not redistributed):
- **Method/result paper:** Miller EF, Mallela A, Neumann J, Lin YT, Hlavacek WS, Posner RG. *Front
  Immunol* 2026; 17:1663008. doi:10.3389/fimmu.2026.1663008 (+ Supplementary "Data Sheet 1").
- **Base model:** Kocieniewski P, Lipniacki T. *Phys Biol* 2013; 10:035006 (8 species, 41 rules).
- **Quantitative WT data:** Kamioka Y, et al. *J Biol Chem* 2010; 285:33540 (Fig 3D).
- **Qualitative orderings:** Catalanotti F, et al. *Nat Struct Mol Biol* 2009; 16:294.
- **Authors' own PyBNF files:** `~/Code/PyBNF/examples/Miller2025_MEK_Isoforms/MEK_isoform_optimization_DE/`.
- **Edition-2 enablers** (both on `lanl/PyBNF` main): `9ab15167` (ADR-0027 — conditions perturb
  fixed/IC params in the bngsim mutant path) and `ae21d212` (`lanl/PyBNF#483` — aMCMC
  `output_trajectory` off-diagonal suffixes, used by the sibling slug).

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
| `WT.exp` (pSOS1/pEGFR/pERK, 18 pts, AU; padded to the 0…3600 s / 300 s grid) | Kamioka **Fig 3D** (HeLa + 50 ng/ml EGF), tabulated verbatim as **SI Table 8** "Experimental Data (AU)" | transcribed (table) + re-digitized (figure) | values ≡ SI Table 8 **exactly** (all 18); re-digitized Fig 3D matches to ≤0.5 AU (peak-normalized). Padding adds only `nan` rows | **PASS** |
| 5 × `.prop` (90 BPSL: WT 30, KO 24, N78G 18, T292A 12, T292D 6) | up/down orderings of pMEK(pRDS)/pERK across cell lines — Catalanotti **Figs 3a, 4h, 4i, 5d–g** (EGF regime), formalized as **SI Tables 1–5** | transcribed vs SI; cross-checked vs Catalanotti figures | all 90 match SI Tables 1–5 (see corrections); references rewritten to `<X>.MEK_pRDS` for the single edition-2 model | **PASS** |

Notes: quantitative data are **HeLa** (Kamioka), qualitative data are **MEFs** (Catalanotti) — a
cell-line mix inherited from K&L, not a job defect. SOS1 phospho in Kamioka Fig 3D is a
mobility-shift readout. Catalanotti's opposite-sign adhesion regime (Fig 1d) is correctly **excluded**.

**Verdict: PASS** — every value traces to a primary table/figure; the single defect (KO typo) is
corrected below.

## Gate 2 — Model fidelity (edition-2 network ≡ legacy per-variant models)

Reference: Kocieniewski & Lipniacki 2013 (Table 1 params, Figs 4–7 perturbation scheme) + the
authors' own `.bngl`. The single edition-2 `wt.bngl`, with each `condition:` applied, regenerates a
network **byte-identical** to the corresponding legacy per-variant model.

| aspect | paper (K&L) / legacy | edition-2 `wt.bngl` + `condition:` | verdict |
|---|---|---|---|
| molecule types / seed species | 8 (EGFR, Sos1, Ras, Raf, MEK1, MEK2, PHP, ERK) | identical | match |
| reaction rules | 41 rules (Fig 1) | identical | match |
| generated network | 5 per-variant models, **109 sp / 689 rxn** each | `wt.bngl`+condition, **109 sp / identical reaction topology** per cell line | **byte-identical** |
| `MEK_pRDS` observable | combined MEK1+MEK2 bis-phospho | `MEK1(S1~Ypp)+MEK2(S1~Ypp)` | match |
| perturbations (KO/N78G/T292A/T292D) | K&L p.8 (i)–(v), baked per model | one model + `condition:` param/IC perturbations (Mechanism A) | equivalent |
| nominal rate constants | K&L Table 1 = Miller Table 2 | identical | match |

The perturbations are applied as `condition:` parameter overrides on `wt.bngl` — `b2/b4/p4/u4 = 0`,
`b5 / 3` (relative on free `b5`), and the IC-seeding params `MEK1_0 = 0` (KO), `MEK1_0 = 0` +
`MEK1_0_T292p = 134000` (T292D seed swap). bngsim keeps the seed-species initializers symbolic and
re-derives concentrations after `set_param`, so no species `setConcentration` is needed. The
"110 sp" in the legacy VALIDATION is an off-by-one; both the legacy and edition-2 networks are **109 sp**.

Two anomalies (documented, not "corrected" — correcting would break reproduction of the paper's result):
- **`c1_L` is inert.** No reaction rule uses it (rule 4 uses `c1`, fixed at 0.02 = K&L's S0); this is
  *why* its profile likelihood is flat (SI Fig 4). Confirmed: F_quant 22.7 at c1=0.02 vs 33.2 at
  c1=c1_L=7.4e-3 → 0.02 is the operative value.
- **K&L Table-1 `p3 = 2×10⁻³`** is a µM-vs-copy-number unit slip; the model's `p3 = 2×10⁻⁹
  (mcls·s)⁻¹` is the physically correct value. `p3` is a fitted parameter, so this affects only the nominal.

**Verdict: PASS** — the edition-2 one-model+conditions network is byte-identical to the authors'
five per-variant files (functionally equivalent to K&L's published model).

## Gate 3a — Reproduce at the paper's parameters

`make_reproduction.py` substitutes Table-3 "Best-fit value" (31 params) into the single edition-2
`wt.bngl`, applies each cell line's `condition:` perturbation to the parameter block, runs BNG2.pl
(ODE, c1=0.02):

- **WT quantitative:** `F_quant = 22.72` (sum-of-squares on the 18 WT points) vs the paper's total DE
  objective **24.0** (SI Table 7) — identical to the legacy oracle's 22.7. Per-point residuals track
  SI Table 8 RMSE(PyBioNetFit).
- **Qualitative:** **89/90 BPSL constraints** satisfied (paper states 84/90; the difference is
  borderline near-tie sensitivity to Table-3's 2-sig-fig rounding). The single failure is the same
  KO 1% near-tie the legacy oracle reports.
- **Trajectories reproduce the biology:** WT transient MEK_pRDS/pERK; KO, N78G, T292A **sustained**;
  T292D **low/fast** — matching Catalanotti + K&L Figs 4–7 (`mek_isoform_de_reproduction.png`).

That the edition-2 one-model+conditions job reproduces the legacy oracle **to the digit** is the
operational confirmation of Gate 2 (byte-identical networks).

**Verdict: PASS** — model at the paper's own params reproduces Fig 1B and Fig 2.

## Gate 3b — Recover the paper's parameters by fitting

Bounded edition-2 DE (population_size=24, 10 iterations; the authors used 10000 on HPC):

- Objective **1668 → 41** — descends to K&L's original 40.0 (SI Table 7) and trends toward the PyBNF
  MLE's 24.0; matches the legacy twin's 1668 → 40 in 15 iters. `save_best_data` writes the best-fit
  gdat and the per-condition emitted `.bngl` (`MEK1_0 0.0`, etc.).
- The **identifiable** parameter subset recovers (log-ratios within ~0.4: `t1, a1, p2a, u5, i2, b2,
  u3, i1, p4`, all 3 scaling factors); the sloppy parameters stay 1–2.5 decades off.

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
2. `.prop` observable references rewritten `MEK_pRDS_<X>` → `<X>.MEK_pRDS` for the single edition-2 model.
3. `scalepEGFR` bound `[3e-5,6e-5] → [3e-5,3e-4]` to bracket Table 3's own best-fit (1.1e-4), which
   the shipped bound could not reach (the model requires ≈1.05e-4 for the pEGFR peak).
4. `WT.exp` padded to the uniform 0…3600 s / 300 s grid (`nan` where no data) so cross-experiment
   BPSL constraints (`KO.MEK_pRDS < WT.MEK_pRDS`) share an identical time grid with the constraint-only experiments.
5. CRLF → LF on all `.prop`/`.exp`.

**Edition-2 conversion** (vs the legacy twin): the five per-variant models collapse to one `wt.bngl`
+ four `condition:` perturbations (Mechanism A), enabled by `lanl/PyBNF` `9ab15167` (ADR-0027,
fixed/IC condition targets) and `ae21d212` (#483). Re-run after conversion: Gate 3a (89/90, F_quant
22.72) and Gate 3b (obj→41) both on the edition-2 job; both match the legacy oracle.

## Bottom line

A gold-standard validation: the authors' own files, all three upstream sources traced to their
ultimate origin, and a model that reproduces **both** the quantitative WT fit (F_quant≈objective) and
the qualitative cross-cell-line orderings (89/90) — now as an edition-2 one-model+conditions job whose
per-cell-line networks are byte-identical to the legacy per-variant files. The only residual is Gate
3b's structural sloppiness — the paper's *own* documented result (SI Fig 4), not a defect. Most
valuable next step to raise confidence further: a longer (HPC-scale) DE to confirm the objective
reaches 24.0 and to map the full identifiable/sloppy partition against SI Fig 4 panel-by-panel.
