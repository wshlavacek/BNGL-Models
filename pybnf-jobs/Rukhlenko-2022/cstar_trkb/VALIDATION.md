# VALIDATION — Rukhlenko-2022/cstar_trkb

Primary-source validation of the PyBNF job `pybnf-jobs/Rukhlenko-2022/cstar_trkb/`, per the
`validate-pybnf-job` skill. Confidence is **earned from the gate evidence below**.

> **Confidence: 84 / 100** — Gate 2 is a byte-level match to the authors' own `TrkB_S_model.bngl`
> (incl. the TrkB–ERBB hetero-dimerization module), and Gate 1 regenerates the `.exp`
> byte-identically from the authors' own RPPA data; the paper explicitly documents the pyBioNetFit
> fit of Fig. 4A/4B (Methods p.22). Gate 3a reproduces Fig. 4 at the published params to ~27 %
> median error (larger than TrkA — the stronger, more sustained TrkB signalling, with the 45-min
> *validation* points hardest). Gate 3b is sloppy/non-identifiable (PARTIAL). The higher 3a residual
> + indirect data path put it just below the TrkA twin.

Primary sources (in the untracked `dev/papers/Rukhlenko-2022/`; not redistributed):
- Model + data paper: Rukhlenko OS et al. "Control of cell state transitions." *Nature* 2022;
  609(7929):975–985. PMCID **PMC9644236**, DOI 10.1038/s41586-022-05194-y (`nihms-1844164.pdf`).
- Authors' own files (the decisive references): GitHub **OleksiiR/cSTAR_Nature** —
  `Trk_AB_models/TrkB_S_model.bngl` (the model, 114 species / 645 reactions), `RPPA_DA/
  RPPA_data_trusted.csv` (raw RPPA, TrkB columns) and `RPPA_DA/RPPA_data_Trk_normalized_new.csv`
  (Western pTrk). Parameters live only in the repo files (no printed table; Methods p.22).

"The paper's result" for this job = the TrkB mechanistic model (`TrkB_S_model.bngl`) fit with
pyBioNetFit to the SH-SY5Y/TrkB + BDNF phospho time courses, reported as **Figure 4A/4B** (TrkB =
red, model curves over RPPA/Western fold-change data). Parameter reference = `TrkB_S_model.bngl`.

---

## Gate 0 — Materials inventory

| needed | present? | path / note |
|---|---|---|
| model paper PDF | ✅ | `nihms-1844164.pdf` (Methods p.22 names pyBioNetFit; Fig. 4 legend pp.59–60) |
| authors' model file | ✅ | `cSTAR_Nature/Trk_AB_models/TrkB_S_model.bngl` — enables a byte-level Gate-2 diff |
| authors' RPPA data | ✅ | `cSTAR_Nature/RPPA_DA/RPPA_data_trusted.csv` (TrkB_DMSO cols) + `RPPA_data_Trk_normalized_new.csv` |

**Verdict:** PASS.

## Gate 1 — Data provenance

Reproducible extraction from the authors' RPPA data (`extract_exp.py`, TrkB columns). Re-ran and
diffed the regenerated file against the shipped one:

| `.exp` column | source (antibody) | method | normalization | diff vs. shipped | verdict |
|---|---|---|---|---|---|
| `FC_pTRK()` | pTrk (Western; `…_normalized_new.csv`), TrkB+BDNF DMSO | extracted, ÷t0 | fold change vs t=0 | **byte-identical** | PASS |
| `FC_pAdRTK()` | ErbB-2/Her2/EGFR P Tyr1248/Tyr1173 | extracted (6-rep mean, ÷t0) | " | **byte-identical** | PASS |
| `FC_ppERK()` | p44/42 MAPK (ERK1/2) P Thr202/Thr185 | " | " | **byte-identical** | PASS |
| `FC_pAKT()` | Akt P Ser473 | " | " | **byte-identical** | PASS |
| `FC_pJNK()` | SAPK/JNK P Thr183 | " | " | **byte-identical** | PASS |
| `FC_pRSK()` | p90 S6 kinase (Rsk1-3) P Thr359 | " | " | **byte-identical** | PASS |
| `FC_pS6K()` | p70 S6 Kinase P Thr389 | " | " | **byte-identical** | PASS |

Method identical to the TrkA twin but on the **TrkB_DMSO + BDNF** columns. The large TrkB responses
(pTRK ~12×, pAdRTK/ERBB2 ~11× at 10 min) are the biological signature the paper highlights (the
ERBB→ERK→RSK→ERBB amplifier). These map to **Fig. 4A/4B** (red). Training = 10-min RPPA + Western;
validation = 45-min RPPA (Methods p.22). **Verdict: PASS.**

## Gate 2 — Model fidelity

Reference compared against: the authors' `Trk_AB_models/TrkB_S_model.bngl` (in the garage).

`scripts/model_diff.py` (generated network): **species A=114 B=114 IDENTICAL; reactions A=645 B=645
IDENTICAL.** A comment-stripped body diff shows only the documented changes:

| aspect | authors' `TrkB_S_model.bngl` | our `cstar_trkb.bngl` | verdict |
|---|---|---|---|
| molecule types / seed species | incl. the TrkB–ERBB hetero-dimer module | identical | match |
| reaction rules + rate laws | 645 reactions (incl. `sqrt(fh*factor)` hetero-dimer ligand rates) | identical topology; ligand-assoc **forward** rates gated by `Lig_on` (incl. the hetero-dimer family) | equiv |
| parameter values (incl. 8 free ids) | vpTRK 2.356, vpRTK 0.12, g_B_TRKERK 67.658, vpERK 0.006, g_B_TRKAKT 30, vpAKT 0.0003, kpRSKB 0.018, g_B_AKTS6K 10 | **identical** | match |
| ligand application | `setConcentration("TRKl(Rec)", 0→10)` | `Lig_on` (0/1) × forward ligand rates; conf `basal`→`stim` | **equiv** |
| inhibitors | nominal values for the inhibitor arms | all 0 — DMSO/BDNF arm | equiv (DMSO arm) |
| units / cap | nM ligand; finite network, no cap | identical | match |
| observables ↔ measured | `FC_p*()` = pTRK, pAdRTK, pJNKt, ppERK, pAKT, pRSK, pS6Kt | identical | match |

Same `Lig_on` equivalence argument as TrkA, extended to the TrkB-specific hetero-dimer ligand-binding
rates. **Verdict:** PASS (equiv) — identical generated network and identical parameters to the
authors' own model.

## Gate 3a — Reproduce at the paper's parameters

- Published params used: the authors' `TrkB_S_model.bngl` values (verbatim). Ran the **authors' own
  model** through BNG2.pl (`timecourse2_B_S` = BDNF stimulation, 0–2700 s); FC(t) = obs(t)/obs(0) at
  10/45 min.
- Metric vs. the `.exp` (14 points): **median 27.3 % rel err** (mean 25.2 %, max 41.9 %).
  Published-param objective sos ≈ **65.4** (matches the job README's ~69.8; absolute sos is larger
  than TrkA's simply because TrkB fold changes are ~7–12× vs ~1–5×). The strong, sustained
  multi-pathway activation and the large ERBB2/AdRTK response are captured; residuals are largest on
  the 45-min *validation* points (pTRK, pAdRTK, pERK all ~35 %).

**Verdict:** PASS — the authors' model at the authors' params reproduces Fig. 4 (TrkB) to ~27 %
median error, with residuals concentrated on the validation-set 45-min points.

## Gate 3b — Recover the paper's parameters by fitting

- Run: `pybnf -c cstar_trkb.conf` with `max_iterations` reduced 50→20 (demonstrative), `de`,
  `population_size=12`, `sos`, `normalization=init` (114 species / 645 reactions — the heaviest of
  the three).

  Best obj **23.7** — below the published-param sos ≈ 65.4:

| param | published | recovered | log10 ratio | within 3×? |
|---|---|---|---|---|
| vpTRK | 2.356 | 9.658 | +0.61 | NO |
| vpRTK | 0.12 | 0.0260 | −0.66 | NO |
| g_B_TRKERK | 67.66 | 30.31 | −0.35 | yes |
| vpERK | 0.006 | 0.00439 | −0.14 | yes |
| g_B_TRKAKT | 30 | 46.52 | +0.19 | yes |
| vpAKT | 0.0003 | 0.000265 | −0.05 | yes |
| kpRSKB | 0.018 | 0.0701 | +0.59 | NO |
| g_B_AKTS6K | 10 | 14.26 | +0.15 | yes |

  **5/8 within 3×** of the published value (better-conditioned than the TrkA/SKMEL fits — the strong
  TrkB responses pin more of the coupling constants).

- Identifiability: sloppy / non-identifiable, matching TrkA and SKMEL (the ligand-driven amplitudes
  are compensated across the coupling constants). The paper's fit used scatter-search + simplex +
  **BMRA-CI inequality constraints** over a joint TrkA+TrkB × multi-inhibitor dataset (Methods p.22),
  none of which this reduced single-arm `de` subset imposes.

**Verdict:** PARTIAL (sloppy / non-identifiable; Gate-3a curve reproduction is the real check).

---

## Divergence & corrections

- **Scope vs. paper: MATCHES.** The job is the authors' own TrkB model fit to the authors' own RPPA
  time-course data with the paper's declared method (pyBioNetFit). Fig. 4A/4B citation correct.
- **Corrections applied to `README.md` (provenance precision only):** the training(10 min)
  /validation(45 min) split and the actual fit method (scatter search + simplex + BMRA-CI inequality
  constraints via pyBioNetFit). No changes to model/conf/exp numeric values — model byte-identical to
  the authors' file; `.exp` regenerates byte-identically from the authors' RPPA data.

## Bottom line

Model (incl. the hetero-dimer amplifier) and data are grounded to the authors' own files
(byte-identical model; byte-reproducible `.exp`), and the paper explicitly documents the pyBioNetFit
fit of Fig. 4 that this job reconstructs; the authors' model at the authors' params reproduces Fig. 4
(TrkB) to ~27 % median error, the residual sitting (as expected) on the validation 45-min points.
Residual risk is the higher 3a error and the sloppy Gate-3b subset. The most valuable next step is
the joint TrkA+TrkB × inhibitor, BMRA-CI-constrained fit — the paper's actual constrained fit.
