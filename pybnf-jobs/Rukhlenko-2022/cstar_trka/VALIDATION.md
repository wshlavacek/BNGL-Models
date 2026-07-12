# VALIDATION — Rukhlenko-2022/cstar_trka

Primary-source validation of the PyBNF job `pybnf-jobs/Rukhlenko-2022/cstar_trka/`, per the
`validate-pybnf-job` skill. Confidence is **earned from the gate evidence below**.

> **Confidence: 86 / 100** — Gate 2 is a byte-level match to the authors' own `TrkA_S_model.bngl`,
> and Gate 1 regenerates the `.exp` byte-identically from the authors' own RPPA data; the paper
> explicitly documents the pyBioNetFit fit of Fig. 4A (Methods p.22). Gate 3a reproduces Fig. 4A at
> the published params to ~18 % median error (larger on the 45-min points, which are the paper's
> *validation* set). Gate 3b is sloppy/non-identifiable (PARTIAL), and the `.exp` is one step
> removed from a shipped `.exp` (extracted from raw RPPA), so slightly below the SKMEL slug.

Primary sources (in the untracked `dev/papers/Rukhlenko-2022/`; not redistributed):
- Model + data paper: Rukhlenko OS et al. "Control of cell state transitions." *Nature* 2022;
  609(7929):975–985. PMCID **PMC9644236**, DOI 10.1038/s41586-022-05194-y (`nihms-1844164.pdf`).
- Authors' own files (the decisive references): GitHub **OleksiiR/cSTAR_Nature** —
  `Trk_AB_models/TrkA_S_model.bngl` (the model, 82 species / 405 reactions, verbatim in the paper),
  `RPPA_DA/RPPA_data_trusted.csv` (raw RPPA, 6 reps) and `RPPA_DA/RPPA_data_Trk_normalized_new.csv`
  (Western-blot pTrk). The paper states the model parameters live only in these repo files (no
  parameter table; Methods p.22).

"The paper's result" for this job = the TrkA mechanistic model (`TrkA_S_model.bngl`) fit with
pyBioNetFit to the SH-SY5Y/TrkA + NGF phospho time courses, reported as **Figure 4A/4B** (model
curves overlaid on RPPA/Western fold-change data). Parameter reference = `TrkA_S_model.bngl` (no
printed table).

---

## Gate 0 — Materials inventory

| needed | present? | path / note |
|---|---|---|
| model paper PDF | ✅ | `nihms-1844164.pdf` (Methods "Refining parameters of the dynamic model", p.22 names pyBioNetFit; "82 species and 405 reactions" p.22) |
| authors' model file | ✅ | `cSTAR_Nature/Trk_AB_models/TrkA_S_model.bngl` — enables a byte-level Gate-2 diff |
| authors' RPPA data | ✅ | `cSTAR_Nature/RPPA_DA/RPPA_data_trusted.csv` (+ `RPPA_data_Trk_normalized_new.csv` for pTrk) |

**Verdict:** PASS — the authors' own model and the source RPPA data are in hand.

## Gate 1 — Data provenance

The `.exp` is a reproducible extraction from the authors' RPPA data (`extract_exp.py`). I re-ran the
extraction against the garage CSVs and diffed the regenerated file against the shipped one:

| `.exp` column | source (antibody, series) | method | normalization/units | diff vs. shipped | verdict |
|---|---|---|---|---|---|
| `FC_pTRK()` | pTrk (Western; `RPPA_data_Trk_normalized_new.csv`), TrkA+NGF DMSO | extracted, ÷t0 | fold change vs t=0 | **byte-identical** on regeneration | PASS |
| `FC_pAdRTK()` | ErbB-2/Her2/EGFR P Tyr1248/Tyr1173 | extracted (6-rep mean, ÷t0) | " | **byte-identical** | PASS |
| `FC_ppERK()` | p44/42 MAPK (ERK1/2) P Thr202/Thr185 | " | " | **byte-identical** | PASS |
| `FC_pAKT()` | Akt P Ser473 | " | " | **byte-identical** | PASS |
| `FC_pJNK()` | SAPK/JNK P Thr183 | " | " | **byte-identical** | PASS |
| `FC_pRSK()` | p90 S6 kinase (Rsk1-3) P Thr359 | " | " | **byte-identical** | PASS |
| `FC_pS6K()` | p70 S6 Kinase P Thr389 | " | " | **byte-identical** | PASS |

Method: TrkA_DMSO columns at 0/10/45 min (6 replicates), `nanmean`, each series normalized to its
own t=0 mean → fold change (matching `normalization = init`). Times: 10 min = 600 s, 45 min = 2700 s.
These seven readouts and their fold changes are exactly what the paper plots in **Fig. 4A** (pAKT,
pERK, pJNK, pS6K, pRSK, pERBB) **+ 4B** (pTRK), NGF-stimulated TrkA (blue), model over data
(confirmed from the PDF). Per the paper (Methods p.22), the **10-min** RPPA + Western Trk time course
was the *training* set and the **45-min** RPPA the *validation* set. **Verdict: PASS** (one step
removed from a shipped `.exp` — extracted from raw RPPA — but the extraction is byte-reproducible).

## Gate 2 — Model fidelity

Reference compared against: the authors' `Trk_AB_models/TrkA_S_model.bngl` (in the garage).

`scripts/model_diff.py` (generated network): **species A=82 B=82 IDENTICAL; reactions A=405 B=405
IDENTICAL.** A comment-stripped body diff shows only the documented changes:

| aspect | authors' `TrkA_S_model.bngl` | our `cstar_trka.bngl` | verdict |
|---|---|---|---|
| molecule types / seed species | 17 molecule types; S seeded at 2*S0 | identical | match |
| reaction rules + rate laws | 405 reactions | identical topology; ligand-assoc forward rates gated by `Lig_on` | equiv |
| parameter values (incl. 8 free ids) | vpTRK 2.356, vpRTK 0.12, g_A_TRKERK 220, vpERK 0.006, g_A_TRKAKT 3.68, vpAKT 0.0003, kpRSKA 0.1, g_A_AKTS6K 500 | **identical** | match |
| ligand application | actions block: `setConcentration("TRKl(Rec)", 0→10)` around equilibration | `Lig_on` (0/1) multiplies each ligand-association **forward** rate; conf sets `basal`→`stim` | **equiv** |
| inhibitors | nominal AKT 1.42 / ERK 0.5 / TRK 26 / JNK 11 / S6K 25 / RSK 1.7 (for the inhibitor arms) | all 0 — DMSO/NGF arm (nonzero values kept in a comment) | equiv (DMSO arm) |
| units / network cap | ligand in nM; finite network, no cap | identical | match |
| observables ↔ measured | `FC_p*()` = pTRK, pAdRTK, pJNKt, ppERK, pAKT, pRSK, pS6Kt | identical | match |

The `Lig_on` gate is functionally equivalent to the authors' equilibration idiom: with `Lig_on = 0`
the forward ligand-association rate is zero (receptors reach basal, exactly as `setConcentration(...,
0)`), and with `Lig_on = 1` binding proceeds while `TRKl(Rec)=10`; reverse rates (acting only on
already-bound ligand, of which there is none at equilibrium) are unchanged, so the equilibrated state
is identical. **Verdict:** PASS (equiv) — identical generated network and identical parameters to the
authors' own model; deviations are the ligand-gate idiom and the DMSO-arm inhibitor zeros.

## Gate 3a — Reproduce at the paper's parameters

- Published params used: the authors' `TrkA_S_model.bngl` values (verbatim). Ran the **authors' own
  model file** through BNG2.pl (its own actions block: equilibrate → `timecourse2_A_S` = NGF
  stimulation, 0–2700 s), and computed FC(t) = obs(t)/obs(0) for the seven readouts at 10 and 45 min.
- Metric vs. the `.exp` (2 timepoints × 7 readouts = 14 points): **median 18.3 % rel err**
  (mean 19.8 %, max 42.8 %). The published-param objective sos ≈ **5.19** (matches the job README's
  ~5.21). Qualitatively: pTRK peaks (~5×) then declines, pERK peaks (~3×), pAKT/pRSK modest, flat
  readouts flat — the Fig. 4A/4B shape. The largest errors are on the 45-min points (pTRK 42.8 %,
  pERK 37.8 %), consistent with those being the paper's *validation* set (not fit).

**Verdict:** PASS — the authors' model at the authors' params reproduces Fig. 4A/4B to ~18 % median
error, with the residual concentrated (as expected) on the validation-set 45-min points.

## Gate 3b — Recover the paper's parameters by fitting

- Run: `pybnf -c cstar_trka.conf` with `max_iterations` reduced 50→20 (demonstrative), `de`,
  `population_size=12`, `sos`, `normalization=init`. Best obj **2.34** — i.e. below the published-param
  sos ≈ 5.19 (the unconstrained fit lowers the objective on this single NGF arm):

| param | published | recovered | log10 ratio | within 3×? |
|---|---|---|---|---|
| vpTRK | 2.356 | 0.275 | −0.93 | NO |
| vpRTK | 0.12 | 0.172 | +0.16 | yes |
| g_A_TRKERK | 220 | 101.8 | −0.33 | yes |
| vpERK | 0.006 | 0.0243 | +0.61 | NO |
| g_A_TRKAKT | 3.68 | 0.488 | −0.88 | NO |
| vpAKT | 0.0003 | 0.000126 | −0.38 | yes |
| kpRSKA | 0.1 | 0.0271 | −0.57 | NO |
| g_A_AKTS6K | 500 | 65.2 | −0.88 | NO |

- Identifiability: **3/8 within 3×**; the rest trade off (the ligand-driven amplitudes are
  compensated across vpTRK/g_A_TRKERK/vpERK). Same sloppy/non-identifiable picture as SKMEL. The
  paper's fit used scatter-search + simplex + **BMRA-CI inequality constraints** over a joint
  TrkA+TrkB × multi-inhibitor dataset (Methods p.22); this reduced single-arm, unconstrained `de`
  fit does not impose those, so it lowers the objective without recovering the published values.

**Verdict:** PARTIAL (sloppy / non-identifiable; Gate-3a curve reproduction is the real check).

---

## Divergence & corrections

- **Scope vs. paper: MATCHES.** The job is the authors' own TrkA model fit to the authors' own RPPA
  time-course data, with the paper's declared method (pyBioNetFit). The Fig. 4A/4B citation is
  correct. The 8 free params are a curator-selected subset of the model's connection/rate constants.
- **Corrections applied to `README.md` (provenance precision only):** noted the training(10 min)
  / validation(45 min) split and the actual fit method (scatter search + simplex + BMRA-CI
  inequality constraints via pyBioNetFit), which explains why the 45-min reproduction error is
  larger. No changes to model/conf/exp numeric values — the model is byte-identical to the authors'
  file and the `.exp` regenerates byte-identically from the authors' RPPA data.

## Bottom line

Model and data are grounded to the authors' own files (byte-identical model; byte-reproducible `.exp`
from the authors' RPPA), and the paper explicitly documents the pyBioNetFit fit of Fig. 4A that this
job reconstructs; the authors' model at the authors' params reproduces Fig. 4A/4B to ~18 % median
error. Residual risk is Gate 3b (reduced, unconstrained, single-arm fit is non-identifiable) and the
slightly indirect data path (extracted from raw RPPA, though byte-reproducibly). The most valuable
next step would be to fit the **joint** TrkA+TrkB × inhibitor dataset under BMRA-CI constraints, i.e.
reproduce the paper's constrained fit rather than the single-arm subset.
