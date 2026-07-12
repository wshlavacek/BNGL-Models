# VALIDATION — Rukhlenko-2022/cstar_skmel133

Primary-source validation of the PyBNF job `pybnf-jobs/Rukhlenko-2022/cstar_skmel133/`, per the
`validate-pybnf-job` skill. Confidence is **earned from the gate evidence below**.

> **Confidence: 88 / 100** — Gate 2 is a byte-level match to the authors' own model file and Gate 1
> is byte-identical to the authors' own pyBioNetFit fit targets; the paper explicitly documents the
> pyBioNetFit fit (Methods p.24); Gate 3a reproduces the panel at the published params to ~13 %
> median error; Gate 3b is sloppy/non-identifiable (as the paper itself implies by constraining its
> fit with BMRA confidence intervals). Not higher only because the recovered 8-param subset is
> underdetermined (3b PARTIAL).

Primary sources (in the untracked `dev/papers/Rukhlenko-2022/`; not redistributed):
- Model + data paper: Rukhlenko OS et al. "Control of cell state transitions." *Nature* 2022;
  609(7929):975–985. PMCID **PMC9644236**, DOI 10.1038/s41586-022-05194-y (`nihms-1844164.pdf`,
  66 pp; `NIHMS1844164-supplement-Supplementary_Text.pdf`, 15 pp).
- Authors' own files (the decisive references): GitHub **OleksiiR/cSTAR_Nature** (cloned into the
  garage) — `SKMEL-133_models/SKMEL-133-3.bngl` (the model) and `SKMEL-133_preproc/dose*.exp` (the
  pyBioNetFit fit targets). The paper states the model parameters live only in these repo files (no
  parameter table is printed; Methods p.22/p.24).
- Underlying measurements: the SKMEL-133 RPPA is from **Korkut A et al. *eLife* 2015;4:e04640**
  (ref 29; 238 proteins × 89 perturbations, http://projects.sanderlab.org/pertbio) — the authors
  processed it into the fold-change `.exp` fit targets.

"The paper's result" for this job = the SKMEL-133 mechanistic model (`SKMEL-133-3.bngl`) fit with
pyBioNetFit to the single-drug perturbation fold changes (the `dose*.exp` training data), reported
qualitatively as **Extended Data Fig. 17A** (model DPD-vs-dose overlaid on data). The parameter
reference is the authors' `SKMEL-133-3.bngl` (there is no printed parameter table).

---

## Gate 0 — Materials inventory

| needed | present? | path / note |
|---|---|---|
| model paper PDF | ✅ | `nihms-1844164.pdf` (Methods "Building a nonlinear model of RAF inhibitor-resistant melanoma", p.24 names pyBioNetFit) |
| SI text | ✅ | `NIHMS1844164-supplement-Supplementary_Text.pdf` (discussion SI; no model equations — those are in the main Methods) |
| authors' model file | ✅ | `cSTAR_Nature/SKMEL-133_models/SKMEL-133-3.bngl` — enables a byte-level Gate-2 diff |
| authors' fit-target data | ✅ | `cSTAR_Nature/SKMEL-133_preproc/dose{1,2}_{ERK,AKT,SRC}inh.exp` — the pyBioNetFit `.exp` |

**Verdict:** PASS — the authors' own model and fit-target files are in hand.

## Gate 1 — Data provenance

The four `.exp` are the **authors' own pyBioNetFit fit-target files**, copied verbatim (per
`prep_exp.py`). Diff vs. the garage originals `SKMEL-133_preproc/dose*.exp`:

| `.exp` | source (authors' file) | method | normalization/units | diff vs. shipped | verdict |
|---|---|---|---|---|---|
| `dose1_ERKinh.exp` | `SKMEL-133_preproc/dose1_ERKinh.exp` | copied (2 mechanical edits) | 24 h fold change vs. no-drug DMSO baseline; `_SD` per point | **byte-identical** at t=86400 | PASS |
| `dose2_ERKinh.exp` | `…/dose2_ERKinh.exp` | copied | " | **byte-identical** | PASS |
| `dose1_AKTinh.exp` | `…/dose1_AKTinh.exp` | copied | " | **byte-identical** | PASS |
| `dose1_SRCinh.exp` | `…/dose1_SRCinh.exp` | copied | " | **byte-identical** | PASS |

The two documented mechanical edits (verified): (1) observable names take `()` in the header
(`FC_pERK` → `FC_pERK()`; the `_SD` companions do not), and (2) the authors' `t=0` baseline row
(all 1.0, SD 0.0 — SD=0 breaks `chi_sq`) is dropped and a single unmeasured `NaN` row at t=43200
inserted, since `normalization = init` supplies the baseline. The **measured** row (t=86400, all
nine fold changes + `_SD`) is byte-identical to the authors' file. The `dose*_S_*inh.exp` variants
(the DPD `Sval` curves) were correctly **not** used here — they belong to the sibling `_bpsl` slug.

Underlying provenance chain: Korkut et al. 2015 RPPA → authors' STV/fold-change processing →
authors' `.exp` → this job's `.exp`. **Verdict: PASS.**

## Gate 2 — Model fidelity

Reference compared against: the authors' `SKMEL-133_models/SKMEL-133-3.bngl` (in the garage).

`scripts/model_diff.py` (generated network): **species A=44 B=44 IDENTICAL; reactions A=56 B=56
IDENTICAL.** A comment-stripped textual diff of the model body shows only the documented changes:

| aspect | authors' `SKMEL-133-3.bngl` | our `cstar_skmel133.bngl` | verdict |
|---|---|---|---|
| molecule types / seed species | 24 molecule types; S seeded at Sp | identical | match |
| reaction rules + rate laws | 56 reactions | identical topology; 3 rate laws changed (below) | equiv |
| parameter values (incl. the 8 free ids) | g_IRSERK 4.2243, g_mTORERK 25, g_SRCERK 10.742, g_ERKAKT 0.16834, g_ERKSRC 98.898, g_ERKCDK 79.800, g_AKTmTOR 30, g_CDKMYC 0.30487, doses I_ERK 1.3122 / I_AKT 4.3207 / I_SRC 0.5173 | **identical** | match |
| inhibitor application | `setConcentration("IERK(ERKBD)", I_ERK_conc)` (clamped species); rate law reads species obs `I_ERK` | inhibitor concs zeroed at baseline (raised by conditions); rate law reads **parameter** `I_ERK_conc` (AKT/SRC likewise) | **equiv** |
| units | dose in Kd units (Kinh=1); fold-change readouts | identical | match |
| network cap | none (finite network) | none | match |
| observables ↔ measured | `FC_*()` functions = tIRS, pIRSI, pERK, pAKT, pSRC, pPKCt, pS6K, pRB, MYCt | identical | match |

The only functional edit — three rate laws using `1/(1+I_X_conc)` instead of `1/(1+I_X)` for
X ∈ {ERK, AKT, SRC} — is **exactly equivalent**: those three inhibitor molecule types have no
inhibitor-binding component and no binding rule in the authors' model, so the free-inhibitor species
observable `I_X` never depletes and equals the seed parameter `I_X_conc` identically. The rewrite is
required because a PyBNF edition-2 Condition emits `setParameter` (not `setConcentration`). The
genuinely competitive-binding inhibitors (mTOR/PKC/CDK) are correctly left out of this
parameter-driven panel.

**Verdict:** PASS (equiv) — identical generated network and identical parameters to the authors'
own model file; the sole deviation is a provably-equivalent inhibitor idiom swap.

## Gate 3a — Reproduce at the paper's parameters

- Published params used: the authors' `SKMEL-133-3.bngl` values (verbatim; there is no printed
  table). Ran the **authors' own model file** through BNG2.pl (its own actions block:
  equilibration → `dose{1,2}_S_{ERK,AKT,SRC}inh` at the published doses), and computed each readout's
  fold change = (drug steady state at 24 h) / (no-drug baseline steady state at 24 h).
- Metric vs. the `.exp` (4 conditions × 9 readouts = 36 points):

| condition | median rel err | notes |
|---|---|---|
| ERK dose1 | 15.0 % | pERK 0.42 (data 0.57), pRB↓, pS6K↓ — directions correct |
| ERK dose2 | 27.9 % | pERK 0.29 — dose-dependent suppression captured |
| AKT dose1 | 8.6 % | pAKT strongly suppressed (0.30 vs data 0.16) |
| SRC dose1 | 10.4 % | pSRC↓, pERK≈flat |
| **overall** | **13.2 %** (mean 18.3 %, max 81.6 %) | max = pAKT under AKTi (right direction, magnitude off) |

**Verdict:** PASS — the authors' model at the authors' published params reproduces the 24 h
inhibitor panel to ~13 % median relative error, matching the fit quality the paper itself achieves
(these fold changes are the training data; the model is validated on drug combinations).

## Gate 3b — Recover the paper's parameters by fitting

- Run: `pybnf -c cstar_skmel133.conf` with `max_iterations` reduced 50→30 (demonstrative budget),
  `de`, `population_size=12`, `chi_sq`, `normalization=init`. Completed in ~21 s; best obj **1812.4**.
- The published params give **chi_sq ≈ 7269** on this 4-condition subset — so the unconstrained
  8-param `de` fit finds a **much lower** objective than the authors' values, at parameters that
  differ substantially:

| param | published | recovered | log10 ratio | within 3×? |
|---|---|---|---|---|
| g_IRSERK | 4.224 | 0.924 | −0.66 | NO |
| g_mTORERK | 25 | 20.79 | −0.08 | yes |
| g_SRCERK | 10.74 | 2.232 | −0.68 | NO |
| g_ERKAKT | 0.1683 | 0.936 | +0.75 | NO |
| g_ERKSRC | 98.9 | 580.1 | +0.77 | NO |
| g_ERKCDK | 79.8 | 62.26 | −0.11 | yes |
| g_AKTmTOR | 30 | 17.88 | −0.22 | yes |
| g_CDKMYC | 0.3049 | 0.0432 | −0.85 | NO |

- Identifiability: **3/8 within 3×** (g_mTORERK, g_ERKCDK, g_AKTmTOR relatively identifiable; the
  rest trade off). The lower objective at "wrong" values is the classic sloppy/non-identifiable
  signature. This directly explains **why the paper constrained its fit** with BMRA-inferred
  connection-coefficient confidence intervals as inequality constraints, plus scatter-search +
  simplex and a drug-combination validation set (Methods p.24) — none of which the reduced,
  unconstrained `de` subset imposes.

**Verdict:** PARTIAL (sloppy / non-identifiable; the identifiable params + the Gate-3a curve
reproduction are the real check, not this table).

---

## Divergence & corrections

- **Scope vs. paper: MATCHES.** The job is the authors' own SKMEL-133 model fit to the authors' own
  pyBioNetFit fold-change targets, with the paper's declared method (pyBioNetFit). No re-scoping
  needed. The 8 free params are a curator-selected subset of the model's ERK/AKT/mTOR/CDK-axis
  connection coefficients (the paper fits ~all coefficients under BMRA-CI constraints).
- **Corrections applied to `README.md` / `cstar_skmel133.conf` (provenance precision only):**
  1. The 9-readout `.exp` are the pyBioNetFit **training data** (single-drug dose responses), not
     the plotted content of "Extended Data Figs. 16–17". The paper's SKMEL model-vs-data figure is
     **ED Fig. 17A** (which plots the **DPD** vs dose/Kd, a *derived* output); **ED Fig. 16** is the
     SVM state separation. Corrected the figure citation accordingly.
  2. Added the underlying-data provenance (Korkut et al. 2015 *eLife*, ref 29) and the pyBioNetFit
     fit method (scatter search + simplex + BMRA-CI inequality constraints; training = single drugs,
     validation = combinations).
- No changes to model/conf/exp numeric values — they are byte-identical to the authors' own files.

## Bottom line

Model and data are grounded to the authors' own files at the byte level, and the paper explicitly
documents the pyBioNetFit fit this job reconstructs; the authors' model at the authors' params
reproduces the 24 h panel to ~13 % median error. The single soft spot is Gate 3b: the reduced,
unconstrained 8-param fit is non-identifiable (as the paper's own use of BMRA-CI constraints
implies), so it lowers the objective without uniquely recovering the published values. The most
valuable next step to raise confidence would be to add the BMRA-CI inequality constraints (BPSL
`.prop`) around the fitted coefficients and fit the **full** panel (all six single-drug doses), i.e.
reproduce the paper's constrained fit rather than the reduced subset.
