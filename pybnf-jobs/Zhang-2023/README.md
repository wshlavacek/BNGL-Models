# Zhang-2023 — endothelial VEGFR2/CD47 signaling and angiogenesis-driven tumor growth (PyBNF fitting jobs)

PyBNF edition-2 parameter-fitting jobs derived from the coupled signaling/tumor-growth model of:

> **Zhang Y, Popel AS, Bazzazi H.**
> **"Combining Multikinase Tyrosine Kinase Inhibitors Targeting the Vascular Endothelial Growth
> Factor and Cluster of Differentiation 47 Signaling Pathways Is Predicted to Increase the
> Efficacy of Antiangiogenic Combination Therapies."** *ACS Pharmacol Transl Sci* 2023;
> **6**:710–726.
> DOI: [10.1021/acsptsci.3c00008](https://doi.org/10.1021/acsptsci.3c00008).
> — the 870-species endothelial signaling model (Supporting Information File S2), the
> two-state tumor growth model (File S4 / File S3), the parameter values (Table S1, Table S2),
> and the calibration and validation figures (Fig. 4B–4D).

Built with the `curate-pybnf-job` skill. The **library-model** siblings — the same two models
with actions blocks, for reference simulation rather than fitting — live in
[`models/endothelial_vegfr2_and_cd47_signaling_zhang2023/`](../../models/endothelial_vegfr2_and_cd47_signaling_zhang2023/),
together with the verification notebook that documents and re-runs every digitization used
here.

## The models

The paper couples two models.

- **Endothelial signaling** (File S2): VEGF165a binding to VEGFR2, VEGFR1 and neuropilin-1;
  ligand-dependent and ligand-independent receptor dimerization; CD47 preassociation with
  VEGFR2 and its capture by thrombospondin-1; trafficking through two signaling endosomes and
  a recycling endosome; and downstream PLCγ/calcium, PKC–Raf–MEK–ERK with a
  sphingosine-1-phosphate positive feedback, Src–Axl–PI3K–AKT, and eNOS. Network generation
  gives **870 species and 5602 reactions** — the 870-ODE system of the paper.
- **Angiogenesis-driven tumor growth** (File S4): the normalized peak endothelial pERK1/2 and
  ppAKT drive a Hill-type activation of ribosomal protein S6, a delayed driver `TD` relaxes
  toward it, and tumor volume grows exponentially, softly capped at a linear rate, gated by
  `TD`, and shrunk by a first-order kill term.

Only the second is fit here. The signaling model has no free parameters against any data in
this paper: its 180 constants come from the authors' earlier publications, and the only
signaling data reproduced (Fig. 4B) is a four-bar normalized comparison used for validation,
not calibration.

## The jobs

| slug | run cost | fits | flavor | data | status |
|---|---|---|---|---|---|
| [`tumor_growth`](tumor_growth/) | `minutes` | the six growth parameters the Methods name as fitted — `w_OR`, `kTD`, `EC50TD`, `kg`, `klinear`, `kkill` | quantitative, `sos`, **PEtab-exportable** | 34 points digitized from Fig. 4D (Bridgeman et al. 2016 xenograft), four arms | ✅ tier-1 + PEtab round-trip + fit + reproduction · **75/100** ([VALIDATION](tumor_growth/VALIDATION.md)) |

## What this job pins down

> ⚠️ **No published parameter set reproduces Fig. 4D.** The paper reports the tumor growth
> parameters twice and the two reports disagree on five of them. File S4's gate
> (`kTD ≈ 3`, `EC50TD ≈ 2.8e-4`) is saturated in every arm, so the per-capita growth rate has
> the same sign in every arm and no `kkill` can grow the untreated tumor while shrinking the
> combination-treated one. Table S1's gate (`kTD = 29.28`, `EC50TD = 0.584`) does separate the
> arms but Table S1 reports no `kkill` at all, and File S4 sets it to 0. Measured against the
> digitized data: File S4 SSE **93.0**, Table S1 with `kkill = 0` SSE **2622.6**, this fit SSE
> **5.39** (RMSE 0.40 fold, all 34 points).

> ⚠️ **The fit lands on Table S1's gate**, which is the strongest available evidence for which
> of the two published sets is the operative one: `kTD` **30.5** vs. the reported 29.28,
> `EC50TD` **0.573** vs. 0.584, `w_OR` **0.34** vs. 0.3, with `kg` **0.139** matching File S4's
> 0.146. The one parameter neither source reports, `kkill`, comes out at **0.0264 /day**.

> ⚠️ **The paper's own three-arm calibration protocol does not work on these data.** Fitting
> the untreated, sunitinib and trametinib arms and holding out the combination arm — what the
> Methods describe — fits its own arms slightly better (`sos` 1.91 vs. 2.70) and then predicts
> the held-out arm **growing 2.15-fold** where the data shrink to 0.35-fold. All three
> calibration arms grow, so all three sit on the same side of the `TD` gate and none of them
> locates it; only the shrinking arm does. The holdout run is reported in
> [VALIDATION.md](tumor_growth/VALIDATION.md) Gate 4, and this job fits all four arms.

> ⚠️ **BioNetGen's rate-law semantics are load-bearing here.** The authors write the growth
> step as the rule `TD() -> Cells + TD()`, so BioNetGen multiplies the growth law by the
> reactant `TD` on top of the law's own `TD` gate. The Methods equation does not show that
> factor; the authors' own SBML export does (File S3, `rateLaw2 * S2`). Dropping it makes even
> the untreated arm miss.

## Source materials

- **Primary paper:** `dev/papers/Zhang2023/ACSPharmTranslSci.pdf`.
- **Authors' models:** `dev/papers/Zhang2023/pt3c00008_si_002/` — `S1.xml`/`S2.bngl` (signaling,
  SBML and BNGL) and `S3.xml`/`S4.bngl` (tumor growth, SBML and BNGL).
- **Authors' parameter tables:** `dev/papers/Zhang2023/pt3c00008_si_001.pdf` — Table S1 (all
  model parameters) and Table S2 (initial values).
- **Fit data:** Bridgeman VL et al., *Mol Cancer Ther* 2016; 15:172–183,
  DOI [10.1158/1535-7163.MCT-15-0170](https://doi.org/10.1158/1535-7163.MCT-15-0170) — mouse
  xenograft renal cell carcinoma volumes under sunitinib, trametinib and their combination,
  replotted as the open circles of Fig. 4D. Digitized here; the extraction is committed in the
  library-model folder's `verify_zhang2023.ipynb`.

Not built (optional future slugs): a fit of the endothelial signaling model itself, which would
need the upstream data of the authors' earlier papers rather than this one; a bootstrap
resampling job reproducing the practical-identifiability analysis of Fig. S1.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl

cd pybnf-jobs/Zhang-2023/tumor_growth
pybnf -c tumor_growth.conf        # differential evolution over 6 parameters, ~2 min
python make_reproduction.py       # figure + metrics: File S4 vs Table S1 vs this fit
```

`make_reproduction.py` simulates through BNG2.pl, so the figure reproduces without the fitting
toolchain.
