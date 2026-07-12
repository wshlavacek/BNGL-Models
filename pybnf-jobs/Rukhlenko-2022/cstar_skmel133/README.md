# cstar_skmel133 — cSTAR SKMEL-133 melanoma inhibitor-panel fit (PyBNF edition-2 job)

Third cSTAR example, and a **different flavor** from the two Trk time-course fits: a
**steady-state inhibitor-perturbation panel** on a RAF-inhibitor-resistant melanoma line.

> Rukhlenko OS, Halasz M, Rauch N, Zhernovkov V, Prince T, Wynne K, Maher S,
> Kashdan E, MacLeod K, Carragher NO, Kolch W, Kholodenko BN.
> **"Control of cell state transitions."** *Nature* 2022; **609**(7929):975–985.
> PMCID: [PMC9644236](https://pmc.ncbi.nlm.nih.gov/articles/PMC9644236/) ·
> DOI: [10.1038/s41586-022-05194-y](https://doi.org/10.1038/s41586-022-05194-y)

> ⚠️ **Scope: reduced *demonstration* fit — NOT a reconstruction of the paper's fit.**
> This job frees **8 hand-picked parameters** over ±1-decade bounds against **4 of the panel's
> inhibitor conditions**, as a small, runnable PyBNF example built on the paper's **authentic
> model and data**. The paper's *actual* fit freed the full parameter set under **BMRA-derived
> confidence-interval inequality constraints**, across all six single-drug dose responses, with
> scatter search + simplex (Methods p.24). This demo therefore does **not** recover the published
> parameters (`VALIDATION.md` Gate 3b, PARTIAL/sloppy). A BMRA-constrained **real-world** fit is
> planned as a separate slug (see the paper-level README).

## What is fit

**SKMEL-133** is a BRAF-V600E / PTEN-null, **RAF-inhibitor-resistant** melanoma line. Unlike
the Trk models (ligand time courses), the design is a **steady-state inhibitor perturbation
panel**: equilibrate the network to its proliferative steady state (no drug), then apply a
targeted-kinase inhibitor and read out the new steady state at **24 h** (86400 s). Nine
phospho / total-protein **fold changes vs. the no-drug baseline** are fit per perturbation.
The DPD force is **bistable** (Sd = arrest/death, Sp = proliferation); combined inhibition
drives the switch Sp → Sd.

These nine-readout fold changes are the authors' **pyBioNetFit training data** (single-drug
dose responses; Methods p.24: "the model-generated dose responses were fitted to these
training set data" with the objective = sum of squares, using scatter search + simplex and
BMRA-inferred connection-coefficient confidence intervals as inequality constraints). The
paper's *published* SKMEL-133 model-vs-data figure is **Extended Data Fig. 17A**, which plots
the resulting **DPD** (S) vs. inhibitor dose in units of Kd — a derived output — not the nine
individual readouts; **Extended Data Fig. 16** is the SVM state separation. The underlying RPPA
measurements are from **Korkut A et al., *eLife* 2015;4:e04640** (238 proteins × 89
perturbations), which the authors processed into the fold-change targets.

**Perturbations:** ERK inhibitor (dose1 + dose2), AKT inhibitor (dose1), SRC inhibitor
(dose1). dose1 = the published `I_<x>_conc`; dose2 = 2×.

The fit data are the **authors' own pyBioNetFit `.exp` files** (from
`SKMEL-133_preproc/`), so no digitization was needed. Each carries per-point `_SD`
→ **`chi_sq`** objective. See **[`VALIDATION.md`](VALIDATION.md)** for the primary-source
audit (confidence 88/100; model + data byte-identical to the authors' files).

## Files

| file | role |
|---|---|
| `cstar_skmel133.bngl` | edition-2, fitting-ready model (no `actions` block); adapted from `SKMEL-133-3.bngl` |
| `cstar_skmel133.conf` | the edition-2 job setup (4 inhibitor experiments) |
| `dose{1,2}_ERKinh.exp`, `dose1_AKTinh.exp`, `dose1_SRCinh.exp` | authors' 24 h fold-change targets (9 readouts + `_SD`) |
| `prep_exp.py` | documents the two mechanical edits that adapt the authors' files to edition-2 |
| `cstar_skmel133_reproduction.png` | verification: model at published params vs. the RPPA data |

Source model + data: [github.com/OleksiiR/cSTAR_Nature](https://github.com/OleksiiR/cSTAR_Nature).

## Fit observables (9 model functions)

`FC_tIRS`, `FC_IRSI`, `FC_pERK`, `FC_pAKT`, `FC_pSRC`, `FC_pPKC`, `FC_pS6K`, `FC_pRB`,
`FC_MYC` — total IRS, inhibitory-phospho IRS, and the phospho forms of ERK, AKT, SRC, PKC,
S6K, RB, plus total MYC. In the `.exp` header the functions take parentheses
(`FC_pERK()`), while their noise companions do **not** (`FC_pERK_SD`) — PyBNF forms the
noise column as `<entity>_SD`.

## Changes from the published model (all documented in `cstar_skmel133.bngl`)

- **No `begin actions` block** — the equilibrate → perturb → measure protocol is
  synthesized from `preequilibrate: dmso` → `condition: <inh>`.
- **All inhibitor concentrations = 0** (the DMSO baseline). Each perturbation condition
  raises one inhibitor to its published dose.
- **Inhibitors applied as parameters, not clamped species.** ERK/AKT/SRC inhibitors act
  *only* through the rate-law factor `1/(1+I_<x>)` with **no** sequestering binding, so the
  free-inhibitor observable equals the seed parameter `I_<x>_conc` **exactly**. The rate
  laws were rewritten to use `I_<x>_conc` directly, so a PyBNF Condition (which emits
  `setParameter`) actually applies the drug. (mTOR/PKC/CDK inhibitors act by competitive
  binding, so they are left out of this parameter-driven panel.)
- Fitted constants are bare `id nominal`; the conf's `*_var` free params bind by name
  (ADR-0034, no `__FREE`).

## Free parameters (8)

ERK/AKT/mTOR/CDK-axis connection coefficients (bare-id declarations), each `loguniform`
over ≈±1 decade around the published value: `g_IRSERK` (4.22), `g_mTORERK` (25),
`g_SRCERK` (10.74), `g_ERKAKT` (0.168), `g_ERKSRC` (98.9), `g_ERKCDK` (79.8),
`g_AKTmTOR` (30), `g_CDKMYC` (0.305).

## ⚠️ Native-only (not PEtab-v2-exportable)

Fold-change data → `normalization = init` → not PEtab v2-expressible (like the shipped
`igf1r.conf`). Verified **without** the PEtab round-trip.

## Verification (all pass)

- **Tier-1** (`check_conf.py`): edition 2, `job_type=de` resolves, data bound, 8 free params
  bind by id, no `__FREE`. **PASS.**
- **Model build** (BNG2.pl): generates **44 species, 56 reactions** (small, fast).
- **Bounded bngsim fit** (`de`, `max_iterations=2`, `population_size=6`, `chi_sq`): finite
  objective (`≈ 2945` at random params) — the simulate→score→propose loop runs across all
  four inhibitor experiments (not `heavy`).
- **Paper reproduction** (`cstar_skmel133_reproduction.png`): at the published parameters the
  model reproduces the authors' 24 h fold changes across the panel — dose-dependent pERK
  suppression (ERK dose1→dose2: 0.42 → 0.29), strong pAKT loss under AKT inhibition, and
  proliferation markers (pRB, pS6K) down — at **≈ 13 % median relative error** across
  4 inhibitors × 9 readouts.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Rukhlenko-2022/cstar_skmel133
pybnf -c cstar_skmel133.conf
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='cstar_skmel133', conf='cstar_skmel133.conf', simulator='ode',
    observables=('FC_tIRS', 'FC_IRSI', 'FC_pERK', 'FC_pAKT', 'FC_pSRC',
                 'FC_pPKC', 'FC_pS6K', 'FC_pRB', 'FC_MYC'),
    system='cSTAR core signalling, SKMEL-133 RAF-inhibitor-resistant melanoma; '
           'steady-state kinase-inhibitor perturbation panel (Rukhlenko 2022, PMC9644236); '
           'ODE, pre-equilibrate -> inhibitor, 24 h, chi_sq; normalization=init '
           '-> NATIVE-ONLY (not PEtab-exportable)'),
    # native-only: assert export_job raises NotImplementedError instead of a PEtab-lint test.
```
