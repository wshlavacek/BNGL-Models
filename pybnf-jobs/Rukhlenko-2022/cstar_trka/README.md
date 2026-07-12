# cstar_trka — cSTAR TrkA/NGF signalling fit (PyBNF edition-2 job)

A PyBNF edition-2, parameter-fitting job setup derived from:

> Rukhlenko OS, Halasz M, Rauch N, Zhernovkov V, Prince T, Wynne K, Maher S,
> Kashdan E, MacLeod K, Carragher NO, Kolch W, Kholodenko BN.
> **"Control of cell state transitions."** *Nature* 2022; **609**(7929):975–985.
> PMCID: [PMC9644236](https://pmc.ncbi.nlm.nih.gov/articles/PMC9644236/) ·
> DOI: [10.1038/s41586-022-05194-y](https://doi.org/10.1038/s41586-022-05194-y)

Built with the `curate-pybnf-job` skill. The paper's core signalling model was itself
fit with **pyBioNetFit**, so it is a natural real-world example.

> ⚠️ **Scope: reduced *demonstration* fit — NOT a reconstruction of the paper's fit.**
> This job frees **8 hand-picked parameters** over ±1-decade bounds against a **subset** of the
> data (the single DMSO/NGF arm), as a small, runnable PyBNF example built on the paper's
> **authentic model and data**. The paper's *actual* fit freed the full parameter set under
> **BMRA-derived confidence-interval inequality constraints**, across the full joint TrkA+TrkB ×
> multi-inhibitor dataset, with scatter search + simplex (Methods p.22). This demo therefore does
> **not** recover the published parameters (`VALIDATION.md` Gate 3b, PARTIAL/sloppy). A
> BMRA-constrained **real-world** fit is planned as a separate slug (see the paper-level README).

## What is fit

SH-SY5Y/TrkA human neuroblastoma cells stimulated with **NGF** (DMSO, no-inhibitor arm).
Seven phospho-protein **fold changes vs. t=0** at **0 / 10 / 45 min** are fit by the cSTAR
core signalling network (TRK → ERBB → ERK/AKT/JNK/S6K/RSK). The paper overlays this model on
data in **Fig. 4A** (pAKT, pERK, pJNK, pS6K, pRSK, pERBB) + **Fig. 4B** (pTRK), TrkA in blue.

**Fit method (paper).** The authors fit the model with **pyBioNetFit** (scatter search +
simplex, sum-of-squares objective) under **BMRA-inferred connection-coefficient confidence
intervals as inequality constraints** (Methods p.22). The **10-min** RPPA + Western Trk time
course was the *training* set; the **45-min** RPPA the *validation* set — so at the published
parameters the 45-min points fit less tightly than the 10-min ones (see `VALIDATION.md`).

**Experimental design.** Pre-equilibrate with no NGF (`Lig_on = 0`, receptors reach basal
phosphorylation), then add NGF (`Lig_on = 1`) and measure over 0–2700 s (= 45 min). PyBNF
synthesizes the two-phase protocol from `preequilibrate: basal` → `condition: stim`.

See **[`VALIDATION.md`](VALIDATION.md)** for the primary-source audit (confidence 86/100; model
byte-identical to the authors' `TrkA_S_model.bngl`, `.exp` byte-reproducible from their RPPA).

## Files

| file | role |
|---|---|
| `cstar_trka.bngl` | edition-2, fitting-ready model (no `actions` block); adapted from `TrkA_S_model.bngl` |
| `cstar_trka.conf` | the edition-2 job setup |
| `cstar_trka.exp` | fit target: 7 phospho fold changes at 0/10/45 min |
| `extract_exp.py` | reproducible extraction of `cstar_trka.exp` from the authors' RPPA CSVs |
| `cstar_trka_reproduction.png` | verification: model at published params vs. the RPPA data |

Source model + data: [github.com/OleksiiR/cSTAR_Nature](https://github.com/OleksiiR/cSTAR_Nature)
(`Trk_AB_models/TrkA_S_model.bngl`, `RPPA_DA/RPPA_data_trusted.csv`,
`RPPA_DA/RPPA_data_Trk_normalized_new.csv`).

## Model observable → RPPA antibody mapping

`FC_p*()` are model `functions` (they appear **with parentheses** in the `.exp` header).
Data are the TrkA + NGF, DMSO columns, replicate-averaged (r1–r6) and normalized to t=0.

| model fn | reads model observable | RPPA antibody (DMSO/NGF) |
|---|---|---|
| `FC_pTRK()`   | `pTRK`   | pTrk (Western; `RPPA_data_Trk_normalized_new.csv`) |
| `FC_pAdRTK()` | `pAdRTK` | ErbB-2/Her2/EGFR P Tyr1248/Tyr1173 |
| `FC_ppERK()`  | `ppERK`  | p44/42 MAPK (ERK1/2) P Thr202/Thr185 |
| `FC_pAKT()`   | `pAKT`   | Akt P Ser473 |
| `FC_pJNK()`   | `pJNKt`  | SAPK/JNK P Thr183 |
| `FC_pRSK()`   | `pRSK`   | p90 S6 kinase (Rsk1-3) P Thr359 |
| `FC_pS6K()`   | `pS6Kt`  | p70 S6 Kinase P Thr389 |

## Changes from the published model (all documented in `cstar_trka.bngl`)

- **No `begin actions` block** — the protocol is synthesized from the conf.
- **`Lig_on` ligand gate** (0/1) multiplies every ligand-association forward rate, replacing
  the published `setConcentration("TRKl(Rec)", 0 → 10)` equilibration idiom
  (same trick as `examples/real-world/receptor`).
- **All inhibitor concentrations = 0** — this is the DMSO (no-inhibitor) arm.
- Fitted constants are bare `id nominal` declarations; the conf's `*_var` free params bind
  to them **by name** (ADR-0034, no `__FREE`).

## Free parameters (8)

Core rate/coupling constants shaping the seven readouts, each `loguniform` over ≈±1 decade
around the published TrkA value:

| id | published | role |
|---|---|---|
| `vpTRK` | 2.356 | Trk trans-phosphorylation (→ pTRK) |
| `vpRTK` | 0.12 | ERBB/AdRTK phosphorylation (→ pAdRTK) |
| `g_A_TRKERK` | 220 | Trk → ERK (→ ppERK) |
| `vpERK` | 0.006 | ERK phosphorylation |
| `g_A_TRKAKT` | 3.68 | Trk → AKT (→ pAKT) |
| `vpAKT` | 0.0003 | AKT phosphorylation |
| `kpRSKA` | 0.1 | ERK → RSK (→ pRSK) |
| `g_A_AKTS6K` | 500 | AKT → S6K (→ pS6K) |

## ⚠️ Native-only (not PEtab-v2-exportable)

The data are **fold changes vs. t=0**, so the conf uses **`normalization = init`** — which
PyBNF's `export_job` cannot express in PEtab v2 (`NotImplementedError`), exactly like the
shipped `igf1r.conf`. This job is therefore verified **without** the PEtab round-trip.

## Verification (all pass)

- **Tier-1** (`skills/curate-pybnf-job/scripts/check_conf.py`): edition 2, `job_type=de`
  resolves, data bound, 8 free params bind by id, no `__FREE`. **PASS.**
- **Model build** (BNG2.pl): generates **82 species, 405 reactions** — matches the paper.
- **Bounded bngsim fit** (`de`, `max_iterations=2`, `population_size=6`): finite objective
  (best `sos ≈ 2.50`) in ~12 s — the simulate→score→propose loop runs (not `heavy`).
- **Paper reproduction** (`cstar_trka_reproduction.png`): at the published parameters the
  model reproduces Fig. 4A's dominant features — pTRK (≈5×) and pERK (≈3×) peak at ~10 min
  then decline, flat readouts stay flat. Objective at published params `sos ≈ 5.21`,
  median relative error ≈ 18 % (max ≈ 43 %) across the 10- and 45-min points; the bounded
  fit already improves on this single-arm subset (the published values were fit jointly to
  a much larger TrkA+TrkB × multi-inhibitor × 10/45-min dataset).

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Rukhlenko-2022/cstar_trka
pybnf -c cstar_trka.conf
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='cstar_trka', conf='cstar_trka.conf', simulator='ode',
    observables=('FC_pTRK', 'FC_pAdRTK', 'FC_ppERK', 'FC_pAKT',
                 'FC_pJNK', 'FC_pRSK', 'FC_pS6K'),
    system='cSTAR core signalling, TrkA/NGF phospho fold changes (Rukhlenko 2022, '
           'PMC9644236); ODE, pre-equilibrate -> NGF, 45 min; normalization=init '
           '-> NATIVE-ONLY (not PEtab-exportable)'),
    # native-only: assert export_job raises NotImplementedError instead of a PEtab-lint test.
```
