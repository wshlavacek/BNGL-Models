# cstar_trkb — cSTAR TrkB/BDNF signalling fit (PyBNF edition-2 job)

**Run cost: `minutes`** — 600 evaluations (12 × 50 `de`), 8 free parameters.

Twin of [`../cstar_trka`](../cstar_trka), for the **TrkB / BDNF** variant of the cSTAR
core signalling network. Same paper, same method, same verification discipline.

> Rukhlenko OS, Halasz M, Rauch N, Zhernovkov V, Prince T, Wynne K, Maher S,
> Kashdan E, MacLeod K, Carragher NO, Kolch W, Kholodenko BN.
> **"Control of cell state transitions."** *Nature* 2022; **609**(7929):975–985.
> PMCID: [PMC9644236](https://pmc.ncbi.nlm.nih.gov/articles/PMC9644236/) ·
> DOI: [10.1038/s41586-022-05194-y](https://doi.org/10.1038/s41586-022-05194-y)

> ⚠️ **Scope: reduced *demonstration* fit — NOT a reconstruction of the paper's fit.**
> This job frees **8 hand-picked parameters** over ±1-decade bounds against a **subset** of the
> data (the single DMSO/BDNF arm), as a small, runnable PyBNF example built on the paper's
> **authentic model and data**. The paper's *actual* fit freed the full parameter set under
> **BMRA-derived confidence-interval inequality constraints**, across the full joint TrkA+TrkB ×
> multi-inhibitor dataset, with scatter search + simplex (Methods p.22). This demo therefore does
> **not** recover the published parameters (`VALIDATION.md` Gate 3b, PARTIAL/sloppy). A
> BMRA-constrained **real-world** fit is planned as a separate slug (see the paper-level README).

## What is fit

SH-SY5Y/TrkB cells stimulated with **BDNF** (DMSO, no-inhibitor arm). Seven phospho
**fold changes vs. t=0** at **0 / 10 / 45 min**, fit by the TrkB cSTAR model and overlaid on
data in the paper's **Fig. 4A** (pAKT/pERK/pJNK/pS6K/pRSK/pERBB) + **4B** (pTRK), TrkB in red.
Design: pre-equilibrate (no BDNF, `Lig_on = 0`) → add BDNF (`Lig_on = 1`) → measure
over 0–2700 s (= 45 min), synthesized from `preequilibrate: basal` → `condition: stim`.

**Fit method (paper).** Same as the TrkA twin: **pyBioNetFit** (scatter search + simplex,
sum-of-squares) under **BMRA-CI inequality constraints** (Methods p.22); **10-min** RPPA +
Western = training, **45-min** RPPA = validation. See **[`VALIDATION.md`](VALIDATION.md)** for
the primary-source audit (confidence 84/100; the generated network is identical to that of the
authors' `TrkB_S_model.bngl`, checked by [`compare_upstream.py`](compare_upstream.py) below;
`.exp` byte-reproducible from their RPPA).

## How TrkB differs from TrkA (`../cstar_trka`)

The published `TrkB_S_model.bngl` differs from the TrkA model by:

- **TrkB-specific connection coefficients** (`g_B_*` / `K_B_*` instead of `g_A_*`),
- an added **TrkB–ErbB hetero-dimerization** module — the **ERBB→ERK→RSK→ERBB
  autocatalytic amplifier** the paper highlights for TrkB — which enlarges the network to
  **114 species / 645 reactions** (vs. 82 / 405 for TrkA), and
- a **proliferation-only** phenotype force (`beta_B_ERK`, `beta_B_S6K`, both positive; no
  JNK/RSK negative terms).

Biologically this shows up in the data and the fit: TrkB signalling is **stronger and more
sustained** than TrkA, with a large **pAdRTK/ERBB2** response (≈11× at 10 min vs. ≈1.3× for
TrkA) and pERK/pAKT still elevated at 45 min.

## Files

| file | role |
|---|---|
| `cstar_trkb.bngl` | edition-2, fitting-ready model (no `actions` block); adapted from `TrkB_S_model.bngl` |
| `cstar_trkb.conf` | the edition-2 job setup |
| `cstar_trkb.exp` | fit target: 7 phospho fold changes at 0/10/45 min |
| `extract_exp.py` | reproducible extraction of `cstar_trkb.exp` from the authors' RPPA CSVs (TrkB columns) |
| `cstar_trkb_reproduction.png` | verification: model at published params vs. the RPPA data |
| `upstream.json` | pin of the authors' `TrkB_S_model.bngl` (commit, path, sha256); not vendored |
| `compare_upstream.py` | BioNetGen-only comparison of `cstar_trkb.bngl` with the authors' file |
| `upstream_comparison.json` | what that comparison measured, against which tolerance |

## Model observable → RPPA antibody mapping

Same mapping as the TrkA twin; data are the **TrkB + BDNF, DMSO** columns, replicate-averaged
(r1–r6) and normalized to t=0. `FC_p*()` are model `functions` (parentheses in the header).

| model fn | reads model observable | RPPA antibody (DMSO/BDNF) |
|---|---|---|
| `FC_pTRK()`   | `pTRK`   | pTrk (Western; `RPPA_data_Trk_normalized_new.csv`) |
| `FC_pAdRTK()` | `pAdRTK` | ErbB-2/Her2/EGFR P Tyr1248/Tyr1173 |
| `FC_ppERK()`  | `ppERK`  | p44/42 MAPK (ERK1/2) P Thr202/Thr185 |
| `FC_pAKT()`   | `pAKT`   | Akt P Ser473 |
| `FC_pJNK()`   | `pJNKt`  | SAPK/JNK P Thr183 |
| `FC_pRSK()`   | `pRSK`   | p90 S6 kinase (Rsk1-3) P Thr359 |
| `FC_pS6K()`   | `pS6Kt`  | p70 S6 Kinase P Thr389 |

## Changes from the published model (all documented in `cstar_trkb.bngl`)

- **No `begin actions` block** — the protocol is synthesized from the conf.
- **`Lig_on` ligand gate** (0/1) multiplies every ligand-association forward rate,
  **including the hetero-dimer `sqrt(fh*factor)` family unique to TrkB** (13 forward rates
  gated in total).
- **All inhibitor concentrations = 0** — DMSO (no-inhibitor) arm.
- Fitted constants are bare `id nominal`; the conf's `*_var` free params bind by name
  (ADR-0034, no `__FREE`).

## Free parameters (8)

TrkB-specific ids, each `loguniform` over ≈±1 decade around the published value:
`vpTRK` (2.356), `vpRTK` (0.12), `g_B_TRKERK` (67.658), `vpERK` (0.006),
`g_B_TRKAKT` (30), `vpAKT` (0.0003), `kpRSKB` (0.018), `g_B_AKTS6K` (10).

## ⚠️ Native-only (not PEtab-v2-exportable)

Fold-change data → `normalization = init` → not PEtab v2-expressible (like the shipped
`igf1r.conf`). Verified **without** the PEtab round-trip.

## Verification (all pass)

- **Tier-1** (`check_conf.py`): edition 2, `job_type=de` resolves, data bound, 8 free params
  bind by id, no `__FREE`. **PASS.**
- **Model build** (BNG2.pl): generates **114 species, 645 reactions**.
- **Bounded bngsim fit** (`de`, `max_iterations=2`, `population_size=6`): finite objective
  (best `sos ≈ 41.7`) in ~17 s — the simulate→score→propose loop runs (not `heavy`).
- **Paper reproduction** (`cstar_trkb_reproduction.png`): at the published parameters the
  model reproduces Fig. 4B — strong, sustained multi-pathway activation, with pAKT nearly
  exact and the large ERBB2/AdRTK response captured. Objective at published params
  `sos ≈ 69.8`, median relative error ≈ 29 % (max ≈ 42 %) across the 10-/45-min points; the
  bounded fit already improves on this single-arm subset. (Absolute `sos` is larger than
  TrkA's simply because TrkB fold changes are ≈7–12× vs. ≈1–5×.)

## Checked against the authors' file (BioNetGen only)

`compare_upstream.py` compares `cstar_trkb.bngl` as committed with the authors'
`Trk_AB_models/TrkB_S_model.bngl` at
[OleksiiR/cSTAR_Nature@bfb092f](https://github.com/OleksiiR/cSTAR_Nature/tree/bfb092f2e4bffcfc43b3de01c301303fc6d4fd36)
(sha256 `6086bdef...`, pinned in `upstream.json`). That repository grants no licence ("All
rights reserved"), so the file is not copied here: `--fetch` downloads it from the pinned commit
into the git-ignored `dev/upstream/` and checks its sha256, or `--upstream DIR` reads a local
checkout. Nothing uses bngsim; the shared code is `../cstar_upstream.py`. Every number below is
in `upstream_comparison.json` (BioNetGen 2.9.3):

| check | measured | asserted |
|---|---|---|
| species (canonical labels), reactions (canonical reactants and products) | 114 = 114, 645 = 645, both bijections, same order | bijection |
| rate coefficients, dy/dt (to each species' gross flux), functions, observables at 32 sampled states, `Lig_on = 1` | 0 | ≤ 1e-12 |
| dy/dt with the authors' free ligand at 0 and this file's at 10 with `Lig_on = 0` | 0 | ≤ 1e-12 |
| differences in parameters | only `Lig_on`, the six inhibitor doses (0 here) and the thirteen ligand-association rate laws that carry `Lig_on` (four of them the hetero-dimer family) | exactly these |
| trajectories under the published protocol: 2 equilibrations + 7 stimulated arms (DMSO and six inhibitors), every species, observable and function | identical to the printed digits | ≤ 1e-9 of peak |

The published protocol is the authors' own actions block, run verbatim on their file. The slug's
file runs the same actions, translated mechanically: `TRKl(Rec)` 0/10 becomes `Lig_on` 0/1 (the
free ligand then stays at 10 during the equilibrations, checked separately), and an inhibitor
dose the slug sets to 0 is written out as the authors' value.

**The published pre-equilibration is not a steady state.** It is 50000 s without ligand, a reset
of the phenotype coordinate `S` to `2*S0`, and another 50000 s. Run to steady state instead
(reported, not asserted), the AdRTK species end up to 1.8 % away, and without the `S` reset
`Sval` settles at +4.374, in the proliferative `Sp` well, instead of the authors' -4.140 in the
`S0` well. The fitted fold changes move by at most 5.6e-3, but with the reset the stimulated
`Seq`/`S_FSforce` differ by up to 123 % (S6K-inhibitor arm), and without it `Sval` by up to
191 %. A simulation meant to match the paper should therefore replay the finite protocol, reset
included. BioNetGen's default tolerances (atol = rtol = 1e-8) are themselves accurate only to
~2e-3 of peak on species whose peak is below 1e-4, and to 1.5e-2 on one of peak 1e-6 in the
TRK-inhibitor arm (`solver_error`, against 1e-12).

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
uv run python pybnf-jobs/Rukhlenko-2022/cstar_trkb/compare_upstream.py --fetch
```

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Rukhlenko-2022/cstar_trkb
pybnf -c cstar_trkb.conf
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='cstar_trkb', conf='cstar_trkb.conf', simulator='ode',
    observables=('FC_pTRK', 'FC_pAdRTK', 'FC_ppERK', 'FC_pAKT',
                 'FC_pJNK', 'FC_pRSK', 'FC_pS6K'),
    system='cSTAR core signalling, TrkB/BDNF phospho fold changes (Rukhlenko 2022, '
           'PMC9644236); ODE + TrkB-ErbB hetero-dimerization, pre-equilibrate -> BDNF, '
           '45 min; normalization=init -> NATIVE-ONLY (not PEtab-exportable)'),
    # native-only: assert export_job raises NotImplementedError instead of a PEtab-lint test.
```
