# fig7b_timecourse — ATF2 phosphoswitch anisomycin time course (PyBNF edition-2 job)

A PyBNF edition-2, parameter-fitting job setup derived from:

> Kirsch K, Zeke A, Tőke O, Sok P, Sethi A, Sebő A, Kumar GS, Egri P, Póti ÁL,
> Gooley P, Peti W, Bento I, Reményi A, Alexa A.
> **"Co-regulation of the transcription controlling ATF2 phosphoswitch by JNK and p38."**
> *Nat Commun* 2020; **11**(1):5769.
> PMCID: [PMC7666158](https://pmc.ncbi.nlm.nih.gov/articles/PMC7666158/) ·
> DOI: [10.1038/s41467-020-19582-3](https://doi.org/10.1038/s41467-020-19582-3)

Built with the `curate-pybnf-job` skill. The authors provide their rule-based model as
**Supplementary Software 1** (`Bionetgen_JNK-p38-ATF2_model.bngl`); its actions block
generates the paper's **Fig. 7b**, so it is a natural real-world example.

## The model

JNK and p38 dock the ATF2 transactivation domain at distinct sites — JNK at the D-site;
p38 bipartitely at the D-site **and** the F-site (SPFENEF / S90 region) — and
distributively phosphorylate the T69/T71 "phosphoswitch". JNK additionally phosphorylates
the F-site **S90**, and **S90-P sterically blocks p38's F-site (FRS) binding** — the
paper's central mechanism. MKK6/MKK7 activate p38/JNK and phosphatases reverse everything.
Deterministic ODE; **48 species, 152 reactions**.

## What is fit

WT **pp-ATF2(T69/T71)** western-blot time course after **anisomycin** stimulation of
HEK293T cells, **0–41 min** (19 points, mean ± SD over n=3). Fit by the model observable
`pT69pT71` (doubly-phosphorylated ATF2).

**Experimental design.** Pre-equilibrate at basal MAPK activity (`Stim = 0` → `keq6/keq7`),
then apply anisomycin (`Stim = 1` → `kstim6/kstim7`) and measure over 0–2466 s. PyBNF
synthesizes the two-phase protocol from `preequilibrate: basal` → `condition: stim`.

## Files

| file | role |
|---|---|
| `fig7b_timecourse.bngl` | edition-2, fitting-ready model (no `actions` block); reconstructed from Supplementary Software 1 |
| `fig7b_timecourse.conf` | the edition-2 job setup |
| `fig7b_timecourse.exp` | fit target: WT pp-ATF2(T69/T71) treated/untreated, 0–41 min, with SD |
| `extract_exp.py` | reproducible extraction of the `.exp` from the paper's Source Data workbook |
| `make_reproduction.py` | regenerates the reproduction PNG (published vs fitted model vs data) |
| `fig7b_timecourse_reproduction.png` | verification figure |

**Data source.** The paper's **Source Data** workbook
([MOESM6](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19582-3/MediaObjects/41467_2020_19582_MOESM6_ESM.xlsx)),
sheet `Figure_4` / Fig. 4a, block *"Anti-ppATF2(T69/T71), treated/untreated"*, WT column
(mean + SD). Model provided as
[Supplementary Software 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19582-3/MediaObjects/41467_2020_19582_MOESM4_ESM.zip).
No public GitHub / BioModels deposition exists. Published best-fit parameters and their
data provenance are in **Supplementary Table 2**.

## Model observable → data mapping

`pT69pT71` is a `Molecules` observable (plain name in the `.exp` header). The data are
**treated/untreated** ratios (anisomycin-stimulated ÷ unstimulated control), so the conf
compares them to `pT69pT71` under `normalization = init` (each simulated series ÷ its own
basal/t=0 value).

## Changes from the published model (all documented in `fig7b_timecourse.bngl`)

- **No `begin actions` block** — the protocol is synthesized from the conf.
- **`Stim` gate** (0/1) selects basal (`keq*`) vs stimulated (`kstim*`) MAPK-activation
  rates through two `functions`, replacing the authors' `setParameter(kstim7, keq7)`
  equilibration idiom (same trick as `examples/real-world/receptor`). Verified to
  reproduce the authors' trajectories **exactly** (max relative difference 0.0).
- Fitted constants are bare `id nominal` declarations at the published values; the conf's
  `*_var` free params bind to them **by name** (ADR-0034, no `__FREE`).

## Free parameters (8)

The cell-based rate constants the paper fit to Fig. 7b / Fig. 4a (Supplementary Table 2,
"fitted from WB / NanoBit"), each `loguniform` over ≈±1 decade around the published value.
The docking/catalytic rates (`kon*`, `koff*`, `k1`–`k4`) were fixed in vitro (Fig. S7).

| id | published | role |
|---|---|---|
| `keq7` | 4.695e-6 | basal MKK7 → JNK activation |
| `kstim7` | 9.074e-5 | stimulated MKK7 → JNK activation |
| `keq6` | 1.74e-5 | basal MKK6 → p38 activation |
| `kstim6` | 1.16e-4 | stimulated MKK6 → p38 activation |
| `dp1` | 6.26e-4 | phosphatase on pp-JNK |
| `dp2` | 1.76e-3 | phosphatase on pp-p38 |
| `dp3` | 9.54e-3 | phosphatase on ATF2 pT69 / pT71 |
| `dp4` | 4.50e-3 | phosphatase on ATF2 pS90 |

## ⚠️ Native-only (not PEtab-v2-exportable)

The data are treated/untreated ratios, so the conf uses **`normalization = init`** — which
PyBNF's `export_job` cannot express in PEtab v2 (`NotImplementedError`), exactly like the
shipped `igf1r.conf` and the `Rukhlenko-2022/cstar_trka` example. Verified **without** the
PEtab round-trip.

## Verification (all pass)

- **Tier-1** (`skills/curate-pybnf-job/scripts/check_conf.py`): edition 2, `job_type=de`
  resolves, data bound, 8 free params bind by id, no `__FREE`. **PASS.**
- **Native-only assertion**: `pybnf.petab.export_job(...)` raises `NotImplementedError`
  (whole-fit `normalization`). **PASS.**
- **Model faithfulness** (BNG2.pl): the reconstructed BNGL reproduces the authors'
  published-parameter trajectories (`pT69pT71`, `p38ATF2all`, `ppp38`, `ppJNK`)
  **exactly** — max relative difference **0.0** at every checkpoint.
- **Bounded bngsim fit** (`de`, `max_iterations=2`, `population_size=6`): finite objective
  in ~6 s — the simulate→score→propose loop runs (not `heavy`).
- **Full fit + reproduction** (`fig7b_timecourse_reproduction.png`): a `de` + simplex-refine
  run reaches **chi_sq ≈ 64.9** and reproduces the WT data to **median 8.4 % relative error**
  (max 17.8 %), capturing the activation rise and amplitude. The model uses *constant*
  (step) MAPK activation, so it settles to a plateau and does **not** capture the modest
  post-9-min adaptation seen in the data — that shape residual dominates the remaining
  error and is the model's known limitation, not a fit failure.
  - At the **published** parameters the model overshoots the *treated/untreated* data
    (median ≈ 125 % error): those values were fit to the differently-normalized Fig. 7b
    (the model reaching ~80 % ATF2 phosphorylation in absolute terms), not to the Fig. 4a
    fold-change quantification used here. The fit correctly re-scales the basal/stimulated
    balance (chiefly a higher `keq7`) to match the fold amplitude.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Kirsch-2020/fig7b_timecourse
pybnf -c fig7b_timecourse.conf
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='fig7b_timecourse', conf='fig7b_timecourse.conf', simulator='ode',
    observables=('pT69pT71',),
    system='JNK/p38/ATF2 TAD phosphoswitch, WT pp-ATF2(T69/T71) anisomycin time course '
           '(Kirsch 2020, PMC7666158); ODE, pre-equilibrate -> stimulate, 41 min; '
           'normalization=init -> NATIVE-ONLY (not PEtab-exportable)'),
    # native-only: assert export_job raises NotImplementedError instead of a PEtab-lint test.
```
