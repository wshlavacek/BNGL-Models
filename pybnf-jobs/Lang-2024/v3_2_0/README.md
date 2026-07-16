# v3_2_0 — RPE-1 cell-cycle model, full fused v3.2.0 parameter fit (PyBNF edition-2 job)

A PyBNF edition-2 parameter-fitting job that captures the **headline parameter-estimation
problem** of Lang et al. (2024): fit **177 of 205 parameters** of the full fused, rule-based
RPE-1 cell-cycle model to **8 nuclear single-cell observables** over ~2 cell cycles, derived from:

> Lang PF, Penas DR, Banga JR, Weindl D, Novak B. **"Reusable rule-based cell cycle model
> explains compartment-resolved dynamics of 16 observables in RPE-1 cells."**
> *PLoS Comput Biol* 2024; **20**(1):e1011151.
> DOI: [10.1371/journal.pcbi.1011151](https://doi.org/10.1371/journal.pcbi.1011151)

Built with the `curate-pybnf-job` skill by transcribing the authors' PEtab problem
(`paulflang/cell_cycle_petab` @ `versions/v3.2.0`, commit `d562d3d`) into a native edition-2
`.conf` + `.bngl` + `.exp`. This is **Phase 2** of the Lang-2024 curation and uses a **different
model** from the Phase-1 Fig-1 transition submodels in `models/cell_cycle_oscillator_lang2024/`
— it is the *full fused* model those submodels were assembled into.

> ⚠️ **A setup transcription, not a converged-fit reproduction — read [`VALIDATION.md`](VALIDATION.md).**
> The PEtab problem ships `nominalValue`s that are a **reference point, not the cluster-optimized
> (saCeSS) optimum**, so at nominals the model reproduces the oscillation's **period and phase**
> but not calibrated **amplitudes** (pRB is flat). Converging the 177-parameter fit is
> **cluster-scale** (`heavy=True`); this job is scoped to *set it up*, not run it to convergence.

## The model

The full fused human RPE-1 cell-cycle regulatory network — Rb/E2F restriction point, cyclin
E/A/B–CDK activity, the APC/C–FZR1(Cdh1)/CDC20 and FBXO5(Emi1) switches, the Greatwall–ENSA–PP2A
(B55) mitotic phosphatase loop, Wee1/Cdc25, FOXM1 transcription, and the p21(CDKN1A)/p27(CDKN1B)/
Skp2 CDK-inhibitor module. It **has a `cell` compartment → ODE only (NFsim is out)**. The network
is **finite and small — 73 species, 332 reactions, builds in ~0.4 s** — so the model is cheap;
what makes the job **heavy** is the **177-parameter search over a single condition**, which is not
locally convergent (the paper used cooperative scatter search on a cluster).

## What is fit

8 nuclear observables over one wildtype time course (`v3_2_0.exp`, 600 points, t = 0..1.36e5,
the Stallaert et al. *Cell Systems* 2021 single-cell medians duplicated to 2 cell cycles):

| observable (`.exp` column) | measurement model | species sum |
|---|---|---|
| `obs_cycA__nuc_median()` | `oCCNA + sCCNA·tCCNA` | total cyclin A (`CCNA()`) |
| `obs_cycB1__nuc_median()` | `oCCNB + sCCNB·tCCNB` | total cyclin B1 (`CCNB()`) |
| `obs_cycE__nuc_median()` | `oCCNE + sCCNE·tCCNE` | total cyclin E (`CCNE()`) |
| `obs_E2F1__nuc_median()` | `oE2F + sE2F·(tE2F − E2F_pSer332_prom)` | E2F, minus Ser332-p promoter-bound forms |
| `obs_pRB__nuc_median()` | `oRB1_pSer807Ser811 + sRB1_pSer807Ser811·pRB` | free phospho-RB (`RB1(E2F,Ser807_Ser811~p)`) |
| `obs_Skp2__nuc_median()` | `oSKP2 + sSKP2·tSKP2` | total Skp2 |
| `obs_p21__nuc_median()` | `oCDKN1A + sCDKN1A·tCDKN1A` | total p21 (`CDKN1A()`) |
| `obs_p27__nuc_median()` | `oCDKN1B + sCDKN1B·tCDKN1B` | total p27 (`CDKN1B()`) |

Each observable is a model **function** (parens in the `.exp` header; no parens in `_manifest.py`).
The 8 offsets `o*` and 8 scales `s*` are **fitted observation-model parameters**. Noise is a fixed
unit Gaussian σ (`noiseParameter1 + noiseParameter2 = 1`) → the objective is plain **`sos`**.

## Files

| file | role |
|---|---|
| `v3_2_0.bngl` | edition-2, fitting-ready ODE model; 114 kinetic + 16 observation + 73 ic params; all 73 species seeded (incl. the 12 fitted initials generated at 0 in the staged model); no actions block |
| `v3_2_0.conf` | the edition-2 `de` job setup (`objective = sos`; 177 `*_var` free params) |
| `v3_2_0.exp` | fit target — 8 nuclear observables, t = 0..1.36e5 (Stallaert 2-cycle medians) |
| `make_reproduction.py` | regenerates the reproduction PNG (model at nominals, or `--params-file` a PyBNF best fit) |
| `v3_2_0_reproduction.png` | verification figure (8 panels, model vs data) |
| `VALIDATION.md` | provenance / exact model-fidelity mapping / verification record |

## Free parameters (177 estimated of 205)

| block | count | family | note |
|---|---|---|---|
| kinetic rate constants | 114 | `loguniform_var` (111) / `uniform_var` (3) | 3 log-scaled rates with PEtab lower bound 0 (`kDpApc_1`, `kDpE2f1`, `kPhC25A`) can't be log-searched → linear `[0, upper]` |
| initial concentrations | 47 (of 73) | `uniform_var` | the `wt` condition; 26 ic are fixed (promoters, inputs, zeroed complexes) |
| observation model | 16 | `uniform_var` | 8 offsets `o*` over `[0, upper]`, 8 scales `s*` over `[0.1, 10]` |

Bounds bracket the PEtab `parameters_v3.2.0.tsv` lower/upper. Nominals are the PEtab
`nominalValue`s (kinetic values taken exact from the staged `.bngl`, which the rounded TSV column
loses). The PEtab Laplace regularization priors on the offsets/scales are **simplified to uniform**
here to stay inside the PEtab-exportable subset (see `VALIDATION.md`).

## ✅ PEtab-exportable · 🔶 heavy (cluster-scale fit)

Stays inside the PEtab v2 exportable subset (`sos`, function/arithmetic observables, constant
Gaussian noise — no `normalization`/`cumulative`/`neg_bin`/constraints), so it **round-trips
through PEtab v2** (export → lint clean → import). Marked **`heavy=True`**: the model builds fast,
but a converged 177-parameter fit is cluster-scale, so it stays in the backend-free test tier.

## Verification (see `VALIDATION.md`)

- **Tier-1** (`scripts/check_conf.py`): edition 2, `job_type=de` resolves, `v3_2_0.exp` bound,
  177 free params bind by id, no `__FREE`. **PASS.**
- **PEtab v2 round-trip** (`scripts/petab_roundtrip.py`): export (8 observables, 177 params, 4800
  measurements) → `petab.v2` lint clean → import. **PASS.**
- **Real bngsim fit** (`pybnf -c`, `population_size=6, max_iterations=2`): finite objective
  (SOS ≈ 3.0e4, k=177, n=4800) in ~14 s. **PASS.**
- **Reproduction** (`v3_2_0_reproduction.png`, model at PEtab nominals): the 2-cycle oscillation's
  **period and phase** are reproduced for cyclin A/B1/E, E2F1, Skp2; amplitudes are uncalibrated
  (overall median rel err 51 %; pRB flat). A converged reproduction needs the cluster fit.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Lang-2024/v3_2_0
pybnf -c v3_2_0.conf                       # CLUSTER-SCALE (177 params); parallelize with islands
python make_reproduction.py                # model at PEtab nominals vs data
# after a real fit:
python make_reproduction.py --params-file output/.../best_fit.txt   # converged reproduction
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='v3_2_0', conf='v3_2_0.conf', simulator='ode',
    observables=('obs_cycA__nuc_median', 'obs_cycB1__nuc_median', 'obs_cycE__nuc_median',
                 'obs_E2F1__nuc_median', 'obs_pRB__nuc_median', 'obs_Skp2__nuc_median',
                 'obs_p21__nuc_median', 'obs_p27__nuc_median'),
    system='RPE-1 human cell cycle, full fused rule-based model (Lang 2024, '
           'DOI 10.1371/journal.pcbi.1011151; PEtab paulflang/cell_cycle_petab v3.2.0); ODE, '
           '73 species / 332 reactions, 8 nuclear single-cell observables (Stallaert 2021, 2 '
           'cycles), single wt condition, 177/205 params estimated (114 kinetic + 47 initial '
           'conditions + 16 observation), sos objective. PEtab-exportable. Cluster-scale fit '
           '(saCeSS in the paper); nominals are a reference point, not the converged optimum. '
           'See VALIDATION.md.',
    heavy=True,
)
# recover left empty: the PEtab problem does not ship a converged best-fit to recover to;
# a short fit reaches a finite objective (SOS ~3e4). See README for the reproduction table.
```
