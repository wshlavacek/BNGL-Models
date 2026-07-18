# mek_isoform_de — MEK1/MEK2 isoform ERK cascade, global MLE by differential evolution (PyBNF, native-only)

The **data-fusion** flagship of the Miller-2026 study: one global fit that leverages **quantitative**
WT time-series **and** 90 **qualitative** up/down orderings (BPSL) to parameterize a MEK-isoform
model of the ERK cascade, yielding maximum-likelihood estimates for 31 parameters.

This is a PyBNF **edition-2 (new-era)** job (`mek_isoform_de.conf` / `wt.bngl`): **one** model plus
four `condition:` perturbations for the KO/N78G/T292A/T292D cell lines. The SI-faithful **edition-1**
twin — five per-variant models with the perturbations baked in — is kept as
[`../mek_isoform_de_legacy`](../mek_isoform_de_legacy) for provenance and as the independent
cross-engine (BNG2.pl) reproduction oracle.

> Miller EF, Mallela A, Neumann J, Lin YT, Hlavacek WS, Posner RG.
> **"Using PyBioNetFit to leverage qualitative and quantitative data in biological model
> parameterization and uncertainty quantification."**
> *Frontiers in Immunology* 2026; **17**:1663008.
> DOI: [10.3389/fimmu.2026.1663008](https://doi.org/10.3389/fimmu.2026.1663008)

**Sibling job:** [`../mek_isoform_amcmc`](../mek_isoform_amcmc) — Bayesian UQ (adaptive MCMC) on the
same model/data. **Provenance & earned confidence:** see [`VALIDATION.md`](VALIDATION.md).

## The model

An ODE model of the Raf/MEK/ERK cascade that distinguishes the two MEK isoforms — **base model:
Kocieniewski P, Lipniacki T, *Phys Biol* 2013; 10:035006** (8 seed species: EGFR, Sos1, Ras, Raf,
MEK1, MEK2, PHP, ERK; 41 rules; **109 species / 689 reactions**). MEK1 (but not MEK2) carries the
ERK→Thr292 negative-feedback site; MEK1 relays feedback to MEK2 via heterodimerization. Measured
outputs:

- `MEK_pRDS` = combined bis-phosphorylated (Raf-dependent-site) MEK1+MEK2 = `MEK1(S1~Ypp) + MEK2(S1~Ypp)`
- `pERK1_2_wt` = bis-phosphorylated ERK; plus `pEGFR_wt`, `pSos1_wt` (via `scaled_*` functions) for the WT quantitative fit.

Five cell lines. In edition-2 they are **one model (`wt.bngl`) + four named `condition:`
perturbations** (K&L Fig 4–7 scheme):

| experiment | cell line | `condition:` perturbation |
|---|---|---|
| `WT` | wild type | none (the base `wt.bngl`) |
| `KO` | MEK1 knockout | `MEK1_0 = 0` |
| `N78G` | dimerization-dead | `b2 = 0, b4 = 0` (no MEK1 homo/heterodimerization) |
| `T292A` | feedback-dead | `p4 = 0` (no Thr292 feedback phosphorylation) |
| `T292D` | phosphomimetic | `MEK1_0 = 0, MEK1_0_T292p = 134000` (seed swap), `u4 = 0`, `b5 / 3` |

## What is fit

| data | source | how it enters the objective |
|---|---|---|
| **`WT.exp`** — pEGFR, pSOS1, pERK in AU at 0/300/600/900/1800/3600 s (18 points; padded to the uniform 0…3600 s / 300 s grid, `nan` elsewhere) | Kamioka et al. *J Biol Chem* 2010; 285:33540, **Fig 3D** (HeLa + 50 ng/ml EGF); transcribed verbatim as SI Table 8 | `sos` (sum of squares), `F_quant` |
| **90 BPSL constraints** across the 5 `.prop` files (WT 30, KO 24, N78G 18, T292A 12, T292D 6) | qualitative up/down orderings of pMEK/pERK across cell lines, Catalanotti et al. *Nat Struct Mol Biol* 2009; 16:294 (Figs 3a/4h/4i/5d–g); formalized as SI Tables 1–5 | soft weighted penalties (`weight 0.001`), `F_qual` |

Objective: `F(θ) = F_quant + F_qual` (paper Eqs 5 + 11; Mitra et al. 2018). "The paper's result" =
**Table 3 "Best-fit value"** (31 params) + **Fig 1B** (WT quant) + **Fig 2 B/D/F/H** (qualitative
satisfaction); reported DE objective **24.0** (SI Table 7).

## Free parameters (31)

28 rate constants + 3 output scaling factors (copy number → published AU), **bound by id**
(ADR-0034, no `__FREE` marker); bounds = the authors' verbatim optimization bounds (Table 3). The
model uses `p2b_eff = p2b*5`; the free id `p2b` binds the fit base. Twelve rate constants are
practically **identifiable**; the rest are sloppy/non-identifiable (SI Fig 4). See the
nominal↔best-fit table in `VALIDATION.md`.

## How edition-2 expresses the mutants (Mechanism A) — one model + `condition:` perturbations

The four mutants perturb parameters (`b2/b4/p4/u4/b5` — some **also fit**) and initial conditions
(`KO`: `MEK1_0 = 0`; `T292D`: seed all MEK1 at Thr292p via `MEK1_0 = 0` + `MEK1_0_T292p = 134000`).
Edition-2 now applies all of these as `condition:` parameter perturbations on the single `wt.bngl`:

- **Fixed / IC parameters can be perturbed.** A `condition:` may now target a **fixed** model
  parameter or an IC-seeding parameter (`MEK1_0`, `MEK1_0_T292p`), not only a free one — ADR-0027,
  realized in the bngsim mutant path (`lanl/PyBNF` `9ab15167`: seed a non-free condition target from
  the model's nominal before mutating). Enabler #2 (`ae21d212`, `lanl/PyBNF#483`) fixes the aMCMC
  `output_trajectory` path used by the sibling slug; both are on `lanl/PyBNF` main.
- **No species `setConcentration` needed.** bngsim keeps the `.net` seed-species initializers
  symbolic (`… MEK1_wt(...) MEK1_0`) and re-derives the seed concentrations after every `set_param`
  (`_sync_species_initial_concentrations`), so perturbing `MEK1_0` / `MEK1_0_T292p` re-seeds the
  species — the KO knockout and the T292D seed-swap are just parameter perturbations.
- **The regenerated network is byte-identical to the legacy per-variant model** for every cell line
  (Gate 2), so the edition-2 job reproduces the legacy oracle exactly (F_quant 22.72, 89/90 BPSL).

`b2/b4/p4/u4` are set to `0`, and `b5 / 3` is a relative op on the free `b5`. (In the sibling aMCMC
slug these parameters are *fixed*, so its T292D condition writes `b5` as the absolute
`1.3333333333333335e-9` (= `4e-9/3`) — the BNGL-emit path refuses a relative op on a fixed parameter.)

## Corrections applied vs the authors' shipped setup (all documented in `VALIDATION.md`)

1. **`KO.prop` typo fixed:** `...MEK_pRDS at time=1600` → `time=3600` (SI Table 2, KO-C9; the models
   only output at 300/1800/3600 s, so the constraint could not evaluate as written).
2. **`.prop` references rewritten to the edition-2 observable** `<X>.MEK_pRDS` (the legacy per-model
   suffix names `MEK_pRDS_<X>` collapse to the single model's `MEK_pRDS`, keyed per experiment).
3. **`scalepEGFR` bound widened** `[3e-5,6e-5] → [3e-5,3e-4]` to bracket Table 3's own reported
   best-fit (1.1e-4); the model *requires* ≈1.05e-4 to match the pEGFR peak (the shipped bound could
   not reach the paper's value).
4. **`WT.exp` padded** to the uniform 0…3600 s / 300 s grid (`nan` where no data) so it shares an
   identical time grid with the constraint-only experiments — required because BPSL constraints
   compare observables *across* experiments (e.g. `KO.MEK_pRDS < WT.MEK_pRDS`).
5. **CRLF → LF** on all `.prop`/`.exp`.

Two **findings** are documented but *not* "corrected" (correcting them would stop reproducing the
paper's published result): **`c1_L` is inert** (no rule uses it; the effective ligand signal is
`c1 = 0.02`, K&L's S0 — this is why its SI-Fig-4 profile is flat), and **K&L's Table-1 `p3 = 2e-3` is
a µM-vs-copy-number unit slip** (the model's `p3 = 2e-9 (mcls·s)⁻¹` is the physically correct value).

## ⚠️ Native-only (not PEtab-v2-exportable)

Attaching BPSL `.prop` constraints makes the job native-only (`export_job` raises
`NotImplementedError`). Verify with `make_reproduction.py` + a bounded fit, **not**
`petab_roundtrip.py`.

## Verification (all pass — see `VALIDATION.md`)

- **Tier-1** (`check_conf.py`): `fit_type=de` resolves, 5 BPSL constraint sets bound, 31 free params,
  correctly flagged native-only.
- **Gate 2 — identical network:** the edition-2 `wt.bngl` + each `condition:` generates a network
  **byte-identical** (109 species / identical reaction topology) to the corresponding legacy
  per-variant model (`MEK1_{WT,KO,N78G,T292A,T292D}.bngl`). (The legacy VALIDATION's "110 sp" is an
  off-by-one; the true count is 109.)
- **Gate 3a — reproduce at Table-3 best-fit** (`make_reproduction.py`, c1=0.02): `F_quant = 22.72`
  (≈ paper objective 24.0; identical to the legacy oracle 22.7), **89/90 BPSL constraints satisfied**
  (paper states 84/90; the difference is borderline near-tie sensitivity to Table-3's 2-sig-fig
  rounding), WT quant + all 5 qualitative trajectories reproduced (`mek_isoform_de_reproduction.png`).
- **Gate 3b — recover by fitting** (bounded edition-2 DE, 10 iterations): objective **1668 → 41**
  (≈ K&L's original 40.0, trending to PyBNF's 24.0; matches the legacy 1668 → 40 in 15 iters);
  identifiable parameters recover, sloppy ones do not — the paper's own SI Fig 4 finding.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl (network generation)
cd pybnf-jobs/Miller-2026/mek_isoform_de
pybnf -c mek_isoform_de.conf          # edition-2 global MLE fit (bngsim; HEAVY: authors used 10000 iters on HPC)
python make_reproduction.py           # Gate-3a reproduction at Table-3 best-fit -> the PNG
```

The edition-2 fit runs on the **bngsim** backend (`edition >= 2 ⇒ bngsim`); BNG2.pl is still used
once to expand the rules into a reaction network.

## `_manifest.py` note (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='mek_isoform_de', conf='mek_isoform_de.conf', simulator='ode',
    observables=('MEK_pRDS','pERK1_2_wt','pEGFR_wt','pSos1_wt'),
    system='MEK1/MEK2 isoform ERK cascade (Kocieniewski & Lipniacki 2013); Miller 2026 global MLE '
           'fusing WT quantitative time-series (Kamioka 2010) + 90 BPSL up/down orderings '
           '(Catalanotti 2009); EDITION-2 one model + 4 condition: perturbations (fixed/IC targets, '
           'ADR-0027); BPSL -> NATIVE-ONLY.',
    # BPSL: assert export_job raises NotImplementedError; verify via make_reproduction.py, NOT PEtab.
    # Legacy edition-1 twin (5 per-variant models) at ../mek_isoform_de_legacy for cross-engine repro.
```
