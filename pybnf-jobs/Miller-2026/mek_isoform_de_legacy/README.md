# mek_isoform_de — MEK1/MEK2 isoform ERK cascade, global MLE by differential evolution (PyBNF, native-only)

The **data-fusion** flagship of the Miller-2026 study: one global fit that leverages **quantitative**
WT time-series **and** 90 **qualitative** up/down orderings (BPSL) to parameterize a MEK-isoform
model of the ERK cascade, yielding maximum-likelihood estimates for 31 parameters.

> Miller EF, Mallela A, Neumann J, Lin YT, Hlavacek WS, Posner RG.
> **"Using PyBioNetFit to leverage qualitative and quantitative data in biological model
> parameterization and uncertainty quantification."**
> *Frontiers in Immunology* 2026; **17**:1663008.
> DOI: [10.3389/fimmu.2026.1663008](https://doi.org/10.3389/fimmu.2026.1663008)

**Sibling job:** [`../mek_isoform_amcmc`](../mek_isoform_amcmc) — Bayesian UQ (adaptive MCMC) on the
same models/data. **Provenance & earned confidence:** see [`VALIDATION.md`](VALIDATION.md).

## The model

An ODE model of the Raf/MEK/ERK cascade that distinguishes the two MEK isoforms — **base model:
Kocieniewski P, Lipniacki T, *Phys Biol* 2013; 10:035006** (8 seed species: EGFR, Sos1, Ras, Raf,
MEK1, MEK2, PHP, ERK; 41 rules; 110 species/ODEs). MEK1 (but not MEK2) carries the ERK→Thr292
negative-feedback site; MEK1 relays feedback to MEK2 via heterodimerization. The measured outputs:

- `MEK_pRDS_<line>` = combined bis-phosphorylated (Raf-dependent-site) MEK1+MEK2 = `MEK1(S1~Ypp) + MEK2(S1~Ypp)`
- `pERK1_2_wt` = bis-phosphorylated ERK; plus `pEGFR_wt`, `pSos1_wt` for the WT quantitative fit.

Five cell lines, each a **separate model** with its perturbation baked in (K&L Fig 4–7 scheme):

| model | cell line | perturbation |
|---|---|---|
| `MEK1_WT.bngl` | wild type | none |
| `MEK1_KO.bngl` | MEK1 knockout | `MEK1_0 = 0` |
| `MEK1_N78G.bngl` | dimerization-dead | `b2 = b4 = 0` (no MEK1 homo/heterodimerization) |
| `MEK1_T292A.bngl` | feedback-dead | `p4 = 0` (no Thr292 feedback phosphorylation) |
| `MEK1_T292D.bngl` | phosphomimetic | MEK1 seeded `T292~Yp`; `u4 = 0`; `b5 → b5/3` |

## What is fit

| data | source | how it enters the objective |
|---|---|---|
| **`WT.exp`** — pEGFR, pSOS1, pERK in AU at 0/300/600/900/1800/3600 s (18 points) | Kamioka et al. *J Biol Chem* 2010; 285:33540, **Fig 3D** (HeLa + 50 ng/ml EGF); transcribed verbatim as SI Table 8 | `sos` (sum of squares), `F_quant` |
| **90 BPSL constraints** across the 5 `.prop` files (WT 30, KO 24, N78G 18, T292A 12, T292D 6) | qualitative up/down orderings of pMEK/pERK across cell lines, Catalanotti et al. *Nat Struct Mol Biol* 2009; 16:294 (Figs 3a/4h/4i/5d–g); formalized as SI Tables 1–5 | soft weighted penalties (`weight 0.001`), `F_qual` |

Objective: `F(θ) = F_quant + F_qual` (paper Eqs 5 + 11; Mitra et al. 2018). "The paper's result" =
**Table 3 "Best-fit value"** (31 params) + **Fig 1B** (WT quant) + **Fig 2 B/D/F/H** (qualitative
satisfaction); reported DE objective **24.0** (SI Table 7).

## Free parameters (31)

28 rate constants + 3 output scaling factors (copy number → published AU), bounds = the authors'
verbatim optimization bounds (Table 3). Twelve rate constants are practically **identifiable**; the
rest are sloppy/non-identifiable (SI Fig 4). See the nominal↔best-fit table in `VALIDATION.md`.

## Corrections applied vs the authors' shipped setup (all documented in `VALIDATION.md`)

1. **`KO.prop` typo fixed:** `...T292A.MEK_pRDS_T292A at time=1600` → `time=3600` (SI Table 2, KO-C9;
   the models only output at 300/1800/3600 s, so the constraint could not evaluate as written).
2. **Model action blocks standardized:** the authors' variant files carried inconsistent trailing
   action junk (loose `setParameter`/`simulate` lines, a capitalized `Begin actions … writeXML()`
   block) that breaks BNG2 dispatch ("Not a CODE reference"). Each model now ends `end model` + one
   clean `setParameter("c1",0.02); simulate(suffix, ode, print_functions)`.
3. **`scalepEGFR__FREE` bound widened** `[3e-5,6e-5] → [3e-5,3e-4]` to bracket Table 3's own reported
   best-fit (1.1e-4); the model *requires* ≈1.05e-4 to match the pEGFR peak (the shipped bound could
   not reach the paper's value).
4. **CRLF → LF** on all `.prop`/`.exp`.

Two **findings** are documented but *not* "corrected" (correcting them would stop reproducing the
paper's published result): **`c1_L__FREE` is inert** (no rule uses it; the effective ligand signal is
`c1 = 0.02`, K&L's S0 — this is why its SI-Fig-4 profile is flat), and **K&L's Table-1 `p3 = 2e-3` is
a µM-vs-copy-number unit slip** (the model's `p3 = 2e-9 (mcls·s)⁻¹` is the physically correct value).

## ⚠️ Native-only (not PEtab-v2-exportable) & why edition-1 (legacy)

Attaching BPSL `.prop` constraints makes the job native-only (`export_job` raises
`NotImplementedError`). It also runs in **edition-1 (legacy)** mode (no `edition`/`bngl_backend`
line): the mutant models zero out parameters that are *also fit* (b2/b4/p4/u4/b5), which the legacy
`X__FREE` token idiom expresses per-model (the token is simply absent where perturbed), whereas
edition-2 bind-by-id overrides a parameter globally by name and would wipe out a variant's `b2 0`.
Verify with `make_reproduction.py` + a bounded fit, **not** `petab_roundtrip.py`.

## Verification (all pass — see `VALIDATION.md`)

- **Tier-1** (`check_conf.py`): `fit_type=de` resolves, 5 BPSL constraint sets bound, 31 free params,
  correctly flagged native-only. (The `edition==2` / `no __FREE` FAILs are edition-2-only checks,
  N/A to a legacy job — cf. `../../Erickson-2019/igf1r/igf1r_legacy.conf`.)
- **Model build:** each of the 5 models generates **110 species / 689 reactions** through BNG2.pl.
- **PyBNF end-to-end** (bounded legacy DE): the multi-model BPSL+quant fit loop runs to a **finite
  objective** (0 failed sims).
- **Gate 3a — reproduce at Table-3 best-fit** (`make_reproduction.py`, c1=0.02): `F_quant = 22.7`
  (≈ paper objective 24.0), **89/90 BPSL constraints satisfied** (paper states 84/90; the difference
  is borderline near-tie sensitivity to Table-3's 2-sig-fig rounding), WT quant + all 5 qualitative
  trajectories reproduced (`mek_isoform_de_reproduction.png`).
- **Gate 3b — recover by fitting** (bounded DE, 15 iterations): objective **1668 → 40.3** (≈ K&L's
  original 40.0, trending to PyBNF's 24.0); **12/31 params within 3×** of Table 3 — precisely the
  identifiable set (SI Fig 4). Sloppy/non-identifiable by construction, as the paper reports.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Miller-2026/mek_isoform_de
pybnf -c mek_isoform_de.conf          # the global MLE fit (HEAVY: authors used 10000 iters on HPC)
python make_reproduction.py           # Gate-3a reproduction at Table-3 best-fit -> the PNG
```

## `_manifest.py` note (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='mek_isoform_de', conf='mek_isoform_de.conf', simulator='ode',
    observables=('MEK_pRDS_WT','pERK1_2_wt','pEGFR_wt','pSos1_wt'),
    system='MEK1/MEK2 isoform ERK cascade (Kocieniewski & Lipniacki 2013); Miller 2026 global MLE '
           'fusing WT quantitative time-series (Kamioka 2010) + 90 BPSL up/down orderings '
           '(Catalanotti 2009); 5 cell-line models; edition-1 legacy; BPSL -> NATIVE-ONLY.',
    # BPSL: assert export_job raises NotImplementedError; verify via make_reproduction.py, NOT PEtab.
```
