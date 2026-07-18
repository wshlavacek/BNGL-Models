# Miller 2026 — leveraging qualitative + quantitative data for MEK-isoform ERK-cascade parameterization & UQ

> Miller EF, Mallela A, Neumann J, Lin YT, Hlavacek WS, Posner RG. Using PyBioNetFit to leverage
> qualitative and quantitative data in biological model parameterization and uncertainty
> quantification. *Frontiers in Immunology*. 2026;17:1663008. doi:10.3389/fimmu.2026.1663008

A methods study that re-parameterizes the **Kocieniewski & Lipniacki (2013)** MEK1/MEK2-isoform model
of the Raf/MEK/ERK cascade using **PyBioNetFit**, replacing the original manual trial-and-error tuning
with a reproducible, automated pipeline that **fuses quantitative time-series with qualitative up/down
orderings** (formalized as BPSL constraints) — and adds the uncertainty quantification the original
study lacked. This folder mirrors the authors' own PyBioNetFit job setups, grounded in primary sources
and audited by the `validate-pybnf-job` skill.

## Shared model & framework

- **Model** (Kocieniewski & Lipniacki 2013, *Phys Biol* 10:035006): 8 seed species (EGFR, Sos1, Ras,
  Raf, MEK1, MEK2, PHP, ERK), 41 rules; generated reaction network **109 species / 689 reactions**
  (K&L cite ~110 species/ODEs). MEK1 (not MEK2) carries the ERK→Thr292
  negative-feedback site and relays feedback to MEK2 via heterodimerization. Five cell lines — WT, KO,
  N78G (dimerization-dead), T292A (feedback-dead), T292D (phosphomimetic), each with its own
  parameter/IC perturbation (K&L Figs 4–7 scheme; see the edition-2 note below for how they are expressed).
- **Quantitative data:** WT pEGFR/pSOS1/pERK time-series (AU), Kamioka et al. *J Biol Chem* 2010;
  285:33540, Fig 3D (HeLa + 50 ng/ml EGF); tabulated verbatim as the paper's SI Table 8.
- **Qualitative data:** 90 up/down orderings of pMEK(pRDS)/pERK across the five cell lines,
  Catalanotti et al. *Nat Struct Mol Biol* 2009; 16:294 (EGF-stimulation figures); formalized as
  BPSL (`.prop`) statements, SI Tables 1–5.
- **Edition-2 (primary) + edition-1 legacy twins.** Both jobs are now **edition-2** (`edition = 2`,
  bngsim): **one** model (`wt.bngl`) plus four named `condition:` perturbations for the KO/N78G/T292A/
  T292D cell lines (**Mechanism A**). This is enabled by two `lanl/PyBNF` fixes on main — `9ab15167`
  (ADR-0027: a `condition:` may perturb a **fixed** or **IC-seeding** parameter, `MEK1_0`/
  `MEK1_0_T292p`, not only a free one) and `ae21d212` (`#483`: aMCMC `output_trajectory` under the
  one-model+conditions layout). Each experiment's regenerated network is **byte-identical** (109 sp)
  to the corresponding legacy per-variant model. The SI-faithful edition-1 files (five per-variant
  models with the perturbation baked in) are kept as `*_legacy/` twins for provenance and as the
  independent BNG2.pl reproduction oracle.
- **Native-only:** attaching BPSL constraints makes both jobs non-PEtab-exportable (`export_job`
  raises `NotImplementedError`); verify via each slug's `make_reproduction.py`, not PEtab. See
  `skills/curate-pybnf-job/references/bpsl-constraints.md`.

## Jobs

Each job is an **edition-2 primary** slug plus its **edition-1 legacy twin** (`*_legacy/`):

| slug (edition-2 primary) | legacy twin | analysis | free params | objective | reproduces | confidence |
|---|---|---|---|---|---|---|
| [`mek_isoform_de`](mek_isoform_de/) | [`mek_isoform_de_legacy`](mek_isoform_de_legacy/) | global MLE by differential evolution (data fusion) | 31 (28 rates + 3 scales) | `sos` = F_quant + F_qual | Table 3 best-fit; Fig 1B; Fig 2 (obj 24.0) | **92** |
| [`mek_isoform_amcmc`](mek_isoform_amcmc/) | [`mek_isoform_amcmc_legacy`](mek_isoform_amcmc_legacy/) | Bayesian UQ by adaptive MCMC on d3, u3 | 6 (d3, u3, σ, 3 scales) | `chi_sq_dynamic` | (d3,u3) posterior; Fig 5–7; Table 4 | **85** |

The two jobs are **distinct analyses** (the DE fits all 31 parameters; the aMCMC samples 2 rate
constants around K&L's original baseline), so they are separate slugs — see each slug's `VALIDATION.md`.
Within a slug, the edition-2 primary and its `*_legacy` twin generate byte-identical networks and
reproduce the same Gate-3a numbers (DE: F_quant 22.72, 89/90 BPSL; aMCMC: F_quant 50.05, 87/90 BPSL).

## Source materials (untracked `dev/papers/Miller-2026/`; not redistributed)

Miller 2026 paper + Supplementary Data Sheet 1; Kocieniewski & Lipniacki 2013; Kamioka 2010;
Catalanotti 2009; and the authors' own PyBioNetFit files
(`~/Code/PyBNF/examples/Miller2025_MEK_Isoforms/`).

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Miller-2026/<slug>          # e.g. mek_isoform_de (edition-2) or mek_isoform_de_legacy
pybnf -c <slug>.conf          # the fit / sampling job (both HEAVY — authors used an HPC cluster)
python make_reproduction.py   # reproduce the paper's figures at the reported parameters (via BNG2.pl)
```

The `mek_isoform_*` slugs are edition-2 (bngsim); the `mek_isoform_*_legacy` twins are edition-1
(BNG2.pl). Both reproduce the same figures at the paper's parameters.
