# mek_isoform_amcmc — MEK-isoform ERK cascade, Bayesian UQ by adaptive MCMC (PyBNF, native-only)

The **uncertainty-quantification** half of the Miller-2026 study: a focused Bayesian analysis that
samples the posterior of the two profile-likelihood-identifiable feedback rate constants (**d3, u3**),
propagating parametric uncertainty into prediction credible bands.

This is a PyBNF **edition-2 (new-era)** job (`mek_isoform_amcmc.conf` / `wt.bngl`): **one** model plus
four `condition:` perturbations for the KO/N78G/T292A/T292D cell lines. The SI-faithful **edition-1**
twin (five per-variant models) is kept as [`../mek_isoform_amcmc_legacy`](../mek_isoform_amcmc_legacy)
for provenance and as the independent BNG2.pl reproduction oracle.

> Miller EF, Mallela A, Neumann J, Lin YT, Hlavacek WS, Posner RG.
> **"Using PyBioNetFit to leverage qualitative and quantitative data in biological model
> parameterization and uncertainty quantification."** *Frontiers in Immunology* 2026; **17**:1663008.
> DOI: [10.3389/fimmu.2026.1663008](https://doi.org/10.3389/fimmu.2026.1663008)
> aMCMC sampler: Neumann J, et al. *Bioinformatics* 2022; **38**:1770.

**Sibling / primary job:** [`../mek_isoform_de`](../mek_isoform_de) — the 31-parameter MLE fit. See
its README + [`VALIDATION.md`](../mek_isoform_de/VALIDATION.md) for the **shared** model, data
provenance (Kocieniewski & Lipniacki 2013 base model; Kamioka 2010 WT quant; Catalanotti 2009
qualitative orderings), and the Gate 0–2 evidence. This job's own audit: [`VALIDATION.md`](VALIDATION.md).

## Scope — a separate analysis from the MLE

The aMCMC baseline is Kocieniewski & Lipniacki's **original parameterization**: 25 rate constants are
held FIXED at their Table-2 defaults (the `wt.bngl` nominal), and only **six** parameters are sampled —
`d3` (EGFR-dimer degradation), `u3` (Sos1 feedback dephosphorylation), the noise hyperparameter
`sigma`, and the three output scaling factors. Paper §2.4: *"for simplicity, we considered only two
adjustable model parameters."* This is **not** a UQ around the MLE; it is a distinct, deliberately
minimal Bayesian problem → its own slug.

## What differs from the DE slug

| aspect | `mek_isoform_de` | `mek_isoform_amcmc` |
|---|---|---|
| algorithm | differential evolution (`de`) | adaptive MCMC (`am`) |
| objective | `sos` (`F_quant + F_qual`) | `chi_sq_dynamic` (dynamic noise model) |
| free params | 31 (all rate constants + scales) | 6 (d3, u3, sigma, 3 scales) |
| fixed params | none | 25 rate constants at K&L Table-2 defaults |
| BPSL penalty family | hinge (`weight 0.001`) | **probit** (`confidence 0.98 tolerance 1`) |
| result | Table 3 best-fit; Fig 1B, Fig 2 | (d3,u3) posterior: Fig 5–7, Table 4 (ESS, R̂) |

Sampling schedule (Table 4 / §2.4): **5 independent chains**, each 25,000 burn-in + 25,000 adaptation
+ 250,000 production iterations (~1.25M objective evaluations, ~5.15 days on 25 CPUs). The conf defines
**one** chain; run it 5× (different seeds) for the ensemble.

## How edition-2 expresses the mutants (Mechanism A) — and the fixed-parameter subtlety

Same one-model + `condition:` scheme as the DE slug, but with one twist: **here the perturbed
parameters `b2/b4/p4/u4/b5` are FIXED** (only d3, u3 + noise + scales are sampled). Edition-2
conditions may now perturb a **fixed** model parameter, or an IC-seeding parameter (`MEK1_0`,
`MEK1_0_T292p`), not only a free one — ADR-0027, realized in the bngsim mutant path (`lanl/PyBNF`
`9ab15167`). Two edition-2 details this slug exercises:

- **T292D's `b5 / 3` is written as the absolute value `1.3333333333333335e-9`.** Because `b5` is
  fixed at K&L's `4e-9` here (not fit), the division is a constant — and the BNGL-emit path refuses a
  *relative* op on a fixed parameter, so the condition uses the pre-divided absolute value, `4e-9 / 3`
  to full double precision (bit-identical to the authors' `b5 4e-9/3`).
- **`output_trajectory`** saves full observable trajectories for every sampled parameter set (→
  posterior-predictive bands). Under the one-model+conditions layout the base model runs every action
  suffix under every mutant; keying only the scored data-keys used to `KeyError` — fixed upstream by
  `lanl/PyBNF#483` (`ae21d212`, on main), so the chain now samples end-to-end with `output_trajectory`.

Each experiment's regenerated network is byte-identical to the legacy per-variant model (Gate 2), and
bngsim re-derives the seed-species concentrations from the perturbed `MEK1_0` / `MEK1_0_T292p`, so no
species `setConcentration` is needed.

## Corrections & findings

Same audit corrections as the DE slug where applicable: CRLF→LF, `KO.prop` `time=1600`→`3600` typo
fixed, `.prop` references rewritten to the edition-2 `<X>.MEK_pRDS`, `WT.exp` padded to the uniform
0…3600 s / 300 s grid (cross-experiment BPSL needs an identical time grid). The aMCMC scale-factor
bounds need **no** widening — with the K&L-baseline model (`c2=2e-7`) the pEGFR copy number is ~2×
higher than at the MLE, so `scalepEGFR∈[3e-5,6e-5]` reaches the data. `c1_L`/`p3` findings: see the DE
slug's VALIDATION.

## ⚠️ Native-only (not PEtab-v2-exportable)

Native-only (BPSL probit `.prop` + `chi_sq_dynamic` → `export_job` raises `NotImplementedError`).
Verify with `make_reproduction.py` + a short chain, **not** `petab_roundtrip.py`.

## Verification (see `VALIDATION.md`)

- **Tier-1** (`check_conf.py`): `fit_type=am` resolves, 5 BPSL constraint sets bound, 6 free params,
  native-only flagged.
- **Gate 2 — identical network:** the edition-2 `wt.bngl` + each `condition:` generates a network
  byte-identical (109 species) to the legacy per-variant model — the same shared network as the DE slug.
- **Gate 3a — posterior-median predictive** (`make_reproduction.py`, d3≈1.0e-3, u3≈0.4e-3 from Fig 7,
  on the K&L baseline): `F_quant = 50.05` (looser than the MLE's 22.7 — expected; the original
  parameterization scores 40 vs the MLE's 24, SI Table 7), **87/90 BPSL constraints satisfied** (both
  identical to the legacy oracle), WT quant (Fig 6) + all 5 MEK_pRDS trajectories (Fig 5) reproduced
  (`mek_isoform_amcmc_reproduction.png`).
- **Gate 3b — sample end-to-end:** a short chain (50 burn-in + 50 adapt + production over all 5 cell
  lines, 90 probit constraints, **with `output_trajectory`**) reports `Fitting complete`, exit 0, and
  writes `samples.txt` **and** `constraint_samples.txt` (100 samples each) — i.e. the cross-model BPSL
  constraint-satisfaction tracker runs cleanly. The full 5-chain × 300k posterior (Fig-7 marginals,
  Table-4 ESS/R̂) is heavy (~5 days) but no longer blocked.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"
cd pybnf-jobs/Miller-2026/mek_isoform_amcmc
pybnf -c mek_isoform_amcmc.conf     # edition-2 one chain (bngsim; HEAVY — cut max_iterations for a workstation)
python make_reproduction.py         # posterior-median predictive (Fig 5/6) via BNG2.pl
```

## `_manifest.py` note

```python
RealWorldExample(
    folder='mek_isoform_amcmc', conf='mek_isoform_amcmc.conf', simulator='ode',
    observables=('MEK_pRDS','pERK1_2_wt','pEGFR_wt','pSos1_wt'),
    system='MEK-isoform ERK cascade Bayesian UQ (Miller 2026); adaptive MCMC over d3,u3,sigma + 3 '
           'scales on the Kocieniewski & Lipniacki 2013 baseline; EDITION-2 one model + 4 condition: '
           'perturbations of FIXED/IC params (ADR-0027); output_trajectory (lanl/PyBNF#483); '
           'probit BPSL -> NATIVE-ONLY.',
    # verify via make_reproduction.py + a short chain; NOT PEtab.
    # Legacy edition-1 twin (5 per-variant models) at ../mek_isoform_amcmc_legacy for cross-engine repro.
```
