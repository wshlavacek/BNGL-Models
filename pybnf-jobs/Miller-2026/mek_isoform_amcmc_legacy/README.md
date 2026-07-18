# mek_isoform_amcmc — MEK-isoform ERK cascade, Bayesian UQ by adaptive MCMC (PyBNF, native-only)

The **uncertainty-quantification** half of the Miller-2026 study: a focused Bayesian analysis that
samples the posterior of the two profile-likelihood-identifiable feedback rate constants (**d3, u3**),
propagating parametric uncertainty into prediction credible bands.

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
held FIXED at their Table-2 defaults (baked into the models), and only **six** parameters are sampled —
`d3` (EGFR-dimer degradation), `u3` (Sos1 feedback dephosphorylation), the noise hyperparameter
`sigma`, and the three output scaling factors. Paper §2.4: *"for simplicity, we considered only two
adjustable model parameters."* This is **not** a UQ around the MLE; it is a distinct, deliberately
minimal Bayesian problem. → its own slug.

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

## Corrections & findings

Same audit corrections as the DE slug where applicable: CRLF→LF, standardized single action block
(`end model` slice), `KO.prop` `time=1600`→`3600` typo fixed. The aMCMC scale-factor bounds need **no**
widening — with the K&L-baseline model (`c2=2e-7`) the pEGFR copy number is ~2× higher than at the MLE,
so `scalepEGFR∈[3e-5,6e-5]` reaches the data (unlike the DE slug). `c1_L`/`p3` findings: see the DE
slug's VALIDATION.

## ⚠️ Native-only + edition-1 (legacy) + a PyBNF-version caveat

Native-only (BPSL `.prop` → `export_job` raises `NotImplementedError`); edition-1 legacy (per-model
perturbation via `X__FREE` tokens + `am` + `chi_sq_dynamic`). **The authors used PyBioNetFit v1.1.9.**
The initially-installed **v1.6.0** hit a regression in the aMCMC *cross-model BPSL
constraint-satisfaction tracker* (`adaptive_mcmc.py:191 → evaluate_constraints`: `res.simdata` was
`None` when a probit constraint references another cell line's observable — `constraint.py:473`),
filed as **[lanl/PyBNF#480](https://github.com/lanl/PyBNF/issues/480)** and since **fixed + verified
upstream on this example**. With the fix the full 5-model chain samples end-to-end (see Verification);
the WT-only isolation chain and the DE slug's identical 5-model + 90-constraint machinery independently
confirm the slug was correct throughout.

## Verification (see `VALIDATION.md`)

- **Tier-1** (`check_conf.py`): `fit_type=am` resolves, 5 BPSL constraint sets bound, 6 free params,
  native-only flagged. (edition-2 FAILs are N/A to a legacy job.)
- **Model build:** each of the 5 models → 110 species / 689 reactions (K&L-baseline parameters).
- **aMCMC sampler works** (isolation test): WT-only, quant-only chain draws 100 posterior samples,
  `Fitting complete`, exit 0.
- **Full 5-model + cross-model BPSL chain:** runs to completion after the
  [lanl/PyBNF#480](https://github.com/lanl/PyBNF/issues/480) fix — a 300-sample smoke chain reports
  `Fitting complete` and writes `samples.txt` + `constraint_samples.txt` (the tracker that had crashed).
  The full 5-chain × 300k posterior is heavy (~5 days) but no longer blocked.
- **Gate 3a — posterior-median predictive** (`make_reproduction.py`, d3≈1.0e-3, u3≈0.4e-3 from Fig 7,
  on the K&L baseline): `F_quant = 50` (looser than the MLE's 22.7 — expected; the original
  parameterization scores 40 vs the MLE's 24, SI Table 7), **87/90 BPSL constraints satisfied**, WT
  quant (Fig 6) + all 5 MEK_pRDS trajectories (Fig 5) reproduced (`mek_isoform_amcmc_reproduction.png`).

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"
cd pybnf-jobs/Miller-2026/mek_isoform_amcmc
pybnf -c mek_isoform_amcmc.conf     # one chain (HEAVY; needs PyBNF ~v1.1.9 for the BPSL tracker)
python make_reproduction.py         # posterior-median predictive (Fig 5/6) via BNG2.pl
```

## `_manifest.py` note

```python
RealWorldExample(
    folder='mek_isoform_amcmc', conf='mek_isoform_amcmc.conf', simulator='ode',
    observables=('MEK_pRDS_WT','pERK1_2_wt','pEGFR_wt','pSos1_wt'),
    system='MEK-isoform ERK cascade Bayesian UQ (Miller 2026); adaptive MCMC over d3,u3,sigma + 3 '
           'scales on the Kocieniewski & Lipniacki 2013 baseline; probit BPSL -> NATIVE-ONLY; '
           'cross-model BPSL tracker needs PyBNF ~v1.1.9 (v1.6.0 regression documented).',
    # verify via make_reproduction.py + WT-only isolation chain; NOT PEtab.
```
