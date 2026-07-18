# VALIDATION — Miller-2026 / mek_isoform_amcmc

Primary-source validation of `pybnf-jobs/Miller-2026/mek_isoform_amcmc/`, per `validate-pybnf-job`.

> **Confidence: 85 / 100** — Gates 0–2 PASS at gold-standard level (**shared** with the DE slug:
> authors' own files + all three upstream sources); Gate 3a PASS (posterior-median predictive on the
> K&L baseline reproduces Fig 5/6 and 87/90 BPSL); Gate 3b PARTIAL — the full 5-model + cross-model
> BPSL chain now **runs to completion** (a 300-sample smoke chain samples + tracks constraints
> cleanly) after the **lanl/PyBNF#480** fix (committed locally, verified upstream on this very
> example); the full 5-chain × 300k-iteration posterior is heavy (~5 days) but no longer blocked.
> Below the DE slug only because the full-scale posterior (Fig-7 marginals, Table-4 ESS/R̂) was not run.

Primary sources & "shared" gates: see [`../mek_isoform_de/VALIDATION.md`](../mek_isoform_de/VALIDATION.md).
Authors' files: `~/Code/PyBNF/examples/Miller2025_MEK_Isoforms/MEK_isoform_aMCMC/`.

"The paper's result" for this job = the **(d3, u3) posterior** — Fig 7 marginals, Fig 5/6
posterior-predictive credible bands, Table 4 (ESS, R̂). UQ bounds = Table 3, last column.

---

## Gate 0 — Materials inventory
**PASS** (identical to the DE slug — all sources present, including the aMCMC sampler paper: Neumann
et al. *Bioinformatics* 2022; 38:1770).

## Gate 1 — Data provenance
**PASS.** Same `WT.exp` (≡ SI Table 8 ≡ Kamioka Fig 3D) as the DE slug. The 90 BPSL constraints are
the same orderings (SI Tables 1–5 ≡ Catalanotti), re-authored in the **probit** family
(`confidence 0.98 tolerance 1`) instead of the DE hinge (`weight 0.001`) — the Bayesian
likelihood form (Mitra & Hlavacek 2020). KO typo `time=1600`→`3600` corrected. CRLF→LF.

## Gate 2 — Model fidelity
**PASS (equiv).** Same 5-cell-line model structure and perturbations as the DE slug (110 sp / 689
rxn). The distinguishing feature: **25 rate constants are FIXED at Kocieniewski & Lipniacki's
Table-2 defaults** (baked in), so only 6 parameters (`d3, u3, sigma`, 3 scales) carry `__FREE`
tokens. Verified: the fixed values equal K&L Table 1 / Miller Table 2 (spot-checked c2=2e-7,
b1=4e-8, a1=1.5e-7, p2a=1e-6, p4=1.2e-9, u4=2e-4, b5=4e-9, u5=20). This is the deliberate
"for simplicity, two adjustable parameters" scope of paper §2.4 — a separate analysis from the MLE,
hence a separate slug (divergence policy: multiple published fits → multiple slugs).

## Gate 3a — Reproduce at the paper's parameters
`make_reproduction.py` at the Fig-7 posterior medians (d3≈1.0e-3, u3≈0.4e-3) on the K&L baseline:
- **WT quantitative:** `F_quant = 50` — looser than the MLE's 22.7, **as expected**: the K&L original
  parameterization scores 40 vs the MLE's 24 (SI Table 7) and global RMSE 1.47 vs 1.076 (SI Table 9).
  The posterior-predictive median captures the data trend (Fig 6), with the K&L baseline's
  characteristic pSos1 overshoot.
- **Qualitative:** **87/90 BPSL constraints** satisfied at the posterior median; MEK_pRDS trajectories
  reproduce Fig 5 (WT transient; KO/N78G/T292A sustained; T292D low). `mek_isoform_amcmc_reproduction.png`.

**Verdict: PASS.**

## Gate 3b — Recover the posterior by sampling
**PARTIAL (runs; full-scale posterior heavy).** The full 5-model + cross-model BPSL chain **now runs
to completion**: a 300-sample smoke chain (100 burn-in + 100 adapt + production over all 5 cell lines,
90 probit constraints) reports `Fitting complete`, exit 0, and writes both `samples.txt` and
`constraint_samples.txt` — i.e. the cross-model constraint-satisfaction bookkeeping that previously
crashed now works.

History (resolved): under the initially-installed **v1.6.0** the run aborted in the aMCMC cross-model
BPSL satisfaction tracker (`adaptive_mcmc.py:191 → constraint.py:473 get_key`: `res.simdata` is
`None`) — a regression vs the authors' **v1.1.9**. Filed as
**[lanl/PyBNF#480](https://github.com/lanl/PyBNF/issues/480)** and since **fixed + verified upstream
on this very example** (`85157141 fix(samplers): feed accepted-move simdata to constraint tracking`;
`1465997f docs: verified on MEK_isoform_aMCMC`). Isolation cross-check (still valid): a WT-only,
quant-only chain also samples cleanly, and the DE slug drives the identical 5-model + 90-constraint
machinery to a finite objective.

What is **not** run: the paper's full schedule — 5 chains × (25k burn-in + 25k adapt + 250k
production) ≈ 1.25M evaluations, ~5.15 days on 25 CPUs — so the Fig-7 marginals and Table-4 ESS/R̂
are not reproduced at scale (heavy, not blocked).

**Verdict: PARTIAL** — the sampler + constraint tracking run end-to-end; the full-scale posterior is
a compute-time matter, no longer a tooling block.

---

## Divergence & corrections
Scope matches the paper's Bayesian analysis. Corrections: KO typo, action-block standardization,
CRLF→LF (as the DE slug). Scale-factor bounds left at the authors' values (no widening needed on the
K&L baseline). **Upstream filing:** the v1.6.0 aMCMC cross-model BPSL `res.simdata=None` regression is
filed as **[lanl/PyBNF#480](https://github.com/lanl/PyBNF/issues/480)** (analogous to #473/#474).

## Bottom line
A faithful port of the paper's Bayesian-UQ setup with gold-standard Gate 0–2 grounding (shared with
the DE slug), a passing posterior-median reproduction, and — after the lanl/PyBNF#480 fix — a chain
that samples end-to-end with cross-model BPSL satisfaction tracking. Confidence sits at 85 (just below
the DE slug) only because the full-scale 5-chain × 300k posterior was not run at HPC scale. Most
valuable next step: run the full schedule to reproduce the Fig-7 marginals and Table-4 ESS/R̂.
