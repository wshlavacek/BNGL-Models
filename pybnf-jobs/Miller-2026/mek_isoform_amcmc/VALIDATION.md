# VALIDATION — Miller-2026 / mek_isoform_amcmc

Primary-source validation of `pybnf-jobs/Miller-2026/mek_isoform_amcmc/`, per `validate-pybnf-job`.

> **Reorganized for edition-2 (Mechanism A).** The job is now PRIMARY edition-2
> (`mek_isoform_amcmc.conf` / `wt.bngl` — one model + four `condition:` perturbations). The SI-faithful
> edition-1 files (five per-variant models) are kept as
> [`../mek_isoform_amcmc_legacy`](../mek_isoform_amcmc_legacy) for provenance and as the BNG2.pl
> reproduction oracle. Gate 0–2 provenance and model chemistry are **shared** with the DE slug and the
> legacy twin; the edition-2 job reproduces the legacy oracle **exactly** (Gate 3a: F_quant 50.05,
> 87/90 BPSL) and its short chain samples end-to-end with `output_trajectory` (Gate 3b).

> **Confidence: 85 / 100** — Gates 0–2 PASS at gold-standard level (**shared** with the DE slug:
> authors' own files + all three upstream sources); the edition-2 network is byte-identical to the
> legacy per-variant models; Gate 3a PASS (posterior-median predictive on the K&L baseline reproduces
> Fig 5/6 and 87/90 BPSL, identical to the legacy oracle); Gate 3b PARTIAL — the full 5-model +
> cross-model BPSL chain **runs to completion** with `output_trajectory` (a short chain samples +
> tracks constraints cleanly, writing `samples.txt` + `constraint_samples.txt`); the full 5-chain ×
> 300k-iteration posterior (Fig-7 marginals, Table-4 ESS/R̂) is heavy (~5 days) but no longer blocked.
> Below the DE slug only because the full-scale posterior was not run.

Primary sources & "shared" gates: see [`../mek_isoform_de/VALIDATION.md`](../mek_isoform_de/VALIDATION.md).
Authors' files: `~/Code/PyBNF/examples/Miller2025_MEK_Isoforms/MEK_isoform_aMCMC/`.
Edition-2 enablers (both on `lanl/PyBNF` main): `9ab15167` (ADR-0027 — conditions perturb fixed/IC
params) and `ae21d212` (`lanl/PyBNF#483` — aMCMC `output_trajectory` off-diagonal suffixes).

"The paper's result" for this job = the **(d3, u3) posterior** — Fig 7 marginals, Fig 5/6
posterior-predictive credible bands, Table 4 (ESS, R̂). UQ bounds = Table 3, last column.

---

## Gate 0 — Materials inventory
**PASS** (identical to the DE slug — all sources present, including the aMCMC sampler paper: Neumann
et al. *Bioinformatics* 2022; 38:1770).

## Gate 1 — Data provenance
**PASS.** Same `WT.exp` (≡ SI Table 8 ≡ Kamioka Fig 3D; padded to the 0…3600 s / 300 s grid) as the
DE slug. The 90 BPSL constraints are the same orderings (SI Tables 1–5 ≡ Catalanotti), re-authored in
the **probit** family (`confidence 0.98 tolerance 1`) instead of the DE hinge (`weight 0.001`) — the
Bayesian likelihood form (Mitra & Hlavacek 2020). References rewritten to the edition-2 `<X>.MEK_pRDS`.
KO typo `time=1600`→`3600` corrected. CRLF→LF.

## Gate 2 — Model fidelity (edition-2 network ≡ legacy per-variant models)
**PASS (equiv).** The single edition-2 `wt.bngl` + each `condition:` regenerates a network
byte-identical (109 sp) to the legacy per-variant model — the same shared 5-cell-line network as the
DE slug. The distinguishing feature: **25 rate constants are FIXED at Kocieniewski & Lipniacki's
Table-2 defaults** (the `wt.bngl` nominal), so only 6 parameters (`d3, u3, sigma`, 3 scales) are free.
Verified: the fixed values equal K&L Table 1 / Miller Table 2 (spot-checked c2=2e-7, b1=4e-8,
a1=1.5e-7, p2a=1e-6, p4=1.2e-9, u4=2e-4, b5=4e-9, u5=20). This is the deliberate "for simplicity, two
adjustable parameters" scope of paper §2.4 — a separate analysis from the MLE, hence a separate slug.

**Fixed-parameter conditions.** Because `b2/b4/p4/u4/b5` are fixed (not fit), the T292D condition
writes `b5 / 3` as the absolute `1.3333333333333335e-9` (= `4e-9 / 3` to full double precision,
bit-identical to the authors' `b5 4e-9/3`) — the BNGL-emit path refuses a relative op on a fixed parameter. The other perturbations (`b2/b4/p4/u4 = 0`, `MEK1_0 = 0`,
`MEK1_0_T292p = 134000`) apply directly; bngsim re-derives the seed concentrations. All are enabled by
ADR-0027 (`9ab15167`).

## Gate 3a — Reproduce at the paper's parameters
`make_reproduction.py` at the Fig-7 posterior medians (d3≈1.0e-3, u3≈0.4e-3) on the K&L baseline
(single edition-2 `wt.bngl` + each `condition:`):
- **WT quantitative:** `F_quant = 50.05` — looser than the MLE's 22.7, **as expected**: the K&L
  original parameterization scores 40 vs the MLE's 24 (SI Table 7) and global RMSE 1.47 vs 1.076 (SI
  Table 9). The posterior-predictive median captures the data trend (Fig 6), with the K&L baseline's
  characteristic pSos1 overshoot.
- **Qualitative:** **87/90 BPSL constraints** satisfied at the posterior median; MEK_pRDS trajectories
  reproduce Fig 5 (WT transient; KO/N78G/T292A sustained; T292D low). `mek_isoform_amcmc_reproduction.png`.
- Both figures (50.05, 87/90) are **identical to the legacy oracle** (50, 87/90) — the operational
  confirmation of Gate 2 (byte-identical networks).

**Verdict: PASS.**

## Gate 3b — Recover the posterior by sampling
**PARTIAL (runs; full-scale posterior heavy).** The full 5-model + cross-model BPSL chain **runs to
completion with `output_trajectory`**: a short chain (50 burn-in + 50 adapt + production over all 5
cell lines, 90 probit constraints) reports `Fitting complete`, exit 0, and writes both `samples.txt`
and `constraint_samples.txt` (100 samples each) — i.e. the cross-model constraint-satisfaction
bookkeeping that had crashed on `output_trajectory` now works.

History (resolved): under the one-model+conditions layout, `output_trajectory` made the base model run
every action suffix under every mutant (`WTn78g`, …); the sampler keyed `output_run_current` only by
the scored data-keys → `KeyError: 'WTn78gMEK_pRDS'`. Filed as
**[lanl/PyBNF#483](https://github.com/lanl/PyBNF/issues/483)** and **fixed upstream** (`ae21d212`, on
main — off-diagonal suffixes are skipped). (A separate v1.6.0 aMCMC cross-model BPSL regression,
`res.simdata=None`, was filed earlier as **[lanl/PyBNF#480](https://github.com/lanl/PyBNF/issues/480)**
and also fixed; the DE slug drives the identical 5-model + 90-constraint machinery to a finite objective.)

What is **not** run: the paper's full schedule — 5 chains × (25k burn-in + 25k adapt + 250k
production) ≈ 1.25M evaluations, ~5.15 days on 25 CPUs — so the Fig-7 marginals and Table-4 ESS/R̂ are
not reproduced at scale (heavy, not blocked).

**Verdict: PARTIAL** — the sampler + constraint tracking + `output_trajectory` run end-to-end; the
full-scale posterior is a compute-time matter, no longer a tooling block.

---

## Divergence & corrections
Scope matches the paper's Bayesian analysis. Corrections: KO typo, `.prop` reference rewrite, `WT.exp`
padding, CRLF→LF (as the DE slug). Scale-factor bounds left at the authors' values (no widening needed
on the K&L baseline). **Edition-2 conversion:** five per-variant models collapse to one `wt.bngl` +
four `condition:` perturbations (Mechanism A); T292D's fixed-`b5`/3 written absolute. **Upstream
filings:** `lanl/PyBNF#483` (`output_trajectory` off-diagonal suffixes) and `#480` (v1.6.0 aMCMC
cross-model `res.simdata=None`) — both fixed on main.

## Bottom line
A faithful edition-2 port of the paper's Bayesian-UQ setup with gold-standard Gate 0–2 grounding
(shared with the DE slug), a passing posterior-median reproduction that matches the legacy oracle to
the digit (F_quant 50.05, 87/90), and — after the lanl/PyBNF#483 fix — a chain that samples
end-to-end with `output_trajectory` and cross-model BPSL satisfaction tracking. Confidence sits at 85
(just below the DE slug) only because the full-scale 5-chain × 300k posterior was not run at HPC scale.
Most valuable next step: run the full schedule to reproduce the Fig-7 marginals and Table-4 ESS/R̂.
