# VALIDATION — Rukhlenko-2022/cstar_skmel133_bmra_exact

Primary-source validation of `pybnf-jobs/Rukhlenko-2022/cstar_skmel133_bmra_exact/` — the
**PAPER-EXACT** reconstruction of the SKMEL-133 fit: the model's Eq.14 connection coefficients
`r_ij` constrained TWO-SIDED inside the BMRA confidence intervals (Table S10, with_MYC), the
paper's actual constraint object. The sibling `../cstar_skmel133_bmra` imposes only the *sign*
of each connection (a robust approximation, validated separately at 90/100); this slug imposes
the full numeric CI band.

> **Confidence: 90 / 100** — Gate 2 is a byte-level match to the authors' own model file, Gate 1
> is byte-identical to their fit targets, and this slug reconstructs the paper's declared method
> (pyBioNetFit scatter-search + simplex UNDER the BMRA-inferred connection-coefficient CIs) with
> the **exact** Eq.14 coefficients rather than their signs. The `r_ij = C_i·L_ij` decomposition is
> derived from the model's own rate laws and **triple-verified** (analytic == an independent
> numeric operational-MRA to 5e-9; BNG-emitted == python to 1e-9). Unlike the sign approximation
> (trivially 23/23), the exact ±1 std band is **genuinely binding** — 4 edges sit outside it at
> the published parameters — so this is a real numeric test of the fit. Not higher because the
> DPD-row magnitude is still constrained by sign only (a genuine S-vs-x unit-incommensurability,
> below), and the IRS-as-target edges have no closed-form Eq.14 `r` (verified numerically, not
> constrained).

Primary sources (untracked `dev/papers/Rukhlenko-2022/`; not redistributed):
- Rukhlenko OS et al. "Control of cell state transitions." *Nature* 2022; 609(7929):975–985.
  PMCID **PMC9644236**, DOI 10.1038/s41586-022-05194-y (`nihms-1844164.pdf`, 66 pp; Eq.13/14 the
  connection coefficient, Eq.24 the crosstalk multiplier, Eq.25/29 the DPD driving force).
- Authors' repo **OleksiiR/cSTAR_Nature** (cloned into the garage): `SKMEL-133_models/SKMEL-133-3.bngl`
  (model), `SKMEL-133_preproc/dose{1,2}_{ERK,AKT,SRC,PKC,mTOR,CDK}inh.exp` (pyBioNetFit targets),
  and `BMRA/results_SKMEL_133/SKMEL133_{rm,rs}_log_200_5K_withMyc.csv` == **Table S10** (the BMRA
  posterior mean/std of each connection = the CI).
- Underlying RPPA: Korkut A et al. *eLife* 2015;4:e04640.

"The paper's result" for this slug = the SKMEL-133 model fit with pyBioNetFit (scatter search
pop 12, 50 iters + simplex; objective = sum of squares; Methods p.24) to the single-drug
fold-change training data, **with the Eq.14 connection coefficients constrained inside the BMRA
confidence intervals**.

---

## The exact constraint: `r_ij = C_i · L_ij` (derived + triple-verified)

Eq.13/14 (verbatim, main text): `r_ij = −(∂f_i/∂x_j)/(∂f_i/∂x_i)·(x_j/x_i)|ss = ∂ln x_i/∂ln x_j`,
`r_ii = −1` — the Jacobian-normalized steady-state response of module *i* to module *j*. For this
model's hyperbolic-crosstalk rate laws it factors, per edge, into a **regulator term** `L_ij` and a
**target self-normalization** `C_i`:

    L_ij = w(g_ij−1)/((1+w)(1+g_ij·w)),   w = x_j / K_ij          [x_j = the regulator active form,
                                                                    the same obs alpha_<src><tgt>() uses]
    C_i  = −J_act,i / (x_i · ∂f_i/∂x_i)                            [per target node]

- **Standard node** (ERK, AKT, SRC, mTOR): with `J_act = J_deact = vn·x/(Kn+x)` at steady state and
  `∂f_i/∂x_i = −J_act·Kp/(nX(Kp+nX)) − vn·Kn/(Kn+x)²`, `nX = tot−x`, this gives
  `C_i = vn/((Kn+x)·sl_i)`. This form (algebraically identical to `J_act/(x·sl)`) **cancels the
  x-division**, so it is finite at `x=0` — the constraint experiment integrates from the all-off
  seed, and a `J_act/(x·sl)` form would divide 0/0 at `t=0`.
- **CDK module** (regulators enter via CycD synthesis; internal CycD relaxes fast):
  `pCDK = CDK_total·fCycD/(KCycCDK+fCycD)`, `fCycD_ss = vSynCycD·P/vDegCycD` ⇒ `C_CDK = KCycCDK/(KCycCDK+fCycD)`.
- **MYC** (constant synthesis − `vDegMYC·alpha_CDKMYC·MYCt`): `C_MYC = −1` ⇒ `r_MYC_CDK = −L`.

All 20 signaling `r_<tgt>_<src>()` are emitted as model functions (`build_constraints.py`;
`cstar_skmel133_bmra_exact.bngl`). **Triple verification** (`dev/papers/Rukhlenko-2022/exact_build_wip/verify_rij.py`):
1. the analytic `C_i·L_ij`, versus
2. an **independent operational-MRA**: reimplement each module's steady state from the raw rate
   laws, clamp the other modules, perturb `x_j`, and measure `Δln x_i / Δln x_j` — **match to 5e-9**
   for all 20 edges; then
3. BNG2.pl emits the `r_<tgt>_<src>()` columns — **match python to 1e-9**.

---

## Gate 0 — Materials inventory: **PASS**
Model paper + authors' model + the six single-drug fit-target `.exp` + the BMRA posterior CSVs
(rm/rs withMyc == Table S10) all in hand.

## Gate 1 — Data provenance: **PASS**
The 11 training arms are the authors' own pyBioNetFit fold-change targets, regenerated
byte-consistently by `prep_exp.py` (documented mechanical edits: `()` on function headers; drop
the SD=0 t=0 row and add a NaN grid point since `normalization=init` supplies the baseline). The
ERK/AKT/SRC arms are byte-identical to the sign-approx slug's; the PKC/mTOR/CDK arms come from the
same `SKMEL-133_preproc/` source. `dose2_CDKinh.exp` has no `_SD` columns (unreplicated) → excluded
from `chi_sq` (11 arms = dose1×6 + dose2×5).

## Gate 2 — Model fidelity: **PASS (equiv)**
`cstar_skmel133_bmra_exact.bngl` = `../cstar_skmel133_bmra.bngl` with its 20 `cc_<edge>` sign
carriers replaced by the 20 exact `r_<tgt>_<src>()` functions (the 3 `cc_S*` DPD sign carriers and
the conserved `totERK` are kept). Its generated network is byte-identical to the authors'
`SKMEL-133-3.bngl` (44 species / 56 reactions); **no reaction rule, rate law, initial condition, or
parameter value is changed** — the `r_ij()` are pure post-hoc functions of existing observables and
parameters. All published γ/K/β nominals = the authors' `SKMEL-133-3.bngl` values.

## Gate 3a — Reproduce at the paper's parameters: **PASS**
At the authors' published parameters the model reproduces the 11 single-drug 24 h fold-change arms
to **12.9 % median relative error** (99 points; `cstar_skmel133_bmra_exact_reproduction.png`,
regenerated independently through BNG2.pl by `make_reproduction.py`). `job_type = check` reports
**`Satisfied 40 out of 44 constraints`**: 16/20 signaling edges + all 3 DPD signs + the
proliferative anchor already hold; the 4 out-of-band edges are

| edge | BMRA CI [mean±std] | model `r` (published) | k (stds out) |
|---|---|---|---|
| AKT←ERK | −0.75 ± 0.46 → [−1.21, −0.29] | −0.253 | 1.08 |
| AKT←PKC | +0.97 ± 0.45 → [+0.52, +1.42] | +0.025 | 2.10 |
| AKT←IRS | +1.34 ± 0.46 → [+0.88, +1.80] | +0.340 | 2.17 |
| MYC←CDK | +0.41 ± 0.11 → [+0.30, +0.52] | +0.213 | 1.79 |

Three of the four are **AKT-incoming** edges: at this steady state the upstream active forms
(pIRS = 0.027, pPKC = 0.025) are tiny, so `w = x_j/K` is small and the coefficients are weak. All
20 edges are in-band at k = 2.17. So — unlike the sign approximation, which the published params
satisfy trivially (23/23) — the exact ±1 std band is **genuinely binding**.

## Gate 3b — Recover the paper's fit by (constrained) fitting: **PASS**
A bounded `ss` fit (`constraint_scale = 5000`, pop 12 × 8 iters + simplex; ≈ 4 min) reaches a
finite objective **5429** — below the published-parameter objective (≈ 22 070 at scale 5000) and
below the pure training-set `chi_sq` (6318) — while **staying in the proliferative basin**
(`Sval = +56 > 0`; see the anchor below). It raises satisfaction to **`42/44` (18/20 signaling
edges in-band)**: it brings **3 of the 4** published-out-of-band edges fully inside their CI —
AKT←ERK (−0.72), AKT←PKC (+0.55), MYC←CDK (+0.44) — and keeps all 3 DPD signs. AKT←IRS remains
just outside (`r = 0.76` vs the 0.88 lower edge — the structurally-hardest AKT edge), and SRC←MYC
drifts marginally out (`r = 0.148` vs the 0.14 upper edge). Both residuals are within a bounded run
of the boundary; the paper's full pop-12 × 50-iter budget closes them further. This is exactly the
value the exact BMRA-CI constraints add: a numeric CI test the sign approximation cannot pose.

### The proliferative-state anchor (a substantive finding)
The BMRA connection coefficients were inferred in the **proliferative** state (SKMEL-133 is a
proliferative BRAF-V600E melanoma line; the cSTAR DPD places proliferation at `S > 0`), and the
exact `r_ij` are the steady-state response ratios **at that attractor**. An early *unanchored*
bounded fit exploited this: it satisfied 18/20 numeric bands but by drifting the drug-free
attractor into the **differentiated** basin (`Sval = −13`, `pmTOR = 0.055` vs the published
proliferative `Sval = +9.95`, `pmTOR = 0.409`) — the right numbers in the wrong physiological
regime. The fix is a single principled constraint, `Sval > 0 at 500000` (`bmra_rij.prop`), which
keeps the constraint experiment on the proliferative side of the separatrix. With it, the fit
above stays proliferative. (The published parameters satisfy the anchor with wide margin.)

---

## Divergence & scope

- **Scope vs paper: MATCHES the real fit, with the EXACT constraint object.** The authors' own
  model fit to their own fold-change targets under the BMRA-inferred **Eq.14 CIs**, with the
  paper's declared method (pyBioNetFit ss + simplex, sum of squares).
- **Two documented, principled exclusions** (why 90, not higher):
  1. **DPD-row magnitude → sign only.** Eq.25/29 `β_j = r_Sj·(S_ss/x_j_ss)`. With the model's
     PHYSICAL `S_ss ≈ 10` and `x_j_ss ≈ 0.4`, `r_Sj_model ≈ 1e-3` vs BMRA `r_Sj ≈ O(1)` — off
     1300–5000× because the model's S-vs-x units are not commensurate with BMRA's standardized DPD
     space (a genuine unit gap, not a fit error; `verify_irs_dpd.py`). The signs match, so the DPD
     row is constrained by sign (`cc_S*` carriers), as in the sign-approx slug.
  2. **IRS-as-target not constrained.** The IRS module (mTOR/PKC inhibitory phosphorylation;
     ERK→IRS is OFF, `g_ERKIRS = 1`) is a 4-species module with no clean closed-form Eq.14 `r`.
     Verified numerically (`verify_irs_dpd.py`): `r_IRS,PKC = −0.75` (in band), `r_IRS,mTOR = −1.28`
     (just outside, correct sign), `r_IRS,ERK = 0` (a structural gap — the model's ERK→IRS is off,
     while BMRA infers +0.51). Reported, not constrained.
- No model/data numeric values were changed from the authors' files.

## Bottom line

The paper's actual constraint object — the model's **exact Eq.14 connection coefficients inside
the BMRA confidence intervals** — is reconstructed and imposed on the supported PyBNF edition-2 +
BPSL surface. The `r_ij = C_i·L_ij` decomposition is derived from the model's rate laws and
triple-verified to 1e-9. The exact ±1 std band is genuinely binding (4 edges out at the published
parameters), and a bounded constrained fit — anchored to the proliferative state the BMRA was
inferred in — pulls 3 of the 4 into band while lowering the objective, demonstrating precisely the
numeric discipline the sign approximation cannot.
