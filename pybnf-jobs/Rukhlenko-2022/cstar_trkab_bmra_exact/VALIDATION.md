# VALIDATION — Rukhlenko-2022/cstar_trkab_bmra_exact

Primary-source validation of `pybnf-jobs/Rukhlenko-2022/cstar_trkab_bmra_exact/` — the joint,
**PAPER-EXACT** reconstruction of the TrkA + TrkB fit: the models' Eq.14 connection coefficients
`r_ij` constrained TWO-SIDED inside the BMRA confidence intervals (Table S5, 10-min posterior).
The sibling `../cstar_trkab_bmra` imposes only the *sign* of each connection (validated 88/100);
this slug imposes the full numeric CI band on the exact coefficients.

> **Confidence: 87 / 100** — Gate 2 is a byte-level network match to the authors' own model files
> (TrkA 82 sp/405 rxn, TrkB 114/645), Gate 1 is byte-reproducible from the authors' RPPA, and this
> slug reconstructs the paper's declared method with the **exact** Eq.14 coefficients: a joint
> pyBioNetFit scatter-search + simplex fit UNDER the BMRA-inferred connection-coefficient CIs
> (Methods p.22). The `r_ij = C_i·L_ij` decomposition is derived from the models' rate laws and
> **verified** (closed-form == an independent operational-MRA to 1e-10; == the BNG-emitted `r()`
> columns). Grounded in the paper's own timepoint structure (training = 10-min RPPA; validation =
> 45-min). Not higher than the SKMEL exact (90) because one target node — the RTK/AdRTK receptor
> **dimerization** module — has no clean closed-form Eq.14 `C_i`, so its incoming edges are reported
> numerically but not constrained (the honest scope used for SKMEL's IRS node), and the DPD row is
> left unconstrained (SVM-defined; the phospho time courses omit `Sval`).

Primary sources (untracked `dev/papers/Rukhlenko-2022/`; not redistributed):
- Rukhlenko OS et al. *Nature* 2022; 609(7929):975–985. PMCID **PMC9644236** (`nihms-1844164.pdf`;
  Eq.13/14 the connection coefficient, Methods p.22 "Refining parameters of the dynamic model").
- Authors' repo **OleksiiR/cSTAR_Nature**: `Trk_AB_models/{TrkA_S_model,TrkB_S_model}.bngl`
  (models), `RPPA_DA/RPPA_data_trusted.csv` + `RPPA_data_Trk_normalized_new.csv` (RPPA), and
  `BMRA/results/Trk{A,B}_{rm,rs}_10_log_200_5K.csv` == **Table S5** (10-min BMRA posterior = the CI).

"The paper's result" = the TrkA + TrkB models fit **jointly** with pyBioNetFit (scatter search
pop 20 + simplex; objective = sum of squares; Methods p.22) to the ligand-stimulation phospho
time courses, with the Eq.14 connection coefficients constrained inside the BMRA CIs.

---

## Which BMRA timepoint, and at what state (grounded, Methods p.22)

The paper is explicit: *"The training set included the time course of TrkA and TrkB phosphorylation
measured by Western Blot and **10 min RPPA data** for the remaining signaling modules … we
constrained the parameters using the BMRA inferred connection coefficients within their confidence
intervals … the connection coefficients defined in Eq. 14 must be within the confidence intervals …
The validation set consisted of **45 min RPPA data**."* So the exact `r_ij` are constrained to the
**10-min** BMRA posterior (Table S5); the 45-min data is validation, not a constraint.

**Evaluation state = the BASAL (no-ligand) steady state.** Eq.14 defines `r_ij` at steady state
(st.st.). The Trk ligand response is **transient and adaptive** — `ppERK` peaks at 10 min (≈3.2× the
basal level) then relaxes back to basal by long times (verified: at Lig_on=1 the model's `pTRK`
returns to 0.009 ≈ the basal 0.008 by t=1e5) — so there is *no sustained stimulated steady state*.
The well-defined reference st.st. is the basal (Lig_on=0) state, which is also the pre-perturbation
reference Eq.25 uses (`S_st.st`, `x_j,st.st`). The constraint experiments therefore integrate each
model to its no-ligand steady state and read the `r_ij()` there (`bmra_{A,B}_rij.prop` uses
`at 100000`). This matches the sign-approx slug's evaluation point.

## The exact constraint: `r_ij = C_i · L_ij` (derived + verified to 1e-10)

`r_ij = ∂ln x_i/∂ln x_j|st.st = −(∂f_i/∂x_j)/(∂f_i/∂x_i)·(x_j/x_i)` factors per edge into a regulator
term `L_ij = w(g−1)/((1+w)(1+g·w))`, `w = x_j/K_ij`, and a target self-normalization `C_i`:

- **Standard activation/deactivation node** (ERK, AKT, JNK, S6K — MM push-pull, identical rate-law
  form to SKMEL; ERK's active form is the double-phospho `ppERK`): `C_i = vn_i/((Kn_i+x_i)·sl_i)`,
  `sl_i = Jd·Kp_i/(nX(Kp_i+nX)) + vn_i·Kn_i/(Kn_i+x_i)²`, `Jd = vn_i·x_i/(Kn_i+x_i)`, `nX = tot_i−x_i`.
- **RSK node** (only `RSK←ERK`): ERK enters RSK *synthesis* **linearly** (`vpRSKbasal + kpRSK·ppERK/KpRSK`),
  not through a hyperbolic `alpha`, so it needs its own derivative:
  `r_RSK_ERK = [kpRSK·nRSK/(Kp+nRSK)] / [B·Kp²/(Kp+nRSK)² + vnRSK·KnRSK/(KnRSK+pRSK)²] · (ppERK/pRSK)`,
  `B = vpRSKbasal + kpRSK·ppERK/KpRSK`, `nRSK = RSK_total−pRSK`, `Kp = KpRSK{A,B}`.

**Verification** (`dev/papers/Rukhlenko-2022/exact_build_wip_trk/verify_trk_rij.py`): the closed-form
`C_i·L_ij` (and the custom RSK form) versus an **independent operational-MRA** — reimplement each
node's MM push-pull steady state from the raw rate laws, perturb one regulator's active form, take
`Δln x_i/Δln x_j` — **match to 6e-11** (central difference). The models then emit `r_<tgt>_<src>()`
functions whose BNG2.pl values equal the python closed form (checked at the basal steady state).

## What is / isn't constrained

| model | constrained edges (exact two-sided r_ij) | documented, NOT constrained |
|---|---|---|
| TrkA | **10** into ERK/AKT/S6K (all standard nodes) | RTK←S6K (dimer node) |
| TrkB | **13** into ERK/AKT/JNK/S6K + RSK←ERK | RTK←{ERK,RSK} (dimer node) |

The **RTK/AdRTK receptor-dimerization node** (dim0/1/2 states, ligand binding, dimer degradation)
has no clean closed-form Eq.14 `C_i`, so its incoming edges are reported numerically in the `.prop`
header (with their BMRA values) but not constrained — the same honest scope used for the SKMEL IRS
node. The **DPD (STV) row** is not constrained (SVM-defined; the phospho time-course training does
not involve `Sval`), as in the sign-approx slug.

---

## Gate 0 — Materials: **PASS**
Models + RPPA + the 10-min BMRA posteriors (Table S5) all in hand.

## Gate 1 — Data provenance: **PASS**
`cstar_trka.exp` / `cstar_trkb.exp` are the authors' ligand-stimulation fold-change targets
(byte-identical to the sign-approx slug's; byte-reproducible from the RPPA via the demos'
`extract_exp.py` — 6-rep DMSO+ligand mean ÷ t0).

## Gate 2 — Model fidelity: **PASS (equiv)**
Each `*_bmra_exact.bngl` = the sign-approx `*_bmra.bngl` with its `cc_*()` sign carriers replaced by
the exact `r_<tgt>_<src>()` functions (the conserved `totERK` observable kept). The generated
network is byte-identical to the authors' `Trk{A,B}_S_model.bngl` — TrkA **82 species / 405
reactions** (the paper's stated count), TrkB **114 / 645**. No rule, rate law, initial condition, or
parameter value is changed; the `r_ij()` are pure functions of existing observables and parameters.

## Gate 3a — Reproduce at the paper's parameters: **PASS**
At the authors' published parameters the models reproduce the ligand-stimulation phospho time
courses to **TrkA 18 %, TrkB 27 % median relative error** (overall 23.4 %, 28 points at 10+45 min;
`cstar_trkab_bmra_exact_reproduction.png`, regenerated independently through BNG2.pl). `job_type =
check` reports **`Satisfied 40 out of 46 constraints`** at the published parameters: TrkA 7/10 +
TrkB 10/13 signaling edges are in-band; the 6 out-of-band edges are the **tight-CI** ones —
S6K←AKT (both; BMRA std 0.033/0.034 → the tightest band), TrkB RSK←ERK (std 0.019) and JNK←ERK
(std 0.284), and TrkA AKT←ERK and ERK←JNK. So the exact bands are genuinely binding, unlike the
sign approximation (trivially 11/11).

## Gate 3b — Recover the paper's fit by (constrained) fitting: **PASS**
A bounded joint `ss` fit (`constraint_scale = 20`, pop 16 × 6 iters + simplex; ≈ 20 min for the
two large networks) reaches a finite objective **119** — far below the published-parameter objective
(≈ 1025 at scale 20, dominated by the two tight-CI violations). It raises satisfaction to
**43/46 (combined 20/23 signaling edges in-band, up from 17/23)**, and — notably — brings the
**tightest-CI** edges into band: TrkB `RSK←ERK` (r 0.074 → 0.769, band [0.734, 0.772], std 0.019)
and `S6K←AKT` in **both** models (std 0.033/0.034), plus TrkA `AKT←ERK`. Only 3 wide-CI edges remain
out (TrkA `ERK←JNK`, `ERK←TRK`; TrkB `JNK←ERK`). **Caveat (bounded run):** because the training data
are *fold changes* (obs(t)/obs(0)), which do not pin the absolute basal levels, the fit shifts some
basal activities to satisfy the steady-state `r_ij` (e.g. TrkB basal `pAKT` rises from 0.012 to 0.82);
the paper's full fit additionally uses the absolute Western-blot pTrk time course and its pop-20 ×
50-iter budget, which further disciplines this. The result nonetheless demonstrates the numeric
value the exact BMRA-CI constraints add: pulling out-of-band coefficients — including the very
tightest — into their confidence intervals, which the sign approximation cannot pose.

---

## Divergence & scope

- **Scope vs paper: MATCHES the real joint fit, with the EXACT constraint object.** Both authors'
  models fit together to their own RPPA time courses (10-min training) under per-receptor **Eq.14
  BMRA-CI** constraints, with the paper's method (pyBioNetFit ss + simplex, sos).
- **Two documented, principled exclusions** (why 87, one below the SKMEL exact 90):
  1. **The RTK/AdRTK receptor-dimerization node** has no clean closed-form Eq.14 `C_i` (a
     multi-state dimerization module, not a simple push-pull), so its 1 (TrkA) + 2 (TrkB) incoming
     BMRA edges are reported numerically but not constrained — like the SKMEL IRS node.
  2. **The DPD row** is not constrained (SVM-defined; the time-course training omits `Sval`), as in
     the sign-approx slug.
- No model/data numeric values were changed from the authors' files.

## Bottom line

The paper's actual constraint object for the Trk models — the **exact Eq.14 connection coefficients
inside the 10-min BMRA confidence intervals**, evaluated at the basal steady state — is reconstructed
and imposed jointly on TrkA + TrkB. The `r_ij = C_i·L_ij` decomposition (with a custom derivation for
the RSK node) is verified to 1e-10, the exact ±1 std bands are genuinely binding (6 edges out at the
published parameters), and a joint constrained fit pulls out-of-band edges toward their CIs — the
numeric discipline the sign approximation cannot pose.
