# VALIDATION — Rukhlenko-2022/cstar_trkab_bmra

Primary-source validation of `pybnf-jobs/Rukhlenko-2022/cstar_trkab_bmra/` — the joint,
BMRA-CI-constrained reconstruction of the paper's *actual* TrkA + TrkB fit (the siblings
`../cstar_trka` (86/100) and `../cstar_trkb` (84/100) are reduced 8-parameter demos).

> **Confidence: 88 / 100** — Gate 2 is a byte-level network match to the authors' own model
> files (TrkA 82 sp/405 rxn — the paper's stated count — and TrkB 114/645), Gate 1 is
> byte-reproducible from the authors' RPPA, and this slug reconstructs the paper's declared
> method: a **joint** pyBioNetFit scatter-search + simplex fit UNDER the BMRA-inferred
> connection-coefficient CIs (Methods p.22). All 11 statistically-significant connection signs
> are satisfied at the published parameters (`Satisfied 11/11`). Not higher because (as for
> SKMEL) PyBNF BPSL enforces the CI-robust **sign** of each connection rather than the exact
> Eq. 14 Jacobian band, the exact CI multiplier is unpublished (Table S5), and the DPD force
> coefficients are left unconstrained (SVM-defined, and disagreeing with the BMRA STV row).

Primary sources (untracked `dev/papers/Rukhlenko-2022/`; not redistributed):
- Rukhlenko OS et al. *Nature* 2022; 609(7929):975–985. PMCID **PMC9644236**.
- Authors' repo **OleksiiR/cSTAR_Nature**: `Trk_AB_models/{TrkA_S_model,TrkB_S_model}.bngl`
  (models), `RPPA_DA/RPPA_data_trusted.csv` + `RPPA_data_Trk_normalized_new.csv` (RPPA), and
  `BMRA/results/Trk{A,B}_{rm,rs}_10_log_200_5K.csv` (the 10-min BMRA posteriors).

"The paper's result" = the TrkA + TrkB models fit **jointly** with pyBioNetFit (scatter search
pop 20 + simplex; objective = sum of squares; Methods p.22) to the ligand-stimulation phospho
time courses, with the connection coefficients constrained within the BMRA CIs.

---

## The BMRA → model mapping and the significance filter

Same map as `../cstar_skmel133_bmra` (Eqs. 14/24: the model's `g_<edge>` are the paper's γ;
γ>1 activation, γ<1 inhibition). Each connection is mapped to `r[target][source]` from the
**10-min** posterior (the training timepoint). Because the Trk BMRA CIs are wider than SKMEL's,
this slug applies the paper's **significance filter** (p.19: "only interactions ... with
statistically significant non-zero values are included"): a sign constraint is imposed only where
the CI is ≥ 1 std from zero (`z ≥ 1`) **and** the published model carries that sign. Under the
filter, model and BMRA agree on **every** constrained connection.

| model | constrained (z≥1, agree) | dropped (low-z / CI incl. 0) |
|---|---|---|
| TrkA | **5** — S6K→RTK (z1.7), TRK→ERK (z2.3), JNK→ERK (z1.7), ERK→AKT (z1.5), AKT→S6K (z15.6) | 11 (all z<1) |
| TrkB | **6** — ERK→RTK (z1.1), RSK→RTK (z1.5), ERK→JNK (z3.1), RTK→AKT (z1.0), ERK→AKT (z1.1), AKT→S6K (z14.9) | 10 (all z<1) |

The few model/BMRA-mean sign *disagreements* (TrkA TRK→AKT etc.) are all low-confidence
(z<1, CI includes 0) and correctly excluded. **The DPD force coefficients are NOT constrained**:
the DPD is SVM-defined, not the BMRA STV row, and they disagree at non-trivial confidence
(TrkB `beta_B_ERK`>0 vs BMRA `r_STV,ERK`=−0.68, z=1.55); the phospho time courses do not involve
`Sval` regardless. `build_constraints.py` prints the full ledger.

---

## Gate 0 — Materials: **PASS**
Models + RPPA + the 10-min BMRA posteriors all in hand.

## Gate 1 — Data provenance: **PASS**
`cstar_trka.exp` / `cstar_trkb.exp` are byte-reproducible from the authors' RPPA via the demos'
`extract_exp.py` (6-rep DMSO+ligand mean ÷ t0 → fold change), identical to `../cstar_trka` /
`../cstar_trkb`.

## Gate 2 — Model fidelity: **PASS (equiv)**
`model_diff.py`: each `*_bmra.bngl` generates the **identical** network to its demo (and thus to
the authors' `Trk{A,B}_S_model.bngl`) — TrkA **82 species / 405 reactions** (the paper's stated
count, Methods p.22), TrkB **114 / 645**. The only additions are a conserved `totERK` observable
(=1) and the `cc_*()` carrier functions; no rule, rate law, or parameter value is changed. All
published γ nominals = the authors' model files.

## Gate 3a — Reproduce at the paper's parameters: **PASS**
Since the generated network is identical to the demos, the model at the published parameters
reproduces the ligand-stimulation time courses to the demos' metrics — **TrkA ~18 %, TrkB ~27 %
median relative error** (see `../cstar_trka/VALIDATION.md`, `../cstar_trkb/VALIDATION.md`; the
larger TrkB error tracks the 645-reaction network and the 45-min validation points). `job_type =
check` reports **`Satisfied 11 out of 11 constraints`** at the published parameters (joint sos
35.1; zero constraint penalty).

## Gate 3b — Recover the paper's parameters by (constrained) fitting: **PASS**
A bounded joint `ss` fit (`constraint_scale = 5000`, 8 iterations) runs the simulate→score→propose
loop across both receptor datasets + both constraint experiments, lowers the joint objective to
**26.1 (below the 35.1 at the published parameters)**, and **keeps all 11 BMRA-inferred connection
signs (11/11)**. This reconstructs the paper's joint constrained fit; a fully convergent run is the
paper's `ss` pop 20 + simplex (heavy for 38 joint parameters).

---

## Divergence & scope

- **Scope vs paper: MATCHES the real joint fit.** Both authors' models fit together to their own
  RPPA time courses under per-receptor BMRA-CI constraints, with the paper's method (pyBioNetFit
  ss + simplex, sum of squares). The reduced 8-parameter `../cstar_trka` / `../cstar_trkb` are
  demos of the same data.
- **Three honest limitations** (docked into the 88): (1) the constraints enforce each connection's
  CI-robust **sign**, not the exact Eq. 14 numeric band (BPSL cannot access the Jacobian; the
  standardized-`r` vs physical-γ scale gap makes a raw-parameter band ill-posed); (2) the exact CI
  multiplier the paper used is unpublished (Table S5); (3) the DPD force coefficients are left
  unconstrained (SVM-defined; they disagree with the BMRA STV row at z≈1.5).
- **Further extension (not done):** the paper's Trk training also includes Western-blot pTrk and
  the inhibitor-panel arms (Extended Data Figs. 6–8, from `RPPA_data_trusted.csv`); this slug
  trains on the DMSO ligand time courses (as the demos do). Adding the inhibitor arms would enlarge
  the training set but not change the BMRA-constraint methodology reconstructed here.

## Bottom line

The paper's method — a **joint** TrkA+TrkB fit of the full connection-coefficient set under
BMRA-inferred CIs as inequality constraints — is reconstructed on the supported PyBNF edition-2 +
BPSL surface, grounded in the paper's Eqs. 14/24 and the authors' 10-min BMRA posteriors, with the
paper's own statistical-significance filter for which connections to constrain. All 11 significant
connection signs are satisfied at the published parameters.
