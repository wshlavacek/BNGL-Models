# cstar_trkab_bmra_exact — the joint Trk fit with the EXACT Eq.14 BMRA constraints

The joint TrkA + TrkB reconstruction of the paper's fit, imposing the paper's real constraint
object: each model's **Eq.14 connection coefficients `r_ij`** constrained TWO-SIDED inside the BMRA
confidence intervals (Table S5, 10-min posterior). The sibling
[`../cstar_trkab_bmra`](../cstar_trkab_bmra) imposes only the **sign** of each connection; this slug
imposes the full numeric CI band on the exact coefficients.

> Rukhlenko OS, Halasz M, Rauch N, ... Kholodenko BN.
> **"Control of cell state transitions."** *Nature* 2022; **609**(7929):975–985.
> PMCID: [PMC9644236](https://pmc.ncbi.nlm.nih.gov/articles/PMC9644236/) ·
> DOI: [10.1038/s41586-022-05194-y](https://doi.org/10.1038/s41586-022-05194-y)

## What the paper did (Methods p.22), and what this slug reconstructs

> *"The training set included the time course of TrkA and TrkB phosphorylation measured by Western
> Blot and **10 min RPPA data** … we constrained the parameters using the BMRA inferred connection
> coefficients within their confidence intervals … the connection coefficients defined in Eq. 14
> must be within the confidence intervals … pyBioNetFit … inequalities … scatter search … simplex
> … The validation set consisted of **45 min RPPA data**."*

Both receptor models are fit **jointly** (shared core rate constants bind by name across the two)
with `job_type = ss`, `refine = 1` (simplex), objective `sos`, to the ligand-stimulation phospho
**time courses** (training), under the models' **exact Eq.14 `r_ij`** constrained two-sided inside
the **10-min** BMRA CIs. The 45-min points are the validation set.

## The exact constraint: `r_ij = C_i · L_ij`, and where it is evaluated

`r_ij = ∂ln x_i/∂ln x_j|st.st` factors per edge into `L_ij = w(g−1)/((1+w)(1+g·w))`, `w = x_j/K_ij`,
and a target self-normalization `C_i`:

- **standard node** (ERK, AKT, JNK, S6K — MM push-pull; ERK's active form is `ppERK`):
  `C_i = vn_i/((Kn_i+x_i)·sl_i)` (same form as the SKMEL exact slug);
- **RSK node** (`RSK←ERK` only): ERK enters RSK synthesis **linearly**, so a custom
  `r_RSK_ERK` is derived (see `build_constraints.py` / `VALIDATION.md`).

Both are **verified** to 1e-10 against an independent operational-MRA (`dev/papers/Rukhlenko-2022/
exact_build_wip_trk/verify_trk_rij.py`), and the BNG-emitted `r_<tgt>_<src>()` columns equal the
python closed form.

**Evaluated at the BASAL (no-ligand) steady state.** Eq.14 is a steady-state quantity, but the Trk
ligand response is transient/adaptive (`ppERK` peaks at 10 min then relaxes back to basal), so there
is no sustained stimulated steady state — the well-defined reference st.st. is the basal (Lig_on=0)
state, the pre-perturbation reference Eq.25 uses. The constraint experiments integrate each model
to its no-ligand steady state and read `r_ij()` there (`.prop` uses `at 100000`).

**BMRA source:** `BMRA/results/Trk{A,B}_{rm,rs}_10_log_200_5K.csv` == Table S5 (10-min posterior),
embedded in `build_constraints.py`.

## What is / isn't constrained

| model | constrained edges (exact two-sided r_ij) | documented, NOT constrained |
|---|---|---|
| TrkA | 10 into ERK/AKT/S6K | RTK←S6K |
| TrkB | 12 into ERK/AKT/JNK/S6K + RSK←ERK | RTK←{ERK,RSK} |

The **RTK/AdRTK receptor-dimerization node** (dim0/1/2 states, ligand binding, dimer degradation)
has no clean closed-form Eq.14 `C_i`, so its incoming edges are reported numerically (in the `.prop`
header) but not constrained — the honest scope used for the SKMEL IRS node. The **DPD row** is not
constrained (SVM-defined; the phospho time courses omit `Sval`).

## Free parameters (51)

The shared core rate constants (`vpTRK`, `vpRTK`, `vpERK`, `vpAKT`, bound by name across both
models) + per **constrained** edge BOTH `g_<tag>_<edge>` and `K_<tag>_<edge>` (the exact `r_ij`
depends on both) + the RSK synthesis params (`kpRSKA`; `kpRSKB`, `KpRSKB` for the constrained TrkB
`RSK←ERK`). Each `loguniform` over ≈±1 decade around the authors' `Trk{A,B}_S_model.bngl` nominal.

## Files

| file | role |
|---|---|
| `cstar_trk{a,b}_bmra_exact.bngl` | the two models (= sign-approx models with the exact `r_<tgt>_<src>()` functions) |
| `cstar_trkab_bmra_exact.conf` | the joint fit: `ss` + simplex, `sos`, per-observable `init`, 51 free params |
| `cstar_trkab_bmra_exact_check.conf` | `job_type = check` — constraint satisfaction at published params (40/46) |
| `bmra_{A,B}_rij.prop` | the two-sided 10-min BMRA CIs per model (+ the documented RTK edges) |
| `build_constraints.py` | generator for the props + the `r_ij()` functions (Table S5 provenance + the C_i·L_ij + RSK derivations) |
| `cstar_trk{a,b}.exp` | ligand-stimulation phospho time-course training data (authors' targets) |
| `make_reproduction.py` · `cstar_trkab_bmra_exact_reproduction.png` | Gate-3a figure (time courses + the exact BMRA-CI test) |
| `VALIDATION.md` | the earned-confidence gate write-up (**87/100**) |

## ⚠️ Native-only (not PEtab-v2-exportable)

`normalization = init` (fold-change data) **and** the BPSL `.prop` make this job native-only. Verify
with tier-1 + `job_type = check` + a bounded fit + the fold-change reproduction, **not** the PEtab
round-trip.

## Verification

- **`job_type = check` at published params** (`cstar_trkab_bmra_exact_check.conf`): **`Satisfied 40
  out of 46`** — TrkA 7/10 + TrkB 10/13 edges in-band; the 6 out are the tight-CI ones (S6K←AKT in
  both, TrkB RSK←ERK / JNK←ERK, TrkA AKT←ERK / ERK←JNK). The exact bands are genuinely binding.
- **Reproduction** (`cstar_trkab_bmra_exact_reproduction.png`): the models at the published
  parameters reproduce the ligand phospho time courses to **TrkA 18 %, TrkB 27 % median** (overall
  23.4 %), with the bottom panels showing each exact `r_ij` vs its 10-min BMRA CI.
- **Bounded joint `ss` fit** (`constraint_scale = 20`, 6 iters + simplex, ≈ 20 min): objective
  **119** (≪ the published-param ≈ 1025). It raises satisfaction to **43/46 (20/23 edges in-band)**,
  bringing the **tightest-CI** edges into band — TrkB `RSK←ERK` (0.074 → 0.769, band [0.734, 0.772])
  and `S6K←AKT` in both models. (Fold-change training does not pin absolute basal levels, so the fit
  shifts some basal activities to satisfy the steady-state `r_ij`; the paper's full fit adds the
  absolute WB pTrk course + a pop-20 × 50-iter budget. See `VALIDATION.md`.)

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Rukhlenko-2022/cstar_trkab_bmra_exact
python build_constraints.py                          # (re)generate the models + bmra_{A,B}_rij.prop
pybnf -c cstar_trkab_bmra_exact_check.conf           # Satisfied 40/46 at published params
pybnf -c cstar_trkab_bmra_exact.conf                 # the joint BMRA-CI-constrained fit (heavy: 51 params)
python make_reproduction.py                          # Gate-3a reproduction figure
```

## `_manifest.py` note (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='cstar_trkab_bmra_exact', conf='cstar_trkab_bmra_exact.conf', simulator='ode',
    observables=('FC_pTRK','FC_pAdRTK','FC_ppERK','FC_pAKT','FC_pJNK','FC_pRSK','FC_pS6K'),
    system='cSTAR SH-SY5Y TrkA+TrkB: JOINT fit of both receptor models to the ligand phospho time '
           'courses UNDER the EXACT Eq.14 connection coefficients constrained inside the 10-min BMRA '
           'CIs (Table S5) as two-sided BPSL inequalities (Rukhlenko 2022, PMC9644236, Methods p.22); '
           'ODE; scatter search + simplex; native-only (normalization=init + BPSL .prop)', heavy=True),
    # BPSL + normalization: assert export_job raises NotImplementedError and/or a `check`
    # run reports "Satisfied 40 out of 46" -- do NOT add a PEtab-lint test.
```
