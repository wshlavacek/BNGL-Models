# cstar_skmel133_bmra_exact — the paper's SKMEL-133 fit with the EXACT Eq.14 BMRA constraints

The PAPER-EXACT reconstruction of the SKMEL-133 fit: the model's **Eq.14 connection coefficients
`r_ij`** constrained TWO-SIDED inside the BMRA confidence intervals (Table S10, with_MYC) — the
paper's actual constraint object. The sibling [`../cstar_skmel133_bmra`](../cstar_skmel133_bmra)
imposes only the **sign** of each connection (a robust approximation); this slug imposes the full
numeric CI band. Both fit the full connection-coefficient set to the single-drug fold-change
training data with pyBioNetFit scatter-search + simplex.

> Rukhlenko OS, Halasz M, Rauch N, Zhernovkov V, Prince T, Wynne K, Maher S,
> Kashdan E, MacLeod K, Carragher NO, Kolch W, Kholodenko BN.
> **"Control of cell state transitions."** *Nature* 2022; **609**(7929):975–985.
> PMCID: [PMC9644236](https://pmc.ncbi.nlm.nih.gov/articles/PMC9644236/) ·
> DOI: [10.1038/s41586-022-05194-y](https://doi.org/10.1038/s41586-022-05194-y)

## What the paper did (Methods p.24), and what this slug reconstructs

> *"we constrained the parameters using the BMRA inferred connection coefficients within their
> confidence intervals, and used the pyBioNetFit software, which allows adding parameter
> constraints in the forms of inequalities to the parameter fitting process."*

The model was fit with **scatter search (pop 12, 50 iters) + simplex**, objective = **sum of
squares**, to the single-drug perturbation fold changes (the training set), with the connection
coefficients held **inside the BMRA-inferred confidence intervals**. This slug reproduces that
with the *exact* coefficients: `job_type = ss`, `refine = 1`, `chi_sq` (the `.exp` carry `_SD`),
the full connection-coefficient set free, and the BMRA CIs as BPSL `.prop` **two-sided inequality**
constraints on the model-computed `r_ij`.

## The exact constraint: `r_ij = C_i · L_ij`

Eq.13/14: `r_ij = −(∂f_i/∂x_j)/(∂f_i/∂x_i)·(x_j/x_i)|ss = ∂ln x_i/∂ln x_j` — the Jacobian-normalized
steady-state response of module *i* to module *j*. For this model's hyperbolic-crosstalk rate laws
it factors, per edge, into a regulator term and a target self-normalization:

```
L_ij = w(g_ij−1)/((1+w)(1+g_ij·w)),   w = x_j/K_ij        # x_j = regulator active form
C_i  = vn_i/((Kn_i+x_i)·sl_i)   (standard node)          # finite at x_i=0
C_CDK = KCycCDK/(KCycCDK+fCycD)   ;   C_MYC = −1
```

The 20 `r_<tgt>_<src>()` are emitted as model functions and each is pinned two-sided inside its
BMRA CI `[mean−std, mean+std]` (`bmra_rij.prop`), weight = 1/std. The decomposition is derived
from the model's own rate laws and **triple-verified** — analytic == an independent numeric
operational-MRA (5e-9) == BNG-emitted `r()` columns (1e-9); see `VALIDATION.md` and
`build_constraints.py`. **This is the whole point of the slug:** unlike the sign approximation
(trivially satisfied), the exact ±1 std band is genuinely binding.

**BMRA source:** `BMRA/results_SKMEL_133/SKMEL133_{rm,rs}_log_200_5K_withMyc.csv` == Table S10 —
the posterior **mean (`rm`) and std (`rs`)** of each connection (the model matches the *withMyc*
network). Embedded in `build_constraints.py`.

## Free parameters (43)

Per signaling edge **both** `g_XY` (the crosstalk strength, Eq.24 γ) **and** `K_XY` (its
half-saturation) — 20 g + 20 K = 40 — plus the 3 DPD force coefficients `beta_{mTOR,PKC,SRC}`.
Both g and K are freed (unlike the sign-approx, which freezes K) because the exact `r_ij` depends
on both (`L_ij = w(g−1)/((1+w)(1+g·w))`, `w = x_j/K`). Each `loguniform` over ≈±1 decade around the
authors' `SKMEL-133-3.bngl` nominal, spanning γ = 1 so each two-sided CI is reachable. The
restoring-force / landscape geometry and module rate constants stay at the published values.

## Training data (11 arms) & the setConcentration caveat

Training = the authors' single-drug fold-change `.exp` (24 h steady-state vs the no-drug baseline,
per-point `_SD`) for **all six inhibitors**: ERK/AKT/SRC (which act through the rate-law factor
`1/(1+I_<x>_conc)`, applied exactly by an edition-2 `setParameter` condition) **and** PKC/mTOR/CDK
(competitive binders, whose dose re-seeds a bound-inhibitor species via a quoted-species
`setConcentration`, per lanl-PyBNF#474). 11 arms = dose1 ×6 + dose2 ×5 (`dose2_CDKinh.exp` is
unreplicated — no `_SD` — so excluded from `chi_sq`). Drug **combinations** are the paper's
validation set (not fit here).

## The proliferative-state anchor

The BMRA was inferred in the **proliferative** state (SKMEL-133 is proliferative; the cSTAR DPD
places proliferation at `S > 0`), and the `r_ij` are the response ratios at that attractor. A free
(g, K) fit can otherwise drift the drug-free attractor into the **differentiated** basin (`S < 0`),
satisfying the numeric bands in the wrong regime. `bmra_rij.prop` therefore adds one principled
constraint, `Sval > 0`, keeping the constraint experiment on the proliferative side of the
separatrix. The published parameters satisfy it with wide margin (`Sval = +9.95`).

## Files

| file | role |
|---|---|
| `cstar_skmel133_bmra_exact.bngl` | model (= `../cstar_skmel133_bmra` with the 20 exact `r_<tgt>_<src>()` + the `cc_S*` DPD sign carriers) |
| `cstar_skmel133_bmra_exact.conf` | the fit: `ss` + simplex, `chi_sq`, per-observable `init`, 43 free params (g + K + β) |
| `cstar_skmel133_bmra_exact_check.conf` | `job_type = check` — constraint satisfaction at published params (40/44) |
| `bmra_rij.prop` | 40 two-sided signaling CIs + 3 DPD signs + the proliferative anchor (BPSL) |
| `build_constraints.py` | generator for `bmra_rij.prop` + the `r_ij()` functions (Table S10 provenance + the C_i·L_ij derivation) |
| `dose{1,2}_{ERK,AKT,SRC,PKC,mTOR}inh.exp`, `dose1_CDKinh.exp` | single-drug fold-change training data (authors' fit targets) |
| `prep_exp.py` | the mechanical `.exp` edits (src = authors' `SKMEL-133_preproc/`) |
| `make_reproduction.py` · `cstar_skmel133_bmra_exact_reproduction.png` | Gate-3a reproduction figure (11 arms + the exact BMRA-CI test) |
| `VALIDATION.md` | the earned-confidence gate write-up (**90/100**) |

## ⚠️ Native-only (not PEtab-v2-exportable)

Both `normalization = init` (fold-change data) **and** the BPSL `.prop` make this job native-only:
`export_job` raises `NotImplementedError`. Verify with tier-1 + `job_type = check` + a bounded fit
+ the fold-change reproduction, **not** the PEtab round-trip.

## Verification

- **`job_type = check` at published params** (`cstar_skmel133_bmra_exact_check.conf`):
  **`Satisfied 40 out of 44`** — 16/20 signaling edges + 3 DPD signs + the anchor hold; 4 edges
  (AKT←ERK k=1.08, AKT←PKC k=2.10, AKT←IRS k=2.17, MYC←CDK k=1.79) sit just outside the ±1 std band
  (all 20 in-band at k = 2.17). The exact band is genuinely binding.
- **Reproduction** (`cstar_skmel133_bmra_exact_reproduction.png`): the model at the published
  parameters reproduces the 11 fold-change arms to **12.9 % median relative error** (99 points),
  and the bottom panel shows each exact `r_ij` against its BMRA CI (16/20 in-band).
- **Bounded `ss` fit** (`constraint_scale = 5000`, 8 iters + simplex, ≈ 4 min): reaches objective
  **5429** (< the published-param objective ≈ 22 070 and < the pure `chi_sq` 6318) while staying
  proliferative (`Sval = +56`). It raises satisfaction to **42/44 (18/20 signaling in-band)**,
  bringing **3 of the 4** published-out-of-band edges into their CI (AKT←ERK, AKT←PKC, MYC←CDK) and
  keeping all 3 DPD signs. AKT←IRS remains just outside (the structurally-hardest AKT edge); a full
  pop-12 × 50-iter run closes the residual further.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Rukhlenko-2022/cstar_skmel133_bmra_exact
python build_constraints.py                          # (re)generate bmra_rij.prop + the r_ij() functions
pybnf -c cstar_skmel133_bmra_exact_check.conf        # Satisfied 40/44 at published params
pybnf -c cstar_skmel133_bmra_exact.conf              # the BMRA-CI-constrained fit (heavy: 43 params)
python make_reproduction.py                          # Gate-3a reproduction figure
```

## `_manifest.py` note (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='cstar_skmel133_bmra_exact', conf='cstar_skmel133_bmra_exact.conf', simulator='ode',
    observables=('FC_tIRS','FC_IRSI','FC_pERK','FC_pAKT','FC_pSRC','FC_pPKC','FC_pS6K','FC_pRB','FC_MYC'),
    system='cSTAR SKMEL-133 melanoma: FULL connection-coefficient fit to single-drug fold changes '
           'UNDER the EXACT Eq.14 connection coefficients constrained inside the BMRA CIs (Table '
           'S10) as two-sided BPSL inequalities (Rukhlenko 2022, PMC9644236, Methods p.24); ODE; '
           'scatter search + simplex; native-only (normalization=init + BPSL .prop)', heavy=True),
    # BPSL + normalization: assert export_job raises NotImplementedError and/or a `check`
    # run reports "Satisfied 40 out of 44" -- do NOT add a PEtab-lint test.
```
