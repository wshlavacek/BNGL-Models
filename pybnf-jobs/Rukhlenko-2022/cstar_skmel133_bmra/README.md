# cstar_skmel133_bmra — the paper's ACTUAL SKMEL-133 fit (BMRA-CI-constrained, full parameter set)

The real-world reconstruction of the fit the paper actually performed for SKMEL-133: the
**full connection-coefficient set** fit with pyBioNetFit to the single-drug fold-change
training data, **under the BMRA-inferred connection-coefficient confidence intervals imposed
as inequality constraints**. Its reduced sibling [`../cstar_skmel133`](../cstar_skmel133) is
an 8-parameter *demo*; this slug extends it to the paper's real fit.

> Rukhlenko OS, Halasz M, Rauch N, Zhernovkov V, Prince T, Wynne K, Maher S,
> Kashdan E, MacLeod K, Carragher NO, Kolch W, Kholodenko BN.
> **"Control of cell state transitions."** *Nature* 2022; **609**(7929):975–985.
> PMCID: [PMC9644236](https://pmc.ncbi.nlm.nih.gov/articles/PMC9644236/) ·
> DOI: [10.1038/s41586-022-05194-y](https://doi.org/10.1038/s41586-022-05194-y)

## What the paper did (Methods p.24), and what this slug reconstructs

> *"we constrained the parameters using the BMRA inferred connection coefficients within
> their confidence intervals, and used the pyBioNetFit software, which allows adding
> parameter constraints in the forms of inequalities to the parameter fitting process."*

The model was fit with **scatter search (population 12, 50 iterations) + simplex**, objective
= **sum of squares**, to the single-drug perturbation fold changes (the training set), while
the connection coefficients were held inside the **BMRA-inferred confidence intervals**. This
slug reproduces that: `job_type = ss`, `refine = 1`, `chi_sq` (the `.exp` carry `_SD`), the
full connection-coefficient set free, and the BMRA CIs as BPSL `.prop` inequality constraints.

## The BMRA → model mapping (the crux)

Three paper equations tie the BMRA network to the model:

- **Eq. 24** (crosstalk multiplier): `α_Y^X = (1 + γ·Y_a/K)/(1 + Y_a/K)`. The paper states
  **γ > 1 ⇒ activation, γ < 1 ⇒ inhibition, γ = 1 ⇒ no interaction**. The model's `g_<edge>`
  ids *are* these γ.
- **Eq. 14** (connection coefficient): `r_ij = −(∂f_i/∂x_j)/(∂f_i/∂x_i)·(x_j/x_i)|_ss`, a
  Jacobian-normalized ratio (`r_ii = −1`). The constraints "connection coefficients ... within
  the CIs" act on these model-computed `r_ij`.
- **Eq. 25** (DPD driving force): `σ(t) = Σ_j r_Sj·(S_ss/x_j_ss)·x_j(t)`, so the force
  coefficient `β_j = r_Sj·(S_ss/x_j_ss)` — the DPD row of the connection matrix.

**BMRA source:** `BMRA/results_SKMEL_133/SKMEL133_{rm,rs}_log_200_5K_withMyc.csv` — the
posterior **mean (`rm`) and std (`rs`)** of each connection coefficient (top-5000-of-200000
Occam posterior networks; the model matches the *withMyc* network). Matrix rows = target,
cols = source, order `ERK,AKT,mTOR,SRC,CDK,PKC,IRS,MYC,DPD`.

**Why sign constraints.** PyBNF BPSL constrains model **observables**, not the Eq. 14 Jacobian
ratio the paper bounds — so we enforce the robust, CI-supported content of each connection:
its **sign** (which side of γ = 1 it sits on, i.e. activation vs inhibition), weighted by the
BMRA confidence `z = |mean|/std`. The magnitudes are not bounded, because the BMRA `r` are
*standardized* (dimensionless) while the model γ/β are in physical units — even the published
parameters do not reproduce the BMRA `r` magnitudes (e.g. PKC→AKT model `r`≈0.02 vs BMRA 0.96),
so a numeric band on the raw parameters is not well-posed without the paper's (unpublished,
Table S10) normalization. **Verified: at the published parameters the model's γ-signs match all
20 BMRA-inferred signalling signs and the 3 DPD force-coefficient signs exactly** — the
`job_type = check` job reports **Satisfied 23/23**.

Each constraint rides a carrier function `cc_<edge>() = (g_<edge> − 1)·totERK` (`totERK ≡ 1`,
a conserved observable — BPSL needs an observable, and a parameter-only function is not emitted
to the gdat), whose sign is the connection's sign. For a connection acting on a **degradation**
term (CDK→MYC) the γ-side is flipped relative to the BMRA `r` sign. See `bmra_signs.prop` and
`build_constraints.py` (the generator, with the full rm/rs provenance and the sign logic).

## Free parameters (23)

The **full connection-coefficient set**: the 20 `g_XY` crosstalk strengths + the 3 DPD force
coefficients `beta_{mTOR,PKC,SRC}`. Each `loguniform` over ≈±1 decade around the published
value, bracketing γ = 1 so the sign constraint is genuinely active. The `K_XY` and the
restoring-force / landscape geometry are held at the authors' published values (the BMRA
constraints act on the connection strengths). Published values = the authors'
`SKMEL-133-3.bngl` nominals.

## Training data (6 arms) & the mTOR/PKC/CDK caveat

Training = the authors' single-drug fold-change `.exp` (24 h steady-state, vs the no-drug DMSO
baseline, per-point `_SD`): **ERK / AKT / SRC × dose1 / dose2** — the three inhibitors that act
purely through the rate-law factor `1/(1+I_<x>_conc)` with no sequestering binding, so an
edition-2 Condition (`setParameter`) applies them exactly. The **mTOR/PKC/CDK** single-drug
arms act by **competitive binding** (their dose must re-initialise a seed-species concentration
via `setConcentration`, which an edition-2 Condition cannot emit — the igf1r / lanl-PyBNF#474
limitation), so they are omitted; the drug-**combination** arms are the paper's *validation* set
either way (not fit here). The 6 `.exp` are byte-identical to the authors' files (2 mechanical
edits documented in `prep_exp.py`; the 4 that overlap `../cstar_skmel133` are byte-identical to
that slug).

## Files

| file | role |
|---|---|
| `cstar_skmel133_bmra.bngl` | model (= `../cstar_skmel133` + `totERK` + 23 `cc_*()` carriers) |
| `cstar_skmel133_bmra.conf` | the fit: `ss` + simplex, `chi_sq`, per-observable `init`, 23 free params |
| `cstar_skmel133_bmra_check.conf` | `job_type = check` — constraint satisfaction at published params |
| `bmra_signs.prop` | the 23 BMRA-CI sign constraints (BPSL), weighted by confidence `z` |
| `build_constraints.py` | generator for `bmra_signs.prop` + the `cc_*()` functions (rm/rs provenance) |
| `dose{1,2}_{ERK,AKT,SRC}inh.exp` | single-drug fold-change training data (authors' fit targets) |
| `prep_exp.py` | the 2 mechanical `.exp` edits (src = authors' `SKMEL-133_preproc/`) |
| `make_reproduction.py` · `cstar_skmel133_bmra_reproduction.png` | Gate-3a reproduction figure |
| `VALIDATION.md` | the earned-confidence gate write-up |

## ⚠️ Native-only (not PEtab-v2-exportable)

Both `normalization = init` (fold-change data) **and** the BPSL `.prop` make this job
native-only: `export_job` raises `NotImplementedError`. Verify with tier-1 + `job_type = check`
+ a bounded fit + the fold-change reproduction, **not** the PEtab round-trip.

## Verification

- **Tier-1** (`check_conf.py`): edition 2, `job_type = ss` resolves, **1 BPSL constraint set**,
  **23 free params** bind by id, no `__FREE`; correctly reported native-only. **PASS.**
- **`job_type = check` at published params** (`cstar_skmel133_bmra_check.conf`):
  **`Satisfied 23 out of 23 constraints`**, training-set objective 12076.5 — the model
  satisfies every BMRA-inferred connection sign. **PASS.**
- **Reproduction** (`cstar_skmel133_bmra_reproduction.png`): the model at the published
  parameters reproduces the 6 single-drug fold-change arms to **≈ 11 % median relative error**
  (54 points), and every BMRA-inferred sign agrees. **PASS.**
- **Bounded `ss` fit** (`constraint_scale = 5000`, 8 iterations): the simulate→score→propose
  loop runs across all 6 training experiments + the constraint experiment and lowers the
  objective (12076.5 → 11881.8). Critically, it **keeps all 23 BMRA-inferred connection signs
  (23/23)** — `g_IRSERK` stays at 3.76 (activating), `g_ERKAKT` at 0.023 (inhibiting),
  `g_CDKMYC` at 0.197 (degradation) — whereas the reduced, *unconstrained* demo fit **flips**
  `g_IRSERK` from 4.22 to 0.92 (activation → inhibition). This is the point of the BMRA-CI
  constraints: they resolve the sloppiness / non-identifiability the demo exposed. (A full
  convergent fit is `ss` population 12 × 50 iterations + simplex — heavy for 23 parameters.)

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Rukhlenko-2022/cstar_skmel133_bmra
python build_constraints.py                          # (re)generate bmra_signs.prop
pybnf -c cstar_skmel133_bmra_check.conf              # Satisfied 23/23 at published params
pybnf -c cstar_skmel133_bmra.conf                    # the BMRA-constrained fit (heavy: 23 params)
python make_reproduction.py                          # Gate-3a reproduction figure
```

## `_manifest.py` note (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='cstar_skmel133_bmra', conf='cstar_skmel133_bmra.conf', simulator='ode',
    observables=('FC_tIRS','FC_IRSI','FC_pERK','FC_pAKT','FC_pSRC','FC_pPKC','FC_pS6K','FC_pRB','FC_MYC'),
    system='cSTAR SKMEL-133 melanoma: FULL connection-coefficient fit to single-drug fold '
           'changes UNDER BMRA-inferred connection-coefficient CIs as BPSL inequality '
           'constraints (Rukhlenko 2022, PMC9644236, Methods p.24); ODE; scatter search + '
           'simplex; native-only (normalization=init + BPSL .prop)', heavy=True),
    # BPSL + normalization: assert export_job raises NotImplementedError and/or a `check`
    # run reports "Satisfied 23 out of 23" -- do NOT add a PEtab-lint test.
```
