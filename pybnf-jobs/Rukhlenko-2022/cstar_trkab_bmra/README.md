# cstar_trkab_bmra — the paper's ACTUAL joint TrkA + TrkB fit (BMRA-CI-constrained)

The real-world reconstruction of the fit the paper actually performed for the Trk models:
the **TrkA and TrkB models fit jointly** with pyBioNetFit to the ligand-stimulation phospho
time courses, **under the BMRA-inferred connection-coefficient confidence intervals imposed as
inequality constraints**. The reduced siblings [`../cstar_trka`](../cstar_trka) /
[`../cstar_trkb`](../cstar_trkb) are 8-parameter *demos*; this slug extends them to the paper's
real, joint, constrained fit.

> Rukhlenko OS, Halász M, Rauch N, Zhernovkov V, Prince T, Wynne K, Maher S,
> Kashdan E, MacLeod K, Carragher NO, Kolch W, Kholodenko BN.
> **"Control of cell state transitions."** *Nature* 2022; **609**(7929):975–985.
> PMCID: [PMC9644236](https://pmc.ncbi.nlm.nih.gov/articles/PMC9644236/) ·
> DOI: [10.1038/s41586-022-05194-y](https://doi.org/10.1038/s41586-022-05194-y)

## What the paper did (Methods p.22), and what this slug reconstructs

> *"in addition to the training dataset, we constrained the parameters using the BMRA inferred
> connection coefficients within their confidence intervals ... we used ... pyBioNetFit ...
> which allows adding parameter constraints in the forms of inequalities to the parameter
> fitting process. ... Scatter search with a population size 20 was used to obtain the initial
> parameter set, and the simplex algorithm was used for the local refinement."*

The TrkA (**82 species, 405 reactions** — the count the paper reports) and TrkB (114 species,
645 reactions) models are fit **together**: they share the core enzymatic rate constants
(`vpTRK`/`vpRTK`/`vpERK`/`vpAKT`, bound by name across both models) but carry **separate
connection strengths** (`g_A_*` vs `g_B_*`) constrained by **separate BMRA posteriors**. The
joint fit couples the two datasets through the shared core while each receptor's connections are
held inside its own BMRA CIs.

## The BMRA → model mapping and the significance filter

Same map as [`../cstar_skmel133_bmra`](../cstar_skmel133_bmra) (Eq. 24 `α=(1+γ·x/K)/(1+x/K)`,
γ>1 activation / γ<1 inhibition; the model's `g_<edge>` ids *are* γ). Each connection is mapped
to the BMRA-inferred coefficient `r[target][source]` from the **10-min posterior** (the training
timepoint; `BMRA/results/Trk{A,B}_{rm,rs}_10_log_200_5K.csv`, order `TRK,ERK,AKT,JNK,S6K,RSK,RTK,STV`).

Unlike SKMEL (where all 20 connection signs agreed with BMRA), the Trk BMRA CIs are wider, so
this slug applies the paper's own **significance filter** — *"only interactions ... with
statistically significant non-zero values are included"* (p.19). A sign constraint is imposed
**only** where the BMRA CI is ≥ 1 std from zero (`z = |mean|/std ≥ 1`) **and** the published
model already carries that sign; low-confidence connections (CI includes 0) are left free. Under
this filter, model and BMRA agree on **every** constrained connection — **11 total** (5 for TrkA,
6 for TrkB). The few model/BMRA-mean sign disagreements (e.g. TrkA TRK→AKT) are all
low-confidence (z < 1, CI includes 0) and correctly excluded. `build_constraints.py` prints the
full kept/dropped ledger.

**The DPD (STV) force coefficients are NOT constrained.** The DPD is defined by the SVM state
signature, not directly by the BMRA STV row, and the two disagree at non-trivial confidence
(e.g. TrkB `beta_B_ERK` > 0 vs BMRA `r_STV,ERK` = −0.68, z = 1.55). The phospho time-course
training data does not involve `Sval` either, so the DPD force coefficients are held fixed.

Each constraint rides a carrier `cc_<edge>() = (g_<edge> − 1)·totERK` (`totERK ≡ 1`, a conserved
observable), evaluated at each model's no-ligand steady state.

## Free parameters (38)

The **full connection-coefficient set of both receptors** (16 `g_A_*` + 16 `g_B_*`) + the
4 shared core rate constants + the two receptor-specific RSK rates (`kpRSKA`, `kpRSKB`). Each
`loguniform` over ≈±1 decade around the published value, spanning γ = 1 where a constraint
applies. `K_XY`, ligand kinetics, and the DPD landscape are held at the authors' values.

## Training data

The two ligand-stimulation phospho fold-change time courses (0/10/45 min, 7 readouts each),
byte-reproducible from the authors' RPPA (via `../cstar_trka/extract_exp.py`,
`../cstar_trkb/extract_exp.py`): `cstar_trka.exp` (NGF), `cstar_trkb.exp` (BDNF). The paper's
validation set (45-min RPPA for the remaining modules) is not part of the training fit.

## Files

| file | role |
|---|---|
| `cstar_trka_bmra.bngl`, `cstar_trkb_bmra.bngl` | models (= demos + `totERK` + `cc_*()` carriers; generated network **identical** to the demos: 82/405 and 114/645) |
| `cstar_trkab_bmra.conf` | the joint fit: 2 models, `ss` + simplex, `sos`, per-observable `init`, 38 free params |
| `cstar_trkab_bmra_check.conf` | `job_type = check` — constraint satisfaction at published params |
| `bmra_A.prop`, `bmra_B.prop` | the 5 + 6 BMRA-CI sign constraints (BPSL), weighted by `z` |
| `build_constraints.py` | generator (embeds the 10-min BMRA posteriors; prints the kept/dropped ledger) |
| `cstar_trka.exp`, `cstar_trkb.exp` | ligand-stimulation phospho time-course training data |

## ⚠️ Native-only (not PEtab-v2-exportable)

Both `normalization = init` and the BPSL `.prop` make this job native-only. Verify with tier-1
+ `job_type = check` + a bounded fit, **not** the PEtab round-trip.

## Verification

- **Tier-1** (`check_conf.py`): edition 2, `ss` resolves, **2 models**, **2 BPSL constraint
  sets**, **38 free params** bind by id, no `__FREE`; native-only. **PASS.**
- **`job_type = check` at published params** (`cstar_trkab_bmra_check.conf`):
  **`Satisfied 11 out of 11 constraints`** (5 TrkA + 6 TrkB), joint objective 35.1 — the models
  satisfy every statistically-significant BMRA connection sign. **PASS.**
- **Model fidelity**: `model_diff.py` — each `*_bmra.bngl` generates the **identical** network to
  its demo (TrkA 82 sp/405 rxn, TrkB 114/645), so the `totERK`/`cc_*()` additions do not change
  the dynamics. Gate-3a reproduction therefore equals the demos': TrkA ~18 %, TrkB ~27 % median
  relative error at the published parameters (see `../cstar_trka/VALIDATION.md`,
  `../cstar_trkb/VALIDATION.md`).
- **Bounded joint `ss` fit** (`constraint_scale = 5000`, 8 iterations): the simulate→score→propose
  loop runs across both receptor datasets + both constraint experiments, lowers the joint objective
  to **26.1 (< 35.1 at the published parameters)**, and **keeps all 11 BMRA-inferred connection
  signs (11/11)**. A fully convergent run is the paper's `ss` pop 20 + simplex (heavy: 38 params).

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Rukhlenko-2022/cstar_trkab_bmra
python build_constraints.py                          # (re)generate bmra_A/B.prop
pybnf -c cstar_trkab_bmra_check.conf                 # Satisfied 11/11 at published params
pybnf -c cstar_trkab_bmra.conf                       # the joint BMRA-constrained fit (heavy: 38 params)
```

## `_manifest.py` note (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='cstar_trkab_bmra', conf='cstar_trkab_bmra.conf', simulator='ode',
    observables=('FC_pTRK','FC_pAdRTK','FC_ppERK','FC_pAKT','FC_pJNK','FC_pRSK','FC_pS6K'),
    system='cSTAR TrkA+TrkB joint fit: full connection-coefficient set fit to NGF/BDNF phospho '
           'time courses UNDER BMRA-inferred connection-coefficient CIs as BPSL inequality '
           'constraints (Rukhlenko 2022, PMC9644236, Methods p.22); two ODE models, scatter '
           'search + simplex; native-only (normalization=init + BPSL .prop)', heavy=True),
    # BPSL + normalization: assert export_job raises NotImplementedError and/or `check` reports
    # "Satisfied 11 out of 11" -- do NOT add a PEtab-lint test.
```
