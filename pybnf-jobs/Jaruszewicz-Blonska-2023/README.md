# Jaruszewicz-Blonska-2023 — identifiable NF-κB model (PyBNF fitting jobs)

PyBNF **edition-2** parameter-fitting jobs derived from the model-reduction / identifiability study
of the canonical NF-κB signaling pathway:

> **Jaruszewicz-Błońska J, Kosiuk I, Prus W, Lipniacki T.** **"A plausible identifiable model of
> the canonical NF-κB signaling pathway."** *PLoS ONE* 2023; **18**(6):e0286416.
> DOI: [10.1371/journal.pone.0286416](https://doi.org/10.1371/journal.pone.0286416).

The paper reduces the 15-variable **Lipniacki et al. 2004** NF-κB model (two negative-feedback loops,
IκBα and A20) to a **6-ODE / 13-parameter** non-dimensional model (Eqs 1–6) that stays structurally
and practically **identifiable** while still reproducing the original model's dynamics. The reduced
model's parameters (**Table 1**) are obtained by **refitting the reduced model to the original one**
with PyBioNetFit. Both jobs here reconstruct that refit as edition-2 jobs.

Built + validated with the `validate-pybnf-job` skill from the paper's own **S1 Code** — the authors
ship their reduced and original models as BNGL (`pone.0286416.s015.zip`), so Gate 2 is a **direct
diff against the authors' files** and Gate 1's target is **byte-reproducible from the authors' own
original model**. Each slug is a self-contained folder with its own model, conf, target data,
reproduction figure, README, and `VALIDATION.md`.

## Synthetic-target / ground-truth-recovery jobs

Unlike a fit to *experimental* data, the paper fits the reduced model to reproduce the **original
model's trajectories** (an in-silico experiment). So "the paper's result" is **Table 1** (the fitted
values), the fit **target** is the original Lipniacki-2004 model's output, and **Gate 3b = recover
Table 1** — analogous to `Gupta-2018/fceri` (synthetic data). Gate 3a checks the reduced model *at*
Table 1 reproduces the original (measured by the paper's Average Multiplicative Distance, AMD).

**Objective.** The paper minimizes a sum of squared **geometric-mean-normalized log-differences**
(Eq 7) with a ρ=0.03·max measurement floor (Methods). Both jobs use this **exact** objective, expressed
with standard PyBNF tokens via [`lanl/PyBNF#479`](https://github.com/lanl/PyBNF/pull/479) (ADR-0066):
`noise_model = lognormal, sigma = fix_at 1` (a fixed-σ Gaussian residual on log₁₀ = the paper's
squared-log objective) + `normalization <obs> = floor 0.03, scale`. `floor` = `x' = x + 0.03·max(x)`
(the paper's ρ floor); `scale` = the per-series analytic optimal scale, which for a log-family
objective is `geomean(data)/geomean(sim)` = the paper's ÷geomean exactly. Both are applied
**symmetrically to the simulated and the experimental column** at scoring time, so the `.exp` ship
**raw** original-model output. Edition-2-native (not PEtab-exportable).

**Edition-2 form (both slugs).** TNF (`TR`, 0/1) is delivered as a time-dependent function `TRfunc()`
gated on an inert simulation clock (`Clock()`, rule `0->Clock() 1`, observable `simtime`) — replacing
the authors' `setParameter`/`simulate_ode` phase block — and there is **no `begin actions` block**
(the conf synthesizes each time course). The seed `(1,0,0,0,0,0)` is the exact TR=0 steady state, so
**no pre-equilibration** is needed. This is the Lin-2021/`nyc_multiphase` idiom.

## The jobs

Two slugs fit the **same reduced model to the same original model**, under the paper's two
in-silico experiments. Both recover the same Table 1.

| slug | protocol | fit target | data | status |
|---|---|---|---|---|
| [`reduced_combination`](reduced_combination/) | **combination experiment** (Fig 2): tonic + 5 pulsatile, WT + A20KO, **914 independent points** (S1 Table) | original model, 6 protocols × 2 genotypes | 12-model joint `de` fit, exact Eq-7 objective | Gate 3a **reproduces the paper's reported AMD\*_WT=1.29 / AMD\*_A20KO=1.16 exactly**; Gate 3b recovers **all 13 within 3×** (ε now ~1.5×; Table 1 = the objective minimum) · **94/100** ([VALIDATION](reduced_combination/VALIDATION.md)) |
| [`reduced_onoff`](reduced_onoff/) | **on-off** (Table 2): 2 h TNF on / 10 h off, WT + A20KO, **41 independent points** — the paper's proposed identifiability-optimal protocol | original model, on-off | `de` fit, exact Eq-7 objective | Gate 3a AMD\*≈1.4; Gate 3b recovers **11/13 within 3×** (the 2 that drift = the paper's least-identifiable δ, ε) · **89/100** ([VALIDATION](reduced_onoff/VALIDATION.md)) |

The **combination** slug is the paper's headline result (Fig 2). Because a `condition:` can perturb
only *free* parameters under the bngsim backend (and the TNF window parameters are not fit), each of
the 6 protocols × 2 genotypes is a **generated per-protocol model file** (windows baked in; A20KO =
A20-synthesis rule removed), joined in one multi-model fit sharing the 13 free parameters. The
`reduced_onoff` slug keeps a **single** model (one on-off window), so it is the cleaner showcase of
the approach.

## Source materials (untracked parking garage `dev/papers/Jaruszewicz-Blonska-2023/`)

- **Paper + SI:** `journal.pone.0286416.pdf`; `pone.0286416.s001–s014` (S1–S8 Figs, S1–S4 Tables,
  S1/S2 Text). Table 1 (fitted params) and Table 2 (on-off protocol) are in the main PDF; S1 Table
  (combination protocol time points) is `s009`.
- **Authors' code (S1 Code, `pone.0286416.s015.zip` → `S1_Codes/`):** `BNGL files/Reduced2023 - BNGL
  file/` (the reduced model at Table-1 values, on-off + combination protocols) and `BNGL files/
  Lipniacki_et_al_2004_model - BNGL file/` (the original 15-variable model = the fit target). Also
  MATLAB implementations of 5 models (Hoffmann-2002, Lipniacki-2004, Ashall-2009, Murakawa-2015,
  Reduced-2023).
- **Validation tooling (untracked):** `dev/papers/Jaruszewicz-Blonska-2023/validation/` —
  `gen_onoff.py`, `gen_combination.py`, `gen_combination_models.py`, `obj_at_table1.py`,
  `published_table1.json`.

Not built (optional future slugs): the 5 perturbed-parameter refits (Fig 3, S3/S4 Tables); the
Monte Carlo practical-identifiability ensemble (Fig 7); the "reinterpreted" Krishna-2006 fit to
Lipniacki-2002 (S8 Fig, AMD≈2.06 — a deliberately poor-fit comparison).

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl

cd pybnf-jobs/Jaruszewicz-Blonska-2023/reduced_combination   # the headline Fig-2 fit
pybnf -c reduced_combination.conf
python make_reproduction.py                                   # Fig-2 reproduction at Table-1 params

cd ../reduced_onoff                                           # the identifiability-optimal on-off fit
pybnf -c reduced_onoff.conf
python make_reproduction.py
```
