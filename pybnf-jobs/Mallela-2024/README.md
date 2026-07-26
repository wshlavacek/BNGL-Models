# Mallela-2024 — COVID-19 with vaccination and the Alpha/Delta variants (PyBNF fitting jobs)

PyBNF edition-2 parameter-fitting jobs derived from the compartmental model that reproduces daily
COVID-19 case counts in four US metropolitan statistical areas (MSAs) through the vaccination
rollout and the Alpha and Delta variant waves:

> **Mallela A, Chen Y, Lin YT, Miller EF, Neumann J, He Z, Nelson KE, Posner RG, Hlavacek WS.**
> **"Impacts of Vaccination and Severe Acute Respiratory Syndrome Coronavirus 2 Variants Alpha and
> Delta on Coronavirus Disease 2019 Transmission Dynamics in Four Metropolitan Areas of the United
> States."** *Bull Math Biol* 2024; **86**:31.
> DOI: [10.1007/s11538-024-01258-4](https://doi.org/10.1007/s11538-024-01258-4).
> — the model (Appendix Eqs. 1–40), the fixed constants (Table 3), and the per-MSA maximum a
> posteriori (MAP) estimates (Table 1) reproduced in Figs. 2–5.

Built with the `curate-pybnf-job` skill. Each slug is a **self-contained folder** — its own model,
conf, data, `reference/` exports (SBML · SED-ML · COMBINE), reproduction figure, README, and
VALIDATION with the gate-by-gate evidence.

The **library-model** siblings (the same four models with an actions block, for reference
simulation rather than fitting) live in
[`models/covid19_vaccination_and_variants_mallela2024/`](../../models/covid19_vaccination_and_variants_mallela2024/).
The **predecessor** model this one extends — Lin et al. 2021, without vaccination or variants — is
fit in [`pybnf-jobs/Lin-2021/`](../Lin-2021/).

## The model

A deterministic ODE compartmental model. The Lin-2021 core — **S → E₁..E₅** (5-stage Erlang
incubation) **→ A** (asymptomatic) / **I** (mildly symptomatic) **→ H** (severe) **→ R** / **D**,
each population carrying a behavioral state **`~M`** (mixing) ↔ **`~P`** (social-distancing) or
**`~Q`** (quarantined) — extended with:

- **vaccination:** susceptible and recovered persons move at the empirical per capita rate μ(t)
  into a six-stage Erlang immune response V₁..V₆; persons leaving V₆ partition into
  S<sub>V,1..4</sub> by how many of the three circulating strains their response neutralizes;
  breakthrough infections run E<sub>V</sub> → A<sub>V</sub>/I<sub>V</sub> → H<sub>V</sub> → R/D
  with a reduced risk m<sub>h</sub> of severe disease;
- **multi-period social distancing:** two coordinated step functions `Ptau()` (setpoint) and
  `Lambdatau()` (approach rate) switch at σ, τ₁..τ<sub>n</sub>;
- **variants:** `Ytheta()` multiplies the transmission rate constant β by y₁ at the Alpha takeoff
  θ₁ and by y₂ at the Delta takeoff θ₂, while `Utheta1()`/`Utheta2()` open S<sub>V</sub>
  reinfection to each new strain.

A `counter()` molecule (`0->counter() 1`, seed 0) integrates simulation time into the observable
`t`, which every `if(t≥…)` rate law reads — so the models are **ODE-only**. Each network is finite
and small (**42 species, 88 reactions**, ~0.2 s to build); the regions differ in S₀, the empirical
vaccination series, and the number of social-distancing periods.

The fit target is **daily new detected cases** — the day-to-day increment of `fDCs_Cum`, the
cumulative detected-case count f<sub>D</sub>·C<sub>S</sub>(t). It is scored with a
**negative-binomial** count likelihood, mean-centered, with the dispersion **fixed at the authors'
published r**. `neg_bin` + `location = mean` + `cumulative` are all outside the PEtab v2 subset, so
every slug here is **native-only** (verified by `export_job` raising `NotImplementedError`).

## The jobs

One slug per MSA, each a **gradient (L-BFGS-B) refinement of the published MAP in two parameters**
— β (transmission rate constant) and f<sub>D</sub> (detected fraction) — with every other
parameter, including **every switch time**, pinned at the paper's value.

| slug | MSA | S₀ | periods | fig | published MAP → fit | status |
|---|---|---|---|---|---|---|
| [`nyc`](nyc/) | New York-Newark-Jersey City, NY-NJ-PA | 19.2 M | 4 (n=3) | 4A | NLL **5313.95 → 5235.80** (−78) | ✅ tier-1 + guard + fit + reproduction · **85/100** ([VALIDATION](nyc/VALIDATION.md)) |
| [`dallas`](dallas/) | Dallas-Fort Worth-Arlington, TX | 7.6 M | 5 (n=4) | 2A | NLL **4875.96 → 4762.00** (−114) | ✅ tier-1 + guard + fit + reproduction · **85/100** ([VALIDATION](dallas/VALIDATION.md)) |
| [`houston`](houston/) | Houston-The Woodlands-Sugar Land, TX | 7.1 M | 4 (n=3) | 3A | NLL **4776.28 → 4650.97** (−125) | ✅ tier-1 + guard + fit + reproduction · **85/100** ([VALIDATION](houston/VALIDATION.md)) |
| [`phoenix`](phoenix/) | Phoenix-Mesa-Scottsdale, AZ | 4.9 M | 5 (n=4) | 5A | NLL **6324.09 → 4523.99** (−1800) | ✅ tier-1 + guard + fit + reproduction · **85/100** ([VALIDATION](phoenix/VALIDATION.md)) |

**Every slug's 2-parameter fit lowers the paper's own objective**, most dramatically for Phoenix
(−28 %, and the median relative error against the 7-day rolling mean halves from 50 % to 26 %).
That is expected rather than a correction to the paper: the published values are the *mode of an
adaptive-MCMC posterior* over 19–22 parameters, not an optimum, and this job optimizes two of them
directly.

| slug | med \|rel err\| vs 7-day mean | peak of 7-day mean vs data | cumulative cases vs data |
|---|---|---|---|
| `nyc` | 23.0 % → 22.4 % | +31.6 % → +43.1 % | +14.7 % → +28.8 % |
| `dallas` | 28.6 % → 26.6 % | +23.6 % → −1.8 % | +18.0 % → +4.0 % |
| `houston` | 34.8 % → 25.7 % | +50.0 % → −3.3 % | +26.1 % → +8.1 % |
| `phoenix` | 50.4 % → 26.3 % | +49.8 % → −17.2 % | +13.3 % → +8.4 % |

(NYC is the one region where the fit trades a taller winter-2020 peak for a better fit through the
long low-count Delta tail — the negative-binomial likelihood weights the many low days above the
few peak days. Its NLL still improves.)

## Two things these jobs pin down

> ⚠️ **A gradient job starts at the BOX CENTER of the priors, not at the `.bngl` nominals**
> (`gradient_base.py`, `local_base.py` `_resolve_start_pset`). The authors' verbatim bounds
> (β 0–20) would start every search at β = 10 — ~30× the MAP — where CVODES stalls and the
> objective comes back `inf`. Each conf here narrows the box to the published MAP ± 50 %, so the
> box center *is* the paper's estimate. This is the single most important choice in the confs.

> ⚠️ **Gradient-fitting this class of model needed two upstream PyBNF fixes**
> ([lanl/PyBNF#522](https://github.com/lanl/PyBNF/pull/522)), each of which made the objective
> non-finite everywhere before it was fixed: (1) PyBNF asked bngsim for an output sensitivity of
> *every* global function, and bngsim refuses any `if()`-bearing body — all 14 of this model's
> functions are `if()` chains, so every simulation died; (2) the negative-binomial mean-slope
> `r(mean−obs)/(mean(r+mean))` is `0/0 → nan` at a prediction of exactly 0, which is what an
> epidemic model predicts before its start time `t0`. Both are documented per slug in
> `VALIDATION.md` under *Toolchain requirement*.

**No switch time is free in any slug.** `t0`, `t_sigma`, `t_tau*`, `theta1`, `t_theta2` appear only
inside `if(t≥…)` conditions, so their in-branch ∂f/∂p is identically zero — a gradient optimizer
cannot move them without the crossing-jump sensitivity of bngsim issue #48 / PR #50. Wiring that
through PyBNF and freeing the switch times is the natural next step for these jobs; the chained
parameterization here (σ = t₀+t_σ, τᵢ = τᵢ₋₁+t_τᵢ, θ₂ = θ₁+t_θ₂, with 5–7 switch times per region)
makes this paper the scaled-up companion to that PR's Lin-2021 test case.

## Source materials

- **Primary paper:** Mallela 2024 *Bull Math Biol* 86:31 (`dev/papers/Mallela2024/BullMathBiol24.pdf`).
- **Authors' own PyBioNetFit setup:** `~/Code/PyBNF/examples/Vax_and_Variants/{NYC,Dallas,Houston,Phoenix}/`
  — `<MSA>.bngl`, `<MSA>.conf` (the legacy `fit_type = am` / `objfunc = neg_bin_dynamic` adaptive-MCMC
  sampler), `<MSA>.exp` (daily NYT case counts, t = 0..648), and
  `Output/adaptive_files/MLE_params.txt` (the published MAP == Table 1).
- **Case data:** *The New York Times* COVID-19 county-level dataset
  (github.com/nytimes/covid-19-data), aggregated to each MSA (model day 0 = 2020-01-21). Copied
  verbatim from the authors' `.exp` files; a few days carry negative counts (a state revising its
  cumulative total downward), which PyBNF's count-domain guard scores as 0.
- **Predecessor model:** Lin YT et al., *Emerg Infect Dis* 2021; 27(3):767–778,
  DOI [10.3201/eid2703.203364](https://doi.org/10.3201/eid2703.203364) — fit in
  [`pybnf-jobs/Lin-2021/`](../Lin-2021/).

Not built (optional future slugs): the full adaptive-MCMC posterior (`fit_type = am`) capturing the
published uncertainty; a switch-time-free fit once bngsim PR #50 is wired through; joint multi-MSA
fits sharing the universal constants.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl

cd pybnf-jobs/Mallela-2024/phoenix   # the biggest improvement over the published MAP
pybnf -c phoenix.conf                # L-BFGS-B gradient fit of beta and fD (a few minutes)
python make_reproduction.py          # figure + metrics: published MAP vs. this fit vs. NYT data
```

Each slug runs the same way (`nyc`, `dallas`, `houston`, `phoenix`). `make_reproduction.py`
simulates through BNG2.pl, so the figure reproduces without the gradient toolchain.
