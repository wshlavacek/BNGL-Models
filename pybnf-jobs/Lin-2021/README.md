# Lin-2021 — regional COVID-19 epidemic forecasting (PyBNF fitting jobs)

PyBNF edition-2 parameter-fitting jobs derived from the compartmental (SEIR-type) forecasting
model that reproduces daily COVID-19 case counts for US metropolitan statistical areas (MSAs):

> **Lin YT, Neumann J, Miller EF, Posner RG, Mallela A, Safta C, Ray J, Thakur G,
> Chinthavali S, Hlavacek WS.** **"Daily Forecasting of Regional Epidemics of Coronavirus
> Disease with Bayesian Uncertainty Quantification, United States."**
> *Emerg Infect Dis* 2021; **27**(3):767–778.
> DOI: [10.3201/eid2703.203364](https://doi.org/10.3201/eid2703.203364).
> — the compartmental model (social distancing, self-isolation, quarantine) with `n`+1
> social-distancing phases (`n` = 0, 1, 2 selected per region), and the NYC ("BigApple") fit.

**Both slugs here fit the NYC MSA with this Lin-2021 model** — the two-phase (n=1) variant
directly, and the single-phase (n=0) variant. A later paper — Mallela A, Lin YT, Hlavacek WS,
*Epidemics* 2023; 45:100718 ([10.1016/j.epidem.2023.100718](https://doi.org/10.1016/j.epidem.2023.100718))
— **reused this same Lin-2021 model** (the n=0 case) to estimate R₀ across 280 MSAs. That R₀
analysis is not what these jobs implement; Mallela-2023 enters here only as the **source of the
single-phase `nyc/` slug's `DataS1` packaging** — which shipped a placeholder `S0` (see below). So
the folder is named for **Lin-2021**, the paper that actually reports the NYC fit.

Built with the `curate-pybnf-job` skill. Each slug is a **self-contained folder** — its own model,
conf, data, reproduction figure, and README with the exact adaptations, verification results, and
a ready-to-paste `_manifest.py` snippet.

## The model

A deterministic ODE compartmental model: **S → E₁..E₅** (a 5-stage non-exponential incubation
period) **→ A** (asymptomatic) / **I** (symptomatic) **→ H** (hospitalized) **→ R** / **D**, with
testing/contact-tracing **quarantine** (`~Q`) and a social-distancing split **`~M`** (mobile) ↔
**`~P`** (protected) that switches on at day `sigma`. A `counter()` molecule (rule `0->counter() 1`)
integrates simulation time into the observable `t`, which the `if(t>…)` rate-law functions read to
gate the time-dependent social distancing and self-isolation. Fixed epidemiological constants
(incubation/recovery/hospitalization rates, asymptomatic and severe fractions) are Mallela 2023
Table 2. The network is **finite and tiny** (31 species, 62 reactions) — a plain ODE, not heavy.

The fit target is **daily new detected cases** — the day-to-day increment of the running detected-
case count `CumNum_detected_cases_Cum()`. It is scored with a **negative-binomial count
likelihood** (over-dispersed integer counts), which is **not PEtab-exportable**, so the job is
**native-only** (like `Erickson-2019/igf1r`).

## The jobs

Two slugs fit the **same NYC daily case data** with two model complexities. Together they tell a
complete story: the single-phase model needs a correction and doesn't fully reproduce, while the
two-phase model reproduces the published fit exactly — and explains why.

| slug | model | fits | flavor | status |
|---|---|---|---|---|
| [`nyc`](nyc/) | **single-phase** (n=0), Mallela-2023 per-MSA, 6 params | Daily detected cases, first wave (t=50..151) | ODE, **native-only** (neg_bin), `de` fit | ✅ tier-1 + guard + fit (nll 839.9) + reproduction (median **16.3 %**) · **78/100** ([VALIDATION](nyc/VALIDATION.md)) |
| [`nyc_multiphase`](nyc_multiphase/) | **two-phase** (n=1), Lin-2021 forecasting ("BigApple"), 10 params | Daily detected cases, t=0..170 | ODE, **native-only** (neg_bin), `de` fit | ✅ tier-1 + guard + fit + **published-MAP reproduction** (median **15.9 %**, peak within 9 %) · gold-standard **90/100** ([VALIDATION](nyc_multiphase/VALIDATION.md)) |

The authors' published setup is **adaptive-MCMC Bayesian sampling** (`fit_type = am`,
`objfunc = neg_bin_dynamic`); both jobs build a plain **`de` optimization** instead — cheaper and
simpler. `neg_bin` is not PEtab-exportable, so both are **native-only** (verified by `export_job`
raising `NotImplementedError`, like `Erickson-2019/igf1r`).

> ⚠️ **The single-phase `nyc/` slug carries a real provenance finding** (see
> [`nyc/VALIDATION.md`](nyc/VALIDATION.md)). The shipped Mallela-2023 **DataS1** PyBioNetFit `.bngl`
> files carry a **placeholder `S0` (4,903,185 = Alabama's population, identical across every MSA —
> not NYC's ~19.2M) and `sigma`**, so the shipped single-phase `MLE_params.txt` do **not** reproduce
> the data (they overshoot ~4×). `S0` there is **corrected to census** and the nominals are an
> independent `de` fit.
>
> ✅ **The two-phase `nyc_multiphase/` slug is the control that proves this.** Its BigApple files
> ship the **correct census S0** and its **published MAP reproduces the data** (median ~16 %,
> capturing the broad plateau) — a gold-standard reproduction (like `Salazar-Cavazos-2019`). So the
> single-phase discrepancy is fully explained: **a placeholder `S0` plus insufficient model
> complexity** (NYC genuinely needs two social-distancing phases), not a data problem.

## Source materials

- **Primary paper (model, method, NYC fit):** Lin 2021 *Emerg Infect Dis* 27:767–778.
- **`nyc_multiphase` (two-phase) setup:** the authors' own PyBioNetFit files at
  `~/Code/PyBNF/examples/COVID19forecasting_aMCMC/BigApple/` — `m1.bngl` (the two-phase model,
  **correct census S0**), `m1.conf` (`am`/`neg_bin_dynamic` sampler), `m1.exp` (daily NYT case data,
  t = 0..170), `adaptive_files/MLE_params.txt` (the published MAP — **reproduces**).
- **`nyc` (single-phase) setup:** Mallela 2023 *Epidemics* 45:100718 `DataS1/PyBioNetFit/New
  York-Newark-Jersey City, NY-NJ-PA/` — the same Lin-2021 model reduced to `n=0`, but the shipped
  `…​.bngl` has a **placeholder S0/σ**, so its `…_output/adaptive_files/MLE_params.txt` does **not**
  reproduce (see the caveat above). Byte-identical copy at
  `~/Code/PyBNF/examples/Mallela2022MSAs/New York-Newark-Jersey City, NY-NJ-PA/`.
- **Case data:** *The New York Times* COVID-19 county-level dataset (github.com/nytimes/covid-19-data),
  aggregated to the NYC MSA (model day 0 = 2020-01-21).

Not built (optional future slugs): other MSAs (Detroit, Miami, Appleton, …); the `n=2` three-phase
variant; the full adaptive-MCMC posterior (`fit_type = am`) capturing the published uncertainty.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl

cd pybnf-jobs/Lin-2021/nyc_multiphase   # the gold-standard two-phase fit (published MAP reproduces)
pybnf -c nyc_multiphase.conf            # de fit (raise population_size/max_iterations for a tighter optimum)
python make_reproduction.py             # reproduction figure (model at published MAP vs. NYT daily cases)

cd ../nyc                               # the single-phase fit (corrected S0; independent de fit)
pybnf -c nyc.conf
python make_reproduction.py
```
