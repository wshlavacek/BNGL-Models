# nyc_multiphase — New York City MSA COVID-19, TWO-PHASE model, ODE (PyBNF edition-2 job)

A PyBNF edition-2 parameter-fitting job that fits **daily new detected COVID-19 case counts** for
the New York City MSA (NY-NJ-PA) over the first wave, using the **two-phase (n=1)** social-
distancing model, derived from:

> Lin YT, Neumann J, Miller EF, Posner RG, Mallela A, Safta C, Ray J, Thakur G, Chinthavali S,
> Hlavacek WS. **"Daily Forecasting of Regional Epidemics of Coronavirus Disease with Bayesian
> Uncertainty Quantification, United States."** *Emerg Infect Dis* 2021; **27**(3):767–778.
> DOI: [10.3201/eid2703.203364](https://doi.org/10.3201/eid2703.203364) (NYC = "BigApple").

Built with the `curate-pybnf-job` skill. The authors' published setup is **adaptive-MCMC Bayesian
sampling** (`fit_type = am`); this job builds a plain **`de` optimization** instead. It ports the
authors' own PyBioNetFit setup `examples/COVID19forecasting_aMCMC/BigApple/` (`m1.bngl` / `m1.conf` /
`m1.exp`) to the edition-2 surface.

> ✅ **Gold-standard reproduction.** UNLIKE the sibling single-phase [`nyc/`](../nyc/) slug — whose
> shipped DataS1 file had a **placeholder `S0`** so the published MLE did *not* reproduce — this
> model ships the **correct NYC census S0 (19,216,182)** and its **published MAP reproduces the
> data** (median ~16 % rel err, peak within 9 %). The **two-phase** social distancing captures the
> broad NYC plateau (days ~68–85) that a single phase cannot. This slug is *why* we know the
> `nyc/` discrepancy was a placeholder-`S0` + model-complexity issue, not a data problem.

## The model

Same SEIR-type compartmental structure as the single-phase slug (S → E₁..E₅ → A/I → H → R/D, with
`~Q` quarantine and `~M`/`~P` social distancing), but with **two** social-distancing phases:

- **phase 1** (setpoint `p0`, timescale `1/lambda0`) switches on at `sigma = t0 + t_delta`;
- **phase 2** (`p1`, `lambda1`) at `tau1 = t0 + t_delta + t_delta2`.

A `counter()` molecule (rule `0->counter() 1`, **seed `counter()=1`**) integrates simulation time
into the observable `t` (so `t = 1 + sim_time` — the authors' convention, ported verbatim), read by
the `if(t≥…)` functions to gate the phase changes and self-isolation. Fixed constants are the
region-independent Lin-2021 values. Network is finite (**~35 species, 61 reactions**), ODE, **not
heavy**.

## What is fit

Daily new detected cases over the first wave (one time course, `nyc_multiphase.exp`, t = 0..170):

| observable | design | source |
|---|---|---|
| `fDCs_Cum` (differenced → daily incidence) | time course, model day t = 0..170 (2020-01-21..07-09) | `nyc_multiphase.exp` (authors' MSA-aggregated NYT case counts) |

`fDCs_Cum` is a **Molecules** observable — the true cumulative detected-case count (integral of the
detection rate `fD·(1-fA)·kL·nE5`). The conf's per-observable **`cumulative`** flag differences it
into daily incidence; counts are over-dispersed integers → a **negative-binomial** count likelihood
with the dispersion estimated as `r_disp` and mean centering.

## Free parameters (9 model + 1 nuisance) — published MAP (which reproduces)

Nominals in the `.bngl` are the authors' published MAP (`BigApple/adaptive_files/MLE_params.txt`);
the conf brackets each on the authors' original uniform bounds.

| id | published MAP | bounds | role |
|---|---|---|---|
| `t0` | 32.832 | 0–60 | epidemic start day |
| `t_delta` | 0.073 | 0–60 | delay to phase-1 onset (σ = t0+t_delta = 32.9 ≈ 2020-02-23) |
| `t_delta2` | 110.223 | 0–200 | delay to phase-2 (τ₁ = 143.1 ≈ 2020-06-12) |
| `beta` | 2.033 | 0–10 | transmission rate constant (/d) |
| `lambda0` | 0.098 | 0–10 | phase-1 social-distancing eigenvalue (/d) |
| `p0` | 0.872 | 0–1 | phase-1 setpoint (fraction protected) |
| `lambda1` | 5.181 | 0–10 | phase-2 social-distancing eigenvalue (/d) |
| `p1` | 0.736 | 0–1 | phase-2 setpoint (fraction protected) |
| `fD` | 0.115 | 0–1 | fraction of mild symptomatic cases detected |
| `r_disp` | 10.27 | 1–30 (log) | negative-binomial dispersion (likelihood only) |

## 🚫 Native-only (NOT PEtab-exportable) · not heavy

`neg_bin` and mean-centering are outside the PEtab v2 subset, so this job is native-only (verified by
`export_job` raising `NotImplementedError`, not by a round-trip). The ~35-species network builds in
~0.5 s.

## Verification (see `VALIDATION.md`)

- **Tier-1** (`scripts/check_conf.py`): edition 2, `job_type=de` resolves, data bound, 10 free params
  bind by id, no `__FREE`. **PASS.**
- **Native-only guard** (`scripts/petab_roundtrip.py`): `export_job` raises `NotImplementedError`.
  **PASS.**
- **Real bngsim fit:** converges to a finite objective (neg_bin nll 1018.99); the published MAP
  scores **1003.5** — *better* than the de run, confirming the published params are near-optimal.
- **Reproduction at the PUBLISHED MAP** (`nyc_multiphase_reproduction.png`): **median 15.9 % rel err**
  over t=50..151 (same window as the single-phase slug), peak 14,435/day @ day 77 vs. data 15,925 @
  day 73 (within 9 %), cumulative 498,528 vs. 487,332 (+2 %). Tracks the data across ~5 orders of
  magnitude and captures the broad plateau. **This is the authors' own reported fit** — gold-standard.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Lin-2021/nyc_multiphase
pybnf -c nyc_multiphase.conf     # raise population_size/max_iterations for a tighter optimum
python make_reproduction.py      # reproduction figure (model at the published MAP vs. data)
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='nyc_multiphase', conf='nyc_multiphase.conf', simulator='ode',
    observables=('fDCs_Cum',),
    system='COVID-19 regional epidemic, New York City MSA, TWO-PHASE (n=1) social distancing '
           '(Lin 2021, DOI 10.3201/eid2703.203364; NYC "BigApple"); ODE SEIR-type compartmental, '
           'daily detected-case counts, neg_bin count likelihood with cumulative->incident '
           'differencing; native-only. de fit; the published MAP reproduces the data (median '
           '~16% rel err) -- gold-standard, correct census S0. See VALIDATION.md.'),
# Native-only: assert export_job raises NotImplementedError (not a PEtab round-trip).
# The published MAP reproduces, so recover={} could target it; the 10-param fit is sloppy, so the
# reproduction (not parameter recovery) is the validation. See README for the table.
```
