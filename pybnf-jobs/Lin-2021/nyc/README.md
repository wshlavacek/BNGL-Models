# nyc — New York-Newark-Jersey City MSA COVID-19 first wave, ODE (PyBNF edition-2 job)

A PyBNF edition-2 parameter-fitting job that fits **daily new detected COVID-19 case counts**
for the New York-Newark-Jersey City MSA (NY-NJ-PA) over the first pandemic wave, derived from:

> Mallela A, Lin YT, Hlavacek WS. **"Differential contagiousness of respiratory disease across
> the United States."** *Epidemics* 2023; **45**:100718.
> DOI: [10.1016/j.epidem.2023.100718](https://doi.org/10.1016/j.epidem.2023.100718)
> · model lineage: Lin YT, et al. *Emerg Infect Dis* 2021; **27**(3):767–778
> ([10.3201/eid2703.203364](https://doi.org/10.3201/eid2703.203364)).

Built with the `curate-pybnf-job` skill. The authors' published setup is **adaptive-MCMC Bayesian
sampling** (`fit_type = am`, `objfunc = neg_bin_dynamic`); this job builds a plain **`de`
optimization** instead — cheaper and simpler. It is the **single-phase (n=0), 6-parameter** model
from the Mallela-2023 supplement `DataS1/`.

> ⚠️ **Not a gold-standard reproduction — read [`VALIDATION.md`](VALIDATION.md).** The DataS1
> shipped model files carry a **placeholder `S0`** (`4,903,185` = Alabama's population, identical
> across every MSA — not NYC's ~19.2M) **and `sigma`**, so the shipped `MLE_params.txt` does
> **not** reproduce the data (it overshoots ~4×). Here `S0` is **corrected to NYC's 2019 census
> MSA population** and the nominals / reproduction use an **independent PyBNF `de` fit**, not the
> non-reproducing published MLE.

## The model

A deterministic ODE compartmental (SEIR-type) model: **S → E₁..E₅** (5-stage non-exponential
incubation) **→ A** (asymptomatic) / **I** (symptomatic) **→ H** (hospitalized) **→ R** / **D**,
with testing/contact-tracing **quarantine** (`~Q`) and a social-distancing split **`~M`** (mobile)
↔ **`~P`** (protected) that switches on at day `sigma = 63` (≈ 2020-03-24, NY's stay-at-home
order). A `counter()` molecule (rule `0->counter() 1`) integrates simulation time into the
observable `t`, which the `if(t>…)` rate-law functions read to gate the time-dependent social
distancing and self-isolation. Fixed constants are Mallela 2023 Table 2. The network is finite and
tiny — **31 species, 62 reactions**, ODE, **not heavy**, no network cap.

## What is fit

Daily new detected cases over the first wave (one time course, `nyc.exp`):

| observable | design | source |
|---|---|---|
| `CumNum_detected_cases_Cum()` (differenced → daily incidence) | time course, model day t = 50..151 (2020-03-11..06-20; day 0 = 2020-01-21) | `nyc.exp` (authors' MSA-aggregated NYT case counts) |

`CumNum_detected_cases_Cum()` is a model **function** (parens in the `.exp` header; no parens in
`_manifest.py`). It is a running (cumulative-like) detected-case count; the conf's per-observable
**`cumulative`** flag differences it row-to-row into daily incidence before scoring. Counts are
over-dispersed integers → a **negative-binomial** count likelihood with the dispersion estimated
as a free nuisance `r_disp` and mean centering (mandatory — the ODE gives the mean count).

## Files

| file | role |
|---|---|
| `nyc.bngl` | edition-2, fitting-ready ODE model; 5 free params (best-fit nominals), `r` deleted (→ `r_disp`), corrected census `S0`, no actions block |
| `nyc.conf` | the edition-2 `de` job setup (neg_bin count likelihood + per-observable `cumulative`) |
| `nyc.exp` | fit target — the authors' daily NYT case data, t = 50..151 |
| `make_reproduction.py` | regenerates the reproduction PNG (model at the `de` best fit vs. data) |
| `nyc_reproduction.png` | verification figure (linear + log) |
| `VALIDATION.md` | provenance / model-fidelity / reproduction record + the placeholder-`S0` finding |

## Free parameters (5 model + 1 nuisance) — independent `de` fit vs. published MLE

Nominals in the `.bngl` are an independent PyBNF `de`+refine fit (population 30 × 60 iters,
neg_bin nll **839.88**); the conf brackets each on the authors' original uniform bounds.

| id | `de` fit (nominal) | published MLE | bounds | role |
|---|---|---|---|---|
| `t0` | 26.633 | 27.367 | 0–60 | epidemic start day (**recovers the published value** once S0 is corrected) |
| `beta` | 0.538 | 0.671 | 0–10 | transmission rate constant (/d) |
| `lambda0` | 9.160 | 0.105 | 0–10 | social-distancing relaxation eigenvalue (/d) |
| `p0` | 0.592 | 0.752 | 0–1 | social-distancing setpoint (fraction protected) |
| `fD` | 0.883 | 0.262 | 0–1 | fraction of mild symptomatic cases detected |
| `r_disp` | 9.887 | 11.21 | 1–30 (log) | negative-binomial dispersion (likelihood only, not a model rate) |

The published MLE does **not** reproduce the data with the shipped placeholder `S0`/`sigma` (nll
6.1M, 4.4× peak overshoot); see `VALIDATION.md` Gate 3b. With `S0` corrected, the `de` fit recovers
`t0` near the published value; `lambda0`/`fD` differ (the fit compensates with faster social
distancing + higher detection to match the wave with a single social-distancing phase).

## 🚫 Native-only (NOT PEtab-exportable) · not heavy

`neg_bin` and mean-centering were removed from PEtab v2, so this job cannot round-trip through
PEtab (like `Erickson-2019/igf1r`). Verified by `export_job` raising `NotImplementedError` — **not**
by a PEtab round-trip. The 31-species network builds in ~0.5 s, so a full fit runs on a workstation.

## Verification (see `VALIDATION.md`)

- **Tier-1** (`scripts/check_conf.py`): edition 2, `job_type=de` resolves, `nyc.exp` bound, 6 free
  params bind by id, no `__FREE`. **PASS.**
- **Native-only guard** (`scripts/petab_roundtrip.py`): `export_job` raises `NotImplementedError`.
  **PASS.**
- **Real bngsim fit:** converges to a finite objective (neg_bin nll 839.88); `counter()` == sim
  time and the day-63 social-distancing switch fires under bngsim.
- **Reproduction** (`nyc_reproduction.png`, model at the `de` fit): **median relative error 16.3 %**
  over the 102 fit points; tracks the data across ~5 orders of magnitude on the log scale. The
  single-phase model fits a sharper/earlier peak than the data's broad plateau (a structural n=0
  limit, not a setup error).

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Lin-2021/nyc
pybnf -c nyc.conf                # raise population_size/max_iterations for a tighter optimum
python make_reproduction.py      # reproduction figure (model vs. NYT daily-case data)
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='nyc', conf='nyc.conf', simulator='ode',
    observables=('CumNum_detected_cases_Cum',),
    system='COVID-19 regional epidemic, New York-Newark-Jersey City MSA (Mallela 2023, '
           'DOI 10.1016/j.epidem.2023.100718; model lineage Lin 2021); ODE SEIR-type '
           'compartmental, daily detected-case counts, neg_bin count likelihood with '
           'cumulative->incident differencing; native-only (neg_bin not PEtab-exportable). '
           'de fit; the shipped DataS1 MLE does not reproduce (placeholder S0/sigma) -- S0 '
           'corrected to census, nominals are an independent de fit. See VALIDATION.md.'),
# Native-only: assert export_job raises NotImplementedError (do NOT assert a PEtab round-trip).
# recover left empty: the shipped MLE does not reproduce, so there is no published point estimate
# to recover to; the de fit reaches nll 839.88 / 16.3% median rel err. See README for the table.
```
