# dallas — Dallas MSA COVID-19 with vaccination and variants, ODE (PyBNF edition-2 job)

A PyBNF edition-2 parameter-fitting job that fits **daily new detected COVID-19 case counts** for
the Dallas MSA (Dallas-Fort Worth-Arlington, TX) over 2020-01-21..2021-10-30, using the vaccination-and-variants
compartmental model of:

> Mallela A, Chen Y, Lin YT, Miller EF, Neumann J, He Z, Nelson KE, Posner RG, Hlavacek WS.
> **"Impacts of Vaccination and Severe Acute Respiratory Syndrome Coronavirus 2 Variants Alpha
> and Delta on Coronavirus Disease 2019 Transmission Dynamics in Four Metropolitan Areas of the
> United States."** *Bull Math Biol* 2024; **86**:31.
> DOI: [10.1007/s11538-024-01258-4](https://doi.org/10.1007/s11538-024-01258-4)
> — model = Appendix Eqs. 1–40; fit target = panel A of **Fig. 2**; published MAP = **Table 1**.

Built with the `curate-pybnf-job` skill. The authors' published setup is **adaptive-MCMC Bayesian
sampling** (`fit_type = am`, `objfunc = neg_bin_dynamic`) over all 22 model parameters plus the
dispersion; this job builds a **gradient (L-BFGS-B) refinement of just two of them** instead. It
ports the authors' own PyBioNetFit setup `~/Code/PyBNF/examples/Vax_and_Variants/Dallas/`
(`Dallas.bngl` / `Dallas.conf` / `Dallas.exp`) to the edition-2 surface.

> ✅ **The 2-parameter fit improves on the published MAP: neg_bin NLL 4875.96 → 4762.00
> (Δ = -113.96, 2.3 % lower).** This is a clean improvement on every summary metric at once. The published MAP is an MCMC-sampled *mode*, not an
> optimum, so this is expected — and it is what the job demonstrates.

## The model

The Lin et al. (2021) SEIR-type compartmental model (S → E₁..E₅ → A/I → H → R/D with `~M`/`~P`
social distancing and `~Q` quarantine) extended with **vaccination** and the **Alpha** and
**Delta** variants:

- susceptible and recovered persons are vaccinated at the empirical per capita rate μ(t) into a
  six-stage Erlang immune response V₁..V₆; persons leaving V₆ partition into S<sub>V,1..4</sub> by
  how many circulating strains their response neutralizes; breakthrough infections run
  E<sub>V</sub> → A<sub>V</sub>/I<sub>V</sub> → H<sub>V</sub> → R/D;
- social distancing is **5 periods** (n = 4): two coordinated step functions `Ptau()`
  (setpoint) and `Lambdatau()` (approach rate) switch at σ, τ₁..τ<sub>4</sub>;
- `Ytheta()` multiplies the transmission rate constant β by y₁ at the Alpha takeoff θ₁ and by y₂ at
  the Delta takeoff θ₂.

A `counter()` molecule (rule `0->counter() 1`, seed `counter()=0`) integrates simulation time into
the observable `t`, which every `if(t≥…)` rate law reads — so the model is **ODE-only**. S₀ =
**7,573,136**. Network is finite (**42 species, 88 reactions**), builds in ~0.2 s, **not heavy**.

`dallas.bngl` is the curated library model
[`models/covid19_vaccination_and_variants_mallela2024`](../../../models/covid19_vaccination_and_variants_mallela2024/)
minus its *simulation* actions (pybnf synthesizes the run from the conf); nothing else differs, and
its parameter nominals are the published MAP.

## What is fit

Daily new detected cases over the paper's full window (one time course, `dallas.exp`, t = 0..648):

| observable | design | source |
|---|---|---|
| `fDCs_Cum` (differenced → daily incidence) | time course, model day t = 0..648 (2020-01-21..2021-10-30) | `dallas.exp` (authors' MSA-aggregated NYT case counts) |

`fDCs_Cum` is a **Molecules** observable — the cumulative detected-case count f<sub>D</sub>·C<sub>S</sub>(t)
(Eqs. 39–40). The conf's per-observable **`cumulative`** flag differences it into daily incidence;
counts are over-dispersed integers → a **negative-binomial** count likelihood, **mean**-centered
(the ODE gives the mean), with the dispersion **fixed at the authors' published r = 3.02166**
rather than fit.

## Free parameters (2) — and why only 2

| id | published MAP | box | L-BFGS-B fit | role |
|---|---|---|---|---|
| `beta` | 0.304173 | 0.15–0.46 | **0.310366** | transmission rate constant (/d) |
| `fD` | 0.481101 | 0.24–0.72 | **0.400519** | fraction of new symptomatic infections detected |

Everything else stays pinned at the published MAP. Two reasons:

1. **No switch time can move under a gradient.** `t0`, `t_sigma`, `t_tau1..t_tau4`, `theta1`,
   `t_theta2` appear only inside `if(t≥…)` conditions, so their in-branch ∂f/∂p is identically
   **zero**. A gradient optimizer can move them only if the simulator supplies the finite jump each
   crossing contributes — bngsim issue #48 / PR #50. Freeing them here would add columns of zeros.
2. **Scope.** The setpoints/eigenvalues (p0..p4, lambda0..lambda4) and variant transmissibilities
   (y₁, y₂) are pinned to keep the job small, fast, and identifiable. β and f<sub>D</sub> are the
   two the daily-case series constrains most directly — an overall transmission scale and an
   overall reporting fraction — and both are smooth everywhere.

**The boxes are narrow on purpose.** A gradient job's start 0 is the **box center** of the priors,
not the `.bngl` nominals (`gradient_base.py` `_n_starts_key = 'population_size'`; `local_base.py`
`_resolve_start_pset`). The authors' verbatim bounds (β 0–20) would start the search at β = 10,
where CVODES stalls and the objective is `inf`. Each box is the published MAP ±50 %, rounded to two
decimals, so the box center *is* the MAP to within 0.5 %.

## 🚫 Native-only (NOT PEtab-exportable) · not heavy

`neg_bin`, mean-centering, and the cumulative→incident differencing are all outside the PEtab v2
subset, so this job is native-only — verified by `export_job` raising `NotImplementedError`, not by
a round-trip.

## Verification (see [`VALIDATION.md`](VALIDATION.md))

- **Tier-1** (`scripts/check_conf.py`): edition 2, `job_type=lbfgs` resolves, data bound, 2 free
  params bind by id, no `__FREE`. **PASS.**
- **Native-only guard:** `export_job` raises `NotImplementedError` (mean centering). **PASS.**
- **Real bngsim gradient fit:** `pybnf -c dallas.conf` converges to a **finite** objective
  **4762.00** (3 L-BFGS-B starts × ≤20 iterations, a few minutes).
- **Reproduction** (`dallas_reproduction.png`, metric = the job's own objective): NLL
  **4875.96 → 4762.00**; median |rel err| vs. the 7-day mean **28.6 % → 26.6 %**;
  peak of the 7-day mean **+23.6 % → -1.8 %** of the data's; cumulative cases
  **+18.0 % → +4.0 %**.

## `reference/` — SBML · SED-ML · COMBINE

`reference/` ships the generated network and verified-faithful exports of this job-form model:
`dallas.net`, `dallas.xml` (SBML L3V2), `dallas.sedml` (SED-ML L1V3 uniform time course over the
fit window, CVODE, atol/rtol 1e-7), and `dallas.omex` (COMBINE archive bundling all of it plus the
source `.bngl` and a conversion report). Produced by `bngsim.convert.net_to_omex(..., gate="full")`;
the **L0–L3 faithfulness ladder passes** (L3: source and conversion agree on all 27,258
species·time cells, max rel |Δ| = 0). L4 (symbolic) is inconclusive by construction — the
67-branch empirical vaccination table `v_rate()` is a piecewise the symbolic checker punts on.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Mallela-2024/dallas
pybnf -c dallas.conf          # L-BFGS-B gradient fit of beta and fD
python make_reproduction.py   # reproduction figure + metrics (published MAP vs. this fit)
```

Requires **PyBNF with the piecewise-gradient fixes** (`fix: gradient-fit a piecewise model with a
count likelihood`) — see [`VALIDATION.md`](VALIDATION.md#toolchain-requirement).

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='dallas', conf='dallas.conf', simulator='ode',
    observables=('fDCs_Cum',),
    system='COVID-19 transmission with vaccination and the Alpha/Delta variants, Dallas MSA '
           '(Mallela 2024, DOI 10.1007/s11538-024-01258-4, Fig. 2A); ODE SEIR-type '
           'compartmental with a six-stage Erlang vaccine response, 5 social-distancing '
           'periods and two variant step changes; daily detected-case counts, neg_bin count '
           'likelihood (dispersion fixed at the published r) with cumulative->incident '
           'differencing; native-only. GRADIENT (lbfgs) refinement of beta and fD from the '
           'published MAP, which it improves (nll 4875.96 -> 4762.00). See VALIDATION.md.'),
# Native-only: assert export_job raises NotImplementedError (not a PEtab round-trip).
# recover={'beta': 0.304173, 'fD': 0.481101} would FAIL by design -- the fit deliberately
# beats the published MAP; the NLL improvement, not parameter recovery, is the validation.
```
