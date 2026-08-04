# `micelle_addition` — on-to-off pathway switching, fatty acid added at 3 h and 24 h

**Run cost: `hours`** — 18,000 evaluations (60 × 300 `ss`), 14 free parameters over an 80-reaction network.

PyBNF edition-2 fitting job derived from:

> Rana P, Bose P, Vaidya A, Rangachari V, Ghosh P.
> **"Global fitting and parameter identifiability for amyloid-β aggregation with competing
> pathways."**
> *2020 IEEE 20th Int. Conf. on BioInformatics and BioEngineering (BIBE)*:73–78.
> DOI: [10.1109/BIBE50027.2020.00020](https://doi.org/10.1109/BIBE50027.2020.00020)

## What it fits, and why this problem

Sec. IV-B of the paper: *"Parameter Estimation and Identifiability with the pseudo-micelle
addition Event (on-to-off pathway switching)."* Rana et al. build a model *"with events at 3h and
24 h, which simulate the addition of pseudo-micelles and monomers"*, fit it to three experiments,
and report **SSE 4.12**. They find every parameter of this fit identifiable except `k_swi` and
`k_con`.

The experiment is the interesting part. An Aβ42 sample in the on-pathway lag phase is given 5 mM
C12 fatty acid, and the ThT signal *"increased ... without a lag phase"* — the monomer that was
slowly filling the nucleation chain is diverted into off-pathway oligomer instead. Doing that at
two different times probes the same switch from two different starting states, which is what makes
the off-pathway constants identifiable from a signal that only reports aggregate.

| | |
|---|---|
| model | `micelle_addition.bngl` — 23 species, 80 reactions, deterministic ODE |
| free parameters | 12 rate constants + 2 ThT mapping constants = 14 |
| data | 322 points from Fig. 1c, digitized |
| objective | `sos` (unweighted sum of squares; PyBNF prints half of it) |
| design | fatty-acid reference + two pre-equilibrate/switch protocols |
| flavor | quantitative, **PEtab.v2-exportable** |

## Files

| file | role |
|---|---|
| `micelle_addition.bngl` | the model, fitting-ready: no simulation actions, `generate_network` retained |
| `micelle_addition.conf` | the edition-2 job setup, banner-commented |
| `fa_control.exp` | Fig. 1c, 5 mM C12 fatty acid throughout — 73 points |
| `addition_3h.exp` | Fig. 1c, micelles added at 3 h — 147 points, clock restarted at the event |
| `addition_24h.exp` | Fig. 1c, micelles added at 24 h — 102 points, clock restarted at the event |
| `make_exp.py` | rebuilds the `.exp` files from the committed digitized CSVs |
| `make_reproduction.py` | scores the committed model through BioNetGen and writes the figure |
| `micelle_addition_reproduction.png` | that figure |

Run it (with `BNGPATH` set, from this folder):

```sh
pybnf -c micelle_addition.conf
python make_reproduction.py          # scores micelle_addition.bngl's nominal values
```

## How the timed event is expressed

There is no "event" directive in an edition-2 conf, and none is needed: an addition experiment is a
fixed-length **pre-equilibration** phase under one condition followed by the measured phase under
the other.

```conf
condition: no_fa_pre, perturbations: "Mic()" = 0
condition: fa_add,    perturbations: "Mic()" = 100
experiment: addition_3h, preequilibrate: no_fa_pre, equil_t_end: 3, condition: fa_add, data: addition_3h.exp
```

Three things make this work, and each is a trap if missed:

1. **`equil_t_end` forces a fixed interval** instead of the ODE default of a steady-state solve.
   A lag-phase sample has no steady state on this timescale, so without it the phase is meaningless.
2. **State carries across the phase boundary** (ADR-0052), so the measured phase starts from the
   oligomer distribution the sample had actually reached at 3 h — which is the whole point of doing
   the experiment at two different times.
3. **The measured phase restarts the clock at zero**, so `make_exp.py` shifts the data by the event
   time and drops the pre-event points. Nothing is lost: those points *are* the on-pathway
   experiment, which the `on_pathway` and `global` jobs fit in their own right.

Two further constraints shaped the conditions. A species `setConcentration` is inline-only inside a
pre-equilibration protocol (ADR-0062), so the no-event control cannot carry `"Mic()" = 100` as a
plain `condition:` — it holds its micelles in the model's seed species instead and takes no
condition. And one condition cannot serve both a pre-equilibrated and a standalone experiment
(ADR-0052), which is why the 0 µM state appears as `no_fa_pre` rather than being shared.

## Adaptations from the published model

`micelle_addition.bngl` is a fitting-ready copy of
[`models/amyloid_beta_competing_aggregation_pathways_rana2020`](../../../models/amyloid_beta_competing_aggregation_pathways_rana2020),
which transcribes Eq. 4, Eq. 4-II and Eq. 5 with the rate laws taken from the paper's own flux
equations. Four differences:

1. **No simulation actions.** Only `generate_network({overwrite=>1})` is kept. The network is
   finite without a cap — the chains terminate at A₁₂ and A′₁₂ — so there is no `max_stoich` to
   retain.
2. **The ThT observable is conf-side.** The curated model carries BNGL functions `ThT_on()`,
   `ThT_off()` and `ThT_total()`; here the signal is `observable: ThT, formula: map_on * Obs_F +
   map_off * Obs_Fp1`, which makes both mapping constants free parameters and keeps the job inside
   the PEtab-exportable subset. `normalization = ...` would not be.
3. **`Mic_0 = 100` in the seed species,** so the fatty-acid control needs no condition (see above).
4. **Nominals are this job's best fit,** not Table I.

`Mic_present = 100` µM is 5 mM C12 fatty acid at an aggregation number of 50. It is **fixed, not
fitted**: while micelles are in excess only the product `k_con·[L]` is constrained by these data,
so the aggregation number and `k_con` are perfectly correlated and `k_con` absorbs the choice.
100 µM is a large excess over the ~6 µM of micelle the off pathway can consume from 25 µM peptide.

## Free parameters: published vs recovered

| parameter | units | Table I, On-Off Switching | recovered here |
|---|---|---|---|
| `k_nuon` | /uM/h | **8.83326** | 13.815 |
| `k_nuon_` | /h | **253.824** | 15.798 |
| `k_fbon` | /uM/h | *not tabulated* | 59.708 |
| `k_fbon_` | /h | *not tabulated* | 4.3497 |
| `k_con` | /uM^4/h | 675.528 *(reported unidentifiable)* | 0.78518 |
| `k_con_` | /h | *not tabulated* | 205.32 |
| `k_nuoff` | /uM/h | **28.4215** | 4427.6 |
| `k_nuoff_` | /h | **1.41607** | 3070.5 |
| `k_fboff` | /uM/h | **426.503** | 0.47638 |
| `k_fboff_` | /h | *not tabulated* | 13.192 |
| `k_swi` | /h | 183.97 *(reported unidentifiable)* | 4887 |
| `k_swi_` | /h | *not tabulated* | 0.083996 |
| `map_on` | a.u./uM | *not tabulated* | 1.3744 |
| `k_off` | a.u./uM | **97.024** | 65.448 |

Bold marks the parameters Rana et al. report as identifiable for this fit. `k_con` is not directly
comparable across the two columns even in principle: the paper never states its pseudo-micelle
concentration, and only `k_con·[L]` is constrained, so the value depends on a convention it does
not give.

The recovered constants differ from the published column by up to three orders of magnitude
even where Rana et al. mark them identifiable, and the fit is nonetheless as good as theirs.
Three things account for that and none of them is reconcilable from the paper alone: the data
here are digitized rather than the authors' own tables; the pseudo-micelle concentration is a
convention this reconstruction had to choose, which moves `k_con` and everything that trades
against it; and five of the thirteen constants (`k_fbon`, `k_fbon_`, `k_con_`, `k_fboff_`,
`k_swi_`) are absent from Table I altogether, so the published vector is not a point in the
same 14-dimensional space this job searches. What the comparison does establish is that the
published *model* reproduces the published *data* to the published accuracy.

## Verification

| gate | result |
|---|---|
| tier-1 (`check_conf.py`) | **PASS** — edition 2, `job_type=ss` resolves, 3 experiments bound, 14 free params bound by id, no `__FREE` |
| PEtab.v2 round-trip (`petab_roundtrip.py`) | **PASS** — export → `petab.v2` lint clean → import |
| real bngsim fit | **PASS** — 4 scatter-search iterations at population 12, 1 min, objective 9.556 (SSE 19.11); the simulate/score/propose loop closes on a finite objective |
| reproduction vs the paper | SSE 4.267 over 322 points (RMSE 0.1151); Rana et al. report SSE 4.12 |

### Reproduction

`make_reproduction.py` re-runs `micelle_addition.bngl` through BioNetGen under the same two-phase
protocol the conf synthesizes, and scores it against the same 322 points:

| experiment | panel | n | SSE | median \|rel. err.\| |
|---|---|---|---|---|
| `addition_24h` | Fig. 1c | 102 | 2.0570 | 9.7% |
| `addition_3h` | Fig. 1c | 147 | 1.5936 | 10.6% |
| `fa_control` | Fig. 1b/c | 73 | 0.6166 | 16.8% |
| **total** | | **322** | **4.2672** | RMSE **0.1151** |

SSE is extensive in the number of data points and Rana et al. do not say how many they used, so the
RMSE is the comparable quantity: their SSE 4.12 over a three-experiment set of this size is an RMSE
near 0.11.

### On the search

The nominals were found by a differential-evolution search run in the SciPy transcription of the
same ODEs (`models/…/independent_rana2020.py`), not inside PyBNF: PyBNF's Dask-backed scatter
search costs tens of seconds per iteration on this 23-species model, which is not a budget under
which a 14-parameter landscape converges. Three implementations were then checked to agree on the
objective at one arbitrary parameter set — PyBNF 7565.14, BioNetGen through
`make_reproduction.py` 7565.46, SciPy 7565.48 (SSE/2) — so the conf, the committed BNGL and the
search all express the same problem. **Search budget: one differential-evolution run, population 140 (10 per parameter), 60 generations, log-scaled, over the brackets in the conf. No local polish: on this model a Nelder-Mead or Powell step costs seconds per evaluation in the fitted region and did not converge inside a usable budget. This is a real global search but not an exhaustive one, and the landscape is sloppy, so the recovered vector should be read as one good parameterization rather than as the optimum.**

## `_manifest.py` entry

```python
RealWorldExample(
    folder="Rana-2020/micelle_addition",
    conf="micelle_addition.conf",
    simulator="ode",
    observables=["ThT"],
    system=(
        "Amyloid-beta 42 aggregation along competing on- and off-pathways with a "
        "pseudo-micelle addition event (Rana et al. 2020, IEEE BIBE, "
        "doi:10.1109/BIBE50027.2020.00020, Sec. IV-B). Fatty-acid pseudo-micelles are "
        "added at 3 h and at 24 h to an Abeta42 sample still in its on-pathway lag "
        "phase, diverting monomer from fibril formation into a kinetically trapped "
        "off-pathway oligomer. Fits all twelve rate constants of Eq. 4, Eq. 4-II and "
        "Eq. 5 plus two ThT mapping constants to 322 points digitized from Fig. 1c. "
        "The timed event is expressed as preequilibrate + equil_t_end + condition, "
        "which is the edition-2 idiom for a mid-protocol stimulus."
    ),
    stochastic=False,
    recover={"k_nuon": 13.815, "k_nuon_": 15.798, "k_fbon": 59.708, "k_fbon_": 4.3497, "k_con": 0.78518, "k_con_": 205.32, "k_nuoff": 4427.6, "k_nuoff_": 3070.5, "k_fboff": 0.47638, "k_fboff_": 13.192, "k_swi": 4887, "k_swi_": 0.083996, "map_on": 1.3744, "map_off": 65.448},
    tol=10.0,   # sloppy landscape: recovery is to within an order of magnitude
)
```
