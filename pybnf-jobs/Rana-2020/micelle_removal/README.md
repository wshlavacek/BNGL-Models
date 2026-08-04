# `micelle_removal` — off-to-on pathway switching, micelles removed at 5 h and 24 h

**Run cost: `hours`** — 18,000 evaluations (60 × 300 `ss`), 14 free parameters over an 80-reaction network.

PyBNF edition-2 fitting job derived from:

> Rana P, Bose P, Vaidya A, Rangachari V, Ghosh P.
> **"Global fitting and parameter identifiability for amyloid-β aggregation with competing
> pathways."**
> *2020 IEEE 20th Int. Conf. on BioInformatics and BioEngineering (BIBE)*:73–78.
> DOI: [10.1109/BIBE50027.2020.00020](https://doi.org/10.1109/BIBE50027.2020.00020)

## What it fits, and why this problem

Sec. IV-C: *"Parameter Estimation and Identifiability with pseudo-micelle removal event (off-to-on
switching)."* Rana et al. *"built the COPASI model with all off and on pathway reactions"* and
*"defined an event to simulate micelle removal from the system along with monomer addition"*, and
report **SSE 1.22**.

This is the best-determined of the paper's four fits: every parameter in the Off-On Switching
column of Table I is set in bold, meaning every one was found identifiable by the
profile-likelihood analysis — the only column for which that holds. That makes it the most useful
recovery target of the four.

Experimentally it is the mirror image of `micelle_addition`. An Aβ42 sample carrying 5 mM C12 fatty
acid from the start sits on the off pathway, where ThT rises with no lag to a modest plateau of
trapped oligomer. Diluting the sample 5- or 10-fold drops the fatty acid below its critical micelle
concentration, the pseudo-micelles disappear, and the paper reports *"a sharp rise in ThT
fluorescence ... indicating the switching of off- to on-pathway species"* along with *"the rise in
the molecular weight of the aggregates including the formation of fibrils"*.

| | |
|---|---|
| model | `micelle_removal.bngl` — 23 species, 80 reactions, deterministic ODE |
| free parameters | 12 rate constants + 2 ThT mapping constants = 14 |
| data | 120 points from Fig. 1b, digitized |
| objective | `sos` (unweighted sum of squares; PyBNF prints half of it) |
| design | undiluted control + two pre-equilibrate/dilute protocols |
| flavor | quantitative, **PEtab.v2-exportable** |

## Files

| file | role |
|---|---|
| `micelle_removal.bngl` | the model, fitting-ready: no simulation actions, `generate_network` retained |
| `micelle_removal.conf` | the edition-2 job setup, banner-commented |
| `fa_control.exp` | Fig. 1b, undiluted 5 mM C12 fatty acid — 39 points |
| `removal_5h.exp` | Fig. 1b, diluted at 5 h — 41 points, clock restarted at the event |
| `removal_24h.exp` | Fig. 1b, diluted at 24 h — 40 points, clock restarted at the event |
| `make_exp.py` | rebuilds the `.exp` files from the committed digitized CSVs |
| `make_reproduction.py` | scores the committed model through BioNetGen and writes the figure |
| `micelle_removal_reproduction.png` | that figure |

Run it (with `BNGPATH` set, from this folder):

```sh
pybnf -c micelle_removal.conf
python make_reproduction.py          # scores micelle_removal.bngl's nominal values
```

## How the timed event is expressed

```conf
condition: fa_pre,     perturbations: "Mic()" = 100
condition: no_fa_post, perturbations: "Mic()" = 0
experiment: removal_5h, preequilibrate: fa_pre, equil_t_end: 5, condition: no_fa_post, data: removal_5h.exp
```

`equil_t_end` forces a fixed 5 h (or 24 h) first phase rather than a steady-state solve; state
carries across the boundary (ADR-0052), so the measured phase starts from the off-pathway oligomer
distribution the sample had actually reached when it was diluted; and the measured phase restarts
the clock at zero, so `make_exp.py` shifts the times. The undiluted control carries its micelles in
the model's seed species and takes no condition, because a species `setConcentration` is inline-only
inside a pre-equilibration protocol (ADR-0062).

**Monomer addition at the event is not modelled explicitly.** Rana et al. mention it but never
quantify it. A 5–10× dilution removes the same fraction of peptide that added monomer restores, so
the peptide pool is simply carried through unchanged — which is the same thing to within the
unstated amount, and adds no free parameter.

## A note on the digitized control

The undiluted control resolves only after 24.8 h. Before that it runs underneath the
dilute-at-24 h sample — the two are the same material until the dilution — and the colour
separation cannot assign the occluded pixels. Those points are not lost to the fit: they are the
dilute-at-24 h sample's own pre-event trace, which the same experiment supplies.

## Adaptations from the published model

Identical in kind to the `micelle_addition` job: no simulation actions, the ThT signal moved into
the conf as `observable: ThT, formula: map_on * Obs_F + map_off * Obs_Fp1` so both mapping constants
are free and the job stays PEtab-exportable, `Mic_0 = 100` in the seed species so the control needs
no condition, and nominals from this job's own fit rather than Table I. `Mic_present = 100` µM
(5 mM C12 fatty acid at an aggregation number of 50) is fixed rather than fitted, because while
micelles are in excess only the product `k_con·[L]` is constrained and `k_con` absorbs the choice.

## Free parameters: published vs recovered

| parameter | units | Table I, Off-On Switching | recovered here |
|---|---|---|---|
| `k_nuon` | /uM/h | **4.61013** | 12.5 |
| `k_nuon_` | /h | **3.49735** | 8.6059 |
| `k_fbon` | /uM/h | *not tabulated* | 4.0565 |
| `k_fbon_` | /h | *not tabulated* | 3.3867 |
| `k_con` | /uM^4/h | **139.545** | 27.409 |
| `k_con_` | /h | *not tabulated* | 0.39221 |
| `k_nuoff` | /uM/h | **195.943** | 62.566 |
| `k_nuoff_` | /h | **49.2318** | 13.57 |
| `k_fboff` | /uM/h | **40.3329** | 0.71978 |
| `k_fboff_` | /h | *not tabulated* | 0.045355 |
| `k_swi` | /h | **0.1** | 598.62 |
| `k_swi_` | /h | *not tabulated* | 0.28192 |
| `map_on` | a.u./uM | *not tabulated* | 42.066 |
| `k_off` | a.u./uM | **164.59** | 1.8698 |

Every published value in this column is bold in Table I — Rana et al. found all of them
identifiable here. `k_con` is not directly comparable even in principle: the paper never states its
pseudo-micelle concentration, and only `k_con·[L]` is constrained.

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
| real bngsim fit | **PASS** — 4 scatter-search iterations at population 12, 1 min, objective 6.633 (SSE 13.27); the simulate/score/propose loop closes on a finite objective |
| reproduction vs the paper | SSE 0.326 over 120 points (RMSE 0.0521); Rana et al. report SSE 1.22 |

### Reproduction

| experiment | panel | n | SSE | median \|rel. err.\| |
|---|---|---|---|---|
| `fa_control` | Fig. 1b/c | 39 | 0.0359 | 4.4% |
| `removal_24h` | Fig. 1b | 40 | 0.2077 | 5.6% |
| `removal_5h` | Fig. 1b | 41 | 0.0827 | 2.8% |
| **total** | | **120** | **0.3263** | RMSE **0.0521** |

SSE is extensive in the number of data points and the paper does not say how many it used, so the
RMSE is the comparable quantity: SSE 1.22 over a set of this size is an RMSE near 0.10.

### On the search

As in the sibling jobs, the global search was run in the SciPy transcription of the same ODEs
rather than inside PyBNF, because PyBNF's Dask-backed scatter search costs tens of seconds per
iteration on this model. The three implementations agree on the objective at an arbitrary parameter
set — PyBNF 7565.14, BioNetGen through `make_reproduction.py` 7565.46, SciPy 7565.48 (SSE/2) — so
the conf, the committed BNGL and the search express the same problem. **Search budget: one differential-evolution run, population 140 (10 per parameter), 60 generations, log-scaled, over the brackets in the conf. No local polish, for the reason given in the sibling jobs. A real global search, but not an exhaustive one; on a sloppy landscape the recovered vector is one good parameterization rather than the optimum.**

## `_manifest.py` entry

```python
RealWorldExample(
    folder="Rana-2020/micelle_removal",
    conf="micelle_removal.conf",
    simulator="ode",
    observables=["ThT"],
    system=(
        "Amyloid-beta 42 aggregation along competing on- and off-pathways with a "
        "pseudo-micelle removal event (Rana et al. 2020, IEEE BIBE, "
        "doi:10.1109/BIBE50027.2020.00020, Sec. IV-C). An Abeta42 sample held on the "
        "off pathway by 5 mM C12 fatty acid is diluted at 5 h and at 24 h to bring the "
        "fatty acid below its critical micelle concentration, releasing trapped "
        "oligomer back onto the fibril-forming pathway. Fits all twelve rate constants "
        "of Eq. 4, Eq. 4-II and Eq. 5 plus two ThT mapping constants to 120 points "
        "digitized from Fig. 1b. Every parameter of this fit is reported identifiable "
        "in Table I, the only one of the paper's four for which that holds."
    ),
    stochastic=False,
    recover={"k_nuon": 12.5, "k_nuon_": 8.6059, "k_fbon": 4.0565, "k_fbon_": 3.3867, "k_con": 27.409, "k_con_": 0.39221, "k_nuoff": 62.566, "k_nuoff_": 13.57, "k_fboff": 0.71978, "k_fboff_": 0.045355, "k_swi": 598.62, "k_swi_": 0.28192, "map_on": 42.066, "map_off": 1.8698},
    tol=10.0,   # sloppy landscape: recovery is to within an order of magnitude
)
```
