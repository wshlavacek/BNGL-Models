# `global` — one parameterization for all five ThT experiments

**Run cost: `hours`** — 28,800 evaluations (72 × 400 `ss`) fitting 15 parameters across all protocols at once.

PyBNF edition-2 fitting job derived from:

> Rana P, Bose P, Vaidya A, Rangachari V, Ghosh P.
> **"Global fitting and parameter identifiability for amyloid-β aggregation with competing
> pathways."**
> *2020 IEEE 20th Int. Conf. on BioInformatics and BioEngineering (BIBE)*:73–78.
> DOI: [10.1109/BIBE50027.2020.00020](https://doi.org/10.1109/BIBE50027.2020.00020)

## What it fits, and why this problem

The paper's headline result, Sec. IV-D: *"We also performed a global fit of all the five curves: on
pathway data, micelle addition at 3h and 24 hour, micelle removal at 5h and 24 hour... The achieved
fit had a combined SSE of 8.2."*

The reason it matters is visible in Table I. The two individual switching fits, of the *same* model
to the *same* system, disagree by up to three orders of magnitude on the same constants —
`k_nuoff_` is 1.42 /h in one column and 49.2 /h in the other, `k_swi` is 184 /h against 0.1 /h.
Either the mechanism is wrong or the individual fits are underdetermined, and a global fit is the
experiment that distinguishes those: one parameter vector now has to satisfy five protocols that
start the system from different states and perturb it in opposite directions.

| | |
|---|---|
| model | `global.bngl` — 23 species, 80 reactions, deterministic ODE |
| free parameters | 12 rate constants + 3 ThT mapping constants = 15 |
| data | 427 points from Figs. 1a, 1b and 1c, digitized |
| objective | `sos` (unweighted sum of squares; PyBNF prints half of it) |
| design | one plain time course + four pre-equilibrate/switch protocols |
| flavor | quantitative, **PEtab.v2-exportable** |

## Files

| file | role |
|---|---|
| `global.bngl` | the model, fitting-ready: no simulation actions, `generate_network` retained |
| `global.conf` | the edition-2 job setup, banner-commented |
| `on_pathway.exp` | Fig. 1a, no fatty acid — 97 points |
| `addition_3h.exp` / `addition_24h.exp` | Fig. 1c — 147 and 102 points |
| `removal_5h.exp` / `removal_24h.exp` | Fig. 1b — 41 and 40 points |
| `make_exp.py` | rebuilds the `.exp` files from the committed digitized CSVs |
| `make_reproduction.py` | scores the committed model through BioNetGen and writes the figure |
| `global_reproduction.png` | that figure |

Run it (with `BNGPATH` set, from this folder):

```sh
pybnf -c global.conf
python make_reproduction.py          # scores global.bngl's nominal values
```

## Two ThT mapping constants, not one

Table I gives the global fit **`k_off1` and `k_off2`** where each individual fit has a single
`k_off`. That is not an accident of notation: the addition series and the removal series are
separate ThT experiments with their own arbitrary normalization, and one scale cannot serve both —
Fig. 1b plateaus near 1.5 a.u. while Fig. 1c plateaus near 1.0. The conf mirrors it with two
measurement models:

```conf
observable: ThT_add, formula: map_on * Obs_F + map_off1 * Obs_Fp1
observable: ThT_rem, formula: map_on * Obs_F + map_off2 * Obs_Fp1
```

`map_on` is shared, because the on-pathway fibril signal is the same signal in every panel. Writing
these as measurement models rather than `normalization = ...` is also what keeps the job
PEtab-exportable.

## The five experiments, and only those

`make_exp.py` binds exactly the five curves the paper names. The two fatty-acid controls are *not*
among them and are left to the individual jobs. The on-pathway experiment needs no event and takes
no condition — the micelle pool is empty in the model's seed species — but it still runs the
complete competing-pathway network, so the off pathway remains reachable through the switching
reaction of Eq. 5, exactly as in the paper's own global model.

The four switching experiments use the `preequilibrate` + `equil_t_end` + `condition` idiom
described in the sibling job READMEs. Because a condition cannot serve both a pre-equilibrated and
a standalone experiment (ADR-0052), and a species `setConcentration` is inline-only inside a
pre-equilibration protocol (ADR-0062), the two micelle states appear under four names, one per role.

## Search brackets

Sec. IV-D: *"To avoid overfitting, we considered the minimum reaction set of on and off-pathway. We
chose the range of parameters carefully using the values estimated from the previous fit."* The
conf does that literally — each rate constant is bracketed by the span of the values the
`micelle_addition` and `micelle_removal` jobs recovered, widened one decade either side, and the
two off-pathway mapping constants are free over the paper's own 10⁰–10⁵ window.

## Free parameters: published vs recovered

| parameter | units | Table I, Global | recovered here |
|---|---|---|---|
| `k_nuon` | /uM/h | **0.141644** | 1.3125 |
| `k_nuon_` | /h | 37.8865 *(reported unidentifiable)* | 157.19 |
| `k_fbon` | /uM/h | *not tabulated* | 337.39 |
| `k_fbon_` | /h | *not tabulated* | 43.182 |
| `k_con` | /uM^4/h | 1221.69 *(reported unidentifiable)* | 39.491 |
| `k_con_` | /h | *not tabulated* | 0.58468 |
| `k_nuoff` | /uM/h | **40.8795** | 322.01 |
| `k_nuoff_` | /h | **8.69222** | 43.61 |
| `k_fboff` | /uM/h | **97.3176** | 0.19102 |
| `k_fboff_` | /h | *not tabulated* | 0.015403 |
| `k_swi` | /h | **2.93977** | 48816 |
| `k_swi_` | /h | *not tabulated* | 1.1412 |
| `map_on` | a.u./uM | *not tabulated* | 0.64786 |
| `k_off1` | a.u./uM | **161.832** | 1.9046 |
| `k_off2` | a.u./uM | **48.8039** | 2.0943 |

Bold marks the parameters Rana et al. report as identifiable for the global fit; they find all of
them so except `k_con` and `k_nuon_`. `k_con` is not directly comparable across the two columns even
in principle, because the paper never states its pseudo-micelle concentration and only `k_con·[L]`
is constrained.

The recovered constants differ from the published column by orders of magnitude even where Rana
et al. mark them identifiable, and the fit is within a factor of two of theirs. Three things
account for that and none is reconcilable from the paper alone: the data here are digitized rather
than the authors' own tables; the pseudo-micelle concentration is a convention this reconstruction
had to choose, which moves `k_con` and everything that trades against it; and five of the thirteen
constants are absent from Table I altogether, so the published vector is not a point in the same
15-dimensional space this job searches. What the comparison establishes is that the published
*model* reproduces the published *data* to roughly the published accuracy — not that the published
*numbers* are recoverable, which the paper's own identifiability analysis already says they are
not for `k_con` and `k_nuon_`.

## Verification

| gate | result |
|---|---|
| tier-1 (`check_conf.py`) | **PASS** — edition 2, `job_type=ss` resolves, 5 experiments bound, 15 free params bound by id, no `__FREE` |
| PEtab.v2 round-trip (`petab_roundtrip.py`) | **PASS** — export → `petab.v2` lint clean → import |
| real bngsim fit | **PASS** — 3 scatter-search iterations at population 12, 1 min, objective 33.75 (SSE 67.49); the simulate/score/propose loop closes on a finite objective |
| reproduction vs the paper | SSE 15.653 over 427 points (RMSE 0.1915); Rana et al. report SSE 8.2 |

### Reproduction

| experiment | panel | n | SSE | median \|rel. err.\| |
|---|---|---|---|---|
| `addition_24h` | Fig. 1c | 102 | 6.0163 | 24.0% |
| `addition_3h` | Fig. 1c | 147 | 2.6758 | 8.3% |
| `on_pathway` | Fig. 1a | 97 | 5.9761 | 228.4% |
| `removal_24h` | Fig. 1b | 40 | 0.8186 | 7.3% |
| `removal_5h` | Fig. 1b | 41 | 0.1662 | 5.1% |
| **total** | | **427** | **15.6529** | RMSE **0.1915** |

Rana et al.'s SSE 8.2 over a five-curve set of this size is an RMSE near 0.14 — the global fit is
visibly the loosest of their four, and Fig. 1d shows it: several of their model curves sit well
away from their data. SSE is extensive in the number of points and the paper does not say how many
it used, so the RMSE is the comparable quantity.

**Where the compromise is paid is worth reading off that table.** The two removal curves score
0.17 and 0.82 — better than the dedicated `micelle_removal` job manages on the same data — and
`addition_3h` scores 2.68. The cost falls on the fatty acid-free curve (5.98) and on
`addition_24h` (6.02), which together are three quarters of the total. Reproducing the ~28 h lag
of Fig. 1a needs a nucleation chain that is barely uphill, and the four fatty-acid protocols pull
the same constants somewhere else. That tension is the whole point of the global fit, and the
paper reports it too: its Global column puts `k_nuon` at 0.14 /µM/h against 22.04 in the
On-Pathway column, a 150-fold move in the same constant between the two fits.

### On the optimizer

`job_type = ss`. This is the one methodological choice the paper makes for us: Table II benchmarks
fourteen COPASI optimizers on this exact problem and reports that for the global fit *"only one
algorithm (scatter search) was able to fit the data"* — 8.99 against 26–301 for every other method.
Differential evolution at PyBNF's defaults behaves the same way here as their DE did there; on the
simpler `on_pathway` problem it left an entire population of 48 between objective 0.4145 and 0.4207
after 178 generations while the optimum sits at 0.0108.

The global search behind the nominals was nonetheless run in the SciPy transcription of the same
ODEs rather than inside PyBNF, because PyBNF's Dask-backed scatter search costs tens of seconds per
iteration on this 23-species model. The three implementations were checked to agree on the
objective at an arbitrary parameter set — PyBNF 7565.14, BioNetGen through `make_reproduction.py`
7565.46, SciPy 7565.48 (SSE/2) — so the conf, the committed BNGL and the search all express the
same problem. **Search budget: one differential-evolution run, population 150 (10 per parameter), 70 generations, log-scaled, over the brackets in the conf. No local polish: a Nelder-Mead or Powell step costs seconds per evaluation on this fifteen-parameter, five-protocol problem and did not converge inside a usable budget. A real global search, but not an exhaustive one, and the paper's own analysis says the landscape has unidentifiable directions in it (k_con and k_nuon_), so the recovered vector is one good parameterization rather than the optimum.**

## `_manifest.py` entry

```python
RealWorldExample(
    folder="Rana-2020/global",
    conf="global.conf",
    simulator="ode",
    observables=["ThT_add", "ThT_rem"],
    system=(
        "Amyloid-beta 42 aggregation along competing on- and off-pathways, fit "
        "globally to all five published thioflavin T experiments at once (Rana et al. "
        "2020, IEEE BIBE, doi:10.1109/BIBE50027.2020.00020, Sec. IV-D): the "
        "fatty acid-free time course, pseudo-micelle addition at 3 h and 24 h, and "
        "pseudo-micelle removal at 5 h and 24 h. Twelve rate constants of Eq. 4, "
        "Eq. 4-II and Eq. 5 plus three ThT mapping constants against 427 digitized "
        "points. The two off-pathway mapping constants mirror Table I's k_off1 and "
        "k_off2, one per measurement series. Four of the five experiments are timed "
        "events expressed as preequilibrate + equil_t_end + condition, so the job "
        "exercises that idiom in both directions from one model."
    ),
    stochastic=False,
    recover={"k_nuon": 1.3125, "k_nuon_": 157.19, "k_fbon": 337.39, "k_fbon_": 43.182, "k_con": 39.491, "k_con_": 0.58468, "k_nuoff": 322.01, "k_nuoff_": 43.61, "k_fboff": 0.19102, "k_fboff_": 0.015403, "k_swi": 48816, "k_swi_": 1.1412, "map_on": 0.64786, "map_off1": 1.9046, "map_off2": 2.0943},
    tol=10.0,   # sloppy landscape: recovery is to within an order of magnitude
)
```
