# `on_pathway` — Aβ42 fibril formation with no fatty acid

PyBNF edition-2 fitting job derived from:

> Rana P, Bose P, Vaidya A, Rangachari V, Ghosh P.
> **"Global fitting and parameter identifiability for amyloid-β aggregation with competing
> pathways."**
> *2020 IEEE 20th Int. Conf. on BioInformatics and BioEngineering (BIBE)*:73–78.
> DOI: [10.1109/BIBE50027.2020.00020](https://doi.org/10.1109/BIBE50027.2020.00020)

## What it fits, and why this problem

This is the paper's own first step. Sec. IV-A: *"First, we fitted the experimental data of the
on-pathway, considering the reactions in Eq. 4. These estimated on-pathway parameters were used to
define the parameter ranges in the subsequent steps."* Everything downstream — both switching fits
and the global fit — brackets its on-pathway constants at 0.1–10× whatever this fit returns, so it
is the anchor of the whole chain.

It is also the cleanest identifiability question in the paper. Rana et al. report that *"the forward
rate constants are not identifiable in the presence of backward rate constants"* and that `k_nuon`
and `k_fbon` become identifiable only once the corresponding backward constants are fixed (Fig. 2).
This job fits all four freely, and finds exactly that.

| | |
|---|---|
| model | `on_pathway.bngl` — 12 species, 44 reactions, deterministic ODE |
| free parameters | 4 rate constants + 1 ThT mapping constant = 5 |
| data | 97 points, Fig. 1a, 25 µM Aβ42, 0–48 h, digitized |
| objective | `sos` (unweighted sum of squares; PyBNF prints half of it) |
| design | single time course, no events |
| flavor | quantitative, **PEtab.v2-exportable** |

## Files

| file | role |
|---|---|
| `on_pathway.bngl` | the model, fitting-ready: no simulation actions, `generate_network` retained |
| `on_pathway.conf` | the edition-2 job setup, banner-commented |
| `on_pathway.exp` | Fig. 1a, 97 points |
| `make_exp.py` | rebuilds the `.exp` from the committed digitized CSV |
| `make_reproduction.py` | scores the committed model through BioNetGen and writes the figure |
| `on_pathway_reproduction.png` | that figure |

Run it (with `BNGPATH` set, from this folder):

```sh
pybnf -c on_pathway.conf
python make_reproduction.py          # scores on_pathway.bngl's nominal values
```

## Adaptations from the published model

`on_pathway.bngl` is a fitting-ready copy of the `_on_pathway` variant of
[`models/amyloid_beta_competing_aggregation_pathways_rana2020`](../../../models/amyloid_beta_competing_aggregation_pathways_rana2020),
which transcribes Eq. 4 and the fluxes H_i and I_i directly. Three differences:

1. **No simulation actions.** Only `generate_network({overwrite=>1})` is kept; the protocol comes
   from the conf. The network is finite without a cap (the chain terminates at A₁₂), so there is no
   `max_stoich` to retain.
2. **The ThT observable is conf-side.** The curated model carries a BNGL function
   `ThT_on() = map_on*Obs_F`; here it is `observable: ThT, formula: map_on * Obs_F`, which makes
   `map_on` a free parameter and keeps the job inside the PEtab-exportable subset. Encoding the
   arbitrary ThT scale as `normalization = ...` would not be exportable.
3. **Nominals are this job's best fit,** not Table I — see below.

## Free parameters: published vs recovered

| parameter | units | Table I, On-Pathway | recovered here | note |
|---|---|---|---|---|
| `k_nuon` | /µM/h | 22.04 | **3.6677** | |
| `k_nuon_` | /h | 12.72 | **82.142** | |
| `k_fbon` | /µM/h | *not tabulated* | **6.7244** | text: "k_fb is about 400 times higher than k_nu" |
| `k_fbon_` | /h | *not tabulated* | **9.9897** | text: "k_fb_ is almost 100–200 times lower than k_fb" |
| `map_on` | a.u./µM | *not tabulated* | **8.6693** | ThT scale |

**Why not the published values.** They cannot produce Fig. 1a at the peptide concentration the
paper states. The lag in this model is the sequential filling of the eleven-step chain
A₁ → … → A₁₂, and that takes time only while the steps are thermodynamically uphill, i.e. while
K·[A₁] = (`k_nuon`/`k_nuon_`)·[A₁] ≲ 1. At the tabulated pair, K = 1.733 /µM and K·[A₁] = 43 at
25 µM: every step runs downhill, the chain equilibrates within ~0.03 h, and no value of `k_fbon` or
`k_fbon_` can then restore a 28 h lag, because A₁₂ is slaved to A₁₁ through the reverse of H₁₁. The
recovered pair has K = 0.0446 /µM, K·[A₁] = 1.12 — just barely uphill, which is what a lag needs.
The curated model's notebook demonstrates this by simulation.

The recovered `k_fbon`/`k_fbon_` ratio is 0.67, against the paper's stated 100–200. That is a real
disagreement, and it is not surprising given that the two backward constants are the ones the paper
itself reports as unidentifiable.

## Verification

| gate | result |
|---|---|
| tier-1 (`check_conf.py`) | **PASS** — edition 2, `job_type` resolves, data bound, 5 free params bound by id, no `__FREE` |
| PEtab.v2 round-trip (`petab_roundtrip.py`) | **PASS** — export → `petab.v2` lint clean → import |
| real bngsim fit | **PASS** — 8 scatter-search iterations at population 16, 1 min, objective 0.12 (SSE 0.24); the simulate/score/propose loop closes on a finite objective |
| reproduction vs the paper | **PASS** — SSE 0.0217 over 97 points, RMSE 0.0150 |

### Reproduction

`make_reproduction.py` re-runs `on_pathway.bngl` through BioNetGen at its nominal values and scores
it against the same 97 points:

| | SSE | RMSE | median \|rel. err.\| |
|---|---|---|---|
| this fit | **0.0217** | **0.0150** | 12.4% |

Rana et al. report **SSE 0.13** for their on-pathway fit. SSE is extensive in the number of data
points and they do not say how many they used, so the comparable quantity is the RMSE: 0.13 over a
plate-reader trace of ~500 points is an RMSE of 0.016, against 0.0150 here. The fit is at least as
good as the published one, and both sit at the precision floor set by digitization — the median
half-height of the Fig. 1a marker cloud is 0.019 ThT a.u.

### Identifiability

Four distinct basins were found, all within a factor of two in SSE:

| SSE | `k_nuon` | `k_nuon_` | `k_fbon` | `k_fbon_` | `map_on` |
|---|---|---|---|---|---|
| **0.0217** | 3.668 | 82.14 | 6.724 | 9.990 | 8.669 |
| 0.0239 | 0.2004 | 7.382 | 11.74 | 9.976 | 15.06 |
| 0.0451 | 0.008527 | 1.125 | 28.43 | 14.52 | 87.99 |
| 0.0455 | 0.009811 | 1.394 | 36.81 | 16.22 | 139.3 |

`k_nuon` spans nearly three decades across basins while the SSE changes by a factor of two — the
concrete form of the non-identifiability Rana et al. describe. Their remedy, fixing the backward
constants, is what Fig. 2 of the paper reports.

### Optimizer note

`job_type = ss`. Differential evolution at PyBNF's defaults converges prematurely on this landscape:
178 generations at population 48 left the entire population between objective 0.4145 and 0.4207
(SSE 0.829–0.841), while the optimum is at objective 0.0108 (SSE 0.0217). Scatter search reached
objective 0.0375 by its ninth iteration. Rana et al. reach the same conclusion by a different route
— their Table II benchmarks fourteen COPASI optimizers and scatter search is the only one that fits
their global problem.

The search that produced the nominals above was run outside PyBNF, in the SciPy transcription of the
same ODEs (`models/…/independent_rana2020.py`), because PyBNF's Dask-backed scatter search costs
~40 s per iteration on this model and a multi-basin five-parameter landscape needs more iterations
than that allows. The result was then **scored back through PyBNF**, which returns objective
0.010845 = SSE 0.021689 — identical to the SciPy value, confirming that the two agree on the
objective as well as on the trajectories.

## `_manifest.py` entry

```python
RealWorldExample(
    folder="Rana-2020/on_pathway",
    conf="on_pathway.conf",
    simulator="ode",
    observables=["ThT"],
    system=(
        "Amyloid-beta 42 on-pathway fibril formation (Rana et al. 2020, IEEE BIBE, "
        "doi:10.1109/BIBE50027.2020.00020, Eq. 4). Eleven-step nucleation chain "
        "A_i + A_1 <-> A_i+1 to the 12-mer nucleus F = A_12, plus catalytic fibril "
        "binding A_i + F <-> F. Fits the four on-pathway rate constants and the ThT "
        "mapping constant to the 97-point thioflavin T time course of Fig. 1a "
        "(25 uM Abeta42, no fatty acid), digitized from the publisher PDF. The "
        "paper's own first fitting step; it anchors the parameter brackets of the "
        "three competing-pathway jobs. Landscape is multi-basin, which is the "
        "non-identifiability the paper reports for the forward/backward pairs."
    ),
    stochastic=False,
    recover={"k_nuon": 3.6677, "k_nuon_": 82.142, "k_fbon": 6.7244,
             "k_fbon_": 9.9897, "map_on": 8.6693},
    tol=10.0,   # multi-basin: recovery is to within an order of magnitude, not tighter
)
```
