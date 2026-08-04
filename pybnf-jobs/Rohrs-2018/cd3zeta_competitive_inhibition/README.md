# cd3zeta_competitive_inhibition — CD3ζ site-specific phosphorylation by LCK (PyBNF edition-2 job)

**Run cost: `hours`** — 12,000 evaluations (40 × 300 `de`), each integrating a 672-reaction network over nine time courses.

A PyBNF edition-2 parameter-fitting job that reproduces the **headline parameter-estimation
problem** of Rohrs et al. (2018): fit the competitive-inhibition mechanism of CD3ζ ITAM
phosphorylation to site-resolved phosphoproteomic time courses across wild-type and ITAM-mutant
constructs, derived from:

> Rohrs JA, Zheng D, Graham NA, Wang P, Finley SD. **"Computational model of chimeric antigen
> receptors explains site-specific phosphorylation kinetics."**
> *Biophys J* 2018; **115**(7):1116–1129.
> DOI: [10.1016/j.bpj.2018.08.018](https://doi.org/10.1016/j.bpj.2018.08.018) · PMCID: PMC6199440

Built with the `curate-pybnf-job` skill. **Phase 1** of the Rohrs-2018 curation
(`models/car_cd3zeta_phosphorylation_rohrs2018/`) holds the same mechanism, plus the three
rejected alternatives and the phosphatase prediction, as static reference models at the published
parameters; this job recovers those parameters from the data.

> ℹ️ **Read [`VALIDATION.md`](VALIDATION.md).** The setup is verified end-to-end — tier-1 parse,
> PEtab v2 export → lint → import, and a full `bngsim` fit — and the model **at the authors'
> published parameters** fits the digitized data as well as their own published curve does
> (SSE 1.253e4 against 1.266e4 over the same points). A from-scratch search recovers the **site
> ranking of the Michaelis constants exactly** but lands 18% above the published objective, with
> every parameter inside a factor of 2.7 — the correlated valley the paper itself flags, which the
> authors handled by running PSO 100 times and reporting the 50 best sets as ranges.

## The model

`cd3zeta_competitive_inhibition.bngl` — one CD3ζ molecule with six independently phosphorylatable
ITAM tyrosines (A1, A2, B1, B2, C1, C2) and a constant pool of constitutively active LCK. Every
site is a Michaelis-Menten substrate, and **all twelve site states share one denominator** (Eq. 4 of
the paper): unphosphorylated sites inhibit through K<sub>M,i</sub>, phosphorylated sites through
K<sub>I,i</sub> = K<sub>M,i</sub>·ξ. Finite network: **65 species, 192 reactions**, ODE.

Two things distinguish it from the reference model in `models/`:

1. **No simulation actions** — the nine time courses are synthesized from the conf. The network is
   bounded by the rules themselves, so no `generate_network` directive has to be retained.
2. **`live_<site>` gate parameters.** A tyrosine-to-phenylalanine point mutation is expressed by
   setting a site's gate to 0, which removes it from *both* the substrate and the inhibitor sums —
   exactly what the `~F` state does structurally in
   `models/car_cd3zeta_phosphorylation_rohrs2018_mutants.bngl`. All seven constructs therefore
   share one network and each mutant becomes a PEtab **Condition** (a parameter perturbation)
   rather than a separate model.

## What is fit

Percent phosphorylation at the six ITAM tyrosines, over **nine conditions and 312 scored points**
(Fig. S5, A–J; the six mutated sites are `NaN` in their own condition and are skipped):

| `.exp` | condition | times (min) |
|---|---|---|
| `wt_10pops_rep2.exp` | wild type, 10% POPS liposomes (2nd biological replicate) | 0.1, 1, 5, 10, 30, 60, 180 |
| `wt_0pops.exp` | wild type, 0% POPS | 0.1, 1.5, 10, 60, 360 |
| `wt_45pops.exp` | wild type, 45% POPS | 0.1, 1.5, 10, 60, 360 |
| `mut_A1.exp` … `mut_C2.exp` | the six Y→F ITAM point mutants, 10% POPS | 0.1, 1, 5, 10, 30, 60, 180 |

Column headers are the model **functions** `pct_A1()` … `pct_C2()`, so they carry parentheses (a
manifest entry would list them without). There are no `_SD` columns: the paper scores with the sum
of squared error, so the objective is plain **`sos`** — the same quantity, and inside the
PEtab-exportable subset.

Rohrs et al. also fit a first wild-type replicate (Fig. 1C). That panel is raster, was not
digitized, and is the one of their ten data sets this job does not include.

## Free parameters

Seven, following the paper's protocol exactly (Materials and Methods, *Comparison of model
structures*). k<sub>cat</sub> is held at 360 /min because it is correlated with the LCK density,
and K<sub>M,B1</sub> at 270 molecules/µm² because the six Michaelis constants vary together; both
are the Hui & Vale (2014) values. Ranges are the paper's own: "two orders of magnitude up and down
from their baseline values (LCK = 60 molecules/µm², K<sub>M,i</sub> = 270 molecules/µm², and
X<sub>I</sub> = 1)".

| id | meaning | search range (log10) | published best fit (Data S4) |
|---|---|---|---|
| `LCK_T` | total LCK density (molecules/µm²) | 0.6 – 6000 | 19.703 |
| `KmA1` | K<sub>M</sub> at A1 (molecules/µm²) | 2.7 – 27000 | 96.548 |
| `KmA2` | K<sub>M</sub> at A2 | 2.7 – 27000 | 331.68 |
| `KmB2` | K<sub>M</sub> at B2 | 2.7 – 27000 | 153.82 |
| `KmC1` | K<sub>M</sub> at C1 | 2.7 – 27000 | 484.9 |
| `KmC2` | K<sub>M</sub> at C2 | 2.7 – 27000 | 405.43 |
| `Xi` | inhibition scale ξ, K<sub>I,i</sub> = K<sub>M,i</sub>·ξ | 0.01 – 100 | 0.1153 |

## Files

| file | role |
|---|---|
| `cd3zeta_competitive_inhibition.bngl` | edition-2, fitting-ready ODE model; 65 species, 192 reactions; no actions block |
| `cd3zeta_competitive_inhibition.conf` | the edition-2 `de` job setup (`objective = sos`; 6 conditions, 9 experiments, 7 `*_var` free params) |
| `wt_*.exp`, `mut_*.exp` | fit targets — nine time courses digitized from Fig. S5 |
| `petab/` | **generated** PEtab v2 bundle (`make_petab.py`); do not hand-edit |
| `make_petab.py` | regenerates `petab/` from the conf and lints it |
| `nominal_check.conf` | pins the seven parameters at the published Data S4 values and evaluates the objective once — the reference point every fit is measured against (`sos = 6253.58`, i.e. SSE 1.2507e4; PyBNF's `sos` is half the sum of squared error) |
| `make_reproduction.py` | overlays the model on the data for all nine conditions and prints the SSE; `--params-file output/Results/sorted_params.txt` adds a fitted set |
| `cd3zeta_competitive_inhibition_reproduction.png` | that figure at the published parameters |
| `VALIDATION.md` | primary-source validation and the verification gates |

## Run it

```bash
cd pybnf-jobs/Rohrs-2018/cd3zeta_competitive_inhibition && BNGPATH=$HOME/Simulations/BioNetGen-2.9.3 pybnf -c cd3zeta_competitive_inhibition.conf
```

## `_manifest.py` entry (for `lanl/pybnf` `examples/real-world/`)

```python
RealWorldExample(
    folder='Rohrs-2018/cd3zeta_competitive_inhibition',
    conf='cd3zeta_competitive_inhibition.conf', simulator='ode',
    observables=('pct_A1', 'pct_A2', 'pct_B1', 'pct_B2', 'pct_C1', 'pct_C2'),
    system='LCK phosphorylation of the six CD3zeta ITAM tyrosines of a CD28-CD3zeta CAR on '
           'liposomes, competitive inhibition mechanism (Rohrs 2018, PMC6199440, Eqs. 4-5 and '
           'Data S4); ODE, nine phosphoproteomic time courses from Fig. S5 A-J (wild type on '
           '0/10/45% POPS plus six Y->F ITAM point mutants, the mutants as Conditions through '
           'live_<site> gates); 7 free params, k_cat and K_M,B1 held as in the paper',
    ),
```

No `recover` dict: the seven parameters are correlated (the paper says so explicitly), so a full
`de` + simplex run lands within a factor of 2.7 of the published values but at a different point
along the same valley. A `tol` loose enough to pass would assert nothing; the meaningful
assertions are the PEtab round-trip and a finite objective. See VALIDATION Gate 4b.

## Fit result

One run of the committed conf on a workstation (~80 min): DE stopped on tolerance at generation
163 (`sos` 8748.0), the simplex polish ran its 300 iterations to `sos` **7368.0 = SSE 1.4736e4**,
against **6253.6 = SSE 1.2507e4** at the published parameters. The Michaelis-constant ranking is
reproduced exactly, `KmB1` included:

```
published:  A1  97 < B2 154 < B1 270 < A2 332 < C2 405 < C1 485
fitted:     A1 113 < B2 184 < B1 270 < A2 568 < C2 590 < C1 720
```

`cd3zeta_competitive_inhibition_reproduction.png` overlays both parameter sets on all nine
conditions (solid = published, dashed = fitted, circles = data).
