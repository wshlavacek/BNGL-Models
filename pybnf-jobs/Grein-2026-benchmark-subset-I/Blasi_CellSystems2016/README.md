# Blasi_CellSystems2016

**Run cost: `minutes`** — 10,000 evaluations (20 × 500 `gntr`) on a 32-reaction ODE model.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**SOLVED** — a from-scratch `gntr` fit reaches `OG = -4.3e-07 < 1.92`, and it recovers the
published maximum-likelihood parameters of Blasi et al. (2016) to four digits.

## The problem

Per the SBML `<notes>`: *"PEtab implementation of the model from Blasi et al. (2016), Cell
Systems; Volume 2, Issue 1, P49-58"* — combinatorial acetylation of the histone H4 N-terminal
tail in *Drosophila melanogaster*. Lysines K5, K8, K12 and K16 give 2⁴ = 16 acetylation
states, which are the model's 16 species (`x_0ac`, four singles, six doubles, four triples,
`x_4ac`); the 32 reactions are the edges of that 4-cube. Nine estimated parameters: `a_basal`,
seven motif-specific multipliers, and one shared `sigma` (deacetylation `d` is fixed at 1).
The scientific question is site-independent versus motif-specific acetylation.

The data are the **stationary** abundances of the states (they sum to 1 — the model carries
that conservation law), so all 252 measurement rows sit at `time = inf` under the single
condition `control`: no time axis and no dose axis. Observables are log-transformed with
normal noise. That combination is why Blasi needed two separate PyBNF fixes — lanl/PyBNF#509
(`lnnormal`, ADR-0084) to import, then lanl/PyBNF#521 (steady state, ADR-0086) to load and fit.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over 380 optimizer runs) | `-1090.5618245715154` |
| paper-scale NLL at the PyBNF best fit | `-1090.561825` |
| optimality gap of the fit | `-4.3e-07` |
| paper-scale NLL at the PEtab nominal point | `-642.8268965824433` |
| optimality gap at nominal | `447.7349279890721` |
| scored data points `n` | 252 |
| free parameters `k` | 9 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`data/best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

**The nominal point is not this problem's optimum**, unlike the ten 🟢 slugs. Its eight
acetylation rate constants *are* the published MLEs of Table S1/S2, but `nominalValue` for
the noise scale is a placeholder `sigma = 0.1` where the published fit profiles `sigma` to
`0.2532`. That one parameter is the whole of `OG_nominal = 447.7`. Holding the eight
published rates and profiling `sigma` gives `OG = 7.5e-04` — so `J*` *is* the paper's own
motif-specific fit. See `nominal_check.json`.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through
trf's Coleman–Li core, ADR-0068) — handles this problem's estimated noise scale, which plain
`trf` refuses. `sbml_backend = bngsim` is required: the gradient path needs bngsim's forward
sensitivities. The steady-state relaxation is differentiable, so `gntr` is a valid choice
here. 20 starts × 500 iterations converge in about 25 s on a laptop.

## Recovered parameters

The fit was run from the box, with no knowledge of the published values:

| parameter | fitted | published (Table S2, rank 1) |
|---|---|---|
| `a_basal` | 0.066758 | 0.0668 |
| `a_k8` | 0.027271 | 0.0273 |
| `a_k5_k5k12` | 2.062033 | 2.0620 |
| `a_k12_k5k12` | 0.551941 | 0.5519 |
| `a_k16_k12k16` | 0.695939 | 0.6959 |
| `a_k5k12_k5k8k12` | 0.325297 | 0.3253 |
| `a_k12k16_k8k12k16` | 2.205530 | 2.2055 |
| `a_k8k12k16_4ac` | 3.591695 | 3.5917 |
| `sigma` | 0.253210 | (profiled; not tabulated) |

## Contents

- `Blasi_CellSystems2016.conf` — the PyBNF job
- `model_Blasi_CellSystems2016.xml` — SBML model (byte-identical to the upstream PEtab problem)
- `experiment1.exp`, `experiment1_rep2.exp` … `experiment1_rep18.exp` — the 18 biological
  replicates, dealt one row per file by the importer
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `best_fit_params.txt`, `information_criteria.txt` — from the shipped fit
- `VALIDATION.md` — the validation write-up
- `score.py` — scores a run against `J*`

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job` from
upstream commit `4d20850`. The converter preserves both `parameterScale` (lanl/PyBNF#491)
and `observableTransformation` (lanl/PyBNF#499), which plain `petab.v2.petab1to2` drops. The
run recipe (`job_type`, `sbml_backend = bngsim`, `wall_time_sim`) is supplied, not recovered —
PEtab specifies a problem, not a method.

One hand edit to the emitted conf: the importer also writes `noise_model x_k5k16 = lnnormal,
sigma = read_exp_file _SD`. That line is dead — `observable_k5k16` is declared in the PEtab
observables table but carries zero measurement rows (the K5K16 motif is below the LC-MS
quantification limit), so no `.exp` file has an `x_k5k16` column and the line is never
consulted. It is dropped so the conf does not imply a 16th fitted observable. Worth filing
upstream as an importer nit.

## Related

`models/combinatorial_histone_h4_acetylation_blasi2016/` in this repository curates the same
publication as house-style BNGL — the same 16-motif network written as context-resolved rules,
in the three scenarios (unspecific / site-specific / motif-specific) that Fig. 2 compares.

## Running

```bash
pybnf -c Blasi_CellSystems2016.conf -o
python score.py output
```
