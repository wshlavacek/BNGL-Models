# Giordano_Nature2020

**Run cost: `hours`** — 100,000 evaluations (100 × 1,000 `gntr`), 50 free parameters on a 13-reaction model.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**SOLVED** — `OG = 0.135007` from a from-scratch 100-start `gntr` fit in 2 h 26 m, inside the
threshold `OG < 1.92`. This is the collection's largest problem on both axes (k = 50, n = 313), and
its acceptance window is the tightest here: with all 7 σ estimated over n = 313, clearing 1.92
requires matching the reference residual norm to 0.61%, and the fit matches it to 0.043%.

The PEtab nominal point is **not** this problem's optimum (`OG_nominal = 3776`), so this result comes
entirely from the optimizer. Four-fifths of that nominal distance was an artefact of placeholder σ:
with the σ profiled at the nominal dynamics the honest starting distance is 743.7, and that analysis —
run *before* the fit — predicted the fit would need to shrink the residual norm 10.8× below the
reference trajectory. It achieved 10.6×. The objective is independently corroborated against
upstream's own `simulatedData` (`obj ✓`). See `VALIDATION.md`.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `-3488.3414981097503` |
| paper-scale NLL at the PEtab nominal point | `287.62776140983027` |
| optimality gap at nominal | `3775.9692595195806` |
| scored data points `n` | 313 |
| free parameters `k` | 50 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through trf's
Coleman–Li core, ADR-0068) — handles this problem's estimated noise scales, which plain `trf`
refuses. `population_size = 100`, `max_iterations = 1000` — the collection's documented working
default, not the shipped 20 × 500 placeholder.

This slug is where **lanl/PyBNF#546 / ADR-0103** was found. Its assembled gradient disagreed with
central differences on 41 of 50 parameters, by up to 26%; the cause was the absolute tolerance, not
the 13 time gates the disagreement appeared to partition along. Giordano is a population-*fraction*
model whose species sit at `1.7e−8 … 1`, so bngsim's default `atol = 1e-8` buried the relative term
across the whole early trajectory. Deriving `atol` from the model's own scale moved the worst
relative error from 7.7e−02 to 4.5e−04, and is what made this problem tractable at all.

## Contents

- `Giordano_Nature2020.conf` — the PyBNF job
- `model_Giordano_Nature2020.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment____pred1.exp` — experimental data
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `score.py` — scores a run against `J*`
- `best_fit_params.txt`, `information_criteria.txt` — the shipped fit's provenance
- `VALIDATION.md` — the full validation against `J*`

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`. The
converter preserves both `parameterScale` (lanl/PyBNF#491) and `observableTransformation`
(lanl/PyBNF#499), which plain `petab.v2.petab1to2` drops. The run recipe (`job_type`,
`sbml_backend = bngsim`, `wall_time_sim`) is supplied, not recovered — PEtab specifies a
problem, not a method. `wall_time_sim = 10` caps pathological parameter points; raise it
if valid simulations on your machine are being marked as failures.

## Running

```bash
pybnf -c Giordano_Nature2020.conf -o
python score.py output
```
