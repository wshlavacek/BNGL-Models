# Laske_PLOSComputBiol2019

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**Setup only** — the job imports, simulates, and scores, but the PEtab nominal point
is **not** the published optimum (OG = 96.7 ≫ 1.92), so nothing about optimality is
claimed. This is a ready-to-run job, not a result.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `276.05406127180015` |
| paper-scale NLL at the PEtab nominal point | `372.75319082302` |
| optimality gap at nominal | `96.69912955121987` |
| scored data points `n` | 42 |
| free parameters `k` | 13 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (ADR-0068),
which handles this problem's estimated noise scales. This is a **default recipe, not
a tuned one**; expect to tune `population_size` / `max_iterations` (or switch to
`cmaes` with IPOP restarts) before treating a run as a statement about PyBNF's
optimizers.

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`.
This problem exercises **lanl/PyBNF#509**: several observables use
`observableTransformation = log` (**natural** log), which imports to PyBNF's
`lnnormal` noise family — distinct from the log10 `lognormal` family. It is mixed
with linear/Gaussian observables in the same problem. Before #509 the natural-log
family was unimplemented and import refused; it now loads, simulates, and scores. The
run recipe (`job_type`, `sbml_backend = bngsim`, `wall_time_sim`) is supplied, not
recovered.

## Contents

- `Laske_PLOSComputBiol2019.conf` — the PyBNF job
- `model_Laske_PLOSComputBiol2019.xml` — SBML model (emitted by the importer)
- `experiment____virus_infection*.exp` — experimental data (3 replicate files)
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `score.py` — scores a run against `J*`

## Running

```bash
pybnf -c Laske_PLOSComputBiol2019.conf -o
python score.py output
```
