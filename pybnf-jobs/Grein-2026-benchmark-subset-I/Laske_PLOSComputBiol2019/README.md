# Laske_PLOSComputBiol2019

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**SOLVED** — `OG = -1e-06` from a from-scratch 100-start `gntr` fit: the fit reaches the
benchmark reference optimum `J*` itself, to six decimal places. Unlike most solved slugs
here the PEtab nominal point is **not** this problem's optimum (`OG = 39.9`), so this
result comes entirely from the optimizer, not from the imported parameter values. See
`VALIDATION.md`.

> **The nominal check was recomputed on 2026-08-02** after **lanl/PyBNF#531**, and the
> number moved from `OG = 96.7` to `OG = 39.9`. This model is a COPASI export: every rate
> law reads a `ModelValue_*` alias that an SBML `initialAssignment` derives from the real
> name (`ModelValue_79 = k_syn_R_M`, `ModelValue_80 = k_syn_P`, 27 in all), and *none* of
> the source names appears in a rate law directly. PyBNF's fast simulation path never
> recomputed such a derived parameter, so the fitted `k_syn_R_M` and the condition target
> `k_syn_P` were both **inert** — which also meant the `k_syn_P = 0` condition
> (`experiment____virus_infection_3`, the no-protein-synthesis experiment) simulated
> identically to `experiment____virus_infection`. Any number recorded here before that
> date is not comparable.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `276.05406127180015` |
| paper-scale NLL at the fitted optimum | `276.0540604` |
| **optimality gap from the fit** | **`-1e-06`** |
| paper-scale NLL at the PEtab nominal point | `315.90591054673496` |
| optimality gap at nominal | `39.85184927493481` |
| scored data points `n` | 42 |
| free parameters `k` | 13 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (ADR-0068),
which handles this problem's estimated noise scales — with **100 starts × 1000
iterations**. That is a deliberately larger budget than the collection's usual
20 × 500, and this slug is the reason the distinction is worth stating: at 20 × 500 the
same method reaches only `OG = 6.76`. The extra starts are what find the basin.

Getting here took two PyBNF fixes and one budget increase, in that order — see
`VALIDATION.md` for the full account:

| | `OG` |
|---|---|
| 20 × 500, before lanl/PyBNF#531 (stale derived parameters) | 96.7 at nominal, no usable fit |
| 20 × 500, after #531 | 6.76 |
| 100 × 1000, after #531, gradient still missing a column (#534) | 0.104 |
| 100 × 1000, after #531 **and** #534 | **-1e-06** |

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
- `best_fit_params.txt`, `information_criteria.txt` — the shipped fit's provenance
- `VALIDATION.md` — the full validation against `J*`
- `score.py` — scores a run against `J*`

## Running

```bash
pybnf -c Laske_PLOSComputBiol2019.conf -o
python score.py output          # or: python score.py   (scores the shipped provenance)
```
