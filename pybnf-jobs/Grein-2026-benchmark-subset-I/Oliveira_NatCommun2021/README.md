# Oliveira_NatCommun2021

**Run cost: `minutes`** — 100,000 evaluations (100 × 1,000 `gntr`) on a 12-reaction ODE model.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**SOLVED** — `OG = 0.011342` from a from-scratch 100-start `gntr` fit in 28 min, well inside the
threshold `OG < 1.92`. The PEtab nominal point is **not** this problem's optimum
(`OG_nominal = 9.6e+06`), so the optimizer crossed roughly seven orders of magnitude from unbiased
box starts — nothing about the answer was known before the fit ran. This was the collection's first
⚪ → ✅ conversion.

> **This slug has no independent oracle.** Upstream ships no `simulatedData_*.tsv` for it, so unlike
> the `obj ✓` rows its objective has never been cross-checked against anything but `J*` itself. That
> matters because a large `OG_nominal` is also what `Brannmark` and `Weber` looked like while
> lanl/PyBNF#547 had their objective wrong. The difference is that this fit *reaches* `J*` to 0.011,
> which a corrupted objective would not do — and Oliveira does not pre-equilibrate, so it was never
> exposed to #547. See `VALIDATION.md`.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `7904.93431737276` |
| paper-scale NLL at the PEtab nominal point | `9623530.12841406` |
| optimality gap at nominal | `9615625.194096688` |
| scored data points `n` | 120 |
| free parameters `k` | 12 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through trf's
Coleman–Li core, ADR-0068). `population_size = 100`, `max_iterations = 1000` — the collection's
documented working default, not the shipped 20 × 500 placeholder.

Noise is σ ≡ 1 upstream, so the conf ships a plain `objective = sos` with no `noise_model` line; the
collection README records why that is faithful rather than a fidelity break. This fit also sets the
upper end of where the 100 × 1000 default is known to work: it solves at k = 6, 9, **12**, 13 and
falls short of the reference basin at k = 22 (`Fiedler_BMCSystBiol2016`).

## Contents

- `Oliveira_NatCommun2021.conf` — the PyBNF job
- `Oliveira_NatCommun2021_model.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment____interior.exp` — experimental data
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
pybnf -c Oliveira_NatCommun2021.conf -o
python score.py output
```
