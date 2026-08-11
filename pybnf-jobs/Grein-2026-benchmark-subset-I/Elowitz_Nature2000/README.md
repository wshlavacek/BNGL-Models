# Elowitz_Nature2000

**Run cost: `hours`** — the repressilator, 21 free parameters. The solve appeared **~34 minutes** in
and the run was carried to **54 minutes** (1,100 `gntr` starts in eleven batches across three
concurrent slots); the tier stays `hours` because the target basin is reached by only ~1 start in
~1,000, so a single short run is not expected to find it.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**✅ Solved — `OG = 0.000175`, from unbiased starts.** A `gntr` multistart over the full PEtab box
(eleven independent 100-start batches, 1,100 starts, 54 minutes, no seeding) reproduces the benchmark's
`J*` to `1.75e−04` and beats this problem's own nominal point (`OG_nominal = 2.4324`). See
[`VALIDATION.md`](VALIDATION.md) for Gates A/B/C.

The problem is genuinely multimodal and the target basin is rare — reached by roughly **1 start in
~1,000**, against a dominant attractor at `OG ≈ 5.81` that takes ~55% of batches. That rarity, not any
defect, is why two earlier attempts missed it (a 20 h `cmaes` run at `OG = 5.1051`, a 100 × 1000
`gntr` run at `5.8200`).

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `-65.63512012927485` |
| paper-scale NLL at the PEtab nominal point | `-63.20275066664393` |
| optimality gap at nominal | `2.4323694626309234` |
| scored data points `n` | 58 |
| free parameters `k` | 21 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

The shipped conf carries the importer's `job_type = cmaes` recipe (CMA-ES with IPOP restarts,
ADR-0070/0082), which is **not** what solved this problem. `cmaes` reached only `OG = 5.1051` in a
completed 20-hour run.

**Use `gntr`.** The reference region is a box *corner* — 6 of 21 nominal parameters sit at or
essentially on a bound — and a Gaussian sampler is adversarial to corners (to put a coordinate *at* a
bound its sampling mean must sit on or past it) while a clamping gradient method reaches them
naturally. The recipe that solved it:

```
job_type = gntr
population_size = 100        # gradient optimizers use this as the START COUNT
max_iterations = 500
parallel_count = 3
```

run as independent seeded batches across **three concurrent `pybnf` processes**, best kept. The
concurrency is the point: one process plateaus at ~310% CPU here regardless of `parallel_count`, so
three of them roughly triple throughput on a 10-core box. Keep `max_iterations` modest — 500 rather
than 1000 costs nothing on this problem (`5.81` vs the historical `5.82`) and buys twice the
independent draws, which is the only currency that matters on a basin-selection failure.

`population_size` itself does **not** affect steady-state throughput — measured 28.6
start-iterations/s at 100 against 34.9 at 1500. An earlier revision of this file claimed a 12x
collapse; that was an artefact of measuring during iteration 0 and is retracted. Batch for the
incremental scored checkpoints, not for speed.

## Contents

- `Elowitz_Nature2000.conf` — the PyBNF job
- `model_Elowitz_Nature2000.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment1.exp` — experimental data
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `best_fit_params.txt` — the solved fit's full sorted parameter table
- `information_criteria.txt` — AIC/BIC/AICc and the full normalized log-likelihood at that fit
- `VALIDATION.md` — Gates A/B/C for the solved row
- `score.py` — scores a run against `J*`

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`. The
converter preserves both `parameterScale` (lanl/PyBNF#491) and `observableTransformation`
(lanl/PyBNF#499), which plain `petab.v2.petab1to2` drops. The run recipe (`job_type`,
`sbml_backend = bngsim`, `wall_time_sim`) is supplied, not recovered — PEtab specifies a
problem, not a method. `wall_time_sim = 10` caps pathological parameter points; raise it
if valid simulations on your machine are being marked as failures.

## Running

The shipped conf runs and scores, but reaches only `OG ≈ 5.1` — see **Optimizer** above and use
`gntr` batches to reproduce the solved row.

```bash
pybnf -c Elowitz_Nature2000.conf -o
python score.py output
```

`python score.py` with no argument scores the shipped `best_fit_params.txt` /
`information_criteria.txt` and prints `SOLVED`.

Reproduced on bngsim `ffbf015` (local editable install), PyBNF `095a5a14`, numpy 2.5.2, against
upstream PEtab `4d20850`. `nominal_check.json`'s recorded values predate both — they were produced
under bngsim **0.11.35**, below PyBNF's own `>=0.12.2` pin; the current build reads
`OG_nominal = 2.432380054142861` against the recorded `2.4323694626309234`, a `1.1e−05` shift that is
immaterial against a 1.92 threshold. Whether the collection's nominal checks get regenerated under
0.12.2 remains an open decision in issue #38 and is **not** resolved here.
