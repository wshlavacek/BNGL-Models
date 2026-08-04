# Smith_BMCSystBiol2013

**Run cost: `hours`** — 64,000 evaluations (32 × 2,000 `cmaes`) over a 367-reaction network — the subset's most expensive product of budget and model size.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**Setup only — not fitted.** The job runs and scores correctly; the PEtab nominal point is not this problem's published optimum, so no optimality claim is made.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `20922.16424399946` |
| paper-scale NLL at the PEtab nominal point | `6.851921116400144e+32` |
| optimality gap at nominal | `6.851921116400144e+32` |
| scored data points `n` | 62 |
| free parameters `k` | 25 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = cmaes` — CMA-ES with IPOP restarts (ADR-0070/0082), a global search, because
**the gradient path genuinely refuses this problem**. The model contains discrete events
(state-dependent jumps); bngsim's forward output sensitivities are not reinitialised across
such a jump, so they go silently stale at and after one fires, and bngsim refuses to supply
them rather than return wrong derivatives. PyBNF reads that as a pre-flight gate
(`_require_differentiable_dynamics`, lanl/PyBNF#461) and refuses `gntr`/`trf`/`lbfgs` up
front. This is the one remaining gradient refusal in subset I, and it needs event-aware
sensitivities — nothing in the ADR-0076 condition-routing line touches it.

The shipped recipe was verified to start and run on this problem.

## Contents

- `Smith_BMCSystBiol2013.conf` — the PyBNF job
- `model_Smith_BMCSystBiol2013.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment____figure2B__0_0__1_3em10.exp`, `experiment____figure2B__0_0__1_5em08.exp`, `experiment____figure2B__0_0__1_5em09.exp`, `experiment____figure2B__0_0__1_5em10.exp`, `experiment____figure2B__0_0__1_5em11.exp`, `experiment____figure2B__0_0__1em06.exp` … — experimental data
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `score.py` — scores a run against `J*`

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`. The
converter preserves both `parameterScale` (lanl/PyBNF#491) and `observableTransformation`
(lanl/PyBNF#499), which plain `petab.v2.petab1to2` drops. The run recipe (`job_type`,
`sbml_backend = bngsim`, `wall_time_sim`) is supplied, not recovered — PEtab specifies a
problem, not a method. `wall_time_sim = 10` caps pathological parameter points; raise it
if valid simulations on your machine are being marked as failures.

## Running

```bash
pybnf -c Smith_BMCSystBiol2013.conf -o
python score.py output
```
