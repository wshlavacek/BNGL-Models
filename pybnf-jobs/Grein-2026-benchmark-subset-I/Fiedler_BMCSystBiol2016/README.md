# Fiedler_BMCSystBiol2016

**Run cost: `minutes`** — 10,000 evaluations (20 × 500 `gntr`), 22 free parameters.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**Objective validated at the PEtab nominal point** (OG = −0.0022, within the solved threshold 1.92). No optimization run has been performed here.

The observables are linear/Gaussian, so the nominal-point evaluation is a clean
end-to-end check: PyBNF's Eq. 6 NLL at the published point reproduces the reference
`J*` to ~2e−3. It makes no claim about PyBNF's optimizer.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `-58.58395532491055` |
| paper-scale NLL at the PEtab nominal point | `-58.58618893016114` |
| optimality gap at nominal | `-0.0022336052505949056` |
| scored data points `n` | 72 |
| free parameters `k` | 22 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian
through trf's Coleman–Li core, ADR-0068) — handles this problem's estimated noise
scale (`sigma_pErk`, `sigma_pMek`), which plain `trf` refuses. This is a **default
recipe, not a tuned one.**

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`.
This problem exercises **lanl/PyBNF#508**: the per-gel scale factors (`s_pErk_*`,
`s_pMek_*`) are replicate-specific `observableParameters` supplied through the
`measurement_params:` sidecars (`experiment____model1_data*_measparams.tsv`). Before
#508 the replicate dimension was dropped and those bindings were silently lost; the
import now loads, simulates, and scores. The run recipe (`job_type`, `sbml_backend =
bngsim`, `wall_time_sim`) is supplied, not recovered — PEtab specifies a problem, not
a method.

## Contents

- `Fiedler_BMCSystBiol2016.conf` — the PyBNF job
- `model_Fiedler_BMCSystBiol2016.xml` — SBML model (emitted by the importer)
- `experiment____model1_data*.exp` (+ `_rep2`) — experimental data, 3 conditions × 2 replicates
- `experiment____model1_data*_measparams.tsv` — per-measurement observable/noise parameter tables (#508)
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `score.py` — scores a run against `J*`

## Running

```bash
pybnf -c Fiedler_BMCSystBiol2016.conf -o
python score.py output
```
