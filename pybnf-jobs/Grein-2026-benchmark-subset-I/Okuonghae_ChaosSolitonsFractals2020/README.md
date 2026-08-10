# Okuonghae_ChaosSolitonsFractals2020

**Run cost: `minutes`** — 64,000 evaluations (32 × 2,000 `gntr`).

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**SOLVED** — `OG = 0.001206` from a from-scratch 32-start `gntr` fit in 12 min 47 s, well inside the
threshold `OG < 1.92`, from a nominal point `4.7×10⁵` away. This slug is the collection's
**reclassified** case: it was carried as strongly multimodal and `cmaes`-only, and it is neither. See
the Optimizer section below and `VALIDATION.md`.

> **The ✅ is on the objective only.** Okuonghae & Omame estimate six quantities and fix the other ten
> from the literature; the PEtab problem frees all sixteen. The best-fit point puts the symptomatic
> death rate at 7.25 /day against the paper's 0.015 — a mean time from symptoms to death of about
> three hours — with `alpha` and `gamma_a` resting on their lower bounds. **It is biologically
> meaningless, and it still matches `J*` to 0.0012**, so the reference optimum has the same character.
> That is a property of the benchmark problem, not of this fit.

> **This slug has no independent oracle.** Upstream ships no `simulatedData_*.tsv` for it, so it
> cannot take the nominal-point recomputation the `obj ✓` rows carry. It does not pre-equilibrate, so
> it was never exposed to lanl/PyBNF#547.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `373.5476580181036` |
| paper-scale NLL at the PEtab nominal point | `470490.51145926415` |
| optimality gap at nominal | `470116.963801246` |
| scored data points `n` | 92 |
| free parameters `k` | 16 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through trf's
Coleman–Li core, ADR-0068) — handles this problem's estimated noise scale, which plain `trf` refuses.
`population_size = 32`, `max_iterations = 2000`: the budget the old `cmaes` conf carried, left
unchanged so the method is the only variable this result speaks to. The collection's 100 × 1000
default has not been tried here and may well be faster.

### This slug shipped `cmaes`, and the reason was inherited rather than measured

Until 2026-08-07 this job carried `job_type = cmaes` with 12 IPOP restarts, and the collection tracked
it as one of three **multimodal** problems where "a local method from a few starts reliably lands in a
wrong basin; they need a global budget, not a better local step."

For this problem that grouping is simply wrong. There was never a capability reason for `cmaes` — the
model has **0 events, 0 `piecewise`, 0 `and`**, 9 species and 10 reactions, so nothing gates the
gradient path — and `tools/fd_check.py` verifies the assembled gradient against central differences at
**6.27e−05**, the second best in the corpus. Run on `gntr` at the same 32 × 2000 budget, it solves in
under 13 minutes.

The other two slugs in that group behave as the flag claimed, which is what makes this a
reclassification rather than a repudiation:

| slug | k | `OG_nominal` | `gntr` result |
|---|---:|---:|---|
| **`Okuonghae`** | 16 | 4.7e+05 | **`OG = 0.0012` — solved, 12 m 47 s** |
| `Borghans_BiophysChem1997` | 23 | 48.7 | `OG = 68.9` — *worse* than its own nominal |
| `Elowitz_Nature2000` | 21 | 2.43 | `OG = 5.82` — *worse* than its own nominal |

`Borghans` and `Elowitz` show the genuine wrong-basin signature and stay on `cmaes`. Only this one
moved.

## Contents

- `Okuonghae_ChaosSolitonsFractals2020.conf` — the PyBNF job
- `model_Okuonghae_ChaosSolitonsFractals2020.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment____cond1.exp` — experimental data
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
pybnf -c Okuonghae_ChaosSolitonsFractals2020.conf -o
python score.py output
```
