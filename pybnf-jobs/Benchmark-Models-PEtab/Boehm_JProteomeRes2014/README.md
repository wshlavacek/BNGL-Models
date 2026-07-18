# Boehm_JProteomeRes2014 — STAT5A/B phosphorylation dynamics (PEtab benchmark, PyBNF gntr)

A PyBNF **edition-2, PEtab-imported SBML** fitting job for the `Boehm_JProteomeRes2014` problem from
the community PEtab Benchmark-Models collection, scored against the **Grein et al. 2026** optimizer
benchmark. See [`../README.md`](../README.md) for the collection, the scoring, and the objective-fidelity
methodology.

> **Result: SOLVED.** OG = **0.0012** (threshold < 1.92). J\* = 138.2219682, PyBNF J_paper = 138.2231909.

> Boehm ME, Adlung L, Schilling M, Roth S, Klingmüller U, Lehmann WD.
> **"Identification of isoform-specific dynamics in phosphorylation-dependent complex formation by
> quantitative mass spectrometry."** *J Proteome Res* 2014; 13(12):5685. doi:10.1021/pr5006923.

## The problem

STAT5A/STAT5B phosphorylation, homo/heterodimerization, and nuclear import/export in BaF3-EpoR cells
after Epo stimulation (SBML: **8 species, 9 reactions, 2 compartments**). Three observables are
**ratio formulas** of the phospho/total STAT5 species (`observableTransformation = lin`, so linear
`gaussian` noise — imported correctly):

- `pSTAT5A_rel`, `pSTAT5B_rel` — relative phospho-STAT5A / -B
- `rSTAT5A_rel` — relative STAT5A

Data: one condition, **48 measurements** (16 time points × 3 observables), in `experiment1.exp`.

## What is fit (9 free parameters, all log10-scaled)

6 kinetic constants (`Epo_degradation_BaF3`, `k_exp_hetero`, `k_exp_homo`, `k_imp_hetero`,
`k_imp_homo`, `k_phos`) + **3 estimated noise σ** (`sd_pSTAT5A_rel`, `sd_pSTAT5B_rel`,
`sd_rSTAT5A_rel`). Two parameters are fixed in the SBML (`ratio = 0.693`, `specC17 = 0.107`). Because σ
is estimated, the objective is the full Gaussian NLL (Eq. 6). Log10 parameter scales are preserved on
import (PR #491).

## Result

Multi-start **`gntr`** (the EFIM trust-region; 20 starts × 500 iterations, box-center + Latin-hypercube,
`random_seed = 1`), from scratch — not seeded at the published values:

| quantity | value |
|---|---|
| PyBNF reduced objective (minimized) | 94.1141413 |
| + restored `N·½log(2π)` (N=48) | 44.1090496 |
| **J_paper = −log_likelihood** | **138.2231909** |
| reference **J\*** (Grein, `best_fx_marvin`) | 138.2219682 |
| **OG = J_paper − J\*** | **0.0012** → **SOLVED** |

`information_criteria.txt`: `lnL = −138.2231909`, AIC = 294.45, BIC = 311.29 (k=9, n=48). See
[`VALIDATION.md`](VALIDATION.md) for the three-way fidelity confirmation.

## Run

```bash
pybnf -c Boehm_JProteomeRes2014.conf -o -L none     # ~4 min, bngsim forward-sensitivity gntr
python score.py output                              # OG from the fresh run
python score.py                                     # OG from the shipped best-fit provenance
```
