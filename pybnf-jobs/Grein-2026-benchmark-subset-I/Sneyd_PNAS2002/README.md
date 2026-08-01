# Sneyd_PNAS2002 — IP₃-receptor Ca²⁺ gating (PEtab benchmark, PyBNF gntr)

A PyBNF **edition-2, PEtab-imported SBML** fitting job for the `Sneyd_PNAS2002` problem, scored against
the **Grein et al. 2026** optimizer benchmark. See [`../README.md`](../README.md) for the collection,
scoring, and fidelity methodology. This is the collection's **multi-condition** exemplar (9 conditions,
one shared σ).

> **Result: SOLVED.** OG = **1.4×10⁻⁵** (threshold < 1.92). J\* = −319.7923458, PyBNF J_paper = −319.7923321.

> Sneyd J, Falcke M, Dufour JF, Fox C.
> **"A comparison of three models of the inositol trisphosphate receptor."** *Prog Biophys Mol Biol*
> 2004 — the model of Sneyd & Dufour, *PNAS* 2002; 99(4):2398 (doi:10.1073/pnas.032281999).

## The problem

The Sneyd–Dufour IP₃-receptor gating model — receptor states driving the channel open probability
(SBML: **6 species, 10 reactions**). One observable, **`open_probability`** (`observableTransformation =
lin` → linear `gaussian`, imported correctly), measured across **9 conditions**:

- 6 **Ca²⁺ dose-response** conditions (`Ca_dose_response__1…6`)
- 3 **IP₃ dose-response** conditions (`IP3_dose_response__1…3`)

**135 measurements** total across the 9 `experiment____*.exp` files. BPSL-free; the noise is a single
shared estimated σ across all conditions.

## What is fit (15 free parameters, all log10-scaled)

14 kinetic constants (`k1…k4`, `k_1…k_4`, `l2/l4/l6`, `l_2/l_4/l_6`) + **1 shared estimated σ**
(`sigma`). Because σ is estimated, the objective is the full Gaussian NLL (Eq. 6).

## Result

Multi-start **`gntr`** (the EFIM trust-region; 10 starts × 300 iterations, box-center + Latin-hypercube,
`random_seed = 1`, `sbml_backend = bngsim`):

| quantity | value |
|---|---|
| PyBNF reduced objective (minimized) | −443.8490341 |
| + restored `N·½log(2π)` (N=135) | 124.0567018 |
| **J_paper = −log_likelihood** | **−319.7923321** |
| reference **J\*** (Grein) | −319.7923458 |
| **OG = J_paper − J\*** | **1.4×10⁻⁵** → **SOLVED** |

(PyBNF's reported `lnL = +319.79` is *positive* here — the open-probability data with a small fitted σ
gives a Gaussian density > 1, so its log-density is positive; `−lnL` is the paper-scale NLL as usual.)
Sneyd is an easy problem (329/380 Marvin runs solved it); PyBNF matches J\* to 5 significant figures.

## Run

```bash
pybnf -c Sneyd_PNAS2002.conf -o -L none      # bngsim forward-sensitivity gntr, 9 conditions
python score.py output                        # OG from the fresh run
python score.py                               # OG from the shipped provenance
```
