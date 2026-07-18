# Benchmark-Models-PEtab — PyBNF on the Grein et al. 2026 optimizer benchmark

PyBNF fitting jobs for problems from the community **PEtab Benchmark-Models** collection, each
scored against the reference objective from the **Grein et al. 2026** optimizer benchmark. These are
the first data points placing PyBNF's optimizers on that leaderboard.

> Grein T, Penas DR, Weindl D, Lakrisenko P, Banga JR, Hasenauer J.
> **"A benchmark collection for optimizer evaluation in systems biology"** (working title).
> *bioRxiv* 2026.07.11.737731 — 33 optimizers × 30 PEtab problems, >1.5M core-hours.
> Data: `Benchmarking-Initiative/Benchmark-Models-PEtab` (the problems, PEtab v1) +
> `ICB-DCM/optimizer-benchmark-2026-suppl-code-and-data` (reference objectives `best_fx_marvin.csv`).

Unlike the hand-built BNGL jobs elsewhere in `pybnf-jobs/`, these are **PEtab-imported SBML** jobs
(`edition = 2`, `sbml_backend = bngsim`): mass-action ODE models converted from the PEtab v1
collection through PyBNF's importer.

## The scoring: optimality gap (OG)

The benchmark scores each fit by its **optimality gap**

    OG = J_pybnf_paperscale − J*        "solved" iff  OG < 1.92     (χ², α = 0.05, 1 dof)

where **J\*** = `min` over all optimizer runs on the *Marvin* cluster of the paper's Eq. 6 Gaussian
negative log-likelihood (`suppl/data/best_fx_marvin.csv`, column `fx_best`), and
`J_pybnf_paperscale` is PyBNF's best objective put on the paper's Eq. 6 scale.

## Objective fidelity (the make-or-break check)

Every subset-I problem estimates its measurement noise as a **free `sigma`/`sd_*` parameter**, so the
objective is the *full* Gaussian NLL (Eq. 6, with the `log(2πσ²)` normalizer), not a bare
sum-of-squares. PyBNF **minimizes a reduced objective** that drops the parameter-independent per-point
constants — `½log(2π)`, and (for a log10-transformed observable) the change-of-variables Jacobian
`Σ log(y_obs·ln10)` — because they do not affect the argmin. It then **reports the full normalized
log-likelihood** at the best fit in `information_criteria.txt` (matching `scipy.stats.norm.logpdf` /
`lognorm.logpdf`), which restores every dropped constant. Therefore

    J_paper  ==  −log_likelihood            (PyBNF's own information_criteria.txt)

exactly, for **both** linear and log10 observables, and

    OG  =  −log_likelihood − J*.

This is what each slug's `score.py` computes. The identity was verified three independent ways on
Boehm: analytically (from the Gaussian per-point kernel), by evaluating the objective at the published
best-fit (`J_pybnf(nominal) + N·½log(2π)` = J\* to 3×10⁻⁵), and against PyBNF's own reported `lnL`.
On the log10 problems it is confirmed by the near-exact recovery of J\* (Perelson OG ≈ 5×10⁻⁷).

## observableTransformation — a fixed import gap

PEtab's `observableTransformation` sets the *scale* the Gaussian noise is additive on:

| `observableTransformation` | noise on | PyBNF `noise_model` family | residual scored |
|---|---|---|---|
| `lin` | linear | `gaussian` | `sim − obs` |
| `log10` | log10 | `lognormal` (`Gaussian(additive_on=LOG10)`) | `log10(sim) − log10(obs)` |

**PyBNF's `import_job` currently reads `noiseDistribution` but ignores `observableTransformation`**, so
a `log10` observable wrongly imports as linear `gaussian`. The three log10 problems here
(Perelson, Borghans, Elowitz) carry a **hand-applied correction**: `gaussian → lognormal` on the
`noise_model` line (documented per slug). This is filed as a PyBNF import gap alongside the earlier
subset-I blockers (#492–#496). The linear problems (Boehm, Sneyd, Okuonghae) import correctly.

## The 6 problems (subset-I "GO" set)

Selected in `dev/subset_I_triage.md` as the gradient-fittable subset (import + simulate cleanly). J\*
re-derived from `best_fx_marvin.csv`.

| slug | J\* | obs scale | params | fit | OG | status |
|---|---|---|---|---|---|---|
| `Boehm_JProteomeRes2014` | 138.2219682 | lin | 9 | gntr | 0.0012 | ✅ solved |
| `Perelson_Science1996` | 222.2807689 | log10 | 3 | cmaes | 5e-7 | ✅ solved |
| `Sneyd_PNAS2002` | −319.7923458 | lin | 15 | gntr | 1.4e-5 | ✅ solved |
| `Borghans_BiophysChem1997` | −132.0084765 | log10 | 23 | cmaes | — | ⛔ blocked (#498) |
| `Elowitz_Nature2000` | −65.6351201 | log10 | 21 | cmaes | — | ⛔ blocked (#498) |
| `Okuonghae_ChaosSolitonsFractals2020` | 373.5476580 | lin | 16 | cmaes | — | ⛔ blocked (#498) |

**Status:** 3 solved (Boehm, Perelson, Sneyd). The three multimodal problems (Borghans, Elowitz,
Okuonghae) need **multistart / IPOP-CMA-ES** — filed as **lanl/PyBNF#498**; a single CMA-ES run reaches
only a local minimum (paper success rates 3/380, 176/380, 249/380). The log10-observable import gap
(hand-corrected in the 3 log10 slugs) is filed as **lanl/PyBNF#499**.

## Import + fit pipeline (reproduce)

```bash
# 1. import (log10 parameter scales preserved via PR #491; needs PyBNF main)
python -c "from pybnf.petab import petab1to2_preserve_scale, import_job; import tempfile; \
  v2 = petab1to2_preserve_scale('<id>/<id>.yaml', tempfile.mkdtemp()); \
  import_job(v2, 'out', job_type='trf')"
# 2. edit the conf: sbml_backend = bngsim (REQUIRED for gradients); for a log10 observable,
#    change noise_model 'gaussian' -> 'lognormal'; pick job_type (see below); set random_seed.
# 3. fit + score
pybnf -c <id>.conf -o -L none
python score.py output          # OG = -log_likelihood - J*
```

**Optimizer choice.** `gntr` (the general-objective EFIM trust-region — PyBNF's fides-analog) is the
faithful gradient method and solves the well-conditioned **linear** problems. On the **log10
(lognormal)** objective it hits an EFIM `SVD did not converge`, and `lbfgs` drives the wide-bounded
`sigma` to `nan`, so those use **`cmaes`** (robust, gradient-free — the paper's #2 method). A stiff
model whose forward-sensitivity integration floods CVODE with `mxstep` failures (Okuonghae; the
Elowitz repressilator) also falls back to `cmaes`. Fallbacks are noted per slug. (The `#492` NoneType
gradient crash on a failed sim still bites the linear gradient path at bad points.)

## Per-slug contents

`<id>.conf` (the runnable fit) · `model_<id>.xml` (SBML, verbatim) · `experiment*.exp` (data) ·
`best_fit_params.txt` (top best-fit rows) · `information_criteria.txt` (PyBNF's lnL/AIC/BIC at best
fit) · `jstar.txt` (the reference J\*) · `score.py` (recompute OG) · `README.md` + `VALIDATION.md`.
