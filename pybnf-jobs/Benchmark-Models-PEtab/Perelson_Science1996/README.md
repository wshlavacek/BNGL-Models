# Perelson_Science1996 — HIV-1 viral decay dynamics (PEtab benchmark, PyBNF cmaes)

A PyBNF **edition-2, PEtab-imported SBML** fitting job for the `Perelson_Science1996` problem, scored
against the **Grein et al. 2026** optimizer benchmark. See [`../README.md`](../README.md) for the
collection, scoring, and fidelity methodology. This is the collection's cleanest **log10-observable**
exemplar.

> **Result: SOLVED.** OG = **5×10⁻⁷** (threshold < 1.92). J\* = 222.2807689, PyBNF J_paper = 222.2807694.

> Perelson AS, Neumann AU, Markowitz M, Leonard JM, Ho DD.
> **"HIV-1 dynamics in vivo: virion clearance rate, infected cell life-span, and viral generation
> time."** *Science* 1996; 271(5255):1582. doi:10.1126/science.271.5255.1582.

## The problem

The Perelson HIV-1 viral-decay model after protease-inhibitor therapy (SBML: **4 species, 5
reactions**). One observable — plasma virus `V` — measured at **16 time points** (one condition,
`experiment1.exp`), spanning V ≈ 9×10⁴ … 3×10⁶ copies.

## ⚠️ log10 observable — hand-applied `lognormal` correction

Perelson's observable has **`observableTransformation = log10`**: the paper scores the Gaussian NLL on
`log10(V)`, and its Eq. 6 objective includes the change-of-variables Jacobian `Σ log(V_obs·ln10)`.
PyBNF's `import_job` ignores `observableTransformation` (it reads only `noiseDistribution`), so it
emitted a linear `gaussian` line. **Corrected by hand** to the log10-additive Gaussian family:

```
noise_model = lognormal, sigma = fit sd_task0_model0_perelson1_V     # was: gaussian
```

`lognormal` = `Gaussian(additive_on=LOG10)`, so the residual scored is `log10(V_sim) − log10(V_obs)`.
Without this, PyBNF's linear objective bottoms out at J_paper = 232.3 (OG ≈ 10, *unsolvable*); the
linear best-fit noise cannot beat σ̂ ≈ 4.9×10⁵ because the model (a monotone decay) cannot follow the
data's early rise, whereas on the log10 scale the same fit is excellent. The importer's dropped
`observableTransformation` is filed as a PyBNF import gap. **(The PEtab nominal `sd = 1e5` is a
linear-scale artifact — ignore it; the fitted log10 σ ≈ 0.12.)**

## What is fit (3 free parameters, all log10-scaled)

`c` (virion clearance rate), `delta` (infected-cell death rate), and the estimated noise
`sd_task0_model0_perelson1_V`. `NN`, `T0`, `K0` are fixed in the SBML. Best fit:
**c = 1.861, delta = 0.547, σ = 0.123** (log10 scale).

## Result

Multi-start **`cmaes`** (24 individuals × ≤500 generations, `random_seed = 1`, `sbml_backend = bngsim`):

| quantity | value |
|---|---|
| PyBNF reduced objective (minimized, log10 residuals) | −25.5509736 |
| + restored constants (`½log2π` + Jacobian, = 247.83) | |
| **J_paper = −log_likelihood** | **222.2807694** |
| reference **J\*** (Grein) | 222.2807689 |
| **OG = J_paper − J\*** | **5×10⁻⁷** → **SOLVED** |

`cmaes` (not `gntr`/`lbfgs`) because on the lognormal objective the EFIM trust-region hits
`SVD did not converge` and `lbfgs` drives the wide-bounded `sigma` (`[1e-10, 1e10]`) to `nan`. All 380
of the paper's Marvin runs solved Perelson (min = median), and PyBNF matches to 5×10⁻⁷.

## Run

```bash
pybnf -c Perelson_Science1996.conf -o -L none
python score.py output        # OG = -log_likelihood - J*  (the -lnL already carries the log10 Jacobian)
python score.py               # OG from the shipped provenance
```
