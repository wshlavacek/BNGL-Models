# Vanhoefer et al. 2026 — measurement-time uncertainty

Vanhoefer, J., Nakonecnij, V., Binder, N., & Hasenauer, J. (2026). *Efficient Bayesian inference for
ordinary differential equation models from experimental data with uncertain measurement times.*
bioRxiv [2026.05.09.724053](https://doi.org/10.64898/2026.05.09.724053).

The paper formulates the latent measurement time as a random variable and **marginalizes it out** of
the likelihood, so a fit accounts for timing uncertainty instead of mistaking it for dynamics. PyBNF
implements this as the `time_error` clause (lanl/PyBNF #587, ADR-0112 — phase 1, quadrature; #588,
ADR-0113 — phase 2, the gradient). This directory holds the paper's real-data application as a PyBNF
job.

## Shared model

The paper's application (Fig 5/6) is the carotenoid-cleavage model of *Arabidopsis thaliana* (Bruno
et al., *J Exp Bot* 67:5993, 2016): 7 state variables, 13 free parameters, 6 experimental conditions,
77 measurements. The SBML is imported verbatim from the Benchmark-Models-PEtab collection (the same
model the `Grein-2026-benchmark-subset-I/Bruno_JExpBot2016` optimizer-benchmark slug uses; see each
slug's `upstream.json`).

## Jobs

| slug | fits | archetype | data | J\* / oracle | status |
|---|---|---|---|---|---|
| [`carotenoid_time_uncertainty`](carotenoid_time_uncertainty/) | timing-uncertainty correction: standard vs. marginal (phase-2 gradient) fit on synthetic timing-perturbed data | B (SBML-imported) | synthetic, from the model at θ\* (`perturb_times.py`) | θ\* recovery + standard-vs-marginal LRT (Fig 5/6) | ✅ Gate B PASS; Gate A PASS on recovery (one documented fixed-grid LRT caveat) |

The oracle here is not a Grein optimality gap but the paper's own paired comparison: on data
corrupted by a known timing error, does the marginal fit correct the bias the standard fit suffers,
recover the timing scale, and reject the standard model? It does (`carotenoid_time_uncertainty/VALIDATION.md`).
