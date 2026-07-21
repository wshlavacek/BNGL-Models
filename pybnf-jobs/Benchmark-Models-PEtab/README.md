# Benchmark-Models-PEtab — PyBNF on the Grein et al. 2026 optimizer benchmark

PyBNF fitting jobs for problems from the community **PEtab Benchmark-Models** collection, each
carrying the reference objective from the **Grein et al. 2026** optimizer benchmark. These are the
first data points placing PyBNF on that leaderboard.

> Grein T, Penas DR, Weindl D, Lakrisenko P, Banga JR, Hasenauer J.
> **"A benchmark collection for optimizer evaluation in systems biology"** (working title).
> *bioRxiv* 2026.07.11.737731 — 33 optimizers × 30 PEtab problems, >1.5M core-hours.
> Data: `Benchmarking-Initiative/Benchmark-Models-PEtab` (the problems, PEtab v1) +
> `ICB-DCM/optimizer-benchmark-2026-suppl-code-and-data` (reference objectives `best_fx_marvin.csv`).

Unlike the hand-built BNGL jobs elsewhere in `pybnf-jobs/`, these are **PEtab-imported SBML** jobs
(`edition = 2`, `sbml_backend = bngsim`): ODE models converted from the PEtab v1 collection through
PyBNF's importer.

## The scoring: optimality gap (OG)

    OG = −log_likelihood − J*        "solved" iff  OG < 1.92     (χ², α = 0.05, 1 dof)

where **J\*** = `min` over all optimizer runs on the *Marvin* cluster of the paper's Eq. 6 Gaussian
negative log-likelihood (`suppl/data/best_fx_marvin.csv`, column `fx_best`). Each slug ships
`score.py`, which computes exactly this.

## Objective fidelity

Every subset-I problem estimates its measurement noise as a **free `sigma`/`sd_*` parameter**, so the
objective is the *full* Gaussian NLL (Eq. 6, with the `log(2πσ²)` normalizer), not a bare
sum-of-squares. PyBNF **minimizes a reduced objective** that drops the parameter-independent per-point
constants — `½log(2π)`, and (for a log10 observable) the change-of-variables Jacobian
`Σ log(y_obs·ln10)` — because they do not affect the argmin. It then **reports the full normalized
log-likelihood** at the best fit in `information_criteria.txt` (matching `scipy.stats.norm.logpdf` /
`lognorm.logpdf`), restoring every dropped constant. Therefore

    J_paper  ==  −log_likelihood

exactly, for **both** linear and log10 observables.

**This identity is now corroborated across the collection, not just argued.** For nine problems the
PEtab `nominalValue` point is the published optimum, and evaluating PyBNF's objective there
reproduces the paper's `J*` to within the solved threshold — six to ~10⁻⁶ or better (see
`nominal_check.json` in each slug). That is an end-to-end check of the whole imported chain: SBML
model → simulation → observable formulas → noise model → objective.

## Coverage: 19 of the 23 subset-I problems

| slug | J\* | scale | k | n | optimizer | OG | status |
|---|---|---|---|---|---|---|---|
| `Armistead_CellDeathDis2024` | −301.9161878 | lin | 14 | 58 | gntr | 5.8e−06 † | 🟢 objective validated |
| `Boehm_JProteomeRes2014` | 138.2219682 | lin | 9 | 48 | gntr | 0.0012 | ✅ **solved** |
| `Bruno_JExpBot2016` | −46.6881979 | lin | 13 | 77 | cmaes | 3.2e−06 † | 🟢 objective validated |
| `Crauste_CellSystems2017` | 190.4570655 | lin | 12 | 21 | gntr | 0.509 † | 🟢 objective validated |
| `Perelson_Science1996` | 222.2807689 | log10 | 3 | 16 | cmaes | 5e−7 | ✅ **solved** |
| `Rahman_MBS2016` | 21.1534861 | lin | 9 | 23 | gntr | 3.9e−06 † | 🟢 objective validated |
| `Raia_CancerResearch2011` | 345.3097673 | lin | 39 | 205 | gntr | 0.78 † | 🟢 objective validated |
| `SalazarCavazos_MBoC2020` | 366.8615730 | lin | 6 | 18 | gntr | 0.326 † | 🟢 objective validated |
| `Sneyd_PNAS2002` | −319.7923458 | lin | 15 | 135 | gntr | 1.4e−5 | ✅ **solved** |
| `Elowitz_Nature2000` | −65.6351201 | log10 | 21 | 58 | cmaes | 2.43 † | ⚪ setup only |
| `Borghans_BiophysChem1997` | −132.0084765 | log10 | 23 | 111 | cmaes | 48.7 † | ⚪ setup only |
| `Zhao_QuantBiol2020` | 501.2270538 | lin | 28 | 82 | gntr | 276 † | ⚪ setup only |
| `Brannmark_JBC2010` | 141.8248543 | lin | 22 | 43 | gntr | 1.5e+03 † | ⚪ setup only |
| `Giordano_Nature2020` | −3488.3414981 | lin | 50 | 313 | gntr | 3.8e+03 † | ⚪ setup only |
| `Weber_BMC2015` | 296.2020025 | lin | 36 | 135 | gntr | 1.4e+04 † | ⚪ setup only |
| `Okuonghae_ChaosSolitonsFractals2020` | 373.5476580 | lin | 16 | 92 | cmaes | 4.7e+05 † | ⚪ setup only |
| `Oliveira_NatCommun2021` | 7904.9343174 | lin | 12 | 120 | gntr | 9.6e+06 † | ⚪ setup only |
| `Bertozzi_PNAS2020` | 158.8642627 | lin | 8 | 22 | cmaes | 1.8e+11 † | ⚪ setup only |
| `Smith_BMCSystBiol2013` | 20922.1642440 | lin | 25 | 62 | cmaes | 6.9e+32 † | ⚪ setup only |

`k` = free parameters, `n` = scored data points.
**† = optimality gap at the PEtab nominal point, not from a fit.** Only the three ✅ rows report an
OG from an actual optimization run.

Three status levels, and the difference matters:

- ✅ **solved** — a real fit was run and reached `OG < 1.92`.
- 🟢 **objective validated** — no fit was run, but the problem's `nominalValue` point *is* its
  published optimum, and PyBNF's objective there lands within the solved threshold of `J*`. This
  validates the import and the objective; it makes no claim about PyBNF's optimizer.
- ⚪ **setup only** — the job imports, simulates, and scores correctly, but its nominal point is not
  the published optimum, so nothing about optimality is claimed. These are ready-to-run jobs.

### The 4 remaining subset-I problems

| problem | blocker |
|---|---|
| `Blasi_CellSystems2016`, `Laske_PLOSComputBiol2019` | `observableTransformation = log` (natural log). PyBNF's `lognormal` family is log10; the Gaussian-on-ln-scale family is unimplemented, so import refuses rather than silently mis-scaling. |
| `Fiedler_BMCSystBiol2016` | `observableParameters` scaling factors (`s_pErk_*`) import as free parameters that bind to no model entity (residual of the #495 placeholder-mapping class). |
| `Schwen_PONE2014` | measurement grid has no positive time point; `TimeCourse` rejects it. |

The five blockers filed during the earlier triage — **lanl/PyBNF#492–#496** — have all **landed**,
which is what took importable coverage from 15/23 to 21/23 (19 of which also run).

## Optimizer choice

`job_type` per slug was chosen by **running each candidate**, not by assumption:

- **`gntr`** — the general-objective Fisher/Gauss-Newton trust region (EFIM Hessian fed through
  `trf`'s Coleman–Li core, ADR-0068). This is the method that handles an **estimated noise scale**;
  plain `trf` refuses it (`trf` needs "an exact least-squares residual … with a fixed noise scale").
  It is PyBNF's fides-analogue and the default here. It works on log10/`lognormal` objectives too —
  the EFIM Fisher block for noise scales landed in ADR-0079/0080/0081.
- **`cmaes`** with IPOP restarts (ADR-0070, restart trigger fixed in ADR-0082) — used where the
  gradient path genuinely refuses the problem (Bertozzi, Bruno: a per-condition estimated initial
  condition the gradient path cannot yet route, ADR-0076; Smith), and for the three strongly
  multimodal problems (Borghans, Elowitz, Okuonghae), where a local method from a few starts lands
  in a local basin.

Every shipped conf was verified to start and complete a tiny run.

**These recipes are reasonable defaults, not tuned ones — and that is a real limitation.** A full
`gntr` run of the shipped `SalazarCavazos_MBoC2020` conf (20 starts × 500) converged to `OG = 10.2`,
i.e. *worse* than that problem's own nominal point (`OG = 0.33`): 20 box-sampled starts were not
enough to find the reference basin. Expect to tune `population_size` / `max_iterations`, or to switch
to `cmaes` with IPOP restarts, before treating any ⚪ or 🟢 row as a statement about PyBNF's
optimizers. The three ✅ rows are the only ones where a tuned fit was actually driven to `OG < 1.92`.

## Import + fit pipeline (reproduce)

```bash
# 1. import — parameterScale (#491) AND observableTransformation (#499) are preserved,
#    both of which plain petab.v2.petab1to2 drops. No hand-editing of noise_model is needed.
python -c "from pybnf.petab import petab1to2_preserve_scale, import_job; import tempfile; \
  v2 = petab1to2_preserve_scale('<id>/<id>.yaml', tempfile.mkdtemp()); \
  import_job(v2, 'out', job_type='gntr')"
# 2. set the run recipe: sbml_backend = bngsim (REQUIRED for gradients), wall_time_sim, random_seed
# 3. fit + score
pybnf -c <id>.conf -o
python score.py output          # OG = -log_likelihood - J*
```

`wall_time_sim` caps pathological parameter points, where CVODE would otherwise grind for seconds
before failing; a timed-out simulation is simply a failed one, so the optimum is unchanged. Raise it
if valid simulations on your machine are being marked as failures.

## Per-slug contents

`<id>.conf` (the runnable fit) · the SBML model (verbatim) · `experiment*.exp` (data) ·
`*_measparams.tsv` (per-measurement observable/noise parameter tables, where the problem uses them) ·
`jstar.txt` (the reference J\*) · `nominal_check.json` (the nominal-point evaluation) · `score.py` ·
`README.md`. The three solved slugs additionally ship `best_fit_params.txt`,
`information_criteria.txt`, and `VALIDATION.md` from their fits.
