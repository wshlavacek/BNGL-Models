# Verifying a stochastic model

How to establish the two `curate-model` verification levels when the model's defining result is
*noise* — a stationary distribution, a switching statistic, a spectrum, a burst — rather than a
trajectory. A deterministic model is verified by integrating the same equations twice; a
stochastic one is not, because the curated run and the independent check never produce the same
numbers. What you compare instead, and how tightly, is what this document fixes.

Roughly a sixth of the collection is in this class. The protocol below is generalized from
`two_state_gene_expression_noise_munsky2012`, `three_stage_stochastic_gene_expression_shahrezaei2008`,
`noise_induced_bistable_futile_cycle_samoilov2005`, `demographic_noise_predator_prey_cycles_mckane2005`,
and `bursty_autoregulated_gene_expression_lin2016`.

All five run a *finite* reaction network, which is what makes an FSP or an independent Gillespie
available as the exact check. A model simulated network-free has neither, and establishes level 1
a different way — see `network-free-verification.md`.

## Contents
1. The independence rule
2. Choosing the independent check
3. Run protocol: seeding, burn-in, sampling, replicates
4. Statistics and tolerances
5. Comparing against a single published realization
6. Spectra
7. Truncation validation for CME/FSP
8. Worked examples

---

## 1. The independence rule

**The independent implementation is transcribed from the paper's equations or reaction list, and
never reads the BioNetGen-generated network.** Say so explicitly in the notebook — every model in
the class does:

> "assembled directly from four transitions … It does not parse or reuse the BNGL-generated
> network." (munsky2012)
> "transcribed directly into a pure-Python direct Gillespie simulator. It shares no code with
> BioNetGen and never reads the generated network." (samoilov2005)

This is the whole point of level 1. A check that starts from the `.net` file verifies the
integrator, not the model, and will agree with a mis-transcribed rule just as happily as a correct
one. If the paper gives reactions, transcribe reactions; if it gives a master equation, transcribe
the generator; if it gives ODEs, transcribe the ODEs.

Transcribe the propensities too, and state any that are not obvious by inspection. A homodimeric
step `N + N -> E+ + N` can be read as `k*N*(N-1)` or half that; Samoilov's convention is unstated
in the paper, and the two readings differ by a factor of two in the driver noise (§8, and
`when-the-paper-is-wrong.md`).

## 2. Choosing the independent check

Use the strongest one the model admits. They are not exclusive — the best notebooks run two or
three.

| check | when | what it proves |
|---|---|---|
| **SciPy ODE** of the paper's mean-field equations | the model has a deterministic limit | rules, rate laws and ICs are right, independent of noise |
| **Pure-Python Gillespie** from the paper's reaction list | a small, explicit reaction set | the stochastic *mechanism*, exactly — the strongest check available |
| **Finite-state CME / FSP** — assemble the generator, solve `Qᵀp = 0` | a low-dimensional state space (one or two counted species) | the exact stationary distribution, with no sampling error at all |
| **The paper's own analytic result** | the paper publishes a closed form or moment formula | that you and the authors are solving the same problem |

The deterministic check is cheap and catches most transcription errors — run it even when the
paper's result is purely stochastic. `samoilov2005` and `mckane2005` both do, and both use the
deterministic agreement (max relative error < 1e-6) to isolate the stochastic comparison to
genuinely stochastic causes.

Prefer an **exact** independent calculation over a sampled one wherever the state space allows:
an FSP has no sampling error, so a disagreement is unambiguously the model's. Reserve the
looser sampled-vs-sampled comparison for models too large to project.

Where the paper hands you a closed form, evaluate it — `shahrezaei2008` implements Eq. 18's
hypergeometric form, `munsky2012` computes the exact steady-state mean and Fano factor from Eq. 2,
`mckane2005` evaluates the linear-noise spectrum of Eq. 7. This turns a plausible-looking curve
into a number the authors would recognize.

**Rescale to simplify when it is free.** `munsky2012` divides every rate by `gamma_R` so
`gamma_R = 1`; the stationary distribution is unchanged and the generator gets much easier to read.
Record that you did it.

## 3. Run protocol: seeding, burn-in, sampling, replicates

- **Fix the seed and commit the output.** The `.bngl` runs a seeded SSA; `reference/` holds its
  `.gdat`/`.cdat`. The notebook re-runs BioNetGen **in a temporary directory** and checks the fresh
  arrays against the committed ones *exactly*. Copying the model into a temp dir is not fussiness:
  it is what stops the check silently reading a stale committed output instead of producing one.
- **Discard the transient, and say how much in physical units.** "Samples before 6,000 s (20
  transcript lifetimes)" (munsky2012); "Samples through 10 million seconds are burn-in; the
  remaining 10,000 uniformly spaced samples estimate the stationary distribution" (shahrezaei2008,
  of a 210-million-second run). A burn-in quoted in *lifetimes* survives a reparameterization; one
  quoted in seconds does not.
- **One long trajectory or many replicates?** A long ergodic trajectory estimates a *stationary*
  quantity; an ensemble estimates a *distribution over realizations* and is required for anything
  time-resolved (a spectrum, a transient, a first-passage statistic). `mckane2005` runs 128
  replicates, `samoilov2005` 64.
- **Justify a replicate count that differs from the paper's.** `mckane2005`: "The paper used 500
  replicates; 128 gives a stable spectrum while keeping verification comfortably below one minute
  on a typical workstation." State the trade, don't just pick a number.

## 4. Statistics and tolerances

Compare *statistics of the distribution*, never trajectories. Pick metrics the figure type
supports, and justify each threshold from a named source of error — sampling, digitization,
figure resolution, or the precision of the printed number.

Thresholds in current use, as calibration:

| comparison | metric | threshold | why |
|---|---|---|---|
| BNG SSA vs exact FSP, distributions | total-variation distance (½·L1) | ≤ 0.03 | finite-trajectory sampling error |
| BNG SSA vs exact moments | mean relative error | ≤ 3% | as above |
| BNG SSA vs exact moments | Fano-factor relative error | ≤ 5% | a variance converges slower than a mean |
| BNG SSA vs independent Gillespie | stationary means | within **3 combined standard errors** | both sides are sampled |
| theory vs digitized curve | RMSE / panel height | ≤ 0.06 | line thickness + digitization |
| sampled BNG vs digitized curve | RMSE / panel height | ≤ 0.08 | the same, plus sampling |
| vs a number printed to 4–5 s.f. | relative error | ≤ 5e-4 | "the resolution of the report itself" |

Two habits worth copying:

- **The sampled side gets the looser bound.** munsky2012 accepts 0.06 for the FSP and 0.08 for the
  SSA against the same digitized curve, because only one of them carries sampling error.
- **When both sides are sampled, use standard errors, not a fixed percentage.** "every stationary
  mean within 3 combined standard errors" (samoilov2005) scales correctly with run length; a flat
  5% does not.

Normalize a distribution comparison by something physical (the panel height, the mean), so the
number means the same thing across models.

## 5. Comparing against a single published realization

Papers often plot **one** trajectory. You cannot match it, and you must not try to.

Compare summary statistics of the published trace against the **distribution of the same statistic
across your ensemble**, and accept if the published value falls inside it. `samoilov2005` does
exactly this for Fig. 3A — a single trace with only about ten switching events — comparing
digitized time averages and low-state occupancy against the spread of the 64-replicate ensemble.

State the ensemble spread you are testing against. "Inside the ensemble range" is a real claim;
"looks similar" is not.

## 6. Spectra

Compare **normalized shape and peak frequency**, not absolute power: "The absolute periodogram
normalization depends on Fourier conventions" (mckane2005). Fit only what you must — an amplitude
normalization, and a frequency scale if the paper's axis disagrees — and report every fitted
factor as a finding rather than absorbing it silently. See `when-the-paper-is-wrong.md` §"a uniform
transformation".

## 7. Truncation validation for CME/FSP

A projection is only exact if the truncation is. Two checks, both cheap, both mandatory:

- **Boundary mass is negligible.** munsky2012 projects 0–250 transcripts and requires the
  probability at the upper boundary to be negligible; shahrezaei2008 uses `m = 0..12`, `n = 0..220`.
- **The stationarity residual is small.** Confirm `Qᵀp ≈ 0` for the solved `p`.

Report both numbers. A silently-truncated FSP produces a smooth, plausible, wrong distribution —
exactly the failure mode a tolerance check will not catch, because it shifts *both* the theory and
your confidence in it.

## 8. Worked examples

- **`two_state_gene_expression_noise_munsky2012`** — the template. Seeded SSA vs. an independently
  assembled FSP per promoter, exact moments from the paper's Eq. 2, and digitized Fig. 2B; three
  metrics with separate thresholds for theory and sample.
- **`three_stage_stochastic_gene_expression_shahrezaei2008`** — finite-state CME over
  `(promoter, mRNA, protein)`, plus the paper's analytic large-γ distribution (Eq. 18); a
  210-million-second run with a 10-million-second burn-in.
- **`noise_induced_bistable_futile_cycle_samoilov2005`** — an independent pure-Python Gillespie
  from the paper's ten elementary reactions; combined-standard-error acceptance; a single
  published realization compared against a 64-replicate ensemble; and an unstated propensity
  convention bracketed by two model files.
- **`demographic_noise_predator_prey_cycles_mckane2005`** — 128-replicate ensemble vs. the
  linear-noise spectrum of Eq. 7; normalized shape and peak frequency; a documented factor-of-two
  time base between the published figures and the model.
- **`bursty_autoregulated_gene_expression_lin2016`** — ODE relaxation carried into a seeded SSA in
  one actions block; independent SciPy ODE at stricter tolerance than the BNG run.
