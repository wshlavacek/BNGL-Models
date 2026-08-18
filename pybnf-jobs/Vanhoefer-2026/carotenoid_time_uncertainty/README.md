# carotenoid_time_uncertainty — measurement-time uncertainty on a real 13-parameter model

Reproduces the **Fig 5/6** experiment of Vanhoefer, Nakonecnij, Binder & Hasenauer, *Efficient
Bayesian inference for ordinary differential equation models from experimental data with uncertain
measurement times* (bioRxiv [2026.05.09.724053](https://doi.org/10.64898/2026.05.09.724053)), on the
carotenoid-cleavage model of *Arabidopsis thaliana* (Bruno et al., *J Exp Bot* 67:5993, 2016).

It is the real-model companion to PyBNF tutorial lesson 49 (a two-point decay toy). Where the lesson
shows the `time_error` marginalization on one rate, this exercises it at the paper's own scale — **13
free parameters, 6 conditions, 77 measurements, 7 state variables** — and, unlike lesson 49, fits it
with the **phase-2 gradient** engine (PyBNF ADR-0113): the marginal-time objective assembles
`dz_k/dθ` by chaining the forward-sensitivity tensor PyBNF already stores, so `job_type = gntr`
(and `lbfgs`) fit it directly.

## The experiment

The published model is fit (θ\*, the reference optimum in `theta_star_source.txt`). Synthetic data
is then generated **from the model at θ\*** with a controlled timing error (`perturb_times.py`): each
datum is drawn at a *latent* time `τ_k = clip(t_k + N(0, σ_t), 0, 200)` but **reported** at `t_k`, so
θ\* is the exact ground truth and the only corruption is the timing. Two fits are compared on that
data:

* **`standard.conf`** — scores each datum at its reported time `t_k` (the classical assumption). It
  is dragged off θ\* as `σ_t` grows.
* **`marginal.conf`** — integrates the latent time out, `z_k = ∫ p(ȳ_k | x(τ,θ)) p(τ | t_k) dτ`
  (the phase-2 gradient fit). It stays at θ\* and **estimates** `σ_t`.

## Archetype & provenance

**B — PEtab/SBML-imported.** `model_carotenoid.xml` is copied **verbatim** from the
Benchmark-Models-PEtab collection (byte-identical, LF-sha256, to the model in the sibling
`Grein-2026-benchmark-subset-I/Bruno_JExpBot2016` slug — see `upstream.json`). Everything else is
ours: the synthetic `sigmaT*/` data and its generator `perturb_times.py`, the two confs, `score.py`,
`make_demo.py`, and the docs. `theta_star_source.txt` is a copy of that slug's `best_fit_params.txt`.

## Files

| file | role |
|---|---|
| `model_carotenoid.xml` · `upstream.json` | the SBML (copied verbatim) + its provenance pin |
| `theta_star_source.txt` | the reference optimum θ\* (13 params, Obj = 32.556) |
| `reference/experiment____model1_data{1..6}.exp` | the published reporting times + per-point σ (`_SD`) + the 6 conditions' structure; the values are replaced by the synthetic generator |
| `perturb_times.py` · `_simulate.conf` | generate the synthetic `sigmaT*/` data from the model at θ\* (`_simulate.conf` is the dense-simulation harness it drives) |
| `sigmaT{0,2,5,10}/` | the generated datasets, one per timing-error level (`σ_t` in minutes; `σ_t = 0` is θ\* + noise, the control) |
| `standard.conf` · `marginal.conf` | the two fits (full-box 20-start `gntr`; `marginal.conf` is the phase-2 gradient fit) |
| `score.py` | the Fig 5/6 gate scorer (parameter MSE vs θ\*, estimated σ_t, the standard-vs-marginal LRT) |
| `make_demo.py` · `demo_results/` | reproduce + the committed θ\*-seeded demonstration fits behind the gate table |

## Free parameters (13)

`init_b10_1`, `init_bcar1`, `init_bcar2`, `init_bcry_1`, `init_ohb10_1`, `init_zea_1` (per-condition
initial amounts), `k5`, `kb1`, `kb2`, `kc1`, `kc2`, `kc4` (rate constants), `szea` (a shared
multiplier), each reaching the model only through a `condition:` parameter reference. The marginal
fit adds one nuisance, `sigma_t__FREE` — the whole dimensionality saving of marginalization: **one**
extra parameter, not one latent time per observation.

## Status

**Gate B (correction, σ_t = 10): PASS.** The standard fit is dragged off θ\* (log₁₀-parameter MSE
0.00059) while the marginal fit stays closer (0.00027, 2.2× better); the marginal **recovers the
timing scale** (`σ_t = 11.3`, injected 10); and the likelihood-ratio test **rejects** the standard
model (2·Δln L = 70, p ≈ 0). **Gate A (no false positive, σ_t = 0): PASS on recovery** — the marginal
matches the standard MSE and drives `σ_t` to its floor. See `VALIDATION.md` for the full table and
the one honest caveat (the fixed-grid LRT at a sub-grid-spacing σ_t).

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"

# 1. (re)generate the synthetic data from the model at theta*
python perturb_times.py                       # writes sigmaT{0,2,5,10}/

# 2. the controlled demonstration behind the gate table (theta*-seeded, ~minutes/level)
python make_demo.py --level 10                # Gate B
python make_demo.py --level 0                 # Gate A

# 3. or the full from-scratch reproduction (20-start, expensive)
pybnf -c standard.conf                        # then edit sigmaT5 -> sigmaT10 etc.
pybnf -c marginal.conf
python score.py --standard output/standard --marginal output/marginal --injected 5
```
