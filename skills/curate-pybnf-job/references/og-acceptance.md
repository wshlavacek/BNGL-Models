# The acceptance bar: reference objective J\* and optimality gap OG

The completion test for a `pybnf-jobs/` slug. A job that runs is not a job that works; this
document defines the number that tells them apart, where that number's anchor comes from, and
what you are allowed to claim once you have it.

## Contents
1. The identity — what OG is
2. Why a non-likelihood objective cannot be scored
3. Where J\* comes from: three provenance tiers
4. The threshold and the status vocabulary
5. The committed provenance set
6. Reading a large OG
7. Worked cases in the corpus

---

## 1. The identity

```
J_paper = -log_likelihood            # from Results/information_criteria.txt
OG      = J_paper - J*               # optimality gap
solved  = OG < 1.92                  # chi^2, alpha = 0.05, 1 dof
```

PyBNF **minimizes a reduced objective**: it drops the parameter-independent per-point constants —
`½log(2π)`, and for a log-transformed observable the change-of-variables Jacobian
`Σ log(y_obs·ln10)` — because they do not move the argmin. It then **reports the full normalized
log-likelihood** at the best fit in `Results/information_criteria.txt`, restoring every dropped
constant (it matches `scipy.stats.norm.logpdf` / `lognorm.logpdf`; ADR-0056,
`algorithms/base.py:1397-1432`). So `−log_likelihood` is the objective on an **absolute,
transform-agnostic scale** — comparable across models, across observable transformations, and
against a value another tool reported. That is what makes OG meaningful rather than a difference
of two arbitrary numbers.

The corollary matters: **do not compare the value in `sorted_params_final.txt` to a published
objective.** That is the reduced objective. `scripts/score.py` prints the restored constant
explicitly so the arithmetic is auditable. For a **fixed**-σ Gaussian the restoration is exactly

```
J_paper  =  reduced  +  sum log(sigma_i)  +  (N/2) log(2*pi)
```

with `Σ log σᵢ` the Gaussian normalizer (`noise/gaussian.py:122`) — zero when σ ≡ 1
(edition-2 `sos`), and generally *negative* for `chi_sq` when the data's `_SD` are below 1.
Verified end to end on `Korwek-2023/nfkb_tnfa`: its recorded residual
`548.245739 − 703.556619 = −155.310880` equals `Σ log σ` summed directly over the 111 `_SD`
entries in its `.exp` files to six decimals, and `N = 111` matches that count. A log-transformed
observable adds its change-of-variables Jacobian on top. `make_jstar.py` records the residual as
`sum_log_sigma`, so a number that does *not* decompose this way is visible immediately.

`1.92` is `½·χ²₁(0.95)` — the likelihood-ratio half-width of a 95% confidence interval on one
parameter. Two fits within 1.92 NLL of each other are statistically indistinguishable at that
level, which is why it is the standard "solved" bar in optimizer benchmarking and why this corpus
adopts it rather than inventing a tolerance.

---

## 2. Which objectives can be scored

What decides scoreability is the **resolved** objective object, not the token in the conf — and
the two editions resolve `sos` to different classes:

| conf | resolved class | value | `supports_pointwise_log_likelihood` | scoreable |
|---|---|---|---|---|
| **edition 2** `objective = sos` | `LikelihoodObjective` + `Gaussian(sigma=1)` | `½Σr²` — the Gaussian reduced term already carries the ½ (`noise/gaussian.py:60`) | **True** | **yes** |
| edition 1 `objfunc = sos` | `SumOfSquaresObjective` (a plain `SummationObjective`) | `Σr²` — no ½ (`objective.py:1628`) | False | no |

So **an edition-2 `sos` job is a Gaussian likelihood with a unit noise scale and needs no conf
change to be scored.** Verified directly: `McMillan-2021/tnfr1_apo` declares `objective = sos`,
resolves to `LikelihoodObjective`/`Gaussian`, and yields

```
reduced objective            = 0.07335653
J_paper = -log_likelihood    = 25.80363546        n = 28
reduced + (N/2)log(2pi)      = 25.80363546        <- exact
```

That identity — `J_paper == reduced + (N/2)log(2π)` — holds exactly for any **fixed**-σ Gaussian,
and is the arithmetic to re-run when a number looks wrong. (With an *estimated* σ or a
log-transformed observable, extra restored constants enter and the identity no longer closes;
`make_jstar.py` records whether it held.)

`likelihood_information_criteria` returns `None` — deliberately, "rather than a misleading
number" — only when the resolved objective is not a per-point likelihood (`objective.py:73-96`);
`_emit_information_criteria` then writes nothing (`algorithms/base.py:1406-1408`). In this corpus
that is just the two **legacy edition-1** slugs, `Miller-2026/mek_isoform_{de,amcmc}_legacy`.
Their edition-2 twins score fine, which is the point of keeping both.

**Do not read a missing `information_criteria.txt` as "unscoreable".** PyBNF writes it at the end
of a *fit*, so a job that has never been run has none regardless of its objective. In the Grein
subset, `Smith_BMCSystBiol2013` (`sos`) and `Weber_BMC2015` (`chi_sq`) both lack the file for the
same reason: both are ⚪ setup-only. The objective is not what distinguishes them.

**Choosing the noise model** is therefore a modeling decision, not a scoreability workaround:

- data carries per-point `_SD` → `objective = chi_sq` (σ read per point, normalizer dropped);
- the paper's fit is unweighted least squares → edition-2 `objective = sos` (σ ≡ 1) is faithful
  *and* scoreable. Write `noise_model = normal, sigma = fix_at 1` only when you want the unit
  scale to be explicit — it resolves to the same thing;
- the noise scale is unknown and should be estimated → `noise_model = normal, sigma = fit
  sd_<obs>`, and use `job_type = gntr`, the one method that handles an estimated noise scale on
  the gradient path.

---

## 3. Where J\* comes from: three provenance tiers

Declare the tier **before fitting**, in the conf banner and in `nominal_check.json`. The tier is
what the status badge means; a J\* without a stated provenance is a number you cannot defend.

### T1 — published or benchmark objective
The paper, its SI, or a benchmark collection reports the optimum's objective value on a scale you
can match. Transcribe it into `jstar.txt` and cite the exact source (table, file, row).

*Example:* the Grein subset takes J\* from `best_fx_marvin.csv` in
`ICB-DCM/optimizer-benchmark-2026-suppl-code-and-data` — the minimum over 33 optimizers × many
runs of the Eq. 6 Gaussian NLL. That is the strongest anchor available: an OG against it is a
statement about PyBNF versus a published field.

**Match the scale before you subtract.** If the paper reports `½Σr²`, or SSR, or a per-point mean,
convert explicitly and show the conversion in `nominal_check.json`. An unconverted J\* produces a
confident, wrong OG.

### T2 — published-parameter objective (the common case)
The paper reports a best-fit **parameter table** but no objective value. Set the model to those
parameters, evaluate PyBNF's objective there, and take that as J\*. This is the ordinary anchor
for a hand-built job from a biology paper, and it is exactly the point that
`validate-pybnf-job`'s Gate 3a reproduces.

Evaluate at the published point — a zero-budget run, or `job_type = check`, or `make_reproduction.py`
extended to emit the objective — and record both the parameter values and their source. OG then
answers a sharp question: **can PyBNF find a fit at least as good as the authors' own?** A
*negative* OG is a real and publishable result (PyBNF found a better optimum than the paper
reported); it is not an error, and it should be stated plainly with the parameter comparison that
explains it.

**Three cases where T2 is structurally undefined**, found by backfilling the corpus. None is a
tooling failure; each needs a different anchor:

- **An estimated noise nuisance.** A free parameter that lives only in the conf — `noise_model =
  neg_bin, dispersion = fit r_disp` (`Lin-2021/nyc`), or `chi_sq_dynamic`'s `sigma__FREE`
  (`Miller-2026/mek_isoform_amcmc`) — has no model declaration, because a published parameter
  table fixes the *model*, not the noise scale. There is no "published point" to evaluate at.
  Use T3.
- **The published point is outside the job's own search space.** `Zhang-2023/tumor_growth`
  declares `kkill 0` in the model and `loguniform_var = kkill 0.001 1` in the conf: zero is
  unrepresentable on a log scale, so the published point cannot even be assigned. (Here the
  paper is the problem — it "reports no kkill at all, and File S4 sets it to 0, which cannot
  shrink a tumor.") Either re-parameterize or use T3, and say which.
- **A constraint-only BPSL job.** `Kirsch-2020/phosphoswitch_bpsl` and
  `Rukhlenko-2022/cstar_skmel133_bpsl` score zero data points, so there is no data-fit
  likelihood and no OG at all. Verify them with `job_type = check` satisfaction
  (`bpsl-constraints.md`) — that *is* their acceptance bar.

### T3 — corpus best-known objective (regression anchor)
No published anchor exists: the paper reports neither an objective nor a usable parameter table.
J\* is then the best `−lnL` this corpus has ever reached on this job, recorded with the run that
produced it.

This is **not an optimality claim** and must never be badged "solved". Its value is as a
**regression detector**: when a bngsim or PyBNF change moves this OG, something changed in the
forward model, the objective, or the optimizer. Badge it `🔁 regression-anchored`, and update
`jstar.txt` (with a dated note) when a better optimum is genuinely found.

---

## 4. The threshold and the status vocabulary

| badge | meaning | requires |
|---|---|---|
| ✅ **solved** | a real fit was run and reached `OG < 1.92` | T1 or T2 J\*, a completed fit |
| 🟢 **objective validated** | no fit was run, but the model's nominal/published point *is* the optimum, and PyBNF's objective there lands within the threshold of J\* | T1 or T2 J\*, an evaluation at nominal |
| 🔁 **regression-anchored** | OG is measured against the corpus's own best, for drift detection | T3 J\* |
| ⚪ **setup only** | imports, simulates and scores correctly, but nothing about optimality is claimed | a scoreable objective |

Mark any OG measured at a nominal point rather than from a fit with `†` in every table it appears
in. The distinction is the whole point: 🟢 validates the import chain end to end — SBML/BNGL model
→ simulation → observable formulas → noise model → objective — and says **nothing** about PyBNF's
optimizer. Conflating the two is the easiest way to make this corpus dishonest.

⚪ is a legitimate deliverable. A ready-to-run, correctly-scoring job that nobody has tuned is
useful. It just must not be described as reproducing the paper's fit.

---

## 5. The committed provenance set

Four files, beside the conf, so that any reader can re-derive the score without re-running:

| file | contents |
|---|---|
| `jstar.txt` | the bare number, one line. A comment line may follow with the source. |
| `best_fit_params.txt` | the parameter set that produced the reported OG (a copy of the run's `sorted_params_final.txt` best row, or the published point for a 🟢) |
| `information_criteria.txt` | `k`, `n`, `log_likelihood`, AIC/BIC/AICc from that same fit |
| `nominal_check.json` | the audit record — see below |

`nominal_check.json` is the field where you are honest in prose:

```json
{
  "problem": "Bertozzi",
  "slug": "Bertozzi_PNAS2020",
  "jstar": 158.86426270904192,
  "jstar_tier": "T1",
  "jstar_source": "best_fx_marvin.csv (ICB-DCM suppl, commit 4d20850)",
  "J_paper": 158.86426780179107,
  "reduced_objective": 138.64762007128826,
  "n_scored": 22,
  "k": 8,
  "OG_nominal": 5.092749148616349e-06,
  "optimizer": "gntr",
  "status": "solved",
  "interpretation": "PEtab nominalValue reproduces the reference optimum within the solved threshold (OG < 1.92); this validates PyBNF's objective against the paper's Eq. 6 NLL. Recomputed 2026-08-02 after lanl/PyBNF#531: the previous value (OG 1.79e+11) was evaluated against a forward model whose derived parameter beta_N = R0_*gamma_/N_ never tracked its condition-set dependencies, so R0_ was inert and the trajectory was wrong by six orders of magnitude."
}
```

That `interpretation` field is the most valuable line in the file. When an OG moves, write down
**why** — the issue number, the mechanism, the date. It is how the corpus turns a benchmark number
into a bug report.

Run the scorer:

```bash
$PY $SKILL/scripts/score.py $JOB            # shipped provenance
$PY $SKILL/scripts/score.py $JOB output     # a fresh run's output/Results/
```

---

## 6. Reading a large OG

A large OG is a *diagnosis prompt*, not a verdict on PyBNF. Work down this list before reporting:

1. **Scale mismatch.** Is J\* on the same scale as `−lnL`? Reduced vs. full, SSR vs. ½SSR,
   log10 vs. natural log. This is the most common cause of an absurd OG.
2. **Wrong forward model.** Does the model at the *published* parameters reproduce the paper's
   figure? If not, the defect is upstream of the optimizer — go back to the model or the data.
   Bertozzi's OG was `1.79e+11` for exactly this reason until PyBNF#531 landed.
3. **Untuned budget.** The 20 × 500 default is not tuned. In this corpus,
   `Laske_PLOSComputBiol2019` reaches only `OG = 6.76` at 20 × 500 and reaches the reference
   optimum at 100 × 1000; `SalazarCavazos_MBoC2020` converged at 20 starts to `OG = 10.2`,
   *worse* than its own nominal point (`0.33`). Raise the budget or switch method before
   concluding anything.
4. **Wrong method for the landscape.** Strongly multimodal problems want `cmaes` with IPOP
   restarts; an estimated noise scale wants `gntr` (plain `trf` refuses it); a model with discrete
   events refuses the gradient path entirely.
5. **Genuine optimizer gap.** Only after 1–4 are excluded. Then it is a finding — file it.

---

## 7. Worked cases in the corpus

- **`Grein-2026-benchmark-subset-I/`** — the reference implementation of everything above: T1
  J\* for all 23 slugs, per-slug `score.py`, and a coverage table carrying J\*, scale, k, n,
  optimizer, OG and badge. Nine ✅, four 🟢, ten ⚪ — a distribution that is informative precisely
  because it is not all green.
- **`Bertozzi_PNAS2020`** — how an OG becomes a bug report: `1.79e+11` → `5.4e−06` across
  PyBNF#530/#531, with the mechanism recorded in `nominal_check.json`.
- **`Blasi_CellSystems2016`** — the one slug whose nominal point is *not* its optimum despite
  carrying the published parameters, because PEtab's `nominalValue` for its noise scale is a
  placeholder. A reminder that "nominal" and "published optimum" are not synonyms; check, don't
  assume.
- **`Smith_BMCSystBiol2013` / `Oliveira_NatCommun2021`** — `objective = sos` under edition 2, so
  a Gaussian likelihood with σ ≡ 1 and perfectly scoreable; they carry no
  `information_criteria.txt` only because no fit has been run (§2). Both are faithful to their
  upstream problems, which specify σ ≡ 1: Smith's `observables_*.tsv` has `noiseFormula = 1.0`,
  and Oliveira's noise parameters have `nominalValue = 1` with `estimate = 0`.
- **`Salazar-Cavazos-2019/egfr_simpull` vs `Grein-.../SalazarCavazos_MBoC2020`** — the same
  published problem reached two independent ways: a hand-built BNGL reconstruction and a
  PEtab/SBML import. Both give `k = 6`, `n = 18`, and J\* within **0.39 NLL** of each other
  (366.47 vs 366.86) — inside the 1.92 threshold. That agreement cross-validates the BNGL
  reconstruction, the SBML import path, and the objective assembly at once, and is the strongest
  evidence in the corpus that these numbers mean what they claim.
