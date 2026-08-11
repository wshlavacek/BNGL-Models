# VALIDATION — Smith_BMCSystBiol2013

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This is the **last slug in subset I** to be fitted, and
the only one that had to clear two independent blockers first — a gradient path that refused it, and
an SBML units defect that made its objective meaningless. The model is

> Smith GR, Shanley DP. **"Computational modelling of the regulation of Insulin signalling by
> oxidative stress."** *BMC Systems Biology* **7**, 41 (2013).
> <https://doi.org/10.1186/1752-0509-7-41> (open access, PMC3668293)

> **Confidence: 84 / 100.** SOLVED with `OG = 0.501586` from a from-scratch 100-start multi-start, no
> seeding, **13 h 11 m** and 1,130,847 simulations. It beats its own nominal point by a factor of
> **1.7 million** (`OG_nominal = 867276.27`), which is the largest such margin in the collection.
> Deductions: the model is imported, not re-derived; the run is a single seed; **5 of 25 parameters
> rest on a box bound**; and — the reason this is not higher — with `σ ≡ 1` the objective is
> dominated by three of the paper's five fitted panels, so a good `OG` here is much narrower evidence
> than the point count suggests (see "What the objective actually measures").

## Gate A — objective fidelity

Linear observables with **σ fixed at 1** upstream (`objective = sos`, no noise model), so the
restored constant is the plain `(N/2)log(2π)` with no Jacobian. From the fit's own
`information_criteria.txt`:

| term | value |
|---|---|
| PyBNF reduced objective | 20865.691640 |
| restored constant (`−lnL − J_reduced`) | 56.974190 |
| `(N/2)·log(2π)`, N = 62 | 56.9741890587 |
| `J_paper = −log_likelihood` | 20922.665830 |

The identity holds to **1.4e−06 absolute**, and that residual is a recording precision, not a
disagreement: `information_criteria.txt` stores `log_likelihood` to five decimals
(`−20922.66583`), while `J_reduced + (N/2)log(2π)` evaluates to `20922.6658285607`. Restoring the
constant from the unrounded reduced objective reproduces the exact `(N/2)log(2π)`. This is the same
identity `Oliveira_NatCommun2021` confirms as the corpus's other unit-σ `sos` slug, at
`110.272624 = (120/2)log(2π)`.

### This slug DOES have a §2c independent oracle — earlier revisions said it did not

The collection README listed `Smith_BMCSystBiol2013` among the slugs whose upstream rows "could not
be joined", and left its `obj ✓` column blank. That was an artefact of the checker's join key, not a
property of the data.

`measurementData` and `simulatedData` join **62 of 62, one-to-one** on the PEtab identity key
`(observableId, simulationConditionId, time)`. The checker also included `datasetId`, which PEtab
defines as a *visualization grouping label* and which the two tables are free to disagree on — and
here they do, on exactly 13 rows: all ten of figure 2C and all three of figure 2D are tagged
`fig2C` / `fig2c` / `fig2D` in `measurementData` and `fig2A` in `simulatedData`. 62 − 13 = 49, which
is the "joined 49 of 62" this slug was written off for.

With the identity key, recomputing the Eq. 6 NLL at the nominal point straight from the upstream
tables — **no PyBNF in the loop** — gives:

| | reduced objective at the PEtab nominal point |
|---|---|
| oracle, from upstream `simulatedData` | 888141.4613667557 |
| PyBNF, same point, same build | 888141.4641019923 |
| difference | 2.74e−03 absolute, **3.08e−09 relative** |

This is **integrator-tolerance agreement, not a digit-for-digit match**, and the distinction should
be kept. Upstream's `simulation` column carries full double precision (16–17 significant digits), so
the residual is not print rounding — it is two different integrators on the same trajectory. The
absolute gap looks large only because this evaluation point sits at `OG_nominal = 8.7e+05`; against
the 1.92 threshold it is 700× smaller. That is weaker evidence than `Giordano`'s exact
reproduction, and stronger than the internal identity alone.

Two fixes to `tools/sigma_profile.py` were needed to make the repo's own tooling agree with the
above, both verified to change **only this slug's** line across all 23:

- `datasetId` dropped from the join key, for the reason above.
- the `('formula', observableId)` σ branch — which `sigma_key` already routed to but never resolved —
  now looks the observable's `noiseFormula` up. Smith's is the literal `1.0` on all **nine**
  observables. (The collection README said "seven"; there are nine.)

σ-profiling itself remains a **no-op** here, correctly: with no estimated σ there is nothing to
profile, and the tool now says so rather than erroring.

**Verdict: PASS (internal identity + independent oracle to 3.1e−09 relative).**

## Gate B — the gradient is the one the objective implies

`tools/fd_check.py` compares the assembled gradient against central differences of the same objective
in the same sampling space, at the **PEtab nominal point**. That choice matters: most slugs here have
`OG_nominal ≈ 0`, so their nominal point *is* the optimum, the true gradient vanishes and both sides
are noise — every such slug flags red for nothing. Smith is the opposite case. At `OG_nominal = 8.7e+05`
the gradient has real magnitude (columns run to `3.6e+06`), which makes this the rare slug where the
nominal point is the *right* place for the test rather than the worst.

Re-measured on the build that produced the fit (bngsim `dff901e`, PyBNF `095a5a14`), 45 m 46 s
single-threaded on an otherwise idle machine:

| | worst relative error |
|---|---|
| all 25 parameters, `h = 3e-04` (default), nominal point | **8.30e−05** (`kminus1`) |

`point-dependent routings: 0 / 35`. **No column is structurally zero** — the two smallest are
`sc_PTP` at `−6.74841e−04` (against a central difference of `−6.74821e−04`, agreeing to 5.7e−12) and
`sc_FOXO1` at `−0.227125` — both small because their observables carry almost none of the objective
(see the leverage table below), not because the gradient path missed them. Sixteen of the 25 columns
agree to better than 1e−05.

This reproduces the figure recorded before the build churn, **8.30e−05, to three significant
figures** — so the local bngsim working copy and the wheel that preceded it assemble the same
gradient here, and the recorded Gate B result carries over rather than being superseded.

**Verdict: PASS.**

## Gate C — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (**100 starts × 1000 iterations**, `random_seed = 1`,
`sbml_backend = bngsim`, no seeding) converges to `J_pybnf = 20865.691640` ⇒
`J_paper = 20922.665830` ⇒ **`OG = 0.501586` < 1.92**, against `J* = 20922.1642440`.

| | |
|---|---|
| wall clock | **13 h 11 m 32 s** |
| simulations | 1,130,847 |
| starts retired on `step is negligible` | 88 |
| starts retired on `reached max_iterations` | 12 |
| failed simulations | 11 |
| CVODES `mxstep` warnings | 1,234 |
| `wall_time_sim` cut-offs | 385 |
| gradient decline / difference-quotient fallbacks | **0** |

**The fit beats its own nominal point by 1.7 million-fold**, which is the claim §1 asks for and by
the widest margin in the corpus: `OG_nominal = 867276.27` against `OG = 0.50`. Unlike the twelve
slugs whose PEtab `nominalValue` *is* the published optimum, this problem could not have been
"converted" by holding the published point — that point is nowhere near the reference basin.

### The descent, and the local optimum it had to escape

The run spent its first three hours in a region no better than 24,000 and its last five hours flat:

| elapsed | best reduced | `OG` |
|---:|---:|---:|
| 0 h 17 m | 2,305,463 | 2.28e+06 |
| 1 h 31 m | 24,742.5 | 3,877.3 |
| 4 h 01 m | 21,174.8 | 309.6 |
| 5 h 01 m | 20,867.4 | 2.26 |
| 6 h 32 m | 20,866.5 | **1.34** ← crosses the threshold |
| 13 h 11 m | 20,865.7 | **0.50** |

**There is a genuine local optimum at ≈24,636 and it attracts independent starts.** Starts 49, 60 and
51 retired there on `step is negligible` at 24,635.96, 24,636.28 and 24,636.68 — agreement to five
significant figures from three unrelated box samples. Any single-start gradient method on this
problem has a real chance of stopping 3,770 NLL units short, which is 1,960× the solved threshold.
This is the argument for the 100-start budget on this slug, and it is worth recording because
`Smith` was the corpus's most expensive run.

**Verdict: PASS (SOLVED).**

### A note on reading `step is negligible` timing

Four starts retired within twelve minutes with printed step norms in a ±1% band
(1.71041e−07 … 1.74162e−07), which is *tighter* than the ±4% band §2 of the kickoff flags as an
infrastructure signature. It is not one, and the reason generalizes: the number in that message is
**not a measured step**, it is the threshold

```
'step is negligible (‖δ‖ ≤ %g)' % (self.step_tol * (point_norm + self.step_tol))
```

— a function of the parameter-vector norm alone (`trf.py:717`). Starts converging to similar points
necessarily print similar values. The diagnostic that actually identifies the infrastructure failure
is the one that accompanied it: `dlopen` failures and `inf` objectives. This run logged **0 of both**,
with the codegen cache intact and its event count pinned at its startup value for all 13 hours.

### Five of twenty-five parameters rest on a box bound

| parameter | value | bound |
|---|---:|---|
| `k32f` | 0.06 | upper |
| `k42r` | 5e−05 | upper |
| `k8` | 0.00026 | upper |
| `kminus7a` | 0.00875 | upper |
| `sc_SOD2` | 0.0001 | lower |

Five of twenty-five is high — the same fraction that `Schwen_PONE2015`'s `VALIDATION.md` calls high
at five of thirty — and the optimum should be read as **constrained rather than interior**. Any
profile-likelihood analysis of this problem should widen those boxes first. `sc_SOD2` on its floor is
the benign one: it multiplies `cytoplasm_SOD2` in the observable, so only the product is identified
and the split between them is degenerate.

## What the objective actually measures — the caveat this slug needs

§2f asks for a comparison against the paper's own figures rather than against `J*` alone. Doing that
here surfaces the most important qualification on this row.

The 62 points are spread evenly across nine observables — no observable holds more than 17.7% of the
count, unlike `Schwen`'s 88% concentration. **But with `σ ≡ 1` the leverage follows data magnitude,
not point count**, and the two do not agree. Per-observable share of the SSE **at the fitted point**:

| observable | n | n share | y range | SSE share |
|---|---:|---:|---|---:|
| `PI3K_activity__2B` | 11 | 17.7% | 10 … 515 | **29.40%** |
| `Glucose_uptake__240__3B` | 4 | 6.5% | 220 … 600 | **23.80%** |
| `IRSYp__2C` | 10 | 16.1% | 0 … 220 | **22.38%** |
| `Glucose_uptake__120__3B` | 4 | 6.5% | 200 … 670 | 9.82% |
| `Cell_Bound_Ins__2B` | 9 | 14.5% | 0 … 142 | 8.28% |
| `Glucose_uptake__2B` | 11 | 17.7% | 45 … 230 | 6.29% |
| `MnSOD_fold_induction__3C` | 5 | 8.1% | 2 … 10.5 | 0.03% |
| `FOXO4__3C` | 5 | 8.1% | 0.3 … 1 | 0.00% |
| `PTP_activ__2D` | 3 | 4.8% | 0.05 … 0.16 | 0.00% |

Grouped by the paper's panels:

| panel | n | SSE share |
|---|---:|---:|
| Figure 2B (dose–response) | 31 | 43.96% |
| Figure 3B (glucose uptake ± H₂O₂) | 8 | 33.62% |
| Figure 2C (IRS tyrosine phosphorylation) | 10 | 22.38% |
| Figure 3C (SOD2 / FOXO4 vs stress) | 10 | 0.04% |
| Figure 2D (PTP1B oxidation) | 3 | 0.00% |

**Figures 2D and 3C together carry 0.04% of the objective — 13 of 62 points that are effectively
invisible to it.** `PTP_activ__2D` ranges over 0.05 … 0.16, so even a 100% relative error there moves
the objective by less than the fourth decimal. A good `OG` on this problem is evidence about
Figures 2B, 3B and 2C, and almost nothing else. It should never be cited as evidence that the model
reproduces the paper's oxidative-stress readouts.

### Against the figures themselves

**Figure 3B is the paper's headline result and the fit reproduces half of it.** The paper's claim is
that insulin *plus* oxidative stress yields *intermediate* glucose uptake — less than insulin alone —
because stress-kinase serine phosphorylation of IRS reduces its availability. The fit gets that
squarely:

| condition | t | data | model |
|---|---:|---:|---:|
| basal | 120 | 200 | 208.5 |
| H₂O₂ alone | 120 | 270 | **208.5** |
| insulin alone | 120 | 670 | 685.6 |
| insulin + H₂O₂ | 120 | 360 | 360.9 |

Insulin alone 685.6 versus insulin + H₂O₂ 360.9 is the published non-monotonicity, and the combined
condition lands within 1 unit of the measurement. **But the fit predicts H₂O₂ alone to be
indistinguishable from basal** — 208.5 against 208.5, identical to the six digits printed, and the
same at t = 240 (293.9 for both) — where the data show a clear rise (270 and 350). The optimizer has
effectively switched off the stress pathway's effect on *basal* uptake while keeping its effect on
*insulin-stimulated* uptake. That is a real disagreement with Figure 3B, and it is invisible in `OG`
because those two points contribute a few percent at most.

Elsewhere: Figure 2B's dose–response is tracked across five decades of insulin, with the one clear
miss at the lowest dose (`PI3K_activity` data 10, model 0.015 — the model has no basal floor).
Figure 2C's transient is reproduced in shape but with the peak displaced (data peak 220 at t = 1,
model 213 at t = 1.4) — that displacement is most of the 22% SSE share. Figure 3C's biphasic SOD2
response is qualitatively right (model 5.8 → 8.5 → 8.5 → 6.6 → 2.5 against data 5 → 8 → 10.5 → 3.5 →
2), rising at weak stress and falling at strong — but as the table above shows, the objective is
indifferent to whether it is.

`tools/sigma_profile.py` is a no-op here: σ is fixed at 1, so there is nothing to profile.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`; **no hand corrections**.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 100`,
  `max_iterations = 1000`, `wall_time_sim = 10`, `random_seed = 1`, `objective = sos`.
- k = 25 free parameters (16 rate constants + 9 observable scale factors), n = 62 scored points
  across 9 observables and 36 experiments, over a 367-reaction network.
- All 133 species are declared `hasOnlySubstanceUnits="true"` and the compartments run
  `8.3e−12 … 1e−13`.

## Provenance — this row is not reproducible on the collection's usual stack

**This is the only fitted row in subset I produced against a local, editable `bngsim` working copy
rather than the released PyPI wheel, and against PyBNF past lanl/PyBNF#553.** Both differences are
load-bearing and they are two differences, not one:

| | this row | every other fitted row |
|---|---|---|
| bngsim | **working copy at `dff901e`** (2026-08-11), editable install | PyPI 0.12.1 / 0.12.2 wheel |
| PyBNF | **`095a5a14`** (contains #553 at `85b36a96`) | pre-#553 |
| numpy | 2.5.2 | — |

- **lanl/PyBNF#553 changes the objective this slug minimizes.** Before it, every
  `hasOnlySubstanceUnits` species reached the observable formulas as `amount / compartment_size`,
  inflating all 62 points by up to 1.6e13 and putting the nominal point at `6.85e+32`. Smith is the
  only slug of the 23 with both preconditions, so no other row moved.
- **bngsim#160/#161 is what makes the gradient viable at all.** Reaction 7 (`R5f`) is
  `per_species_volume_scaling`; without that commit the whole model declines the analytic `df/dp`,
  all 25 columns fall back to CVODES difference quotients, and every start times out to `inf`. The
  capability gate cannot see this: `BNGSIM_HAS_EVENT_SENS` is a **version floor of exactly 0.12.2**,
  and the PyPI 0.12.2 wheel passes it while lacking the commit. **Probe capability, never version** —
  `'_psvs_row_divisor' in inspect.getsource(bngsim._codegen)` and `hasattr(bngsim, 'AUTO')`.

Measured on this build, the gradient path is fully live: **0 decline warnings, 0 difference-quotient
fallbacks**, and the discrete event was rooted and jumped in flight 1,131,253 times across the run's
1,130,847 simulations (`'PI345P3>pip3_basal'`).

`nominal_check.json`'s recorded values were produced on the PyPI wheel. On `dff901e` the same nominal
point evaluates to `888141.4641019923` against the recorded `888141.4642017893` — **1.12e−10
relative**. The file has not been rewritten over a difference four orders below the last digit that
matters.

### The run this replaces

A first attempt was started and discarded. The environment reverted to the PyPI `bngsim` mid-run —
PyBNF's `uv.lock` pins `bngsim` to the registry, so any `uv sync` or `uv run --project` reinstalls
the wheel over the editable install — and a `uv pip install -e` that *reported success* had in fact
reused a stale cached wheel. `uv pip install -e … --reinstall --no-cache` is required. The
distinguishing check used throughout the second run: `site-packages/bngsim/_codegen.py` exists **only**
when the wheel has clobbered the editable install, because the editable install leaves only the
compiled extension there. That file was absent for all 13 hours of the run that produced this row,
and all four capability probes passed both before and after it.

## Bottom line

The last of subset I's 23 problems and the one that needed two upstream fixes to become fittable at
all: solved at `OG = 0.50` from 100 unbiased starts in 13 h 11 m, beating its own nominal point
1.7 million-fold and escaping a local optimum at 24,636 that captured three independent starts. The
caveat is not about the objective but about what the objective can see — with `σ ≡ 1`, Figures 2B, 3B
and 2C carry 99.96% of it, and the fit's flat H₂O₂-alone prediction contradicts Figure 3B at a place
where `OG` is blind.
