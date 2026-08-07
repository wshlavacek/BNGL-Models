# VALIDATION — Okuonghae_ChaosSolitonsFractals2020

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's **reclassified** case: it
was carried as multimodal and `cmaes`-only, and it is neither.

> **Confidence: 76 / 100.** SOLVED with OG = 0.0012 from a from-scratch 32-start multi-start on
> **`gntr`**, no seeding, in 12 min 47 s — from a nominal point 4.7×10⁵ away. Deductions: the model is
> imported (not re-derived from Okuonghae & Omame 2020); the run is a single seed; and — the material
> one — **this slug has no independent oracle**, since upstream ships no `simulatedData_*.tsv` for it.
> See "Unverified" below.

## Gate A — objective fidelity

Linear observables under Gaussian noise with an estimated scalar σ, so the restored constant is
`(N/2)log(2π)` alone:

| quantity | value |
|---|---:|
| PyBNF reduced objective | 289.006519 |
| restored constants | 84.542345 |
| `(N/2)log(2π)`, n = 92 | 84.542345 |
| `J_paper = −log_likelihood` | 373.548864 |
| reference `J*` | 373.547658 |

The restored constant is `(N/2)log(2π)` to every digit.

**Verdict: PASS.**

## Gate B — the fit reaches the benchmark optimum

From-scratch multi-start **`gntr`** (32 starts × 2000 iterations, box-sampled starts seeded by
`random_seed = 1`, `sbml_backend = bngsim`, no seeding from the nominal point):

    J_pybnf = 289.0065  ⇒  J_paper = 373.5489   against  J* = 373.5477
    OPTIMALITY GAP  OG = 0.001206  <<  1.92

Wall time **12 min 47 s** on 10 cores. The nominal point sits at `OG_nominal = 4.7×10⁵`, so the
optimizer crossed roughly five orders of magnitude from unbiased box starts.

**Verdict: PASS (SOLVED).**

## Parameters versus the source paper — the ✅ is on the objective only

Okuonghae & Omame (§4.3) estimate **six** quantities: the transmission rate `β_c`, the two case
detection rates `θ` and `ψ`, and the three initial infected compartments `E(0)`, `A(0)`, `I(0)`.
Everything else in their Table 3 is **fixed from the literature** — `α = 0.5`, `ν = 0.5`,
`σ = 1/5.2`, `γ_i = 1/15`, `γ_a = γ_0 = 0.13978`, `d_0 = d_D = 0.015`.

**The PEtab problem frees all ten of those**, fitting 16 parameters against 92 points. That is the
source of this problem's difficulty, and the fit shows what it costs:

| parameter | our fit | paper | paper status |
|---|---:|---:|---|
| `d_0` | **7.246** | 0.015 | FIXED |
| `sigma` | **7.535** | 0.1923 (1/5.2) | FIXED |
| `transmission_rate_effective` (`β_c`) | **7.989** | 0.4236 | fitted |
| `alpha` | **1.00e−05** — *at its lower bound* | 0.5 | FIXED |
| `gamma_a` | **1.00e−05** — *at its lower bound* | 0.13978 | FIXED |
| `gamma_0` | 6.13e−05 | 0.13978 | FIXED |
| `nu` | 5.12e−03 | 0.5 | FIXED |
| `gamma_i` | 0.0133 | 0.0667 (1/15) | FIXED |
| `d_D` | 0.0142 | 0.015 | FIXED |
| `psi` (`ψ`) | 1.09e−04 | 0.0135 | fitted |
| `theta` (`θ`) | 6.54e−04 | 1.8999e−11 | fitted |
| `exposed_start` (`E(0)`) | 1.49 | 441 | fitted |
| `asymptomatic_start` (`A(0)`) | 8982 | 188 | fitted |
| `symptomatic_start` (`I(0)`) | 1.25 | 212 | fitted |

`d_0 = 7.25 /day` puts the mean time from symptoms to death at about **three hours**, against the
paper's 0.015 /day. Two parameters rest exactly on their lower bounds. **The best-fit point is
biologically meaningless**, and it still matches `J*` to 0.0012 — so the reference optimum has the same
character. We did not find a different answer from the benchmark; we found the benchmark's answer, and
that answer is an artefact of freeing parameters the authors deliberately pinned.

`OG` is defined on the **objective**, not on parameter distance, so the ✅ is correct as scored. It
certifies that PyBNF's optimizer reaches the reference objective. It makes no claim that the fit
recovers the published parameters, and here it plainly does not.

### The paper already shows the identifiability limit

Its Table 4 fits three different data sets and reports `θ` = **1.8999e−11, 4.2719e−11, 2.3752e−04** — a
seven-order-of-magnitude spread — while `β_c` stays put at 0.4236 / 0.4410 / 0.4385. So `θ` is
unidentifiable in the authors' own hands on the *reduced* six-parameter problem. Ours lands at
6.5e−04, within that spread.

The paper also fits with a **genetic algorithm to locate the basin of attraction, then MATLAB's
`fmincon` to refine** (§4.3) — a global-then-local strategy. That is worth holding beside the
reclassification below: the surface is genuinely awkward, and 32 `gntr` starts happen to supply enough
of a global sweep to land in the same basin.

## Why this slug is on `gntr`, and why it was not

The shipped recipe was `cmaes` with 12 IPOP restarts, and the collection tracked this slug as one of
three **multimodal** problems where "a local method from a few starts reliably lands in a wrong basin;
they need a global budget, not a better local step."

That grouping was inherited rather than measured, and for this problem it is wrong. There was never a
capability reason for `cmaes` — the model has **0 events, 0 `piecewise`, 0 `and`**, 9 species and 10
reactions, so nothing gates the gradient path — and `tools/fd_check.py` verifies the assembled
gradient against central differences at an interior, bounds-clear point to **6.27e−05**, the second
best in the corpus. Run on `gntr` at the same 32 × 2000 budget the `cmaes` conf carried, it solves.

The other two slugs in that group behave as the flag claimed, which is what makes this a
reclassification rather than a repudiation:

| slug | k | `OG_nominal` | `gntr` result |
|---|---:|---:|---|
| **`Okuonghae`** | 16 | 4.7e+05 | **`OG = 0.0012` — solved, 12 m 47 s** |
| `Borghans_BiophysChem1997` | 23 | 48.7 | `OG = 68.9` — *worse* than its own nominal |
| `Elowitz_Nature2000` | 21 | 2.43 | `OG = 5.82` — *worse* than its own nominal |

`Borghans` and `Elowitz` show the genuine wrong-basin signature and stay on `cmaes`. Only this one
moved.

## Unverified — the caveat that keeps this at 76

Upstream ships no `simulatedData_*.tsv` for this problem, so it cannot take the independent-oracle
check that the `obj ✓` rows in the collection README carry: recomputing the Eq. 6 NLL at the nominal
point from the upstream PEtab tables with no PyBNF in the loop.

That matters here more than usual, because a very large `OG_nominal` collapsing to ~0 under a fit is
**also** what `Weber_BMC2015` looked like before lanl/PyBNF#547 — where the objective itself was wrong
by 13,740 and the large `OG_nominal` was the artefact, not a statement about the problem. The
difference is that this slug's fit actually reaches `J*` to 0.0012, which a corrupted objective would
not do, and that its gradient is independently FD-verified. The evidence is good; it is simply not the
independent evidence ten other slugs have.

This slug does not pre-equilibrate, so it is not exposed to #547.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`; **no hand corrections**.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 32`,
  `max_iterations = 2000`, `wall_time_sim = 10`, `random_seed = 1`.
- k = 16 free parameters, n = 92 scored points.
- The budget is the one the `cmaes` conf carried (32 × 2000), not the collection's 100 × 1000 default
  — it was left unchanged so the method is the only variable this result speaks to. A 100 × 1000 run
  has not been tried and may well be faster.

## Provenance

Run 2026-08-07 against **PyBNF `749a90c6`** (ADR-0104) and **bngsim `59c6f38`**.

Source model: Okuonghae D, Omame A, *Chaos, Solitons & Fractals* **139** (2020) 110032 —
`doi:10.1016/j.chaos.2020.110032`, a COVID-19 population-dynamics model for Lagos, Nigeria, recorded
in the SBML metadata as `http://identifiers.org/doi/10.1016/j.chaos.2020.110032`.

## Bottom line

16 parameters against 92 points, solved at `OG = 0.0012` in under 13 minutes on a method the
collection had recorded as unsuitable for it. The finding is not that `cmaes` is wrong in general —
`Borghans` and `Elowitz` both go *backwards* on `gntr` — but that "multimodal" had been applied to
this slug by association rather than by measurement, and cost it a working recipe.

The second finding is about the problem rather than the optimizer: the benchmark estimates 16
parameters where the source paper estimated 6, freeing ten it had fixed from the literature. The
result is a solved objective at a biologically impossible point — a symptomatic death rate 480× the
published value, and two parameters against their bounds. That is a property of the benchmark problem,
not of this fit, since `J*` sits in the same place.
