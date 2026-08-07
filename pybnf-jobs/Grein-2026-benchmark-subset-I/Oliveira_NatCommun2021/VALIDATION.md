# VALIDATION — Oliveira_NatCommun2021

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's first **⚪ setup-only**
conversion: unlike the 🟢 rows, its PEtab `nominalValue` point is nowhere near the optimum
(`OG_nominal = 9.6e+06`), so nothing about the answer was known before the fit ran.

> **Confidence: 78 / 100.** SOLVED with OG = 0.0113 from a from-scratch 100-start multi-start, no
> seeding. Deductions: the model is imported (not re-derived from Oliveira et al. 2021); the run is a
> single seed; and — the material one — **this slug has no independent oracle.** Upstream ships no
> `simulatedData_*.tsv` for it, so unlike the ten `obj ✓` rows in the collection README its objective
> has never been cross-checked against anything but `J*` itself. See "Unverified" below.

## Gate A — objective fidelity

Noise is σ ≡ 1 upstream, so the conf ships a plain `objective = sos` with no `noise_model` line. The
collection README records why that is faithful rather than a fidelity break (`Why objective = sos is
faithful for Oliveira and Smith`): with σ ≡ 1 the Gaussian NLL reduces to the sum of squares plus
`(N/2)log(2π)`, and `Σ log σᵢ = 0` exactly.

| quantity | value |
|---|---:|
| PyBNF reduced objective | 7794.673035 |
| restored constants | 110.272624 |
| `(N/2)log(2π)`, n = 120 | 110.272624 |
| `J_paper = −log_likelihood` | 7904.945659 |
| reference `J*` | 7904.934317 |

The restored constant is `(N/2)log(2π)` to every digit, which is the σ ≡ 1 signature.

**Verdict: PASS.**

## Gate B — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (100 starts × 1000 iterations, box-sampled starts seeded by
`random_seed = 1`, `sbml_backend = bngsim`, no seeding from the nominal point) converges to
`J_pybnf = 7794.6730` ⇒ `J_paper = 7904.9457` against `J* = 7904.9343`:

    OPTIMALITY GAP  OG = 0.011342  <<  1.92

Wall time 28 min on 10 cores.

This is the collection's first ⚪ → ✅ conversion, and the informative part is the starting distance:
the nominal point sits at `OG = 9.6e+06`, so the optimizer crossed roughly seven orders of magnitude
from unbiased box starts. That is a genuine optimizer result rather than a refinement of a known
answer, which is what separates a ⚪ conversion from a 🟢 one.

It also fixes the budget boundary. 100 × 1000 now solves at k = 6, 9, **12**, 13 and fails to reach
the reference basin at k = 22 (`Fiedler_BMCSystBiol2016`, `OG = 1.004`), so the ceiling for this
recipe lies between k=12 and k=22.

**Verdict: PASS (SOLVED).**

## Unverified — the caveat that keeps this at 78

Twelve slugs in this collection carry `obj ✓` in the coverage table: their Eq. 6 NLL has been recomputed
at the nominal point directly from the upstream PEtab tables — `simulatedData` joined to
`measurementData`, with the declared transformation and nominal σ — with no PyBNF in the loop, and it
reproduces PyBNF exactly. **Oliveira cannot take that check**: upstream ships no `simulatedData` for
this problem, so there is no reference trajectory to compare against.

That matters because the same audit found two slugs where PyBNF and upstream disagreed —
`Brannmark_JBC2010` and `Weber_BMC2015`, wrong by 1531 and 13,740 through lanl/PyBNF#547 (fixed
2026-08-07 by ADR-0104; both now reproduce) — and in both cases the symptom was a large `OG_nominal`
that had been read for weeks as "the nominal point is not the optimum". Oliveira's
`OG_nominal = 9.6e+06` is exactly that shape. The difference is that Oliveira's *fit* reaches `J*` to
0.011, which a corrupted objective would not do — so the evidence here is good, it is just not the
independent evidence the `obj ✓` rows have.

Oliveira does not pre-equilibrate (no `preequilibrate:` in its conf), so it is not exposed to #547.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`; **no hand corrections**.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 100`,
  `max_iterations = 1000`, `wall_time_sim = 10`, `random_seed = 1`.
- k = 12 free parameters, n = 120 scored points.

## Provenance

Run 2026-08-07 against **PyBNF `47de23cf`** (ADR-0103, `sbml_atol` derived from the model) and
**bngsim `59c6f38`**. Unlike the three slugs solved on 2026-08-06, this fit ran entirely after that
commit landed, so no build boundary is crossed mid-run.

## Bottom line

The first ⚪ → ✅ conversion: 12 parameters against 120 points, solved at `OG = 0.011` in 28 minutes
from a nominal point seven orders of magnitude away. It sets the upper end of where the collection's
100 × 1000 default is known to work (k=12), and it is the clearest case in the corpus of a ✅ that
rests on `J*` alone, with no independent check available.
