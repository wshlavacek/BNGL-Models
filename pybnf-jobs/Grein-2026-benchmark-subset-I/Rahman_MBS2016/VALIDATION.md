# VALIDATION — Rahman_MBS2016

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's **unit-σ** case: every
measurement carries `_SD = 1`, so the fixed per-point σ term vanishes exactly and the restored
constant is `(N/2)log(2π)` alone.

> **Confidence: 93 / 100.** SOLVED with OG = 0.000000 from a from-scratch 100-start multi-start, no
> seeding — the fit reaches `J*` itself, agreeing to every digit reported. The objective identity is
> exact by construction here (see Gate A), and the assembled gradient is verified against central
> differences across all 9 free parameters at 2.3e−05. Deductions: the model is imported (not
> re-derived from Rahman et al. 2016), and the run is a single seed.

## Gate A — objective fidelity

Rahman's noise is **fixed per-point σ** read from the measurement table's `_SD` columns, and every one
of those values is `1`. The restored constant is therefore

    -lnL - J_pybnf  =  Σ log σᵢ + (N/2)·log(2π)  =  0 + 21.135586

and the `Σ log σᵢ` half is exactly zero rather than merely small — 23 × `log(1)`. Measured:

| quantity | value |
|---|---:|
| PyBNF reduced objective | 0.017900 |
| restored constants | 21.135586 |
| `(N/2)log(2π)`, n = 23 | 21.135586 |
| implied `Σ log σᵢ` | −0.000000 |
| `J_paper = −log_likelihood` | 21.153486 |

This makes Rahman the cleanest fidelity check in the collection: the identity `J_paper == −lnL` holds
with no σ bookkeeping to get wrong, so a discrepancy here could only come from the likelihood itself.

**Verdict: PASS.**

## Gate B — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (100 starts × 1000 iterations, box-sampled starts seeded by
`random_seed = 1`, `sbml_backend = bngsim`, no seeding from the nominal point) converges to
`J_pybnf = 0.0179001` ⇒ `J_paper = 21.1534861` against `J* = 21.1534861`:

    OPTIMALITY GAP  OG = 0.000000  <<  1.92

The fit lands **on** the reference optimum, not merely inside the threshold — the two values agree to
all seven reported decimals. Wall time 16 min 01 s on 10 cores.

Worth recording for the collection's open question about start counts: this was run at 100 × 1000
rather than the shipped default of 20 × 500, on the evidence from `Laske_PLOSComputBiol2019` and
`SalazarCavazos_MBoC2020` that 20 box starts is too few. Rahman was not separately tried at 20 × 500,
so it is not itself evidence about the smaller budget.

**Verdict: PASS (SOLVED).**

## Gradient

`tools/fd_check.py` verifies the assembled gradient against central differences across all 9 free
parameters at an interior, bounds-clear point: worst `|Δ|/‖∇‖∞ = 2.3e−05`. This matters because a
wrong gradient column is silent — the objective stays correct and the fit merely stops short — which
is the failure mode that cost `Laske_PLOSComputBiol2019` a solved verdict until lanl/PyBNF#534.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`. Noise imported as fixed per-point σ from the
  measurement table (`_SD` columns); **no hand corrections**.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 100`,
  `max_iterations = 1000`, `random_seed = 1`.
- k = 9 free parameters, n = 23 scored points, one observable (`prevalence`), one experiment.

## Provenance

Run 2026-08-06 against **PyBNF `d014272e`** and **bngsim `59c6f38`**.

PyBNF advanced to `47de23cf` (ADR-0103 — `sbml_atol` derived from the model's own scale rather than
inheriting the backend default `1e-8`) later the same day. Re-evaluating this slug's recorded best-fit
parameter point under that build gives `0.0179001035` against the recorded `0.0179000822` — a shift of
2.1e−08, i.e. ODE integration noise some eight orders below the 1.92 threshold. The recorded OG is
unaffected. Rahman is also not among the four slugs whose nominal check the derived tolerance moved
(Armistead, Bertozzi, Brannmark, Giordano).

## Bottom line

The unit-σ case, and the collection's tightest result: 9 parameters against 23 points, solved at
`OG = 0.000000` in 16 minutes from unbiased starts. Because every `_SD` is exactly 1, it demonstrates
the fixed per-point σ identity in its degenerate form — `Σ log σᵢ = 0` — which is the one variant that
cannot hide an error in the σ bookkeeping.
