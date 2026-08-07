# VALIDATION — SalazarCavazos_MBoC2020

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug is the collection's **budget** case: it is the
problem that a 20-start run drove *backwards*, and the one that settled how many starts this corpus
needs.

> **Confidence: 85 / 100.** SOLVED with OG = 2.9×10⁻⁵ from a from-scratch 100-start multi-start, no
> seeding. The fixed per-point σ identity is verified, and the assembled gradient checks against
> central differences at 4.0e−06 — the cleanest of the eighteen `gntr` slugs. Deductions: the model is
> imported (not re-derived from Salazar-Cavazos et al. 2020); the run is a single seed; and — the main
> one — **the fit reaches `J*` at a parameter point far from the published one, with one parameter
> resting on its box bound**. See "Non-identifiability" below. That is a statement about the problem,
> not a defect, but it means this row certifies the objective, not the parameters.

## Gate A — objective fidelity

Noise is **fixed per-point σ** from the measurement table's `_SD` columns, so the restored constant
carries both terms:

| quantity | value |
|---|---:|
| PyBNF reduced objective | 359.379668 |
| restored constants | 7.481933 |
| of which `(N/2)log(2π)`, n = 18 | 16.540894 |
| of which `Σ log σᵢ` | −9.058961 |
| `J_paper = −log_likelihood` | 366.861602 |
| reference `J*` | 366.861573 |

The `Σ log σᵢ` term is large and negative here, so this is a genuine test of the sign convention rather
than a case where the term nearly vanishes — the complement of `Rahman_MBS2016`, where it is exactly 0.

**Verdict: PASS.**

## Gate B — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (100 starts × 1000 iterations, box-sampled starts seeded by
`random_seed = 1`, `sbml_backend = bngsim`, no seeding from the nominal point) converges to
`J_pybnf = 359.3796685` ⇒ `J_paper = 366.8616020` against `J* = 366.8615730`:

    OPTIMALITY GAP  OG = 0.000029  <<  1.92

Wall time 5 min 49 s on 10 cores.

### This slug is why the collection's budget is 100 × 1000

The shipped default of 20 × 500 does not merely fall short here — it converges to `OG = 10.2`, which is
**worse than doing nothing**: the problem's own PEtab nominal point scores `OG_nominal = 0.326`. Twenty
box-sampled starts did not contain the reference basin at k = 6.

| budget | OG |
|---|---:|
| 20 × 500 | 10.2 |
| **100 × 1000** | **2.9×10⁻⁵** |

Together with `Laske_PLOSComputBiol2019` (k = 13: `6.76` at 20 × 500, reference optimum at 100 × 1000)
this is the second problem showing the same thing at the opposite end of the k range, which is what
moved the corpus default. The *distribution* the starts are drawn from is unchanged and untested —
these are box-sampled (log-uniform in log space, since every free parameter here is a `loguniform_var`
and its prior is its box).

**Verdict: PASS (SOLVED).**

## Non-identifiability: solved objective, different parameters

The fit matches `J*` to five decimals at a parameter point substantially unlike the published one:

| parameter | box | best fit | PEtab nominal (published optimum) |
|---|---|---:|---:|
| `GRB2_total__FREE` | [1e4, 1e6] | 369,948.97 | 169,853 |
| `SHC1_total__FREE` | [1e4, 1e6] | **999,999.98** | 649,426 |
| `kdephosY1068__FREE` | [0.1, 100] | 2.1743 | 1.6588 |
| `kdephosYN__FREE` | [1e-3, 100] | 0.017014 | 0.017182 |
| `ratio_kpkd_Y1068__FREE` | [0.01, 100] | 0.11567 | 0.15755 |
| `ratio_kpkd_YN__FREE` | [0.01, 100] | 0.44406 | 0.44476 |

`SHC1_total__FREE` has gone to its **upper bound** (`1e6`, reached to within 2×10⁻⁸ of it), and
`GRB2_total__FREE` sits at roughly 2.2× its published value. The two dephosphorylation parameters and
`ratio_kpkd_YN` are recovered closely; `ratio_kpkd_Y1068` is not.

Two very different vectors reaching the same likelihood to within 3×10⁻⁵ is a **flat direction**: the
two total-abundance parameters trade off against each other, and 18 scored points do not pin them down.
The optimizer slid along that floor until one hit its box edge.

This is worth stating rather than burying, for two reasons. A fitted value resting on a bound is
normally a signal that the box, not the data, chose the answer. And the benchmark's `OG` is defined on
the **objective**, not on parameter distance — so "solved" here means *found an equally good optimum*,
which is exactly what the metric asks for and is not the same as *recovered the published parameters*.
The ✅ is correct as scored; it should not be read as a parameter-recovery result.

## Gradient

`tools/fd_check.py` verifies the assembled gradient against central differences across all 6 free
parameters at an interior, bounds-clear point: worst `|Δ|/‖∇‖∞ = 4.0e−06`, the best agreement of the
eighteen `gntr` slugs in the corpus.

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`. Noise imported as fixed per-point σ from the
  measurement table (`_SD` columns); **no hand corrections**.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 100`,
  `max_iterations = 1000`, `random_seed = 1`.
- k = 6 free parameters, n = 18 scored points, 5 experiments (baseline plus `EGF25nM` and three doses).

## Provenance

Run 2026-08-06 against **PyBNF `d014272e`** and **bngsim `59c6f38`**.

PyBNF advanced to `47de23cf` (ADR-0103 — `sbml_atol` derived from the model's own scale rather than
inheriting the backend default `1e-8`) later the same day. Re-evaluating the recorded best-fit point
under that build gives `359.3796673` against the recorded `359.3796685` — a shift of 1.2×10⁻⁶, ODE
integration noise far below the 1.92 threshold; the recorded OG is unaffected. This slug is not among
the four whose nominal check the derived tolerance moved (Armistead, Bertozzi, Brannmark, Giordano).

## Bottom line

The budget case: 6 parameters against 18 points, solved at `OG = 2.9×10⁻⁵` in under 6 minutes at
100 × 1000, having gone *backwards* to `OG = 10.2` at 20 × 500. It is the corpus's clearest evidence
that the shipped start count — not the method — was what stood between a 🟢 row and a ✅ one. It is
also its clearest example that a solved objective does not imply recovered parameters: the fit sits on
a flat direction with `SHC1_total__FREE` against its upper bound.
