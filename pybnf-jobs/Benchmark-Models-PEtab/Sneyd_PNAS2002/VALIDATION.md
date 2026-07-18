# VALIDATION — Sneyd_PNAS2002

Validation against the **Grein et al. 2026** reference objective. Oracle = the benchmark's reference
**J\*** (best Eq. 6 NLL over all Marvin runs). This slug demonstrates the fidelity + scoring pipeline on
a **multi-condition** fit (9 conditions, one shared estimated σ).

> **Confidence: 95 / 100.** SOLVED with OG = 1.4×10⁻⁵ from a from-scratch multi-start; the objective is
> the linear Gaussian NLL (same, verified family as Boehm) and PyBNF's `−lnL` is the paper-scale NLL by
> construction. Deduction: the model is imported (not re-derived from Sneyd & Dufour 2002).

## Gate A — objective fidelity

`open_probability` has `observableTransformation = lin`, so the import emits linear `gaussian` noise
(correct — no correction needed). σ is a single **estimated** shared parameter, so the objective is the
full Gaussian NLL. PyBNF minimizes the reduced objective (`Σ (y−m)²/(2σ²) + log σ`, dropping the
per-point `½log(2π)`) and reports the full normalized `lnL` in `information_criteria.txt`; over
N = 135 points the restored constant is `C = 135·½log(2π) = 124.0567018`, and `−lnL = J_pybnf + C =
−319.7923321` to numerical precision. This is the same linear-Gaussian relationship established three
ways on Boehm; here it holds across 9 conditions summed into one objective. `score.py` reads `−lnL`.

Note the sign: `lnL = +319.79` is positive because the open-probability data (values ≲ 1) with a small
fitted σ gives a Gaussian **density > 1**, whose log is positive — so the NLL `−lnL` is negative. This
is why J\* itself is negative; nothing unusual about the fit.

**Verdict: PASS.**

## Gate B — the fit reaches the benchmark optimum

From-scratch multi-start `gntr` (10 starts × 300 iterations, box-center + Latin-hypercube seeded by
`random_seed = 1`, `sbml_backend = bngsim`) converges to `J_pybnf = −443.8490341` ⇒
`J_paper = −319.7923321` ⇒ **OG = 1.4×10⁻⁵ < 1.92**. Sneyd is an easy benchmark problem (329/380 Marvin
runs solved it); PyBNF reaches J\* to 5 significant figures. (A 10-start run suffices and converges
quickly; the earlier 20-start run was identical in objective but slower, unnecessary for this easy
landscape.)

**Verdict: PASS (SOLVED).**

## Configuration

- Import: `petab1to2_preserve_scale` → `import_job`. Emitted `noise_model = gaussian, sigma = fit sigma`
  — correct (`observableTransformation = lin`); **no hand corrections**.
- `edition = 2`, `sbml_backend = bngsim`, `job_type = gntr`, `population_size = 10`,
  `max_iterations = 300`, `random_seed = 1`.
- 9 conditions imported as 9 `experiment____*.exp` files with per-experiment `condition:` overrides.

## Bottom line

The multi-condition linear exemplar: 9 dose-response conditions with one shared estimated σ, summed
into a single Eq. 6 NLL that PyBNF's `gntr` solves to OG = 1.4×10⁻⁵. Confirms the pipeline scales
cleanly from Boehm's single condition to a 9-condition fit.
