# SalazarCavazos_MBoC2020

**Run cost: `hours`** — 100,000 evaluations (100 × 1,000 `gntr`), each integrating a **618-reaction**
network, the largest model in the subset. The `hours` tier is the collection table's a-priori
estimate; the shipped fit measured **5 min 49 s** on ten cores.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**SOLVED** — `OG = 0.000029` from a from-scratch 100-start `gntr` fit in 5 min 49 s, matching `J*` to
five decimals. The PEtab nominal point is *also* this problem's published optimum
(`OG_nominal = 0.326`), so the objective is validated independently of the optimizer.

> **A solved objective is not a parameter-recovery result here.** The fit reaches `J*` at a parameter
> point substantially unlike the published one: `SHC1_total__FREE` has gone to its **upper bound**
> (`1e6`, reached to within 2e−08) and `GRB2_total__FREE` sits at ~2.2× its published value, while the
> two dephosphorylation parameters and `ratio_kpkd_YN` are recovered closely. Two very different
> vectors reaching the same likelihood to within 3e−05 is a **flat direction** — the two
> total-abundance parameters trade off, and 18 scored points do not pin them down. `OG` is defined on
> the objective, so the ✅ is correct as scored; it means *found an equally good optimum*, not
> *recovered the published parameters*. See `VALIDATION.md`.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `366.86157299694054` |
| paper-scale NLL at the PEtab nominal point | `367.18760380827933` |
| optimality gap at nominal | `0.3260308113387964` |
| scored data points `n` | 18 |
| free parameters `k` | 6 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = gntr` — general-objective Fisher/Gauss-Newton trust region (EFIM Hessian through trf's
Coleman–Li core, ADR-0068). `tools/fd_check.py` verifies its assembled gradient against central
differences at `4.0e−06`, the best agreement of the eighteen `gntr` slugs in the corpus.

### This slug is why the collection's budget is 100 × 1000

The shipped default of 20 × 500 does not merely fall short here — it converges to `OG = 10.2`, which
is **worse than doing nothing**, since the problem's own nominal point scores `0.326`. Twenty
box-sampled starts did not contain the reference basin at k = 6.

| budget | OG |
|---|---:|
| 20 × 500 | 10.2 |
| **100 × 1000** | **2.9×10⁻⁵** |

Together with `Laske_PLOSComputBiol2019` (k = 13: `6.76` at 20 × 500, the reference optimum at
100 × 1000) this is the second problem showing the same thing at the opposite end of the `k` range,
and it is what moved the corpus default. The *distribution* the starts are drawn from is unchanged and
untested — these are box-sampled, log-uniform in log space, since every free parameter here is a
`loguniform_var` and its prior is its box.

## Contents

- `SalazarCavazos_MBoC2020.conf` — the PyBNF job
- `model_SalazarCavazos_MBoC2020.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment____EGF25nM.exp`, `experiment____dose1.exp`, `experiment____dose2.exp`, `experiment____dose3.exp` — experimental data
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `score.py` — scores a run against `J*`
- `best_fit_params.txt`, `information_criteria.txt` — the shipped fit's provenance
- `VALIDATION.md` — the full validation against `J*`

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`. The
converter preserves both `parameterScale` (lanl/PyBNF#491) and `observableTransformation`
(lanl/PyBNF#499), which plain `petab.v2.petab1to2` drops. The run recipe (`job_type`,
`sbml_backend = bngsim`, `wall_time_sim`) is supplied, not recovered — PEtab specifies a
problem, not a method. `wall_time_sim = 10` caps pathological parameter points; raise it
if valid simulations on your machine are being marked as failures.

## Running

```bash
pybnf -c SalazarCavazos_MBoC2020.conf -o
python score.py output
```
