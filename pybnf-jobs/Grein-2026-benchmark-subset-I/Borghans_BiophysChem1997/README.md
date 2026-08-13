# Borghans_BiophysChem1997

**Run cost: `hours`** — 64,000 evaluations (32 × 2,000 `cmaes`). Budget is *not* what stands between this job and a solve; see Status.

PyBNF fitting job imported from the [Benchmark-Models-PEtab](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab) collection, as used
in the Grein et al. (2026) optimizer benchmark (bioRxiv 2026.07.11.737731).

## Status

**Setup only — never solved from an uninformed start, and the reason is the transcription rather than
the budget.** The job imports, simulates, and scores correctly, and every gate this collection applies
passes on it. It *has* been solved once, from a privileged start. What has never been demonstrated is a
fit found from box-sampled starts with no seeding, which is what a ✅ row in this collection means.

### The setup is verified

| gate | result |
|---|---|
| gradient (`tools/fd_check.py`) | worst rel err **1.47e-07** over all 23 columns; no zero and no sign-reversed column |
| objective oracle (independent Eq. 6 recompute, no PyBNF in the loop) | **111 of 111** rows join; oracle `J_paper = -83.3237191` against PyBNF's `-83.3236778` (**5.6e-07** relative) |
| search scale | 20 `loguniform_var` + 3 `uniform_var`, matching upstream's 20 `log10` + 3 `lin` |
| integration (`tools/box_probe.py`) | state 26/31, sensitivity 23/31, 12 s — not a `Weber`/`Brannmark` tolerance case |

Its nominal σ is also **not** a placeholder (0.104923 against an MLE of 0.101708), so σ-profiling moves
`OG_nominal` only 48.685 → 48.579. The 48.7 below is an honest distance, unlike `Giordano_Nature2020`'s.

### It has been solved, from a privileged start

`OG = -1.282656`, reduced objective `-248.069154`, verified three independent ways including through
PyBNF's own objective with `gntr` seeded at the vector (`-248.0691541186748`). The fit is a real
oscillator — Pearson `r = 0.864` against the data, three peaks, profiled σ = 0.0649 (~16% relative) —
not a degenerate flat trajectory. It came from the multiple-shooting prototype of lanl/PyBNF#563,
started at a **radius-0.4 perturbation of the PEtab nominal point**: a trajectory that already
oscillates with the period wrong, which is the regime that transcription is for. In that regime
multiple shooting beat `gntr` **9–0–1** over ten paired starts and never lost.

The prototype, the solved vector and its re-score live in `multiple_shooting_prototype/`:
`solved_seed3_r0.4.json` is the parameter vector, and `verify_gntr.conf` re-scores it through PyBNF's
ordinary config surface at `-248.069154` (the prototype's own recomputation: `-248.069166`).
`reproduce.py 3 0.4 600` replays it in about 100 s on one core.

Two things keep this a basin measurement rather than a result:

* **The start is privileged information.** It establishes that a basin at `OG < 0` exists and is
  reachable; it does not establish that a search finds it unaided.
* **The path does not survive a rebuild of the ODE library.** A **2.2e-16** difference in 1 of 111
  residuals (64 ULPs — six orders of magnitude *inside* the requested `rtol = 1e-9`) at evaluation 199
  of ~8,400 changes the answer, and three legitimate builds of bngsim give three different answers from
  an identical start. The basin is a fact; the trajectory to it is not a property of any method.

### From an uninformed start, everything lands on the same flat line

Running tally across every PyBNF configuration tried: **0 successes in 19 CMA-ES runs, 500+ `gntr`
starts, 1 PSO, 1 scatter search**, plus a 24-start box-drawn sweep run through both single and multiple
shooting (0/24 either way). Best `OG` anywhere is **77.6** (`gntr`, 100 × 1000), against a threshold of
1.92.

The completed 15-run BIPOP-CMA-ES campaign (λ₀ = 32, 12 restarts, `cmaes_run_maxgen = 300`, ~33,000
simulations per run) is the sharpest form of it:

| | `OG` | reduced objective |
|---|---:|---:|
| best | 79.0680 | -167.7185 |
| median | 80.6068 | -166.1797 |
| worst | 80.8044 | -165.9821 |

A **1.74-unit** spread across fifteen independent global searches, against a **76.8-unit** gap. These
are not near-misses scattered around a hard basin. Where they stop is analytic: a flat line at the best
constant with σ at the residual RMS scores `J_paper = -51.204092`, `OG = 80.804`. Every one of the 5,000
retained points in each completed run sits in that no-dynamics band — not one oscillating point
survives, because oscillating points score ~25 NLL units worse and are dropped first.

### Why: a wrong-period oscillator scores worse than no dynamics at all

Rescale time by α — multiply the 9 rate constants, leave the 6 concentration constants, so
`Z(t) → Z(αt)`. Only **α ∈ [0.9548, 1.0234]**, a **−4.5% / +2.3% window in period**, beats a horizontal
line:

| α range | best reduced objective | vs. the flat line (-165.98) |
|---|---:|---|
| 0.50 – 0.90 | -144.9 | **worse than flat** |
| **0.955 – 1.023** | **-198.1** | better — the only window that is |
| 1.10 – 2.00 | -141.3 | **worse than flat** |

So under single shooting the flat line is the **ceiling over essentially the whole box**, and a global
search ranking candidates by the objective is *correctly* pushed away from the only region a solve
lives in. The chance that a box-uniform draw lands in a ~3% period window across 20 log dimensions over
8 decades is effectively zero. This is a statement about the **transcription**, not the search — which
is why more starts, a better global method, and a gradient polish all return the same answer.

### What would settle it

Grein et al. solved this slug with CMA-ES in **2 of 10 runs**, at a per-run budget evidently well above
the ~33,000 simulations/run reached here. The outstanding demonstration is therefore roughly **10× the
per-run budget, run ~10 times** — a cluster-scale experiment whose outcome is predictable to about one
binomial draw, not an open question. wshlavacek/BNGL-Models#38 was closed on that basis; the full
measurement record is in that thread and in lanl/PyBNF#563.

### Three claims this README used to make, all disproved

1. **"Strongly multimodal … needs a large multistart budget."** Budget is not the axis: 400 starts and
   100 × 1000 give the same answer, and the reference optimum is 76 reduced-objective units away.
2. **"Multimodal, therefore `cmaes`."** `cmaes` is *worse* than `gntr` here (-165.98 against -169.19),
   and both land on the flat line. Neither picks a worse basin; neither finds a basin at all.
3. **A box-corner effect** adversarial to a Gaussian sampler (`init_A_state` and `init_Y_state` sit on
   the `[0,1]` bound at nominal). The nominal point is a **coordinate-wise minimum in all 23
   directions** across the full box, on 41-point profiles. The corner is real and irrelevant.

## Reference

| quantity | value |
|---|---|
| reference `J*` (Grein et al., best over all optimizer runs) | `-132.00847649739424` |
| paper-scale NLL at the PEtab nominal point | `-83.32367776169257` |
| optimality gap at nominal | `48.684798735701676` |
| scored data points `n` | 111 |
| free parameters `k` | 23 |

`J*` is the minimum Eq. 6 Gaussian NLL over every optimizer run on Marvin
(`best_fx_marvin.csv`). A fit is "solved" iff `OG = -log_likelihood - J* < 1.92`
(chi-square, alpha = 0.05, 1 dof). `score.py` computes this.

## Optimizer

`job_type = cmaes` — CMA-ES with IPOP restarts (ADR-0070/0082). **This is not a measured choice, and
on the evidence above `cmaes` is the weaker of the two methods here** (-165.98 against `gntr`'s
-169.19, both on the flat line). It is kept because CMA-ES is the only method known to have solved this
problem from an uninformed start — Grein et al. report 2 of 10 runs — so the shipped recipe is the one
aligned with the experiment that would settle the row. Two knobs matter if you run it:
`cmaes_run_maxgen` (lanl/PyBNF#507), because `max_iterations` is a **global** generation budget across
all restarts and the early small-λ runs will otherwise starve the large-population restarts that do the
multimodal work; and `wall_time_fit`, which silently downgrades `refine = 1` to no refine at all
(lanl/PyBNF#564, fixed upstream).

## Contents

- `Borghans_BiophysChem1997.conf` — the PyBNF job
- `model_Borghans_BiophysChem1997.xml` — SBML model (emitted by the importer, byte-reproducible)
- `experiment1.exp` — experimental data
- `jstar.txt` — the reference `J*`
- `nominal_check.json` — the nominal-point evaluation recorded above
- `score.py` — scores a run against `J*`

## Provenance

Imported with `pybnf.petab.petab1to2_preserve_scale` then `pybnf.petab.import_job`. The
converter preserves both `parameterScale` (lanl/PyBNF#491) and `observableTransformation`
(lanl/PyBNF#499), which plain `petab.v2.petab1to2` drops. The run recipe (`job_type`,
`sbml_backend = bngsim`, `wall_time_sim`) is supplied, not recovered — PEtab specifies a
problem, not a method. `wall_time_sim = 10` caps pathological parameter points; raise it
if valid simulations on your machine are being marked as failures.

## Running

```bash
pybnf -c Borghans_BiophysChem1997.conf -o
python score.py output
```
