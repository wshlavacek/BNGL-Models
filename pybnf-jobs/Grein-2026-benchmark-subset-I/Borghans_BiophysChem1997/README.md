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
PyBNF's own objective with `gntr` seeded at the vector (`-248.0691541186748`). The fit reproduces the
record — Pearson `r = 0.864` against the data, three spikes inside the record and a fourth just past
it, profiled σ = 0.0649 (~16% relative) — **but it is not an oscillator**: simulated past the record
it never spikes again, and store calcium grows without bound (see *What the fits are, dynamically*).
It came from the multiple-shooting prototype of lanl/PyBNF#563,
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

### What the fits are, dynamically

A fit only ever simulates to the last measurement (t = 8.98), so its score says nothing about what
the model does afterwards. `campaign/is_periodic.py` runs each reference vector to t = 400 and
applies three checks on the second half of that horizon: the peak train (count, spacing and its
coefficient of variation, amplitude trend), the boundedness of the **state** (largest state norm
over successive windows), and the fixed point with its Jacobian eigenvalues. `OG` is from
`campaign/borghans_ode.py`, an independent scipy implementation whose objective differences match
PyBNF's to 0.01.

| vector | `OG` | peaks, spacing, amplitude trend, envelope growth on [200, 400] | verdict |
|---|---:|---|---|
| paper, Fig. 8 values (BioModels BIOMD0000000044) | 70.6* | 212 / 0.39 (CV 0.89) / 1.000 / 1.000 | **periodic bursting**, two spikes per burst, burst period 2.83 |
| PEtab nominal | 48.6 | 72 / 2.76 (CV 0.00) / 1.000 / 1.000 | **periodic**, one spike per period |
| `best_periodic_fit.json` | 37.8 | 71 / 2.81 (CV 0.00) / 1.000 / 1.000 | **periodic**, one spike per period |
| the `OG = -1.28` vector | −1.28 | 0 / – / – / **1.333** | **no attractor within reach** |

\* The paper never fitted these data; 70.6 is its Fig. 8 parameter set scored with the PEtab nominal
initial values and only `scale`/`offset` profiled.

**The `OG = -1.28` vector is four excitable spikes on a calcium-loading ramp.** Its net calcium
inflow is `v0 + beta*v1 = 3.19` per time unit against an efflux rate `K_par = 0.075` and a store leak
`Kf = 0.003`, so the cell takes in calcium it cannot excrete and stores it: `Y` is 29 at t = 9, 634 at
t = 200, 18,890 at t = 6,000, rising 3.2 per time unit, while `A` pins at 0.2236 and cytosolic `Z`
creeps from 0.13 to 0.75. The system has exactly one fixed point, `Z* = 42.6, Y* = 48,711,
A* = 0.2236` (`Ca* = 43`, forty times the data's range), a stable node — eigenvalues −4.28, −0.077,
−0.003, all real — and the trajectory reaches it at t ≈ 20,000. So "aperiodic" and "a fixed point
after a transient" are both literally true and both misleading: the transient is 2,000 times longer
than the record, it is a monotone ramp rather than a decaying oscillation, and the fixed point is
unphysical. The four spikes are the fitted initial values kicking an excitable system once, in the
first 0.05 % of a trajectory that is filling the cell with calcium. This is also why nothing near the
vector converges to it: there is no basin in the usual sense. From ten per-coordinate 0.4-decade
kicks of it, `gntr` ended on the flat line eight times, at `OG` 54–57 twice, and returned zero times.

**What follows for the benchmark.** The PEtab problem is a likelihood over 111 points on
[0.03, 8.98] and nothing else, so a cell that is filling with calcium is a legal optimum of the
problem as posed, and the benchmark cannot tell it from an oscillator. Grein's `J*` is 1.28 units
from this vector; whether their optimum is of the same kind is unknown (they publish objective
values, not parameter vectors). A *meaningful* solution needs a requirement the problem does not
contain — bounded, periodic dynamics — and a fit conditioned on that is a different problem with a
different optimum. Its value is **not known**. The best periodic fit known is `best_periodic_fit.json`
at `OG = 37.8` (a sustained oscillator, period 2.80; see *From an uninformed start* below for where
it came from); nothing says that is the periodic floor.

**What PyBNF (1.8.1) can express.** BPSL has no periodicity primitive. A usable proxy is a
constraint-only experiment run past the record — `experiment: qualitative, data: late.prop,
t_end: 40` with `Ca > 0.7 once between 30, 40` and `Ca < 0.4 once between 30, 40` — which the
loading-ramp vector fails and every periodic vector above passes. But `gntr`, `lbfgs` and `ms` all
refuse a fit that carries property files, which leaves the metaheuristics; a CMA-ES seeded tightly
around each of the vectors above did not stay local and produced nothing usable. So the
periodic-conditioned optimum remains unmeasured. An objective on the Jacobian's eigenvalues (an
unstable fixed point as the condition for oscillation) does not exist in PyBNF either.

Reproduce the table with any interpreter that has numpy and scipy (PyBNF's venv does):

```bash
python3 campaign/is_periodic.py                                  # the four reference vectors
python3 campaign/is_periodic.py output/Results/sorted_params_final.txt   # any PyBNF result
```

### From an uninformed start, everything lands on the same flat line

Running tally across every PyBNF configuration tried: **0 successes in 19 CMA-ES runs, 500+ `gntr`
starts, 1 PSO, 1 scatter search**, plus a 24-start box-drawn sweep run through both single and multiple
shooting (0/24 either way). Best `OG` from a plain search is **77.6** (`gntr`, 100 × 1000), against a
threshold of 1.92.

One start that is *not* a plain search did better: a random box draw that oscillates (about 1 in
8,000 does, at any timescale), its period set to the data's 2.70 by the exact nine-rate time rescale,
its phase set by taking the limit-cycle state 0.79 before a peak as the initial values, brought into
the box by the unit symmetries and clipping, with `scale`/`offset` profiled — no nominal point and no
fitted vector involved. From that start `ms` (8 segments, coarsening 2) reached `OG 38.3` in-box and a
`gntr` polish `37.8`: `best_periodic_fit.json`, a sustained oscillator. Kicks and continuations from
it stall at 38–42. The start-construction scripts are not in `campaign/` yet; the vector is committed
so its score and dynamics are reproducible.

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

Rescale time by α — multiply the 9 rate constants (`v0`, `v1`, `Vm2`, `Vm3`, `Kf`, `K_par`, `Vp`,
`Vd`, `epsilon_par`), leave the 6 concentration constants, so `Z(t) → Z(αt)` exactly. From the nominal
point, only **α ∈ [0.912, 1.048]**, a **−8.8% / +4.8% window in period**, beats a horizontal line
(σ profiled, 0.001 grid):

| α range | best reduced objective | vs. the flat line (-165.98) |
|---|---:|---|
| 0.50 – 0.90 | -161.7 | **worse than flat** |
| **0.912 – 1.048** | **-198.2** | better — the only window that is |
| 1.10 – 2.00 | -157.9 | **worse than flat** |

(Corrected 2026-09-10. The earlier table, and the −4.5% / +2.3% window quoted from ADR-0109, came from a
7-parameter rescale that omitted `Kf` and `epsilon_par`; that set is not a symmetry of the equations,
so those numbers measured a shape change as well as a period change.)

So under single shooting the flat line is the **ceiling over essentially the whole box**, and a global
search ranking candidates by the objective is *correctly* pushed away from the only region a solve
lives in. The chance that a box-uniform draw lands in a ~14% period window across 20 log dimensions over
8 decades is effectively zero — and a box draw that oscillates at all, at any timescale, is ~1 in 8,000. This is a statement about the **transcription**, not the search — which
is why more starts, a better global method, and a gradient polish all return the same answer.

### Starting from a named vector

Three confs at the top level run the job's own problem — same box, same objective — from a
documented parameter vector, with a gradient polish (`gntr`, 200 iterations) as the recipe. `gntr`
scores the start before it moves, so the first objective it reports is the vector's own; set
`max_iterations = 0` to only score it. Each was verified to load and score at these values through
PyBNF on 2026-09-10.

| conf | vector | dynamics | reduced objective at the start | `OG` |
|---|---|---|---:|---:|
| `Borghans_start_best_ever_ramp.conf` | the `OG = -1.28` vector (`multiple_shooting_prototype/verify_best_fit_params.txt`) | excitable transient on a calcium-loading ramp; no attractor within reach | −248.069 | −1.28 |
| `Borghans_start_paper_fig8.conf` | Borghans 1997 Fig. 8 values (BioModels BIOMD0000000044); `scale`/`offset` profiled against the data, PEtab nominal initial values, σ at the residual | periodic bursting, two spikes per burst | −176.148 | 70.6 |
| `Borghans_start_best_periodic.conf` | `best_periodic_fit.json` | periodic, one spike per period | −208.948 | 37.8 |

The paper's values never saw these data, so the middle row is what the published mechanism scores
as published, not a fit. The other two are fits of the same likelihood that differ by 39 `OG` units
and by whether the model oscillates at all; which one is "better" depends on whether periodic
dynamics is part of the question, and the likelihood alone says it is not.

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
- `best_periodic_fit.json` — the best *periodic* fit known (`OG = 37.8`); the `OG = -1.28` vector is not periodic
- `Borghans_start_best_ever_ramp.conf`, `Borghans_start_paper_fig8.conf`, `Borghans_start_best_periodic.conf` — the job started from each named vector (see *Starting from a named vector*)
- `campaign/is_periodic.py`, `campaign/borghans_ode.py` — the dynamics diagnostic and the scipy model it runs on (no PyBNF needed)
- `multiple_shooting_prototype/` — the prototype that reached `OG = -1.282656` (see above)
- `campaign/` — the drivers and conf templates behind the tallies quoted above; its
  runs are not committed (see that directory's `README.md`)

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
