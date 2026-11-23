# reduced_onoff — identifiable NF-κB model, ON-OFF protocol (PyBNF edition-2 job)

**Run cost: `hours`** — 20,000 evaluations (100 × 200 `de`) over the two-genotype on-off protocol.

A PyBNF **edition-2** parameter-fitting job that recovers the **13 parameters of the reduced 2023
NF-κB model** by refitting it to reproduce the original **Lipniacki-2004** model under the paper's
proposed identifiability-optimal **on-off** stimulation protocol (2 h TNF on, 10 h off), derived
from:

> Jaruszewicz-Błońska J, Kosiuk I, Prus W, Lipniacki T. **"A plausible identifiable model of the
> canonical NF-κB signaling pathway."** *PLoS ONE* 2023; **18**(6):e0286416.
> DOI: [10.1371/journal.pone.0286416](https://doi.org/10.1371/journal.pone.0286416).

Built + validated from the paper's own **S1 Code** (the authors ship their reduced and original
models as BNGL). See [`VALIDATION.md`](VALIDATION.md) for the full per-gate audit (confidence 89/100).

## What is fit

The paper reduces the 15-variable Lipniacki-2004 NF-κB model to a **6-ODE / 13-parameter**
non-dimensional model (Eqs 1–6), then **refits** the reduced model to reproduce the original model's
trajectories (with PyBioNetFit). This is a **synthetic-target / ground-truth-recovery** job: the fit
target is the *original* model's output, and "the paper's result" is **Table 1** (the fitted values).

- **Model** ([`reduced_onoff.bngl`](reduced_onoff.bngl)): the authors' reduced model, identical NF-κB
  network (verified by `model_diff`), ported to edition-2. TNF (`TR`, 0/1) is delivered by a
  time-dependent function `TRfunc()=if(simtime<7200,1,0)` gated on an inert simulation clock
  (`Clock()`, rule `0->Clock() 1`) — replacing the authors' `setParameter`/`simulate_ode` phase
  block (the Lin-2021/nyc_multiphase idiom). No pre-equilibration: the seed `(1,0,0,0,0,0)` **is** the
  exact TR=0 steady state.
- **Target** (`reduced_onoff_wt.exp`, `reduced_onoff_a20ko.exp`): the original Lipniacki-2004 model's
  **raw** on-off output at the **Table-2** measurement times, for WT (5 variables) and A20KO (4; no
  A20) — un-normalized (the objective floors and geomean-normalizes at scoring time). 50 measurements
  = 41 independent points (matches the paper). Regenerable with `validation/gen_onoff.py`.
- **A20KO** = `c_deg = 0` (no A20 synthesis; A20 stays 0 → no A20→IKKa feedback), set via the conf
  `a20ko` condition.
- **Objective** ([`reduced_onoff.conf`](reduced_onoff.conf)): the paper's **exact Eq-7** objective —
  `noise_model = lognormal, sigma = fix_at 1` + `normalization <obs> = floor 0.03, scale`
  (`lanl/PyBNF#479`): a fixed-σ log-normal residual on log₁₀ (= the paper's squared-log objective)
  with the ρ=0.03 measurement floor and per-series geometric-mean scaling, applied symmetrically to
  the simulated and experimental trajectories. (Edition-2-native; not PEtab-exportable.)

## Result

- **Gate 3a** — reduced model at Table-1 params reproduces the original on-off dynamics
  (AMD\*_WT≈1.41, AMD\*_A20KO≈1.32; see `reduced_onoff_reproduction.png` / `onoff_gate3a.png`).
- **Gate 3b** — the `de` fit (pop 100 × 200 + refine) under the paper's **exact Eq-7 objective**
  recovers **11/13 parameters within 3×**, the identifiable directions essentially exact. The 2 that
  drift — `δ`, `ε` — are exactly the parameters the paper's own identifiability analysis (Fig 5) flags
  as least identifiable; vs the earlier analog, the exact objective tightens `ε` (~9.5×→~3.7×) and
  brings `a_3` within 3×. Expected sloppy-fit behavior for this deliberately sparse (41-point)
  protocol; the sibling combination slug recovers all 13.

## Author-file check (BNG2.pl, asserted)

[`check_author_file.py`](check_author_file.py) runs this exact file and the authors' own on-off file
through BNG2.pl and asserts that they are the same model on the same protocol. The authors' file is
S1 Code of the article (`fitted_reduced_model_on_off_protocol.bngl`, the reduced model at the Table 1
fitted values), vendored byte-identical as
[`reference/upstream_fitted_reduced_model_on_off_protocol.bngl`](reference/) and pinned by URL and
sha256 in [`upstream.json`](upstream.json); the script refuses a copy whose sha256 differs. Nothing in
it uses bngsim. It writes [`reference/author_file_check.json`](reference/author_file_check.json) and
[`reduced_onoff_author_file_check.png`](reduced_onoff_author_file_check.png), and exits 1 on any
failure. The file checked is `reduced_onoff.bngl` with sha256 `629bd88d…4228`.

| check | result | tolerance |
|---|---|---|
| (a) the authors' TR = 0 hour (0–3600 s) leaves the seed (IKKn = 1, rest 0) unchanged, WT and A20KO | max change 0 | 1e-12 |
| (b) generated networks | authors 6 species / 14 reactions; this file the same plus `Clock()` and `0 -> Clock()` | identical |
| (c) this file at *t* vs the authors' at *t* + 3600 s: 6 species × WT, A20KO × 145 points (0–12 h after onset, the authors' 300 s grid) | worst max\|Δ\| / peak **7.7e-7** (A20KO IκBα); A20 in A20KO is identically 0 in both | 1e-5 of peak |

(a) is what lets this file start at TNF onset with no pre-equilibration. In (c) the authors' file runs
at its own CVODE tolerances (rtol = atol = 1e-16) and this one at BNG2.pl's defaults (1e-8); the
1e-5 tolerance, fixed before the run, is a thousand times the looser of the two and two orders below
the smallest change a parameter edit could make (the Table 1 values carry three significant
digits). The TNF-off step `TRfunc() = if(simtime<7200,1,0)` is a discontinuity inside a single
integration here, where the authors restart CVODE at a phase boundary; the residual jumps at 120 min
in the figure and stays below 1e-6 of peak.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"
cd pybnf-jobs/Jaruszewicz-Blonska-2023/reduced_onoff
pybnf -c reduced_onoff.conf                       # the fit (Gate 3b)
BNGPATH=$BNGPATH python make_reproduction.py      # Gate 3a overlay at Table-1 params
BNGPATH=$BNGPATH uv run python check_author_file.py   # the author-file check (BNG2.pl only)
```

See the sibling [`../reduced_combination/`](../reduced_combination/) slug for the paper's headline
fit (the full 6-protocol combination experiment, Fig 2).
