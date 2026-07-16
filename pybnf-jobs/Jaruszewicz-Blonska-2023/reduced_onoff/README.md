# reduced_onoff — identifiable NF-κB model, ON-OFF protocol (PyBNF edition-2 job)

A PyBNF **edition-2** parameter-fitting job that recovers the **13 parameters of the reduced 2023
NF-κB model** by refitting it to reproduce the original **Lipniacki-2004** model under the paper's
proposed identifiability-optimal **on-off** stimulation protocol (2 h TNF on, 10 h off), derived
from:

> Jaruszewicz-Błońska J, Kosiuk I, Prus W, Lipniacki T. **"A plausible identifiable model of the
> canonical NF-κB signaling pathway."** *PLoS ONE* 2023; **18**(6):e0286416.
> DOI: [10.1371/journal.pone.0286416](https://doi.org/10.1371/journal.pone.0286416).

Built + validated from the paper's own **S1 Code** (the authors ship their reduced and original
models as BNGL). See [`VALIDATION.md`](VALIDATION.md) for the full per-gate audit (confidence 86/100).

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
  on-off output at the **Table-2** measurement times, for WT (5 variables) and A20KO (4; no A20),
  peak-normalized with the paper's ρ=0.03·max floor. 50 measurements = 41 independent points (matches
  the paper). Regenerable with `validation/gen_onoff.py` from the authors' original-model BNGL.
- **A20KO** = `c_deg = 0` (no A20 synthesis; A20 stays 0 → no A20→IKKa feedback), set via the conf
  `a20ko` condition.
- **Objective** ([`reduced_onoff.conf`](reduced_onoff.conf)): `objective = norm_sos` +
  `normalization = peak` — the edition-2-native analog of the paper's Eq-7 objective (sum of squared
  geometric-mean-normalized log-differences). Native-only (not PEtab-exportable).

## Result

- **Gate 3a** — reduced model at Table-1 params reproduces the original on-off dynamics
  (AMD\*_WT≈1.41, AMD\*_A20KO≈1.32; see `reduced_onoff_reproduction.png` / `onoff_gate3a.png`).
- **Gate 3b** — the `de` fit (pop 100 × 150 + refine) recovers **10/13 parameters within 3×**, most
  within 1.3–1.5× (k_1, k_2, i_1a, c_deg, c_3a near-exact). The 3 that drift — `a_3`, `δ`, `ε` — are
  exactly the parameters the paper's own identifiability analysis (Fig 5) flags as least
  identifiable. This is the expected sloppy-fit behavior, consistent with the paper.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"
cd pybnf-jobs/Jaruszewicz-Blonska-2023/reduced_onoff
pybnf -c reduced_onoff.conf                       # the fit (Gate 3b)
BNGPATH=$BNGPATH python make_reproduction.py      # Gate 3a overlay at Table-1 params
```

See the sibling [`../reduced_combination/`](../reduced_combination/) slug for the paper's headline
fit (the full 6-protocol combination experiment, Fig 2).
