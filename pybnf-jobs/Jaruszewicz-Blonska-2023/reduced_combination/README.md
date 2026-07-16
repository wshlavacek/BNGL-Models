# reduced_combination — identifiable NF-κB model, COMBINATION experiment (PyBNF edition-2 job)

The paper's **headline** fit (Fig 2). A PyBNF **edition-2** job that recovers the **13 parameters of
the reduced 2023 NF-κB model** by refitting it to reproduce the original **Lipniacki-2004** model
across the full **combination experiment** — a tonic protocol plus 5 pulsatile protocols, in WT and
A20KO cells, **914 independent data points** (S1 Table). Derived from:

> Jaruszewicz-Błońska J, Kosiuk I, Prus W, Lipniacki T. **"A plausible identifiable model of the
> canonical NF-κB signaling pathway."** *PLoS ONE* 2023; **18**(6):e0286416.
> DOI: [10.1371/journal.pone.0286416](https://doi.org/10.1371/journal.pone.0286416).

Built + validated from the paper's own **S1 Code**. See [`VALIDATION.md`](VALIDATION.md) for the full
per-gate audit (confidence 90/100). See the sibling [`../reduced_onoff/`](../reduced_onoff/) for the
paper's proposed identifiability-optimal on-off protocol.

## What is fit

Synthetic-target / ground-truth-recovery job: the fit target is the **original 15-variable
Lipniacki-2004 model's** output across the 6 protocols; "the paper's result" is **Table 1** (the
reduced model's fitted values); Gate 3b recovers Table 1.

- **6 protocols** (S1 Table): `continuous` (TNF 0–240 min) and 5 pulsatile — `pulse5_60`,
  `pulse5_100`, `pulse5_200` (5-min pulses at 60/100/200-min periods), `pulse22_5` (22.5-min pulses,
  45-min period), `pulse45` (45-min pulses, 90-min period) — each in **WT** (5 observables) and
  **A20KO** (4; no A20).
- **Model** — the authors' reduced model (identical NF-κB network; `model_diff` confirms), edition-2
  form: TNF via a `TRfunc()` sum of ≤3 rectangular ON windows read from an inert `Clock()`; no
  pre-equilibration (seed = the exact TR=0 steady state). Because the bngsim backend perturbs only
  *free* parameters via a `condition:` (and the TNF windows are not fit), each protocol/genotype is a
  **generated per-protocol model file** — `reduced_combination_<protocol>_{wt,a20ko}.bngl`, windows
  baked in, A20KO = A20-synthesis rule removed — joined in **one multi-model fit sharing the 13 free
  parameters** (`gen_combination_models.py` regenerates them from `reduced_combination.bngl`).
- **Target** — `reduced_combination_<protocol>_{wt,a20ko}.exp` (12 files): the original model's output
  at the S1-Table times, peak-normalized with the paper's ρ=0.03·max floor. Regenerable with
  `validation/gen_combination.py`.
- **Objective** — `norm_sos` + `normalization = peak`, the edition-2-native analog of the paper's
  Eq-7 geometric-mean-normalized log objective. Native-only.

## Result

- **Gate 3a** — the reduced model at Table-1 params reproduces the original model's Fig-2 dynamics,
  with **AMD\*_WT = 1.286 and AMD\*_A20KO = 1.158 — matching the paper's reported 1.29 / 1.16 to 3–4
  significant figures** (`reduced_combination_reproduction.png`, `combination_gate3a.png`). The
  point count is **exactly 914 independent** (matches the paper).
- **Gate 3b** — the 12-model joint `de` fit recovers **12/13 parameters within 3×**, most near-exact;
  only **ε** (the paper's single least-identifiable parameter, Fig 5) stays loose. The richer
  combination experiment identifies **δ and a_3** tightly — which the sparser on-off protocol could
  not — exactly as the paper argues.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"
cd pybnf-jobs/Jaruszewicz-Blonska-2023/reduced_combination
pybnf -c reduced_combination.conf                 # 12-model joint fit (Gate 3b)
BNGPATH=$BNGPATH python make_reproduction.py       # Fig-2 reproduction at Table-1 params (Gate 3a)
```
