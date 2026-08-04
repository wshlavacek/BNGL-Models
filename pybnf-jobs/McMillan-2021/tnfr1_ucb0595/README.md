# tnfr1_ucb0595 — UCB-0595-distorted TNF sequential TNFR1 binding fit

**Run cost: `minutes`** — 600 evaluations (12 × 50 `de`) on a 24-reaction equilibrium network.

This PyBNF edition-2 job fits the three sequential microscopic dissociation
constants of soluble human TNFR1 binding to UCB-0595-bound trimeric human TNF.
It reproduces the quantitative analysis behind the right panel of Fig. 1b in
McMillan et al. (2021), where the compound leaves the first two binding events
tight but weakens the third by roughly five orders of magnitude.

The model is a house-style reconstruction of the authors' `version_test.bngl`
from Supplementary Software 1. It differs from the apo job only in the nominal
and fitted dissociation constants. The seven measured rows in
`fig1b_ucb0595.exp` are copied from the publisher's Source Data workbook
(`Fig 1b`, cells B3:C10 and H3:K10); they are not digitized.

## Fitting problem

| item | value |
|---|---|
| simulator | deterministic ODE, 9 species / 24 reactions |
| experiment | TNFR1 scan, 1–23 uM; TNF fixed at 5 uM; read at 10,000 s |
| outputs | fractions of TNF in the 0R, 1R, 2R, and 3R states |
| free parameters | `Kd_1`, `Kd_2`, `Kd_3` (nM) |
| published estimates | 0.04, 0.19, 9612 nM |
| objective | sum of squared occupancy-fraction residuals |
| flavor | quantitative; PEtab v2 exportable |

A seeded bounded differential-evolution verification fit completed with a
finite objective of 0.06284 and recovered `Kd_1 = 0.0210`, `Kd_2 = 0.0901`, and
`Kd_3 = 9433` nM. The first two estimates are within a factor of 2.2 of the
rounded values in Fig. 1b, and the third differs from 9612 nM by 1.9%.

McMillan D, Martinez-Fleites C, Porter J, Fox D 3rd, Davis R, Mori P,
et al. (2021). "Structural insights into the disruption of TNF-TNFR1
signalling by small molecules stabilising a distorted TNF." *Nature
Communications* 12:582. DOI: 10.1038/s41467-020-20828-3.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"
cd pybnf-jobs/McMillan-2021/tnfr1_ucb0595
pybnf -c tnfr1_ucb0595.conf
```
