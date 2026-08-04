# tnfr1_apo — apo-TNF sequential TNFR1 binding fit

**Run cost: `minutes`** — 600 evaluations (12 × 50 `de`) on a 24-reaction equilibrium network.

This PyBNF edition-2 job fits the three sequential microscopic dissociation
constants of soluble human TNFR1 binding to trimeric human TNF. It reproduces
the quantitative analysis behind the left panel of Fig. 1b in McMillan et al.
(2021), where ion-mobility mass spectrometry resolves TNF carrying zero, one,
two, or three receptors.

The model is a house-style reconstruction of the authors' `version_test.bngl`
from Supplementary Software 1. It retains the supplied rule semantics: three
explicit TNF sites, occupancy-dependent association rates, and one common
dissociation rate. The seven measured rows in `fig1b_apo.exp` are copied from
the publisher's Source Data workbook (`Fig 1b`, cells B3:G10); they are not
digitized.

## Fitting problem

| item | value |
|---|---|
| simulator | deterministic ODE, 9 species / 24 reactions |
| experiment | TNFR1 scan, 1–23 uM; TNF fixed at 5 uM; read at 10,000 s |
| outputs | fractions of TNF in the 0R, 1R, 2R, and 3R states |
| free parameters | `Kd_1`, `Kd_2`, `Kd_3` (nM) |
| published estimates | 0.01, 0.02, 0.22 nM |
| objective | sum of squared occupancy-fraction residuals |
| flavor | quantitative; PEtab v2 exportable |

A seeded bounded differential-evolution verification fit completed with a
finite objective of 0.07228 and recovered `Kd_1 = 0.0192`, `Kd_2 = 0.0297`, and
`Kd_3 = 0.345` nM. All three are within a factor of two of the rounded values
reported in Fig. 1b.

McMillan D, Martinez-Fleites C, Porter J, Fox D 3rd, Davis R, Mori P,
et al. (2021). "Structural insights into the disruption of TNF-TNFR1
signalling by small molecules stabilising a distorted TNF." *Nature
Communications* 12:582. DOI: 10.1038/s41467-020-20828-3.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"
cd pybnf-jobs/McMillan-2021/tnfr1_apo
pybnf -c tnfr1_apo.conf
```
