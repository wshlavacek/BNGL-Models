# McMillan-2021 — TNF-TNFR1 sequential binding fits

Two quantitative PyBNF edition-2 jobs reconstruct the BioNetGen/SciPy analysis
used for Fig. 1b of McMillan et al. (2021). Both fit equilibrium IMS-MS
occupancy data for trimeric human TNF binding zero to three soluble human TNFR1
molecules; they differ in whether TNF is apo or stabilized in a distorted
conformation by UCB-0595.

> McMillan D, Martinez-Fleites C, Porter J, Fox D 3rd, Davis R, Mori P,
> et al. (2021). "Structural insights into the disruption of TNF-TNFR1
> signalling by small molecules stabilising a distorted TNF." *Nature
> Communications* 12:582. DOI: 10.1038/s41467-020-20828-3.

The model structure comes from Supplementary Software 1 (`version_test.bngl`),
and the fit targets are copied exactly from the publisher's official Source
Data workbook, sheet `Fig 1b`. No curve digitization is used.

| job | run cost | fits | flavor | data source | verification target |
|---|---|---|---|---|---|
| [`tnfr1_apo`](tnfr1_apo/) | `minutes` | `Kd_1`, `Kd_2`, `Kd_3` for apo TNF | quantitative, PEtab v2 | Fig. 1b left | 0.01, 0.02, 0.22 nM |
| [`tnfr1_ucb0595`](tnfr1_ucb0595/) | `minutes` | `Kd_1`, `Kd_2`, `Kd_3` for distorted TNF | quantitative, PEtab v2 | Fig. 1b right | 0.04, 0.19, 9612 nM |

Both jobs use a 9-species, 24-reaction ODE network and a seven-dose parameter
scan over the initial TNFR1 concentration. Each native `.conf` is intended to
pass parse validation, PEtab v2 export/lint/import, and a bounded bngsim fit.

## Verification

Both native configurations pass the edition-2 configuration check and a PEtab
v2 export, lint, and import round trip. Seeded bounded differential-evolution
fits reached finite objectives and recovered the reported affinities within a
factor of 2.2:

| job | bounded-fit objective | fitted `Kd_1`, `Kd_2`, `Kd_3` (nM) |
|---|---:|---|
| `tnfr1_apo` | 0.07228 | 0.0192, 0.0297, 0.345 |
| `tnfr1_ucb0595` | 0.06284 | 0.0210, 0.0901, 9433 |

The jobs retain the larger search budgets used for this verification; optimizer
outputs are transient and are not part of the example.
