# Cheemalavagu-2024 — JAK-STAT parameter fitting

PyBNF edition-2 job derived from:

> Cheemalavagu N, Baliban RC, Darville T, Lee REC. “Predicting gene level sensitivity
> to JAK-STAT signaling perturbation using a mechanistic-to-machine learning framework.”
> *Cell Systems* 2024; **15**:37-48.e4. DOI: 10.1016/j.cels.2023.12.002.

The authors fit a rule-based IL-6/IL-10–JAK–STAT–SOCS model to pooled pSTAT1 and
pSTAT3 time courses under six cytokine conditions using parallel tempering. The source
BNGL, experimental workbooks, and 2,067 retained parameter sets are available in
[`ncheemalavagu/STAT_models`](https://github.com/ncheemalavagu/STAT_models) at commit
`28874f0801077f479526157300f70fdccb672013`.

| slug | run cost | fits | flavor | data | verification status |
|---|---|---|---|---|---|
| [`jak_stat_fig2b`](jak_stat_fig2b/) | `hours` | 46 biological parameters plus two measurement scales; six IL-6/IL-10 conditions; pSTAT1 and pSTAT3 | ODE, quantitative, `chi_sq`, PEtab-v2-exportable | Figure 2B pooled workbooks; 84 observations | ✅ tier-1 PASS; PEtab round-trip PASS; executed verification notebook PASS; representative published ensemble χ² 74.5, median relative error 8.8%; bounded PyBNF fit PASS |

BioNetGen builds the 53-species/122-reaction network in about 0.1 s. The bounded fit is
an execution smoke test; a converged 48-dimensional fit still warrants parallel compute.
