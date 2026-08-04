# `jak_stat_fig2b` — six-condition JAK-STAT fit

This job reconstructs the parameter-estimation problem behind Figure 2B of
Cheemalavagu et al. (2024): fit a mechanistic IL-6/IL-10 receptor–JAK1/JAK2–STAT1/STAT3
model with SOCS1/SOCS3 feedback to six 0–90 min cytokine time courses.

## Contents

| file | role |
|---|---|
| `jak_stat_fig2b.bngl` | fitting-ready 53-species/122-reaction ODE model; no simulation actions |
| `jak_stat_fig2b.conf` | edition-2 differential-evolution job, `chi_sq`, 48 free parameters |
| `*.exp` | six condition-specific pSTAT3/pSTAT1 targets with the paper's 15% uncertainty estimate |
| `fig2b_pSTAT_pooled.csv` | provenance table used to generate all `.exp` files |
| `author_ensemble1817_six_conditions.csv` | authors' six raw trajectories for the selected representative ensemble member |
| `make_reproduction.py` | regenerates the quantitative comparison and PNG |
| `jak_stat_fig2b_reproduction.png` | Figure 2B comparison |
| `VALIDATION.md` | provenance, adaptations, and gate evidence |

The biological parameter set is the paper's 46 unknowns. Two additional nuisance
parameters, `scale_pSTAT1` and `scale_pSTAT3`, map arbitrary model amounts to normalized
fluorescence. The original pTempEst likelihood recomputed those normalizers from every
proposal's IL-6 10 ng/mL value at 20 min. Because a cross-experiment dynamic normalizer is
not representable in core PEtab v2, this setup fits two explicit positive scales instead.
That preserves the intended channel-wise normalization while keeping the native config
fully PEtab-v2-exportable.

Bounds span the authors' 2,067 retained parameter sets and are expanded tenfold at both
ends. Nominal values use ensemble member 1817, selected as the minimum six-condition
chi-square member from the supplied trajectories—not as a new fit.

## Run

```bash
export BNGPATH=/path/to/folder/containing/BNG2.pl
cd pybnf-jobs/Cheemalavagu-2024/jak_stat_fig2b
pybnf -c jak_stat_fig2b.conf
python make_reproduction.py
```

For a short executable check, use a temporary copy with `population_size = 6`,
`max_iterations = 2`, and `refine = 0`. The model is not computationally heavy; a full
48-dimensional posterior or converged global fit nevertheless warrants parallel compute.
The bounded check completed with PyBNF v1.6.0 and objective `959.9372`; see
`VALIDATION.md` for the complete gate record.

## Proposed PyBNF corpus registration

```python
RealWorldExample(
    folder='Cheemalavagu-2024/jak_stat_fig2b',
    conf='jak_stat_fig2b.conf',
    simulator='ode',
    observables=('pSTAT3_norm', 'pSTAT1_norm'),
    system='Mouse macrophage IL-6/IL-10 JAK-STAT signaling '
           '(Cheemalavagu 2024, DOI 10.1016/j.cels.2023.12.002); ODE, six '
           'conditions, pSTAT1/pSTAT3 Figure 2B fit, 46 biological plus two '
           'measurement-scale parameters, chi_sq, PEtab-v2-exportable.',
)
```
