# `nfkb_tnfa` — NF-κB module of the innate immune response, refit to the TNF-α data

PyBNF edition-2 fitting job derived from:

> Korwek Z, Czerkies M, Jaruszewicz-Błońska J, Prus W, Kosiuk I, Kochańczyk M, Lipniacki T.
> **"Nonself RNA rewires IFN-β signaling: A mathematical model of the innate immune response."**
> *Sci Signal* 2023; **16**(815):eabq1173.
> DOI: [10.1126/scisignal.abq1173](https://doi.org/10.1126/scisignal.abq1173)

## What it fits, and why this problem

Korwek et al. write that TNF-α "is used to calibrate kinetic rate constants in the NF-κB pathway
module that follows Lipniacki et al. (2004) and Tay et al. (2010)." This job reproduces that
calibration step: it refits the **fifteen NF-κB-module rate constants marked FITTED in table S1**
against the TNF-α data of **fig. S12**.

TNF-α reaches the NF-κB module alone. Under it the poly(I:C), IRF3, IFN-β and STAT1/2 modules stay
quiescent, so the protocol isolates one module *inside* the complete 53-species model rather than
requiring a reduced one — the full network still integrates, the other four modules simply sit at
rest.

| | |
|---|---|
| model | `nfkb_tnfa.bngl` — 53 species, 96 reactions, deterministic ODE |
| free parameters | 15 kinetic + 9 per-blot normalization constants = 24 |
| data | 111 points, 4 blots × 2 replicate blots, digitized from figs. S12A and S12B |
| objective | `chi_sq` (per-point σ = 25% of the measurement) |
| design | pre-equilibrate unstimulated → TNF-α bolus → 6 h time course; WT and A20 KO arms |
| flavor | quantitative, **PEtab.v2-exportable** |

## Files

| file | role |
|---|---|
| `nfkb_tnfa.bngl` | the model, fitting-ready: no simulation actions, `generate_network` retained |
| `nfkb_tnfa.conf` | the edition-2 job setup, banner-commented |
| `nfkb_tnfa_wt_nuclear_r{1,2}.exp` | fig. S12A — nuclear NF-κB (RelA), WT, 9 lanes over 6 h |
| `nfkb_tnfa_wt_fine_r{1,2}.exp` | fig. S12B — p-IKK / IκBα / A20, WT, 10 lanes over 3 h |
| `nfkb_tnfa_wt_long_r{1,2}.exp` | fig. S12B — p-IKK / IκBα / A20, WT, 8 lanes over 6 h |
| `nfkb_tnfa_a20ko_r{1,2}.exp` | fig. S12B — p-IKK / IκBα, A20 KO, 8 lanes over 6 h |
| `make_exp.py` | rebuilds the `.exp` files from the committed digitized CSV |
| `nfkb_tnfa_bestfit.txt` | the recovered best fit (PyBNF `sorted_params_refine_final.txt`, top row) |
| `make_reproduction.py` | scores published vs fitted and writes the reproduction figure |
| `nfkb_tnfa_reproduction.png` | that figure |

Run it (with `BNGPATH` set, from this folder):

```sh
pybnf -c nfkb_tnfa.conf
python make_reproduction.py            # uses the committed best fit; --results for a fresh one
```

## Adaptations from the published model

`nfkb_tnfa.bngl` is a fitting-ready copy of
[`models/innate_immune_response_korwek2023`](../../../models/innate_immune_response_korwek2023),
whose primary file reproduces the authors' data file S2 **bit-for-bit** (max |difference| = 0 over
43 observables × 4321 time points). Four differences from that file:

1. **No simulation actions.** Only `generate_network({overwrite=>1})` is retained; the protocol
   comes from the conf.
2. **`h_A20_gene`** gates TNFAIP3 transcription, in the idiom of the authors' own `h_Mavs` /
   `h_Pkr_gene` / `h_Rnasel_gene`. Data file S2 has no such switch, but the paper simulates A20 KO
   cells in fig. S12B; this makes that arm a `condition:`. It is 1, and inert, in the WT arm.
3. **Nine `s_*` normalization constants.** fig. S12 reports a "relative level" — each protein on
   each blot normalized to a loading control and then to an arbitrary reference. That constant is
   a property of the blot, not the biology, so one is estimated per blot per protein through an
   `observable: ..., formula: s * <obs>` measurement model. Their nominal values are the constants
   recovered by fitting the curated model to each panel's *published model curve*, so the model at
   nominal values reproduces the published figure's vertical scale directly. fig. S12B's own two
   wild-type panels justify the per-blot treatment: they are two blots of one experiment and their
   published model curves differ by a near-constant factor.
4. **The stimulus is a species bolus.** `condition: tnfa, perturbations: "TNFa()" = 1` — a target
   containing `(` routes to `setConcentration` rather than `setParameter` (#474). Setting a
   parameter that a seed species references would *not* change the species amount carried over
   from the equilibration phase; the measured phase would run with no TNF-α and every parameter
   set would score identically.

## Data provenance

Korwek et al. tabulate nothing: the blot quantifications and the fitted model's trajectories exist
only as markers in the supplementary figures.
[`digitize_korwek2023.py`](../../../models/innate_immune_response_korwek2023/digitize_korwek2023.py)
extracts each figure's bitmap from the PDF at native resolution (~230 ppi), calibrates every log
panel on its own axis tick marks (≈53 px per decade ⇒ ≈4% per pixel), separates the three series
by colour — the model's orange open circles, the reproduced blot's dark grey dots, the replicate's
light grey dots — and flags markers that sit on the axis floor or are clipped at the top.
`make_exp.py` writes a flagged marker as `NaN`: its value is not readable from the figure, so it
is not a measurement. The A20 panel of the A20 KO blot contributes no column at all — A20 is
absent in those cells and the figure draws a flat placeholder at 1.

**Per-point σ is a flat 25% of the measurement**, which makes the objective a relative-error least
squares — the multiplicative-error criterion Korwek et al. use throughout. The 25% is the paper's
own number: table 1 reports an average multiplicative error between Western-blot replicates of
1.24 over all 2915 points they fit. The two replicates digitized here scatter more (1.55 over the
35 lanes where both resolve separately), but that figure is biased upward and is not used in its
place — markers only separate into two dots when the replicates *disagree*, so the lanes where the
blots agree are exactly the ones that cannot be measured.

## Verification

| gate | result |
|---|---|
| tier-1 (`check_conf.py`) | **PASS** — edition 2, `job_type=de` resolves, data bound, 24 free params, no `__FREE` |
| PEtab.v2 round-trip (`petab_roundtrip.py`) | **PASS** — export → `petab.v2` lint clean → import |
| real bngsim fit | **PASS** — 300 DE generations (population 48) + 300 simplex refinement iterations, 22 min, χ² 164.89; AIC 271.16, BIC 336.19 (k=24, n=111) |
| reproduction vs the published fit | see below |

### Reproduction

`make_reproduction.py` re-runs the model through BioNetGen under the conf's protocol at three
parameter sets and scores each against the same 111 points:

| parameter set | χ² | median \|relative error\| |
|---|---|---|
| published (table S1 values, nominal blot scales) | 1203.1 | 35.5% |
| published kinetics, nine blot scales refit | 437.6 | 35.0% |
| fitted | 329.8 | 29.2% |

The middle row is the honest comparison. The nominal blot scales were fitted to the *published
model curve*, not to the data, so most of the 1203 → 330 gap is normalization rather than kinetics;
with the normalization taken out, the refit improves on the published kinetics by **1.33×**. That
the improvement is small is the point: the published parameterization is already close to the best
this data supports, which is what one should expect, because the authors fit these same TNF-α blots
as part of a much larger dataset while this job sees only them.

**Parameter recovery** — thirteen of the fifteen kinetic constants land within 3.7× of the
published maximum-likelihood values, and eight within a factor of two, from a fit that sees only
figure-digitized data:

| parameter | fitted | published | ratio |
|---|---|---|---|
| `a_Tak1_by_Tnfa` | 2.91e-1 | 8.90e-1 | 0.33× |
| `a_Ikk` | 1.91e-2 | 1.09e-2 | 1.75× |
| `d_Ikk_1` | 8.15e-4 | 8.83e-4 | 0.92× |
| `d_Ikk_2` | 1.33e-1 | 5.98e-2 | 2.23× |
| `d_Ikk_3` | 3.40e-4 | 2.68e-4 | 1.27× |
| `p_Ikba_by_Ikk` | 8.97e-3 | 8.46e-3 | 1.06× |
| `g_Ikba_u_with_Nfkb` | 2.98e-5 | 3.94e-5 | 0.75× |
| `i_Ikba` | 6.97e-5 | 5.64e-4 | **0.12×** |
| `e_Ikba` | 4.51e-2 | 2.05e-2 | 2.20× |
| `tg_Ikba_mrna` | 3.96e-4 | 2.57e-4 | 1.54× |
| `a_Ikba_gene_by_Nfkb__` | 9.15e-3 | 6.85e-2 | **0.13×** |
| `s_Ikba` | 1.15e-2 | 4.26e-3 | 2.71× |
| `tg_A20_mrna` | 6.39e-4 | 9.41e-4 | 0.68× |
| `a_A20_gene_by_Nfkb__` | 1.35e-2 | 2.57e-2 | 0.53× |
| `sg_A20` | 1.83e-5 | 5.01e-6 | 3.65× |

The two outliers, `i_Ikba` and `a_Ikba_gene_by_Nfkb__`, both move by about 8× *in the same
direction* and are a compensating pair by construction: `i_Ikba` sets how much free IκBα reaches
the nucleus, and `a_Ikba_gene_by_Nfkb__` sets how strongly NF-κB competes against that nuclear IκBα
at the NFKBIA promoter — the two appear in the promoter term as a ratio, and only their combination
is identifiable from a single stimulus. Korwek et al. resolved that degeneracy with protocols this
job does not include; the paper's own identifiability analysis (fig. S15) is over all 38 parameters
and all 2915 points, not this subset.

## `_manifest.py` entry

```python
RealWorldExample(
    folder='Korwek-2023/nfkb_tnfa', conf='nfkb_tnfa.conf', simulator='ode',
    observables=('NFkBn_nuc', 'pIKK_fine', 'IkBa_fine', 'A20_fine',
                 'pIKK_long', 'IkBa_long', 'A20_long', 'pIKK_ko', 'IkBa_ko'),
    system='NF-kB module of the innate immune response to nonself RNA, refit to TNF-alpha '
           'time courses in WT and TNFAIP3 (A20) KO A549 cells (Korwek 2023, '
           'doi:10.1126/scisignal.abq1173, figs. S12A/S12B). Deterministic ODE, 53 species; '
           'pre-equilibration + species bolus; nine per-blot normalization constants estimated '
           'through observable formulas.'),
```

Coverage-matrix row for `examples/real-world/README.md`:

```
| [`Korwek-2023/nfkb_tnfa`](Korwek-2023/nfkb_tnfa/) | NF-kB module vs TNF-alpha, WT + A20 KO (Korwek 2023, fig. S12) | **ODE** | pre-equilibration with `equil_t_end:`, species-bolus condition, `observable:` measurement models, replicate `.exp` files, `chi_sq` | ✅ |
```

## Notes for whoever touches this next

- **`equil_t_end:` is required, not optional.** The model is stiff — `k_FAST` is 1 /s while
  `sg_PROTEIN` is 1.5e-5 /s — and the ODE default (a steady-state solve) makes CVODE exhaust its
  step budget, `mxstep steps taken before reaching tout`, scoring the parameter sets it lands on
  as `inf`. A fixed 1.2e6 s interval integrates cleanly. It is long enough: the slowest approach
  to rest is A20's, whose `sg_A20` gives a 2.2-day time constant, and by 1.2e6 s A20 is within
  0.2% of its 30-day value — far inside the data's ≈4% precision.
- **`ave_norm_sos` cannot be used with this data.** It would be the natural column-normalized
  objective and it *is* PEtab-exportable, but it normalizes by `np.average` rather than
  `np.nanmean` (`pybnf/objective.py:1702`), so a single `NaN` in a column makes every parameter
  set score `inf`. These columns necessarily carry `NaN`s. `norm_sos` works but is not exportable
  (its relative sigma has no PEtab v2 shape in PyBNF 1.6.0; `export_job` refuses it —
  ADR-0021/0023, lanl/PyBNF#423). `chi_sq` with the `_SD` columns gets both.
- **A species-pattern perturbation must be quoted** in the conf: `"TNFa()" = 1`. Unquoted, the
  grammar reads it as a parameter identifier and the line fails to parse.
