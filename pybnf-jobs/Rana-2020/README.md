# Rana-2020 — amyloid-β aggregation along competing pathways (PyBNF fitting jobs)

PyBNF edition-2 parameter-fitting jobs derived from one paper:

> Rana P, Bose P, Vaidya A, Rangachari V, Ghosh P.
> **"Global fitting and parameter identifiability for amyloid-β aggregation with competing
> pathways."**
> *2020 IEEE 20th International Conference on BioInformatics and BioEngineering (BIBE)*:73–78.
> DOI: [10.1109/BIBE50027.2020.00020](https://doi.org/10.1109/BIBE50027.2020.00020)

Built with the `curate-pybnf-job` skill. Each job below is a **self-contained folder** — its own
model, conf, data, reproduction figure, and README with the exact adaptations from the published
model, verification results, and a ready-to-paste `_manifest.py` snippet.

This paper is an unusually good fit for a PyBNF example collection because parameter fitting *is*
its subject: it fits the same competing-pathway model four different ways, runs a profile-likelihood
identifiability analysis on each, and benchmarks fourteen optimizers against the problem (Table II).
The four jobs here are those four fits.

## The shared model

Rana et al. specify an ensemble-kinetics (EKS) model of Aβ42 aggregation with two competing
pathways, completely, in equations:

| | reactions | fluxes |
|---|---|---|
| on pathway | `A_i + A_1 <-> A_i+1` (i = 1..11), `A_i + F <-> F` (i = 1..11), F ≡ A₁₂ | Eq. 4, H_i and I_i |
| off pathway | `4 A_1 + L <-> A'_4`, `A'_i + A_1 <-> A'_i+1` (i = 4..11), `A'_12 + A'_i <-> F'_1` (i = 4..11) | Eq. 4-II, G′₁, H′_i, I′_i |
| switching | `A_4 <-> A'_4` | Eq. 5, J |

L is a C12 fatty-acid pseudo-micelle and F′₁ the lumped, kinetically trapped off-pathway 12–23mer.
The thioflavin T signal is `signal_on + signal_off`, which the jobs express as a conf-side
measurement model `map_on*[F] + map_off*[F'_1]`.

The reconstruction each job uses is a fitting-ready copy of
[`models/amyloid_beta_competing_aggregation_pathways_rana2020`](../../models/amyloid_beta_competing_aggregation_pathways_rana2020),
whose notebook checks BioNetGen against an independent NumPy/SciPy transcription of those flux
equations over every state variable. **23 species, 80 reactions** (12 and 44 for the on-pathway
reduction).

### Two things the paper does not supply

1. **Five of the thirteen rate constants are in no column of Table I** — `k_fbon`, `k_fbon_`,
   `k_con_`, `k_fboff_`, `k_swi_`. Nor is the pseudo-micelle concentration.
2. **The tabulated on-pathway pair is incompatible with the reported peptide concentration.**
   At 25 µM Aβ42, `k_nuon = 22.04` /µM/h against `k_nuon_ = 12.72` /h gives
   K·[A₁] = 43, so every step of the nucleation chain runs downhill and A₁…A₁₂ equilibrates within
   minutes. The lag in Fig. 1a is the *sequential filling* of that chain, so once it is at
   equilibrium no choice of the two fibril-binding constants can restore a ~28 h lag. The curated
   model's notebook demonstrates this by simulation.

So Table I is a comparison target here, not an input. Every job reports its recovered values
beside the published column.

## The jobs

| slug | fits | data | paper's SSE | SSE here | status |
|---|---|---|---|---|---|
| [`on_pathway`](on_pathway/) | the 4 on-pathway constants of Eq. 4 + ThT scale (5 free) | Fig. 1a, 97 points | 0.13 | **0.022** | ✅ |
| [`micelle_addition`](micelle_addition/) | all 12 rate constants + 2 ThT scales (14 free), on→off switching | Fig. 1c, 322 points | 4.12 | **4.27** | ✅ |
| [`micelle_removal`](micelle_removal/) | all 12 rate constants + 2 ThT scales (14 free), off→on switching | Fig. 1b, 120 points | 1.22 | **0.33** | ✅ |
| [`global`](global/) | all 12 rate constants + 3 ThT scales (15 free), all five experiments at once | Figs. 1a+1b+1c, 427 points | 8.2 | **15.65** | ✅ |

All four are quantitative `sos` jobs and all four are **PEtab.v2-exportable**; ✅ means tier-1 +
PEtab round-trip + a real bngsim fit + a reproduction figure. "SSE here" is the committed model
scored through BioNetGen against the same digitized points — SSE is extensive in the number of
data points and the paper does not say how many it used, so treat these as same-order rather than
strictly comparable. Three of the four land at or below the published value; the global fit is
about twice the published SSE, and its own README says where the extra residual sits.

The jobs are chained the way the paper chains them: `on_pathway` first, its result bracketing the
on-pathway constants at 0.1–10× in the two switching jobs (Sec. IV-B's own recipe), and those two
bracketing the `global` fit (Sec. IV-D, "we chose the range of parameters carefully using the values
estimated from the previous fit").

## Two gotchas worth knowing

**`objective = sos` reports half the sum of squares.** The number PyBNF prints and stores is the
negative log-likelihood of unit-variance Gaussian noise up to the (n/2)·log(2π) constant — which is
how the run's AIC/BIC line is computed. On a three-point control (y = 1, 2, 3 against a zero
prediction) it reports 7.0, not 14.0. **Double it before comparing with the paper's SSE.** Every
number quoted in these READMEs is an SSE unless it says "PyBNF objective".

**A timed event is `preequilibrate` + `equil_t_end`.** Adding fatty acid at 3 h is expressed as a
fixed-length unmeasured phase under one condition followed by the measured phase under the other.
`equil_t_end` is what makes the first phase run for exactly 3 h rather than solving for a steady
state — this model has none on that timescale. The measured phase restarts the clock at zero, so
`make_exp.py` shifts the data by the event time and drops the pre-event points (which are
themselves another experiment in the set, so nothing is lost).

## Data provenance

Rana et al. tabulate nothing: the ThT traces exist only as plotted curves in Fig. 1. They were
recovered by
[`models/amyloid_beta_competing_aggregation_pathways_rana2020/digitize_rana2020.py`](../../models/amyloid_beta_competing_aggregation_pathways_rana2020/digitize_rana2020.py),
which rasterizes PDF page 4 at 600 dpi, locates each panel's axis box from the grey frame runs
(verified against the printed tick labels to within 4 px), and separates the curves by the direction
of the ink vector 255 − RGB — which is invariant to anti-aliasing coverage, unlike hue or
chromaticity. Values are binned to 0.5 h. The committed CSVs carry a per-point `spread_au`, the
half-height of the marker cloud in that bin, which is the dominant uncertainty (median 0.019–0.054
ThT a.u. depending on the panel). `make_exp.py` turns the CSVs into `.exp` tables.

**Panel assignment.** The published caption labels Fig. 1b "micelle addition event" and Fig. 1c
"micelle removal event", but the body text cites Fig. 1c for the addition fit ("the three
experiments", SSE 4.12, Sec. IV-B) and Fig. 1b for the removal fit (SSE 1.22, Sec. IV-C). The
kinetics settle it in favour of the body text: every trace in panel (b) rises without a lag, which
only a sample containing fatty acid from t = 0 can do, while panel (c) holds two traces that stay
flat until 3 h and 24 h — the signature of fatty acid arriving at a lag-phase sample. The jobs
follow the body text.

## Source files

- Paper: `dev/papers/Rana-2020/` (primary source, not redistributed)
- Original ThT measurements: Ghosh P, Rana P, Rangachari V, Saha J, Steen E, Vaidya A (2020).
  *A game theoretic approach to deciphering the dynamics of amyloid-β aggregation along competing
  pathways.* R Soc Open Sci 7(8):191814. doi:[10.1098/rsos.191814](https://doi.org/10.1098/rsos.191814)
- Curated model collection: [`models/amyloid_beta_competing_aggregation_pathways_rana2020`](../../models/amyloid_beta_competing_aggregation_pathways_rana2020)
