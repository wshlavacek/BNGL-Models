# Suofu-2017 — mitochondrial MT1 receptor control of cAMP (PyBNF fitting job)

PyBNF edition-2 parameter-fitting job derived from the compartmental cAMP model of:

> **Suofu Y, Li W, Jean-Alphonse FG, Jia J, Khattar NK, Li J, Baranov SV, Leronni D,
> Mihalik AC, He Y, Cecon E, Wehbi VL, Kim J, Heath BE, Baranova OV, Wang X, Gable MJ,
> Kretz ES, Di Benedetto G, Lezon TR, Ferrando LM, Larkin TM, Sullivan M, Yablonska S,
> Wang J, Minnigh MB, Guillaumet G, Suzenet F, Richardson RM, Poloyac SM, Stolz DB,
> Jockers R, Witt-Enderby PA, Carlisle DL, Vilardaga J-P, Friedlander RM.**
> **"Dual role of mitochondria in producing melatonin and driving GPCR signaling to block
> cytochrome c release."** *Proc Natl Acad Sci USA* 2017; **114**:E7997–E8006.
> DOI: [10.1073/pnas.1705768114](https://doi.org/10.1073/pnas.1705768114).
> PMCID: [PMC5617291](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5617291/).
> — the compartmental BioNetGen model of the SI Appendix section *Model constructions and
> analysis*, its reactions (Supplementary Data Table 1) and initial quantities
> (Supplementary Data Table 2), the fit data (Fig. 4H, Fig. 4I), the authors' own fit
> (Fig. 4J) and its posteriors (Supplementary Data Fig. S7).

Built with the `curate-pybnf-job` skill. The **library-model** sibling — the same model
with an actions block, for reference simulation rather than fitting — lives in
[`models/mitochondrial_mt1_camp_signaling_suofu2017/`](../../models/mitochondrial_mt1_camp_signaling_suofu2017/),
together with the verification notebook that documents and re-runs every digitization
used here.

## The model

The paper's claim is that the melatonin type-1 receptor (MT1) is a *mitochondrial* GPCR:
it sits in the outer mitochondrial membrane (OMM) with its ligand-binding domain facing
the cytosol, and its Gi/adenylyl-cyclase machinery is inside. To test whether that
placement can explain the compartmental cAMP measurements, the authors built a small
compartmental BioNetGen model — a 1 pL spherical cell holding 100 spherical mitochondria
of radius 0.3 µm, with

- an ATP → cAMP → AMP → ATP cycle running independently in the cytosol and in the
  mitochondrial lumen, driven by adenylyl cyclase on the adjacent membrane and cleared by
  phosphodiesterase;
- nucleotides, and only nucleotides, permeating the OMM;
- agonist-occupied receptor catalytically converting Gi to Gi\*, and Gi\* sequestering
  adenylyl cyclase in an inactive complex;
- luminal amounts set to cytosolic amounts times the volume ratio, so the two compartments
  start with identical concentrations and identical intrinsic kinetics.

Everything that then distinguishes the two compartments is receptor placement plus
diffusive exchange, and because the cytosol is ~88× the luminal volume that exchange is
one-sided. Melatonin (MT1 on both membranes) drives luminal cAMP *below* cytosolic;
DAMGO (µ-opioid receptor, plasma membrane only) reverses the order. Reproducing that
reversal is the point of the model.

Network generation gives **26 species and 30 reactions** — milliseconds per simulation.

## The jobs

| slug | run cost | fits | flavor | data | status |
|---|---|---|---|---|---|
| [`mito_camp`](mito_camp/) | `hours` | the twelve rate constants of Supplementary Data Table 1, the PM/OMM MT1 concentration ratio, and `MT1_OMM_amp` (which of the paper's two contradictory OMM receptor densities is operative) | quantitative, `sos`, **PEtab-exportable** | 77 points digitized from Fig. 4H (melatonin) and Fig. 4I (DAMGO), two compartment-targeted FRET sensors each | ✅ tier-1 + PEtab round-trip + fit + reproduction · **65/100** ([VALIDATION](mito_camp/VALIDATION.md)) |

> ✅ **The fit scores better on the paper's own data than the paper's own plotted fit:
> SSE 4405.6 → 3684.6 over the same 77 points (RMSE 7.56 → 6.92 percentage points), and it
> settles the OMM receptor density in favour of Supplementary Data Table 2.**

## Why this job exists

> ⚠️ **The published fit is not recoverable from the paper.** Suofu et al. fit the model
> by replica-exchange MCMC and report the result only as twelve log10 posterior histograms
> (Supplementary Data Fig. S7). No table of values is given, and no unit convention — so
> the histograms cannot be mapped onto a runnable model. Three of the twelve posteriors
> (`k1f`, `k3f`, `kdiff`) are narrow spikes pinned against the *lower edge* of their
> priors, and four more (`k3r`, `k4`, `k5`, `k6`) are broadly bimodal across an entire
> prior box, so the published fit is not a point estimate either. Panels 2 and 3 of
> Fig. S7 are both labeled `log_k2`; one of them is presumably `log_k1r`.

> ⚠️ **The paper reports the OMM receptor density two ways that differ 20.7-fold.**
> Supplementary Data Table 2 lists 1.13e4 MT1 on the OMM — exactly 1e6 × (V_lumen/V_cell),
> i.e. the same *concentration* as a plasma membrane carrying 1e6 copies. The SI Appendix
> text instead specifies "2068 MT1 per square µm of OMM", which is the same *areal
> density* as that plasma membrane and comes to 2.34e5 copies. Which one is operative
> decides how much more efficacious mitochondrial MT1 is than plasma-membrane µ-OR, which
> is the paper's central quantitative claim, so this job makes it a bounded free parameter
> (`MT1_OMM_amp` ∈ [1, 20.68]) and lets the data choose.

## Source materials

- **Primary paper:** `dev/papers/Soufu2017/PNAS.pdf` (note the folder name transposes the
  first author's name; the correct spelling is Suofu).
- **Supplement:** `dev/papers/Soufu2017/pnas.1705768114.sapp.pdf` — Supplementary Data
  Table 1 (reactions), Table 2 (initial quantities), Fig. S7 (posteriors), and the SI
  Appendix section *Model constructions and analysis*.
- **Authors' model files:** not deposited. The model was reconstructed from the prose
  specification plus the two tables.
- **Fit data:** Fig. 4H and Fig. 4I of the primary paper, digitized here; the extraction
  is committed in the library-model folder's `verify_suofu2017.ipynb` and archived as
  `reference/suofu2017_fig4HI_camp_timecourses_digitized.csv`.
- **Reference for the authors' fit:** Fig. 4J, digitized as
  `reference/suofu2017_fig4J_model_curves_digitized.csv`.

Not built (optional future slugs): a `pt` sampler job that would reproduce the authors'
replica-exchange posteriors rather than a point estimate — worth doing, since the
identifiability of this model is the interesting question, but it needs a prior
specification the paper only draws.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
```

```bash
cd pybnf-jobs/Suofu-2017/mito_camp && pybnf -c mito_camp.conf
```

```bash
cd pybnf-jobs/Suofu-2017/mito_camp && python make_reproduction.py
```

`make_reproduction.py` simulates through BNG2.pl, so the figure reproduces without the
fitting toolchain.
