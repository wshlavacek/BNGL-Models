# mito_camp — mitochondrial vs. plasma-membrane Gi signaling, ODE (PyBNF edition-2 job)

**Run cost: `hours`** — **measured at 28 min on 16 cores** = 7.5 core-hours for 39,000 evaluations (130 × 300 `de`) — the corpus's per-evaluation anchor, 0.69 s.

A PyBNF edition-2 parameter-fitting job that refits the **twelve rate constants of
Supplementary Data Table 1** — plus the PM/OMM MT1 concentration ratio and one parameter
that resolves a contradiction in the paper — of the compartmental cAMP model of:

> Suofu Y, Li W, Jean-Alphonse FG, Jia J, Khattar NK, Li J, *et al.* **"Dual role of
> mitochondria in producing melatonin and driving GPCR signaling to block cytochrome c
> release."** *Proc Natl Acad Sci USA* 2017; **114**:E7997–E8006.
> DOI: [10.1073/pnas.1705768114](https://doi.org/10.1073/pnas.1705768114) ·
> PMCID: [PMC5617291](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5617291/)
> — model = SI Appendix *Model constructions and analysis* + Supplementary Data Table 1
> (reactions) + Table 2 (initial quantities); fit targets = **Fig. 4H** and **Fig. 4I**;
> the authors' own fit of the same model to the same data = **Fig. 4J**; its posteriors =
> **Supplementary Data Fig. S7**.

Built with the `curate-pybnf-job` skill. The **library-model** sibling — the same model
with an actions block, for reference simulation rather than fitting — is
[`models/mitochondrial_mt1_camp_signaling_suofu2017/`](../../../models/mitochondrial_mt1_camp_signaling_suofu2017/),
whose `verify_suofu2017.ipynb` documents and re-runs every digitization used here.

> ✅ **This refit scores slightly better on the paper's own data than the paper's own
> plotted fit: SSE 4405.6 → 3684.6 over the same 77 points (RMSE 7.56 → 6.92 percentage
> points).** It reproduces the qualitative result the paper draws from those data — under
> melatonin the mitochondrial lumen loses more cAMP than the cytosol (21% vs. 41% of basal
> at 600 s), under DAMGO the order reverses (38% vs. 1%) — and it does so with a
> parameter set that is *stated in units*, which the published one is not.

> ✅ **It also settles which of the paper's two OMM receptor numbers is operative:
> Supplementary Data Table 2's.** `MT1_OMM_amp` comes out at **1.13**, with the whole
> 1%-objective band inside [1.00, 1.27]; the SI Appendix text's 20.68× higher areal-density
> reading is excluded by more than an order of magnitude.

## The model

Two aqueous compartments — cytosol (`CY`, 1000 µm³) and the pooled mitochondrial lumen
(`ML`, 11.31 µm³ = 100 spheres of radius 0.3 µm) — each running the same closed cycle:

- `AC + ATP <-> AC.ATP -> AC + cAMP` on the adjacent membrane (`PM`, `OMM`);
- `PDE + cAMP <-> PDE.cAMP -> PDE + AMP`, then `AMP -> ATP`;
- `MT1/µOR + Gi -> R.Gi -> R + Gi*`, gated by the agonist switches `MT1_isActive` /
  `OR_isActive`, and `AC + Gi* <-> AC.Gi*`, which blocks the substrate site;
- free nucleotide crossing the OMM, with the influx constant `phi = V_ml/V_cy` times the
  efflux constant so the compartments equilibrate at equal *concentration*.

ATP, cAMP and AMP are three states of one molecule type `Nuc`, which conserves the
adenine-nucleotide pool by construction. Luminal amounts are cytosolic amounts × `phi`, so
the two compartments have identical concentrations and identical intrinsic kinetics and
everything that distinguishes them is receptor placement plus diffusive exchange.

Network: **26 species, 30 reactions**, milliseconds per simulation. Not heavy.

## What is fit

Two time courses, one per agonist, both pre-equilibrated under the same receptor-silent
condition — which is the paper's own protocol ("the model was executed in the absence of
MT1 and µ-OR until steady state was reached"):

| experiment | preequilibrate | condition | receptor placement | data | points |
|---|---|---|---|---|---|
| `melatonin` | `basal` | `MT1_isActive = 1` | MT1 on PM **and** OMM | `melatonin.exp` (Fig. 4H) | 29 PM + 19 OMM |
| `damgo` | `basal` | `OR_isActive = 1` | µ-OR on PM only | `damgo.exp` (Fig. 4I) | 19 PM + 10 OMM |

The two fit observables are BNGL **functions** (so the `.exp` headers carry parentheses):
`camp_CY_pct()` and `camp_ML_pct()`, free cAMP as a percentage of the Supplementary Data
Table 2 basal quantity. No per-point uncertainty is digitizable, so the objective is plain
**`sos`** (PyBNF reports it as half the residual sum of squares).

**No scaling parameter.** The data are relative — percent of the pre-agonist maximum — but
the observables are absolute. Fixing the reference at Table 2's cAMP level, instead of
adding a free scale or `normalization = init`, keeps the job PEtab-exportable *and* imposes
the constraint Table 2 asserts: the receptor-silent steady state has to land on the
tabulated level, because the data's earliest points are near 100%. The fit gets to
**7.61e6** cytosolic against Table 2's 1.0e7, i.e. 24% low — the one place the model
visibly strains.

**Data provenance.** All 77 points are digitized from the plotted means of Fig. 4H and
Fig. 4I: page 5 of the article rendered at 1200 dpi, axis calibration from the plotted tick
marks, marker centres by matched filtering the colour mask with a 16 px disk and greedy
peak picking at a 22 px exclusion radius. Connected-component centroids do **not** work
here: the markers are ~36 px across and overlap both each other and the mono-exponential
line drawn through them. The extraction is reproducible in the model folder's
`verify_suofu2017.ipynb` — re-running it from the PDF returns the committed CSV exactly —
and is archived as `reference/suofu2017_fig4HI_camp_timecourses_digitized.csv` there.

## Why this job exists

**The published fit cannot be read back into a model.** Suofu et al. fit by
replica-exchange MCMC (eight chains, 10,000 swaps) and report only twelve log10 posterior
histograms (Fig. S7) — no values, no units. And the histograms are not a point estimate:
`k1f`, `k3f` and `kdiff` are narrow spikes pinned against the *lower edge* of their priors,
and `k3r`, `k4`, `k5` and `k6` are broadly bimodal across an entire prior box. (Panels 2
and 3 of Fig. S7 carry the same axis label, `log_k2`; one is presumably `log_k1r`.)

**The paper reports the OMM receptor density two ways that differ 20.7-fold.** Table 2
lists 1.13e4 MT1 on the OMM — exactly 1e6 × `phi`, the same *concentration* as a plasma
membrane carrying 1e6 copies. The SI Appendix text instead specifies "2068 MT1 per square
µm of OMM", the same *areal density* as that plasma membrane, which is 2.34e5 copies. The
difference decides how much more efficacious mitochondrial MT1 is than plasma-membrane
µ-OR, which is the paper's central quantitative claim. `MT1_OMM_amp` boxes the two
readings (1 = Table 2, 20.68 = the text) and lets the data choose; they choose Table 2.

## Free parameters (14)

Fitted values, and the range each covers among the 1268 archived sets within 1% of the best
objective — the cheap identifiability read on a model this sloppy:

| id | box | fit | 1% band | decades | role |
|---|---|---|---|---|---|
| `k1f` | 1e-11 – 1e-4 | **1.144e-05** | 1.0e-5 – 2.6e-5 | 0.4 | ATP + AC (µm³/molecule/s) |
| `k1r` | 1e-4 – 1e3 | **3.048e-03** | 1.0e-4 – 0.22 | 3.3 | AC.ATP dissociation (/s) |
| `k2` | 1e-4 – 1e3 | **0.2514** | 0.243 – 0.265 | 0.04 | cAMP production (/s) |
| `k3f` | 1e-10 – 1e-2 | **4.457e-05** | 4.35e-5 – 4.46e-5 | 0.01 | cAMP + PDE (µm³/molecule/s) |
| `k3r` | 1e-4 – 1e3 | **1.839e-03** | 1.1e-4 – 0.046 | 2.6 | PDE.cAMP dissociation (/s) |
| `k4` | 1e-3 – 1e4 | **3.306** | 2.6 – 12.5 | 0.7 | cAMP hydrolysis (/s) |
| `k5` | 1e-5 – 1e2 | **32.30** | 5.8 – 99.8 | 1.2 | AMP → ATP (/s) |
| `k6` | 1e-11 – 1e-2 | **1.123e-06** | 1.09e-6 – 1.77e-6 | 0.2 | receptor + Gi (µm³/molecule/s) |
| `k7` | 1e-3 – 1e4 | **4148** | 698 – 9991 | 1.2 | Gi activation and release (/s) |
| `k8f` | 1e-9 – 1e0 | **0.02986** | 0.016 – 0.042 | 0.4 | AC + Gi\* (µm³/molecule/s) |
| `k8r` | 1e-4 – 1e4 | **2.209e-04** | 1.0e-4 – 1.7e-3 | 1.2 | AC.Gi\* dissociation (/s) |
| `kdiff` | 1e-5 – 1e1 | **9.636e-03** | 9.2e-3 – 1.0e-2 | 0.05 | nucleotide efflux across the OMM (/s) |
| `ratio_MT1` | 1e-3 – 1 | **0.08222** | 0.058 – 0.087 | 0.2 | [MT1_PM] / (1e6 / V_cy) |
| `MT1_OMM_amp` | 1 – 20.68 | **1.130** | 1.00 – 1.27 | 0.10 | Table 2 (1) vs. SI text (20.68) |

Nine of the fourteen are pinned to better than half a decade; five (`k1r`, `k3r`, `k5`,
`k7`, `k8r`) are not identified at all by these data — every one of them is a step whose
rate only has to be *fast enough*, so the objective is flat above a threshold. The
authors' own posteriors are broad in overlapping places, which is the same message.

**What the fitted receptor densities mean.** `ratio_MT1 = 0.0822` and
`MT1_OMM_amp = 1.130` put 8.2e4 MT1 on the plasma membrane and 1.28e4 on the OMM — a
concentration ratio of 0.073, but areal densities of 170 and 113 per µm², i.e. within 1.5×
of each other. So the fit says what the SI Appendix text says qualitatively (comparable
areal density on the two membranes, and the mitochondrion far denser per unit volume), at
an absolute density ~20× below the "2068 per µm²" the text names.

## ✅ PEtab.v2-compliant · not heavy

`sos` over two plain `.exp` time courses with `condition:` perturbations and
`preequilibrate:` — entirely inside the PEtab-exportable subset. The round-trip (export →
`petab.v2` lint → import) is clean, and the exported bundle is committed under
[`petab/`](petab/): 2 observables, 14 estimated parameters, 77 measurements, 6 condition
rows over 2 experiments. Regenerate it with `make_petab.py`; do not hand-edit it.

## Verification (see [`VALIDATION.md`](VALIDATION.md))

- **Tier-1** (`scripts/check_conf.py`): edition 2, `job_type=de` resolves, data bound, 14
  free params bind by id, no `__FREE`. **PASS.**
- **PEtab.v2 round-trip** (`scripts/petab_roundtrip.py`): export → lint clean → import.
  **PASS.**
- **Real bngsim fit:** `pybnf -c mito_camp.conf` converges to a **finite** objective
  **1842.28** (`sos`, i.e. SSE 3684.6 over 77 points) in about 28 minutes on 16 cores
  (130 × 300 = 39,000 parameter sets, zero failed simulations).
- **Reproduction** (`mito_camp_reproduction.png`): SSE **4405.6 → 3684.6**, RMSE
  **7.56 → 6.92** percentage points against the same 77 points that the authors' Fig. 4J
  curves are scored on.
- **Cross-check of the two simulation paths:** `make_reproduction.py` drives BNG2.pl and
  reaches the basal state by a long unstimulated integration, where PyBNF uses bngsim's
  steady-state solve for `preequilibrate:`. They agree to the digit: SSE 3684.6 vs.
  2 × 1842.28 = 3684.56.
- **Independent check:** the library-model sibling at the same parameters agrees with a
  hand-written SciPy integration of Supplementary Data Table 1 to 3.8e-8
  (`verify_suofu2017.ipynb`).

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

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='mito_camp', conf='mito_camp.conf', simulator='ode',
    observables=('camp_CY_pct', 'camp_ML_pct'),
    system='Compartmental cAMP model of Gi-coupled receptors on the plasma membrane and '
           'the outer mitochondrial membrane (Suofu 2017, DOI 10.1073/pnas.1705768114, '
           'PMC5617291, Figs. 4H/4I/4J); 26-species ODE over cytosol + mitochondrial '
           'lumen with nucleotide permeation of the OMM; two agonists as condition: '
           'switches on receptor activity, both pre-equilibrated receptor-silent; sos on '
           'compartment-targeted FRET cAMP time courses digitized from Fig. 4H/4I. '
           'PEtab-exportable. GLOBAL (de) fit of the twelve rate constants of '
           'Supplementary Data Table 1 plus the PM/OMM MT1 ratio and the OMM receptor '
           'number, which the paper reports two incompatible ways; the published fit is '
           'given only as unit-less MCMC histograms and cannot be recovered. This fit '
           'beats the paper\'s own plotted curves on the paper\'s own data '
           '(SSE 4405.6 -> 3684.6). See VALIDATION.md.'),
    stochastic=False,
)
```
