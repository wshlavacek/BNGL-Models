# Erickson-2019 — IGF1 / IGF1R receptor activation (PyBNF fitting job)

A PyBNF edition-2 parameter-fitting job. The **fitting problem** — a reduced 3-parameter fit of
the IGF1R site affinities to a radioligand-competition curve — is from Erickson et al. (2019);
the underlying **model and data** are from Kiselyov et al. (2009):

> **[fitting problem / job source]** Erickson KE, Rukhlenko OS, Shahinuzzaman M, Slavkova KP,
> Kholodenko BN, Hlavacek WS, et al. **"Modeling cell line-specific recruitment of signaling
> proteins to the insulin-like growth factor 1 receptor."**
> *PLoS Comput Biol* 2019; **15**(1):e1006706.
> PMCID: [PMC6353226](https://pmc.ncbi.nlm.nih.gov/articles/PMC6353226/) ·
> DOI: [10.1371/journal.pcbi.1006706](https://doi.org/10.1371/journal.pcbi.1006706)
>
> **[model + data origin]** Kiselyov VV, Versteyhe S, Gauguin L, De Meyts P.
> **"Harmonic oscillator model of the insulin and IGF1 receptors' allosteric binding and
> activation."** *Mol Syst Biol* 2009; **5**:243.
> PMID: [19225456](https://pubmed.ncbi.nlm.nih.gov/19225456/) ·
> DOI: [10.1038/msb.2008.78](https://doi.org/10.1038/msb.2008.78)

Built with the `curate-pybnf-job` skill. The job below is a **self-contained folder** — its own
model, conf, data, reproduction figure, and README with the exact adaptations, verification
results, and a ready-to-paste `_manifest.py` snippet.

> **Why cite Erickson, not just Kiselyov?** The *model and competition data* originate with
> Kiselyov 2009, but the *fitting problem* carried by PyBNF's `examples/igf1r/` — fit only the
> two site affinities `K1`, `K2` and the cross-site cooperativity `K1prime` to the Fig-5B curve,
> normalized to the no-competitor baseline — was set up by **Erickson et al. 2019** (with
> BioNetFit 1). Citing Kiselyov alone would mis-attribute the job. (Erickson's *full* paper fits
> a larger rate-constant network to competition + dissociation data; the PyBNF 2019 corpus
> problem **15-igf1r**, Mitra et al. *iScience* 2019, uses that fuller 7-parameter version. This
> example is the distilled 3-affinity, F5B-only fit.)

## The biochemistry

A **radioligand competition assay** on the IGF1 receptor. IGF1R is a pre-formed disulfide-bonded
**dimer** with two ligand sites (Site 1, Site 2) that a single IGF1 can **crosslink** *within*
the dimer (avidity). Labelled ("hot") IGF1 is held fixed while unlabelled ("cold") IGF1 is
titrated; the readout is hot-ligand binding vs. cold dose. Three equilibrium parameters are fit.

## The job

| slug | fits | simulator | flavor | status |
|---|---|---|---|---|
| [`igf1r`](igf1r/) | IGF1_hot_bound vs. cold-IGF1 competition dose-response, 18 doses (3 params: K1, K2, K1prime) | **ODE** (finite network, no cap) | quantitative + `normalization=init` → **NATIVE-ONLY** (not PEtab-exportable), `chi_sq` | ✅ tier-1 + export-refused guard + bounded fit (chi_sq 1.08) + reproduction (median 2.5 % rel err) |

**Native-only** (like the `Kirsch-2020` / `Rukhlenko-2022` fold-change jobs, and unlike the
PEtab-exportable `Monine-2010` / `Kozer-2013-2014` jobs) — here because of `normalization = init`,
the canonical non-PEtab-exportable feature. Verified by asserting `export_job` raises
`NotImplementedError`, **not** with a PEtab round-trip. Lightweight ODE — a full fit runs on a
workstation in a few minutes.

## Two caveats worth knowing

1. **The nominals do not fit.** The model's nominal `K1`/`K2` are Kiselyov's Table-1 site KDs,
   but under this reduced, init-normalized setup they miss the F5B curve by ~1–2 orders of
   magnitude — which is *why* Erickson fit the parameters. The reproduction uses the PyBNF
   best-fit.
2. **The fit is sloppy.** The competition curve determines the *apparent* affinity, not the
   individual `K1`/`K2`/`K1prime`; two parameter basins sit at nearly-equal objective. See the
   slug [`README`](igf1r/README.md).

## Source materials

- **Model / data / classic conf:** PyBNF `examples/igf1r/`
  (`IGF1R_Model_receptor_activation_bnf.bngl`, `igf1r.conf`, `F5B.exp`) and its edition-2 twin
  `examples/real-world/igf1r/`.
- **Related (different parameterization):** RuleHub `Published/Mitra2019/15-igf1r/` fits the
  fuller 7-rate-constant model to F5B + F5D_20min + F5D_60min — **not** directly comparable to
  this 3-affinity fit.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Erickson-2019/igf1r
pybnf -c igf1r.conf
python make_reproduction.py
```
