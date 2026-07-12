# Kozer-2013-2014 — EGFR higher-order oligomerisation & phosphorylation (PyBNF fitting jobs)

PyBNF edition-2 parameter-fitting jobs derived from **two companion papers** on the same
EGFR clustering system. The **ODE** model comes from the 2013 study; the **network-free
NFsim** model comes from the 2014 companion:

> **[2013 — ODE model + data]** Kozer N, Barua D, Orchard S, Nice EC, Burgess AW,
> Hlavacek WS, Clayton AHA. **"Exploring higher-order EGFR oligomerisation and
> phosphorylation — a combined experimental and theoretical approach."**
> *Mol BioSyst* 2013; **9**(7):1849–1863.
> PMCID: [PMC3698845](https://pmc.ncbi.nlm.nih.gov/articles/PMC3698845/) ·
> DOI: [10.1039/c3mb70073a](https://doi.org/10.1039/c3mb70073a) · PMID: 23629589
>
> **[2014 — network-free model]** Kozer N, Barua D, Henderson C, Nice EC, Burgess AW,
> Hlavacek WS, Clayton AHA. **"Recruitment of the Adaptor Protein Grb2 to EGFR Tetramers."**
> *Biochemistry* 2014; **53**(16):2594–2604.
> PMCID: [PMC4010257](https://pmc.ncbi.nlm.nih.gov/articles/PMC4010257/) ·
> DOI: [10.1021/bi500182x](https://doi.org/10.1021/bi500182x)

Built with the `curate-pybnf-job` skill. Each job below is a **self-contained folder** — its
own model, conf, data, reproduction figure, and README with the exact adaptations, verification
results, and a ready-to-paste `_manifest.py` snippet. Both are part of the PyBNF 2019 paper
corpus (Mitra et al., *iScience* 2019, 19:1012–1036; BioNetFit 1 problem 2), re-expressed on the
edition-2 surface.

## The shared biochemistry, two model forms

Both jobs fit the **same biology** — EGF binds EGFR; EGFR ectodomains crosslink into
higher-order oligomers (with **ring-closure** reactions); crosslinked cytosolic tails undergo
an activating conformational change and *trans*-phosphorylate. Binding/crosslinking constants
are fixed from the literature (Macdonald & Pike 2008, Elleman 2001, Low-Nam 2011; Table 1 of
Kozer 2013); **five rate constants** (`k_o`, `k_c`, `kaf`, `kar`, `chi_r`) and **four
observable-scale factors** (`alpha1_pre`..`alpha4_pre`) are fit. They are the two model *forms*
the two papers used:

- [`egfr_ode`](egfr_ode/) — **Kozer 2013**, File S1: oligomers capped at **tetramers**, full
  reaction network generated → **deterministic ODE**.
- [`egfr_nf`](egfr_nf/) — **Kozer 2014**, SI `bi500182x_si_001.txt`: **unbounded** oligomers
  → **network-free NFsim**. Fit to the same 2013 data.

## The jobs

| slug | model paper | fits | simulator | flavor | status |
|---|---|---|---|---|---|
| [`egfr_ode`](egfr_ode/) | Kozer **2013** | cluster density + phospho-EGFR, time course (30 nM) + EGF dose-response (9 params) | **ODE** (network-generating) | quantitative, **PEtab-exportable**, `chi_sq` | ✅ tier-1 + PEtab round-trip; 🔶 heavy (cluster-scale network) |
| [`egfr_nf`](egfr_nf/) | Kozer **2014** | same four datasets, network-free (9 params) | **NF** (NFsim) | quantitative, **PEtab-exportable**, `chi_sq` | ✅ tier-1 + PEtab round-trip; 🔶 heavy (NFsim) |

Both are **quantitative and PEtab-v2-exportable** (unlike the native-only `Kirsch-2020` /
`Rukhlenko-2022` fold-change jobs). Both are **heavy** — the ODE network is cluster-scale to
generate and NFsim is cluster-scale to run — so a full fit is a cluster job; the tier-1 parse and
PEtab round-trip are backend-free.

### ⚠️ A network-generation gotcha this pair exposed

`egfr_ode`'s reaction network is finite **only** under a `max_stoich` cap — the rules alone do
not bound oligomer size. The curated model therefore **retains** exactly one actions-block line,
`generate_network({overwrite=>1,max_stoich=>{EGF=>4,EGFR=>4}})`; without it, pybnf's synthesized
default (`generate_network({overwrite=>1})`) drops the cap and generation never terminates. This
is the principle now baked into the `curate-pybnf-job` skill: **network-definition directives stay
in the model; only simulation actions move to the conf.** (`egfr_nf` is network-free and correctly
carries no actions block.)

## Source materials

| file | content | → |
|---|---|---|
| File S1 of Kozer et al. (**2013**) ([DOI](https://doi.org/10.1039/c3mb70073a)) | ODE BioNetGen model + the preprocessed fit data (Figs. 2B/2D/3B/3D) | `egfr_ode.bngl`, both `.exp` |
| SI `bi500182x_si_001.txt` of Kozer et al. (**2014**) ([DOI](https://doi.org/10.1021/bi500182x), [PMC4010257](https://pmc.ncbi.nlm.nih.gov/articles/PMC4010257/)) | network-free BioNetGen model | `egfr_nf.bngl` |

Both models were edited by Brandon R. Thomas and William S. Hlavacek for BioNetFit 1, further
edited by Eshan D. Mitra for PyBioNetFit, and upgraded to the edition-2 surface here.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Kozer-2013-2014/<slug>
pybnf -c <slug>.conf
```

Both are cluster-scale; raise `population_size`/`max_iterations` and run on a cluster to
reproduce the published fits. See each slug's `README.md` for its full write-up.
