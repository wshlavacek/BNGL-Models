# Kozer-2013 — EGFR higher-order oligomerisation & phosphorylation

PyBNF edition-2 parameter-fitting jobs for the EGFR higher-order-clustering system of Kozer et al.
(2013), in **two model forms** fit to the **same 2013 data**:

> Kozer N, Barua D, Orchard S, Nice EC, Burgess AW, Hlavacek WS, Clayton AHA.
> **"Exploring higher-order EGFR oligomerisation and phosphorylation — a combined experimental
> and theoretical approach."** *Mol BioSyst* 2013; **9**(7):1849–1863.
> PMCID: [PMC3698845](https://pmc.ncbi.nlm.nih.gov/articles/PMC3698845/) ·
> DOI: [10.1039/c3mb70073a](https://doi.org/10.1039/c3mb70073a) · PMID: 23629589

Built with the `curate-pybnf-job` skill and validated against the primary sources with
`validate-pybnf-job` (each slug carries a committed `VALIDATION.md`).

## The shared biochemistry, two model forms

Both slugs fit the **same biology** — EGF binds EGFR; EGFR ectodomains crosslink into higher-order
oligomers (with **ring-closure** reactions); crosslinked cytosolic tails undergo an activating
conformational change and *trans*-phosphorylate. Binding/crosslinking constants are fixed from the
literature (Macdonald & Pike 2008, Elleman 2001, Low-Nam 2011; Kozer 2013 Table 1); **five rate
constants** (`k_o`, `k_c`, `kaf`, `kar`, `chi_r`) and **four observable-scale factors**
(`alpha1_pre`..`alpha4_pre`) are fit.

- [`egfr_ode`](egfr_ode/) — **Kozer 2013** model (File S1): oligomers capped at **tetramers**, full
  reaction network generated → **deterministic ODE**.
- [`egfr_nf`](egfr_nf/) — **Mitra 2019** PyBioNetFit "example 2": an **unbounded**, network-free
  reformulation of the same 2013 model → **NFsim**. Fit to the same 2013 data (the ODE↔NFsim
  demonstration twin of `egfr_ode`; the ODE↔NF equivalence is confirmed in its `VALIDATION.md`).

The fit recipes come from the PyBNF 2019 paper corpus (Mitra et al., *iScience* 2019;
19:1012–1036; BioNetFit 1 problem 2), re-expressed on the edition-2 surface.

> A companion study — Kozer et al., *Biochemistry* 2014; 53:2594–2604 (PMC4010257) — extends the
> 2013 model with Grb2 recruitment (SI `bi500182x_si_001.txt`, a capped-ODE + Grb2 model). It is
> **not** reproduced here. Note that `egfr_nf` is *not* that 2014 SI model (a common mislabel — see
> its `VALIDATION.md` Gate 2).

## The jobs

| slug | model | fit data | simulator | fit | status |
|---|---|---|---|---|---|
| [`egfr_ode`](egfr_ode/) | Kozer **2013** ODE, tetramer cap | cluster density + phospho-EGFR (Figs 2B/2D/3B/3D), 9 params | **ODE** | 5 sloppy kinetics + 4 scale factors, `chi_sq` | ✅ validated (89/100) · 🔶 heavy |
| [`egfr_nf`](egfr_nf/) | **Mitra 2019** unbounded network-free | same 2013 data as `egfr_ode`, 9 params | **NFsim** | 5 kinetics + 4 scale factors, `chi_sq` | ✅ validated (86/100) · 🔶 heavy |

Both are **PEtab-v2-exportable** and **heavy** (the ODE network is cluster-scale to generate; NFsim
is cluster-scale to run) — a full fit is a cluster job; the tier-1 parse and PEtab round-trip are
backend-free.

### ⚠️ A network-generation gotcha this pair exposed

`egfr_ode`'s reaction network is finite **only** under a `max_stoich` cap — the rules alone do not
bound oligomer size. The curated model therefore **retains** one actions-block line,
`generate_network({overwrite=>1,max_stoich=>{EGF=>4,EGFR=>4}})`; without it, pybnf's synthesized
default (`generate_network({overwrite=>1})`) drops the cap and generation never terminates.
Network-definition directives stay in the model; only simulation actions move to the conf.
(`egfr_nf` is network-free and correctly carries no actions block.)

## Source materials

| file (untracked `dev/papers/Kozer-2013/`) | content | → |
|---|---|---|
| File S1 of Kozer et al. (**2013**) | ODE EGFR-clustering model + the preprocessed fit data (Figs 2B/2D/3B/3D) | `egfr_ode.bngl`, both `.exp` |
| Mitra 2019 PyBioNetFit `example 2` | unbounded network-free NFsim model | `egfr_nf.bngl` |

The 2013 model was edited by Brandon R. Thomas and William S. Hlavacek for BioNetFit 1, further
edited by Eshan D. Mitra for PyBioNetFit, and upgraded to the edition-2 surface here.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Kozer-2013/<slug>
pybnf -c <slug>.conf
```

Both are cluster-scale; raise `population_size`/`max_iterations` and run on a cluster to reproduce
the published fits. See each slug's `README.md` and `VALIDATION.md`.
