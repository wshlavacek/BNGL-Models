# igf1r — IGF1 / IGF1R competition dose-response, deterministic ODE (PyBNF edition-2 job)

A PyBNF edition-2, parameter-fitting job setup. The **fitting problem** — a 3-parameter fit of
the two site affinities and the cross-site cooperativity to the Fig-5B competition curve — is from:

> Erickson KE, Rukhlenko OS, Shahinuzzaman M, Slavkova KP, Kholodenko BN, Hlavacek WS, et al.
> **"Modeling cell line-specific recruitment of signaling proteins to the insulin-like growth
> factor 1 receptor."**
> *PLoS Comput Biol* 2019; **15**(1):e1006706.
> PMCID: [PMC6353226](https://pmc.ncbi.nlm.nih.gov/articles/PMC6353226/) ·
> DOI: [10.1371/journal.pcbi.1006706](https://doi.org/10.1371/journal.pcbi.1006706)
> *(parameter estimation performed with BioNetFit 1)*

The underlying **model and data** are from:

> Kiselyov VV, Versteyhe S, Gauguin L, De Meyts P.
> **"Harmonic oscillator model of the insulin and IGF1 receptors' allosteric binding and
> activation."** *Mol Syst Biol* 2009; **5**:243.
> PMID: [19225456](https://pubmed.ncbi.nlm.nih.gov/19225456/) ·
> DOI: [10.1038/msb.2008.78](https://doi.org/10.1038/msb.2008.78) (K1/K2 = Table-1 site KDs;
> data = the Fig-5B radioligand competition curve).

> **Citation note.** The *job* (this 3-affinity fit of the normalized F5B curve) is Erickson
> et al. 2019 — citing Kiselyov alone would mis-attribute the fitting problem. The PyBNF 2019
> corpus problem **15-igf1r** (Mitra et al., *iScience* 2019) fits the *fuller* 7-rate-constant
> model to three datasets; **this example is the distilled 3-affinity, F5B-only fit** carried by
> the classic PyBNF `examples/igf1r/`.

## The model

A **radioligand competition assay**. IGF1R is a pre-formed, disulfide-bonded **dimer**; each IGF1
ligand can engage **Site 1** and **Site 2** (with a steric constraint) and can **crosslink** the
two sites *within* a dimer (avidity). Labelled ("hot") IGF1 is held fixed at 7 pM (= 8852
copies/cell, the seed default); unlabelled ("cold") IGF1 is titrated, and the readout is
hot-ligand binding vs. cold dose. Deterministic **ODE**. **Finite network without a cap** —
crosslinking is confined to the dimer, so the classic model used a bare
`generate_network({overwrite=>1})` and generation terminates; no network directive needs to be
retained (contrast `Kozer-2013-2014/egfr_ode`).

## What is fit

| dataset | observable | design | source |
|---|---|---|---|
| hot-ligand binding vs. cold-IGF1 dose (18 doses, 1.3×10⁻¹²–6.6×10⁻⁷ M) | `IGF1_hot_bound` (normalized to no-competitor) | ODE competition **parameter scan** over `IGF1_cold_conc`, each dose to steady state (t_end=14400 s) | Kiselyov 2009 Fig. 5B |

`F5B.exp` carries per-point `_SD` columns → **`chi_sq`**. Three free params: `K1`, `K2`, `K1prime`.

## Files

| file | role |
|---|---|
| `igf1r.bngl` | edition-2, fitting-ready ODE model (finite network; **no** actions block / `generate_network` directive needed) |
| `igf1r.conf` | the edition-2 job setup (`normalization = init`, `chi_sq`, `ss` + Simplex refine) |
| `F5B.exp` | fit target — IGF1_hot_bound vs. cold dose, pre-normalized, with SD |
| `make_reproduction.py` | regenerates the reproduction PNG (ODE `parameter_scan` at the PyBNF best-fit vs. data) |
| `igf1r_reproduction.png` | verification figure |

## ⚠️ NATIVE-ONLY (not PEtab-exportable)

`normalization = init` normalizes the simulated `IGF1_hot_bound` to its no-competitor
(first-scan-row) value — a PyBNF-native prediction transform PEtab v2 cannot express
(`export.py:798-827`). `export_job` raises `NotImplementedError`. This is a **property of the
flavor, not a defect**: the data are reported relative to the no-competitor baseline, so the init
normalization is intrinsic to the problem. (The shipped `igf1r.conf` in the PyBNF corpus is the
canonical cautionary example of a non-exportable `normalization` job — see
`skills/curate-pybnf-job/references/petab-compliance.md`.)

## ⚠️ Nominals are Kiselyov's KDs — they do NOT fit the normalized F5B curve

Unlike the GenFit-starting-point nominals elsewhere in this corpus, the model's nominal `K1`/`K2`
**are** the Kiselyov 2009 Table-1 site KDs. But under this reduced 3-parameter, init-normalized
setup they place the competition midpoint ~1–2 orders of magnitude to the *right* of the data —
which is exactly why Erickson **fit** the parameters. The **PyBNF best-fit** (this repo's Scatter
Search + Simplex run of `igf1r.conf`, `chi_sq = 1.08`) is used for the figure:

| id | nominal (Kiselyov Table 1) | **PyBNF best-fit** | role |
|---|---|---|---|
| `K1` | 9.2×10⁻⁹ M | **5.42×10⁻⁹ M** | Site-1 dissociation constant |
| `K2` | 4.83×10⁻⁷ M | **8.37×10⁻³ M** | Site-2 dissociation constant (→ weak; binding is Site-1 + avidity driven) |
| `K1prime` | 0.1 | **1.44×10⁻⁸** | cross-site cooperativity (→ strong intradimer crosslinking / avidity) |

> **Sloppy / non-identifiable.** The competition curve pins down the *apparent* affinity
> (IC₅₀ ≈ 1.3×10⁻¹⁰ M), not the individual constants: a second basin (`K1≈2.5×10⁻⁹`,
> `K1prime≈1.7×10⁻³`, `K2` at the 10⁻¹⁵ lower bound) sits at nearly-equal objective (1.101 vs
> 1.083). Both describe "tight Site-1 binding + strong avidity." A `recover` assertion should use
> a loose tolerance or be omitted; the *reproduction* (curve through the data) is the meaningful
> check.

## Changes from the source model (documented in `igf1r.bngl`)

- **No actions block** — the 7 verbose commented-out dissociation-experiment protocols of the
  classic model are dropped; the competition scan is synthesized from `igf1r.conf`. A concise
  protocol note is kept at the file's end.
- Fitted constants (`K1`, `K2`, `K1prime`) bind to model ids **by name** (ADR-0034, no `__FREE`).
  The reaction network is otherwise byte-identical to the classic
  `IGF1R_Model_receptor_activation_bnf.bngl`.

## Verification

- **Tier-1** (`scripts/check_conf.py`): edition 2, `job_type=ss` resolves, the experiment binds
  data, `model.stochastic = False` (→ `simulator='ode'`), 3 free params bind by id, no `__FREE`.
  **PASS.**
- **Native-only guard** (instead of PEtab round-trip): `pybnf.petab.export_job` raises
  **`NotImplementedError`** ("This job normalizes the whole fit ('init') … PEtab v2 cannot
  express"). **PASS** (asserted; the export-refused guard keeps the job native-only).
- **Bounded bngsim fit reaches a finite objective:** the full `igf1r.conf` (Scatter Search,
  pop 12 × 50 it + 20-step Simplex refine) converges to **chi_sq = 1.083** — a finite objective;
  the simulate→score→propose loop ran end-to-end.
- **Reproduction** (`igf1r_reproduction.png`, ODE `parameter_scan` at the best-fit,
  init-normalized): the competition curve passes through **every** Fig-5B point within its SD —
  **median 2.5 % / max 14.5 % relative error** over the 8 doses with y > 0.05; **chi_sq ≈ 1.10**
  (BNG2.pl vs. the fit's bngsim; the near-zero high-cold tail dominates the max rel err).

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Erickson-2019/igf1r
pybnf -c igf1r.conf                   # ODE; a full fit runs on a workstation in a few minutes
python make_reproduction.py           # reproduction figure at the PyBNF best-fit
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='igf1r', conf='igf1r.conf', simulator='ode',
    observables=('IGF1_hot_bound',),
    system='IGF1/IGF1R competition dose-response (Erickson 2019 fit, PMC6353226; Kiselyov 2009 '
           'model+data, PMID 19225456, Fig 5B); ODE parameter_scan, chi_sq, normalization=init '
           '-> NATIVE-ONLY (not PEtab-exportable)'),
    # NATIVE-ONLY: assert export_job raises NotImplementedError (do NOT assert PEtab round-trip).
    # Sloppy fit (apparent affinity identifiable, individual K's not) -> no recover assertion,
    # or recover={'K1': 5.42e-9} with a loose tol. PyBNF best-fit chi_sq = 1.08.
```
