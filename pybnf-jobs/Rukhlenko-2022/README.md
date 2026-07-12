# Rukhlenko-2022 — cSTAR cell-state control (PyBNF fitting jobs)

PyBNF edition-2 parameter-fitting jobs derived from one paper:

> Rukhlenko OS, Halász M, Rauch N, Zhernovkov V, Prince T, Wynne K, Maher S,
> Kashdan E, MacLeod K, Carragher NO, Kolch W, Kholodenko BN.
> **"Control of cell state transitions."** *Nature* 2022; **609**(7929):975–985.
> PMCID: [PMC9644236](https://pmc.ncbi.nlm.nih.gov/articles/PMC9644236/) ·
> DOI: [10.1038/s41586-022-05194-y](https://doi.org/10.1038/s41586-022-05194-y)

Built with the `curate-pybnf-job` skill. Each job below is a **self-contained folder** — its
own model, conf, data/constraints, reproduction figure, and README with the exact
adaptations from the published model, verification results, and a ready-to-paste
`_manifest.py` snippet. The paper's core signalling model was itself fit with
**pyBioNetFit**, so these are natural real-world examples.

## The shared framework

All jobs come from the paper's **cSTAR** (cell State Transition Assessment and Regulation)
core signalling network and its **DPD** (Dynamic Phenotype Descriptor `S`) — a signed
cell-state coordinate on a bistable, Waddington-like landscape. Unlike a single-model paper,
each job reconstructs its **own cell-line-specific model** from the authors' repository
([github.com/OleksiiR/cSTAR_Nature](https://github.com/OleksiiR/cSTAR_Nature)); they share
the framework, idioms, and native-only fold-change flavor rather than one `.bngl`.

## The jobs

| slug | fits | flavor | data source | status |
|---|---|---|---|---|
| [`cstar_trka`](cstar_trka/) | TrkA/NGF 7-phospho fold changes at 0/10/45 min (8 params) | quantitative, **native-only** (`normalization=init`) | authors' RPPA (repo), **Fig. 4A/4B** | ✅ validated **86/100** ([`VALIDATION.md`](cstar_trka/VALIDATION.md)); 3a ~18 % median err |
| [`cstar_trkb`](cstar_trkb/) | TrkB/BDNF 7-phospho fold changes at 0/10/45 min (8 params) | quantitative, **native-only** | authors' RPPA (repo), **Fig. 4A/4B** | ✅ validated **84/100** ([`VALIDATION.md`](cstar_trkb/VALIDATION.md)); 3a ~27 % median err |
| [`cstar_skmel133`](cstar_skmel133/) | SKMEL-133 kinase-inhibitor panel, 9 readouts × 4 perturbations at 24 h (8 params) | quantitative, **native-only**, `chi_sq` | authors' pyBNF `.exp` (repo; Korkut 2015 RPPA), fit fig **ED 17A** | ✅ validated **88/100** ([`VALIDATION.md`](cstar_skmel133/VALIDATION.md)); 3a ~13 % median err |
| [`cstar_skmel133_bpsl`](cstar_skmel133_bpsl/) | DPD / Waddington-landscape bistable switch, IRS+AKT synergy (4 params) | **BPSL** constraints, **native-only** | paper's qualitative landscape claims | ✅ tier-1 + **`check` 8/8 satisfied** |

`cstar_trka` and `cstar_trkb` are ligand time-course twins (neuroblastoma, pre-equilibrate →
ligand, 45 min); `cstar_skmel133` is a different flavor — a steady-state inhibitor
perturbation panel on a RAF-inhibitor-resistant melanoma line; `cstar_skmel133_bpsl` is its
constraint-bearing sibling, expressing the paper's headline *qualitative* claim (combined
IRS + AKT inhibition crosses the landscape barrier into the death basin — a switch neither
single agent achieves) as PyBNF BPSL `.prop` constraints on the DPD `Sval`.

All four jobs are **native-only** (not PEtab-v2-exportable): the three quantitative fits
because the data are fold changes vs. baseline (`normalization = init`), and the BPSL job
because BPSL constraints have no PEtab representation. Each is verified accordingly (bounded
fit + figure reproduction; and `job_type = check` satisfaction for the BPSL job) rather than
with a PEtab round-trip.

## Fit scope: demonstration fits (current) vs. real-world constrained fits (planned)

The three quantitative slugs above are **reduced demonstration fits**: each frees ~8 hand-picked
parameters over ±1-decade bounds against a data subset. They are built on the paper's **authentic
model and data** (validated byte-identical to the authors' own files) and reproduce the paper's
figures **at the published parameters** (Gate 3a), but a from-scratch run of the demo `.conf` does
**not** recover the published parameters — the reduced, unconstrained problem is sloppy /
non-identifiable (each slug's `VALIDATION.md` Gate 3b).

The paper's *actual* fit (Methods pp.22, 24) freed the **full** parameter set under **BMRA-derived
connection-coefficient confidence intervals imposed as inequality constraints**, across the **full
joint dataset**, with scatter search + simplex. Reconstructing that — a **real-world constrained
fit** — is planned as **additional slugs** (the BMRA posteriors ship in the authors' repo under
`BMRA/results/` and `BMRA/results_SKMEL_133/`, so the constraints are recoverable; PyBNF expresses
them as BPSL `.prop` constraints, cf. the existing `cstar_skmel133_bpsl` slug). The demo fits are
kept as-is (clearly labelled) — they remain useful small, runnable examples.

## Source materials

All models and data are from the authors' repository
[github.com/OleksiiR/cSTAR_Nature](https://github.com/OleksiiR/cSTAR_Nature); no BioModels
deposition.

| path | used by |
|---|---|
| `Trk_AB_models/TrkA_S_model.bngl`, `TrkB_S_model.bngl` | `cstar_trka`, `cstar_trkb` |
| `RPPA_DA/RPPA_data_trusted.csv`, `RPPA_data_Trk_normalized_new.csv` | Trk fold-change targets |
| `SKMEL-133_models/SKMEL-133-3.bngl`, `SKMEL-133_preproc/*.exp` | `cstar_skmel133`, `cstar_skmel133_bpsl` |

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Rukhlenko-2022/<slug>
pybnf -c <slug>.conf
```

See each slug's `README.md` for its full write-up.
