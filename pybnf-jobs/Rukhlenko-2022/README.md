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

## ⚠️ Reproducibility caveat — the paper's fit *recipe* is not public

The authors' repository ([OleksiiR/cSTAR_Nature](https://github.com/OleksiiR/cSTAR_Nature)) ships the
**models** (`.bngl`, with the final fitted parameter values), the **data** (`.exp`/RPPA), and the
**BMRA MATLAB code plus its outputs** (the connection-coefficient means/stds that define the CIs).
It does **not** ship the **pyBioNetFit fitting setup**: there are no `.conf`/`.prop`/`.con` files, the
`.bngl` carry no free-parameter or constraint blocks, and — critically — **there is no record of how
the Eq. 14 connection-coefficient constraints were expressed** for the fit. The paper (Methods
pp. 22, 24) describes the method in words but publishes no parameter table or fit configuration.

**Consequence:** the `*_bmra*` slugs below are our **best reconstruction** of the paper's *described*
constrained fit — not a byte-exact reproduction of it. We can (and do) reproduce the models at the
authors' published parameters, the data, and the BMRA numbers; what we reconstruct is the fitting
recipe — which parameters were freed, and how each model's Eq. 14 `r_ij` was turned into a
pyBioNetFit-checkable constraint (we derive each `r_ij` in closed form from the model's rate laws;
see `cstar_*_bmra_exact/build_constraints.py` and `VALIDATION.md`). An email to the authors requesting
their `.conf`/constraint files is outstanding; if received, it would replace this reconstruction.
*(This kind of under-reporting is common — most papers don't publish everything needed to re-run a
fit. Among these jobs it is most acute for Rukhlenko-2022; contrast Erickson-2019/igf1r, where the
authors shipped their actual BioNetFit `.conf`.)*

**Canonical vs. superseded.** For each system we keep the **`*_bmra_exact`** slug as the closest
reconstruction (it constrains the actual Eq. 14 coefficients) and mark the earlier **`*_bmra`**
sign-approximation slugs **DEPRECATED** (retained for reference / as the robust sign-only fallback).

## The jobs

| slug | fits | flavor | data source | status |
|---|---|---|---|---|
| [`cstar_trka`](cstar_trka/) | TrkA/NGF 7-phospho fold changes at 0/10/45 min (8 params) | quantitative, **native-only** (`normalization=init`) | authors' RPPA (repo), **Fig. 4A/4B** | ✅ validated **86/100** ([`VALIDATION.md`](cstar_trka/VALIDATION.md)); 3a ~18 % median err |
| [`cstar_trkb`](cstar_trkb/) | TrkB/BDNF 7-phospho fold changes at 0/10/45 min (8 params) | quantitative, **native-only** | authors' RPPA (repo), **Fig. 4A/4B** | ✅ validated **84/100** ([`VALIDATION.md`](cstar_trkb/VALIDATION.md)); 3a ~27 % median err |
| [`cstar_skmel133`](cstar_skmel133/) | SKMEL-133 kinase-inhibitor panel, 9 readouts × 4 perturbations at 24 h (8 params) | quantitative, **native-only**, `chi_sq` | authors' pyBNF `.exp` (repo; Korkut 2015 RPPA), fit fig **ED 17A** | ✅ validated **88/100** ([`VALIDATION.md`](cstar_skmel133/VALIDATION.md)); 3a ~13 % median err |
| [`cstar_skmel133_bpsl`](cstar_skmel133_bpsl/) | DPD / Waddington-landscape bistable switch, IRS+AKT synergy (4 params) | **BPSL** constraints, **native-only** | paper's qualitative landscape claims | ✅ tier-1 + **`check` 8/8 satisfied** |
| [`cstar_skmel133_bmra`](cstar_skmel133_bmra/) | ⚠️ **DEPRECATED** (superseded by `cstar_skmel133_bmra_exact`) — the SKMEL-133 constrained fit with only the **sign** of each connection (23 params); a robust approximation of the paper's Eq. 14 constraint | quantitative **+ BPSL** (data fusion), **native-only**, `chi_sq` | authors' pyBNF `.exp` + **BMRA posteriors** (`BMRA/results_SKMEL_133/*_withMyc.csv`) | ⚠️ superseded; still valid: validated **90/100**, `check` **23/23**, 3a ~11 % |
| [`cstar_skmel133_bmra_exact`](cstar_skmel133_bmra_exact/) | the SKMEL-133 fit with the **EXACT Eq.14 `r_ij`** constrained two-sided inside the BMRA CIs (43 params: g + K per edge + 3 β) — the paper's real constraint object | quantitative **+ BPSL** (two-sided numeric bands), **native-only**, `chi_sq` | authors' 6-inhibitor pyBNF `.exp` + **Table S10** (`*_withMyc.csv`) | ✅ validated **90/100** ([`VALIDATION.md`](cstar_skmel133_bmra_exact/VALIDATION.md)); `check` **40/44** (band genuinely binding), 3a ~13 %, anchored fit → 42/44 |
| [`cstar_trkab_bmra`](cstar_trkab_bmra/) | ⚠️ **DEPRECATED** (superseded by `cstar_trkab_bmra_exact`) — the joint TrkA+TrkB constrained fit with only the **sign** of each connection (38 params); a robust approximation of the paper's Eq. 14 constraint | quantitative **+ BPSL** (2 models, data fusion), **native-only**, `sos` | authors' RPPA time courses + **BMRA posteriors** (`BMRA/results/Trk{A,B}_*_10_*.csv`) | ⚠️ superseded; still valid: validated **88/100**, `check` **11/11**, network 82/405 + 114/645 |
| [`cstar_trkab_bmra_exact`](cstar_trkab_bmra_exact/) | the joint TrkA+TrkB fit with the **EXACT Eq.14 `r_ij`** constrained two-sided inside the 10-min BMRA CIs (51 params: shared core + g/K per constrained edge) | quantitative **+ BPSL** (2 models, two-sided numeric bands), **native-only**, `sos` | authors' Trk time courses + **Table S5** (10-min `*_10_*.csv`) | ✅ validated **87/100** ([`VALIDATION.md`](cstar_trkab_bmra_exact/VALIDATION.md)); `check` **40/46**, 3a ~23 %, bounded fit → 43/46 (brings the tightest-CI edges in) |

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
connection-coefficient confidence intervals imposed as inequality constraints**, with scatter
search + simplex. Reconstructing that — a **real-world constrained fit** — is delivered as
**additional `*_bmra` slugs** built beside the demos (which are kept as-is, clearly labelled):

- **[`cstar_skmel133_bmra`](cstar_skmel133_bmra/) — ⚠️ DEPRECATED (superseded by `_exact`).** The full 23-parameter
  connection-coefficient set fit to the six single-drug fold-change arms under the BMRA-inferred
  CIs (from `BMRA/results_SKMEL_133/SKMEL133_{rm,rs}_log_200_5K_withMyc.csv`) as BPSL sign
  constraints, grounded in the paper's Eqs. 14/24/25. `job_type = check` → **Satisfied 23/23** at
  the published parameters; a constrained fit **keeps** every BMRA sign where the unconstrained
  demo flips `g_IRSERK`.
- **[`cstar_skmel133_bmra_exact`](cstar_skmel133_bmra_exact/) — BUILT (90/100).** The same
  full-parameter SKMEL-133 fit, but imposing the paper's *exact* constraint object rather than the
  sign: the model's **Eq.14 connection coefficients `r_ij`** (emitted as functions via the
  triple-verified `r_ij = C_i·L_ij` decomposition) pinned **two-sided** inside each BMRA CI
  `[mean−std, mean+std]` (Table S10, with_MYC). 43 free params (g **and** K per edge + 3 β) fit to
  all six single-drug inhibitor panels (11 arms; PKC/mTOR/CDK applied by competitive-binding
  `setConcentration`, #474). Unlike the sign approximation (trivially 23/23), the exact ±1 std band
  is **genuinely binding**: `job_type = check` → **Satisfied 40/44** at the published parameters (4
  edges — 3 AKT-incoming + MYC←CDK — sit just outside), and a proliferative-state–anchored bounded
  fit reaches **42/44** while lowering the objective, bringing 3 of the 4 into band. See its
  `VALIDATION.md` for the derivation, the triple verification, and the anchor.
- **[`cstar_trkab_bmra`](cstar_trkab_bmra/) — ⚠️ DEPRECATED (superseded by `_exact`).** The joint TrkA+TrkB fit: both
  models (82/405 and 114/645 reactions) fit together, the full connection-coefficient set of each
  receptor (38 params, shared core bound by name) under its own 10-min BMRA CIs
  (`BMRA/results/Trk{A,B}_{rm,rs}_10_*.csv`) as BPSL sign constraints, applying the paper's own
  statistical-significance filter (11 constrained; low-confidence connections left free).
  `job_type = check` → **Satisfied 11/11**; a joint constrained fit reaches 26.1 < 35.1 keeping
  every sign. (Trains on the ligand time courses like the demos; the inhibitor-panel arms — ED
  Figs. 6–8 — would enlarge the training set without changing the methodology.)

- **[`cstar_trkab_bmra_exact`](cstar_trkab_bmra_exact/) — BUILT (87/100).** The same joint TrkA+TrkB
  fit, imposing the paper's *exact* constraint object rather than the sign: each model's **Eq.14
  connection coefficients `r_ij`** (emitted as functions via the verified `r_ij = C_i·L_ij`
  decomposition, plus a custom derivation for the RSK node) pinned **two-sided** inside each 10-min
  BMRA CI (Table S5). 51 free params (shared core + g and K per constrained edge). Grounded in the
  paper's timepoint structure (Methods p.22: training = 10-min RPPA + WB time course, constraint =
  the Eq.14 CIs, validation = 45-min). Because the Trk ligand response is transient/adaptive (no
  sustained stimulated steady state), the exact `r_ij` are evaluated at the **basal** steady state.
  The RTK receptor-**dimerization** node (no closed-form Eq.14 `C_i`) has its incoming edges
  documented, not constrained — the honest scope used for SKMEL's IRS. `job_type = check` →
  **Satisfied 40/46** at the published parameters (the tight-CI edges bind); a bounded joint fit
  reaches **43/46**, pulling the tightest-CI coefficients (TrkB RSK←ERK, S6K←AKT in both) into band.

The BMRA→model map (the crux): the model's `g_<edge>` connection strengths *are* the paper's
hyperbolic-multiplier γ (Eq. 24: γ>1 activation, γ<1 inhibition), and the BMRA confidence
intervals pin each connection's **sign**; PyBNF expresses these as BPSL `.prop` constraints on
carrier observables (cf. the qualitative `cstar_skmel133_bpsl` slug). See
`cstar_skmel133_bmra/VALIDATION.md` for the full derivation and its two honest limitations.

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
