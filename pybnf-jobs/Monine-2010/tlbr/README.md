# tlbr — trivalent-ligand / bivalent-receptor aggregation, network-free NFsim (PyBNF edition-2 job)

A PyBNF edition-2, parameter-fitting job setup. The **network-free model** is from:

> Monine MI, Posner RG, Savage PB, Faeder JR, Hlavacek WS.
> **"Modeling multivalent ligand-receptor interactions with steric constraints on
> configurations of cell-surface receptor aggregates."**
> *Biophys J* 2010; **98**(1):48–56.
> PMID: [20074516](https://pubmed.ncbi.nlm.nih.gov/20074516/) ·
> DOI: [10.1016/j.bpj.2009.09.043](https://doi.org/10.1016/j.bpj.2009.09.043)

fit to the trivalent-antigen binding data of:

> Posner RG, Geng D, Haymore S, Bogert J, Pecht I, Licht A, Savage PB.
> **"Trivalent antigens for degranulation of mast cells."**
> *Org Lett* 2007; **9**(18):3551–3554.
> PMID: [17691814](https://pubmed.ncbi.nlm.nih.gov/17691814/) ·
> DOI: [10.1021/ol071231z](https://doi.org/10.1021/ol071231z)

This is the same fitting problem as **BioNetFit 1 example 3** (Thomas et al. 2016) and PyBNF
2019 paper corpus problem **11-TLBR** (Mitra et al., *iScience* 2019), re-expressed on the
edition-2 surface.

## The model

A **trivalent ligand** `L(s,s,s)` and a **bivalent receptor** `R(s,s)`. Free ligand binds a
receptor site with association strength `K1`; a receptor-bound ligand arm crosslinks a second
receptor with strength `K2`; both bonds break at rate `koff`. Because each ligand carries three
arms and each receptor two sites, crosslinking builds **aggregates that grow without bound** —
there is no finite reaction network — so the model is simulated **network-free with NFsim**
(`method: nf`), with `complex` tracking and a raised global molecule limit (`gml`). Cell-fraction
units (a subvolume `f = 0.001` of a cell; ~300 receptors, up to ~1.3×10⁵ ligands at the top
dose). Stochastic.

## What is fit

One dose-response dataset — the **bound-ligand fraction FL** across a titration of total ligand:

| dataset | observable | design | source |
|---|---|---|---|
| FL vs. total ligand LTconc (12 doses, 5×10⁻⁴–213 nM) | `FL()` = `lambda`/`alpha` = (L_total−L_free)/(2·R_total·alpha) | NF dose-response scan over `LTconc`, each dose to t_end=5000 s | Posner 2007 (Fig. as digitized for BioNetFit 1 ex 3) |

`FL()` is a model **function** (parenthesized in the `.exp` header). The `.exp` carries **no
`_SD`** columns → the **`sos`** objective.

## Files

| file | role |
|---|---|
| `tlbr.bngl` | edition-2, fitting-ready network-free model (no actions block — NF needs no `generate_network`) |
| `tlbr.conf` | the edition-2 job setup (`method: nf` dose-response scan, `sos`) |
| `tlbr.exp` | fit target — FL() vs. LTconc, 12 doses (no SD) |
| `make_reproduction.py` | regenerates the reproduction PNG (NFsim `parameter_scan` at the PyBNF best-fit vs. data) |
| `tlbr_reproduction.png` | verification figure |

## ⚠️ Shipped nominals are a STARTING point, not the fit

The model's shipped nominals (`alpha=1`, `K1=0.467`, `K2=87.03` — "from Brandon") are the GenFit
**starting values**, not the fit. The **PyBNF best-fit** (Mitra 2019, RuleHub
[`11-TLBR/fit_ss`](https://github.com/RuleWorld/RuleHub/tree/2019Aug21/Published/Mitra2019/11-TLBR/fit_ss)
init7, `sos = 0.00214` — the lowest objective across the paper's four algorithms) is what
reproduces the data and is used for the figure:

| id | starting nominal | **PyBNF best-fit** | role |
|---|---|---|---|
| `alpha` | 1 | **0.74568** | overall scale of the FL readout (linear) |
| `K1` | 0.467 | **0.10915** | free-ligand / receptor-site equilibrium constant (/nM) |
| `K2` | 87.03 | **33.576** | crosslinking equilibrium constant (/nM) |

The four algorithms (`de`/`ade`/`pso`/`ss`) land at similar objectives with `alpha ≈ 0.74–0.79`,
`K1 ≈ 0.1–0.3`, `K2 ≈ 30–55` — the fit is reasonably identifiable here (contrast the sloppy
Kozer kinetics). See the paper-level [`../README.md`](../README.md).

## Changes from the source model (documented in `tlbr.bngl`)

- **No `begin actions` block** — a network-free model needs no `generate_network`; the
  dose-response scan is synthesized from `tlbr.conf`. (This is the correct NF shape — contrast
  the network-generating `Kozer-2013-2014/egfr_ode`, which *keeps* its `generate_network`
  directive.)
- Fitted constants (`alpha`, `K1`, `K2`) bind to model ids **by name** (ADR-0034, no `__FREE`).
  Otherwise the model network is byte-for-byte the Monine 2010 / BioNetFit-1 model.

## ✅ PEtab-v2-exportable · ✅ workstation-runnable

Verified with `scripts/petab_roundtrip.py --job-type de --method nf` (export → lint clean →
import). Unlike the heavy `egfr_nf`, a single dose-response scan here is **fast** (~11 s for all
12 doses via one `parameter_scan`), so the reproduction runs on a workstation in ~1 minute; a
full metaheuristic fit (population × iterations × doses × replicates) is still best on a cluster.

## Verification

- **Tier-1** (`scripts/check_conf.py`): edition 2, `job_type=de` resolves, the experiment binds
  data, `model.stochastic = True` (→ `simulator='nf'`), 3 free params bind by id, no `__FREE`.
  **PASS.**
- **PEtab round-trip** (`--job-type de --method nf`): export → PEtab v2 lint clean → re-import.
  **PASS.**
- **Reproduction** (`tlbr_reproduction.png`, NFsim `parameter_scan` at the PyBNF best-fit, 7
  replicates averaged): reproduces the Posner titration to **median 3.4 % / max 7.1 % relative
  error** (over the 9 doses with FL > 0.05; the 3 lowest doses sit at the near-zero noise floor,
  one datum even negative); **sos ≈ 0.0033** vs. the published **0.00214** (the gap is NFsim
  stochasticity + finite replicates).

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Monine-2010/tlbr
pybnf -c tlbr.conf                    # raise population_size to reproduce the published fit
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='tlbr', conf='tlbr.conf', simulator='nf', stochastic=True,
    observables=('FL',),
    system='Trivalent-ligand/bivalent-receptor aggregation (Monine 2010, PMID 20074516; '
           'Posner 2007 data, PMID 17691814; BioNetFit 1 ex 3); network-free NFsim, '
           'dose-response, sos; PEtab-exportable',
    recover={'alpha': 0.74568, 'K1': 0.10915, 'K2': 33.576}, tol=0.5),
    # PyBNF best-fit: RuleHub Published/Mitra2019/11-TLBR/fit_ss init7 (sos 0.00214). All three
    # params are reasonably identifiable across the paper's four algorithms.
```
