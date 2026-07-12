# fceri_gamma — FcεRI γ-chain signaling, Gillespie SSA parameter recovery (PyBNF edition-2 job)

A PyBNF edition-2, parameter-fitting job setup. The benchmark **model** is from:

> Gupta A, Mendes P.
> **"An Overview of Network-Based and -Free Approaches for Stochastic Simulation of Biochemical
> Systems."** *Computation* 2018; **6**(1):9.
> DOI: [10.3390/computation6010009](https://doi.org/10.3390/computation6010009)

This is a **synthetic-data parameter-recovery** problem: `fceri_gamma2.exp` was generated at
known **ground-truth parameters** (`fceri_gamma2_ground_truth.bngl`), and the fit tries to
**recover** them — the reproduction below confirms the ground-truth model regenerates the data.
Part of the PyBNF 2019 paper corpus (Mitra et al., *iScience* 2019), on the edition-2 surface.

## The model

The high-affinity IgE receptor **FcεRI** γ-chain signaling module: a bivalent ligand `Lig(l,l)`
crosslinks receptors `Rec` (via the single `a` site); the Src-family kinase **Lyn** binds the β
chain and *trans*-phosphorylates β and γ; **Syk** binds phospho-γ and is *trans*-phosphorylated;
phosphatases (`dm`, `dc`) reverse it. Because each receptor has one aggregation site, aggregates
cap at **receptor dimers** — the network is **finite (~58,000 reactions)**, dominated by the
combinatorial phospho/binding states of the dimer. Small molecule counts (Lig 6000, Rec 400,
Lyn 28, Syk 400) → simulated by the **Gillespie SSA**. Stochastic.

## What is fit

One stochastic time-course experiment; the six observables ARE the `.exp` columns:

| observable | meaning |
|---|---|
| `LynFree` | free (unbound) Lyn |
| `RecMon` | monomeric (unaggregated) receptor |
| `RecPbeta` | β-phosphorylated receptor |
| `RecPgamma` | γ-phosphorylated receptor |
| `RecSyk` | Syk bound to phospho-γ |
| `RecSykPS` | activation-loop-phosphorylated receptor-bound Syk |

Synthetic data carry **no `_SD`** columns → **`sos`**; `smoothing = 3` averages replicate SSA
trajectories per parameter set. **20 free** rate/probability constants (`km1`/`km2` fixed at 0).

## Files

| file | role |
|---|---|
| `fceri_gamma2.bngl` | edition-2, fitting-ready SSA model — **KEEPS** `generate_network({overwrite=>1,max_iter=>100})` |
| `fceri_gamma.conf` | the edition-2 job setup (`method: ssa`, `sos`, `smoothing`, `ss` + refine) |
| `fceri_gamma2.exp` | synthetic fit target (6 observables, t=0–100 s) |
| `fceri_gamma2_ground_truth.bngl` | the model at the **known-true** parameters that generated the data (documentation + recovery target; not referenced by the conf) |
| `make_reproduction.py` | regenerates the reproduction PNG (SSA at ground truth vs. data; builds the 58k net once, averages replicates) |
| `fceri_gamma_reproduction.png` | verification figure |

## ⚠️ RETAINED network directive (the key adaptation)

The reaction network is finite but large (~58k reactions), and the classic model expands it with
`generate_network({max_iter=>100})`. Per the corpus convention (lanl/PyBNF#473,
`references/real-world-anatomy.md`), that **network-definition directive is part of the model and
is RETAINED** here as `generate_network({overwrite=>1,max_iter=>100})` — pybnf captures it
(`pset.py:617`) in place of its bare-default `generate_network({overwrite=>1})`. (The PyBNF
`examples/real-world/fceri_gamma/` mirror had *stripped* it; this curated copy restores it.) Only
the *simulation* action (`simulate`) is dropped — that is synthesized from the conf.

## Recovery targets (ground truth)

The fit should recover `fceri_gamma2_ground_truth.bngl`'s values (nominals in `fceri_gamma2.bngl`
are generic `0.01`/`1` starting points):

| param | truth | param | truth | param | truth | param | truth |
|---|---|---|---|---|---|---|---|
| `kp1` | 1.328×10⁻⁶ | `kmLs` | 0.12 | `pLb` | 30 | `pLS` | 30 |
| `kp2` | 0.25 | `kpS` | 0.06 | `pLbs` | 100 | `pLSs` | 100 |
| `kpL` | 0.05 | `kmS` | 0.13 | `pLg` | 1 | `pSS` | 100 |
| `kmL` | 20 | `kpSs` | 0.06 | `pLgs` | 3 | `pSSs` | 200 |
| `kpLs` | 0.05 | `kmSs` | 0.13 | `dm` | 20 | `dc` | 20 |

> Some constants are only weakly identifiable from a single noisy SSA trajectory of 6 observables;
> the recovery is expected to land the well-constrained constants (aggregation `kp2`,
> phosphorylation/dephosphorylation balance) near truth and leave others loose — use a generous
> `tol` on any `recover` assertion.

## ✅ PEtab-v2-exportable · 🔶 Heavy (cluster-scale)

Verified with `scripts/petab_roundtrip.py --job-type ss --method ssa` (export → lint clean →
import) — `sos` + plain-name observables, no `normalization`/BPSL, so it exports cleanly.
**Heavy:** expanding the ~58k-reaction network takes ~3 min, and a full fit (12 × 50 × 3 SSA
runs on that network) is a cluster job — so the manifest marks `heavy=True`. The tier-1 parse and
PEtab round-trip are backend-free; the reproduction builds the network once and runs a few SSA
replicates (workstation-feasible).

## Verification

- **Tier-1** (`scripts/check_conf.py`): edition 2, `job_type=ss` resolves, the experiment binds
  data, `model.stochastic = True` (→ `simulator='ssa'`), 20 free params bind by id, no `__FREE`.
  **PASS.**
- **PEtab round-trip** (`--job-type ss --method ssa`): export → PEtab v2 lint clean → re-import.
  **PASS.**
- **Reproduction** (`fceri_gamma_reproduction.png`, SSA at the ground-truth params, 8 replicates
  averaged; the 58k network built once, each replicate reset to the seed state): the ground-truth
  model reproduces the synthetic data for **all six observables** — median relative error
  **RecMon 2.6 %, RecPbeta 9 %, LynFree 10 %, RecPgamma/RecSyk/RecSykPS ~11 %** (over points above
  the SSA noise floor). Every curve starts at its correct t=0 value (RecMon 400, LynFree 28,
  phospho 0) and tracks the data through t=100 s, validating the recovery target.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Gupta-2018/fceri_gamma
pybnf -c fceri_gamma.conf              # cluster-scale (58k-reaction SSA); run large fits on a cluster
python make_reproduction.py            # SSA at ground truth vs. data (builds the net once)
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='fceri_gamma', conf='fceri_gamma.conf', simulator='ssa', stochastic=True, heavy=True,
    observables=('LynFree', 'RecMon', 'RecPbeta', 'RecPgamma', 'RecSyk', 'RecSykPS'),
    system='FceRI gamma-chain signaling (Gupta & Mendes 2018, Computation 6:9); Gillespie SSA, '
           '~58k-reaction network, synthetic-data parameter recovery; PEtab-exportable',
    recover={'kp2': 0.25, 'kpL': 0.05, 'dm': 20.0, 'dc': 20.0}, tol=0.5),
    # SYNTHETIC data from fceri_gamma2_ground_truth.bngl; recover the well-constrained constants.
    # heavy=True: ~58k-reaction network, cluster-scale full fit.
```
