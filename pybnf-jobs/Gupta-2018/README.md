# Gupta-2018 — FcεRI γ-chain signaling (PyBNF fitting job)

A PyBNF edition-2 parameter-fitting job derived from a stochastic-simulation benchmark:

> Gupta A, Mendes P.
> **"An Overview of Network-Based and -Free Approaches for Stochastic Simulation of Biochemical
> Systems."** *Computation* 2018; **6**(1):9.
> DOI: [10.3390/computation6010009](https://doi.org/10.3390/computation6010009)

Built with the `curate-pybnf-job` skill. The job below is a **self-contained folder** — its own
model, conf, synthetic data, ground-truth model, reproduction figure, and README with the exact
adaptations, verification results, and a ready-to-paste `_manifest.py` snippet. Part of the PyBNF
2019 paper corpus (Mitra et al., *iScience* 2019, 19:1012–1036), re-expressed on the edition-2
surface.

## The biochemistry

The high-affinity IgE receptor **FcεRI** γ-chain signaling module. A bivalent ligand crosslinks
receptors into **dimers**; the Src-family kinase **Lyn** binds the receptor β chain and
*trans*-phosphorylates β and γ; **Syk** docks on phospho-γ and is *trans*-phosphorylated on its
activation loop; phosphatases reverse the modifications. Even capped at receptor dimers, the
combinatorial phospho/binding states give a **~58,000-reaction** network, simulated stochastically
(Gillespie SSA) because the molecule counts are small.

## The job

| slug | fits | simulator | flavor | status |
|---|---|---|---|---|
| [`fceri_gamma`](fceri_gamma/) | 20 rate/probability constants recovered from synthetic data (6 observables, t=0–100 s) | **SSA** (Gillespie, ~58k-rxn network) | quantitative, **PEtab-exportable**, `sos`, synthetic-data recovery | ✅ tier-1 + PEtab round-trip + ground-truth reproduction · 🔶 heavy (cluster-scale) |

**Synthetic-data recovery.** The data were generated at known ground-truth parameters
(`fceri_gamma/fceri_gamma2_ground_truth.bngl`); the fit's job is to recover them, and the
reproduction confirms the ground-truth model regenerates the data. **PEtab-v2-exportable**
(`sos`, plain observables), and **heavy** — the 58k-reaction network takes ~3 min to expand and a
full fit is a cluster job.

### The network-directive lesson this job carries

`fceri_gamma`'s network is finite but large, and the classic model expands it with
`generate_network({max_iter=>100})`. That **network-definition directive stays in the model** —
the curated `fceri_gamma2.bngl` retains `generate_network({overwrite=>1,max_iter=>100})` (the
PyBNF `examples/real-world/` mirror had stripped it). This is the same principle the
`curate-pybnf-job` skill bakes in from `Kozer-2013-2014/egfr_ode`: **network-definition directives
stay in the model; only simulation actions move to the conf** (lanl/PyBNF#473).

## Source materials

- **Model / data / classic conf:** PyBNF `examples/fceri_gamma/` (with the pre-generated
  `fceri_gamma2.net`, ~58k reactions) and its edition-2 twin `examples/real-world/fceri_gamma/`.
- **Ground truth:** `fceri_gamma/fceri_gamma2_ground_truth.bngl` — the parameter set that
  generated the synthetic `fceri_gamma2.exp`; the recovery target.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Gupta-2018/fceri_gamma
pybnf -c fceri_gamma.conf              # cluster-scale (58k-reaction SSA)
python make_reproduction.py            # SSA at ground truth vs. synthetic data
```
