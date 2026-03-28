# BNGL-Models

A collection models drawn from the literature. Models are formulated using the BioNetGen language (BNGL). Provided BNGL files include annotation and simulation protocols. Model contributions include provenance and reference simulation results. Simulation results generated from BNGL files are independently verified.

## Models

| File | Model Description | Reference(s) |
|------|-------------------|---------------|
| `blbr_rings_posner1995.bngl` | A model for interaction of a population of soluble bivalent ligands with a population of bivalent cell-surface receptors when cyclic receptor aggregates (e.g., a cyclic receptor dimer with 2:2 ligand:receptor stoichiometry) can form | Dembo & Goldstein, 1978; Posner et al., 1995 |
| `blbr_cooperativity_posner2004.bngl` | A model for interaction of a population of soluble bivalent ligands with a population of bivalent cell-surface receptor when ligand-receptor binding is (negatively) cooperative | Wofsy & Goldstein, 1987; Posner et al., 2004 |
| `genetic_switch_gardner2000.bngl` | A model for a genetic toggle switch | Gardner et al., 2000 |
| `blbr_heterogeneity_goldstein1980.bngl` | A model for interaction of a population of soluble bivalent ligands with multiple populations of bivalent cell-surface receptors (i.e., receptors with distinct ligand-binding properties) | Goldstein & Wofsy, 1980 |
| `antigen_pulses_harmon2017.bngl` | A model for interaction of a population of soluble (monovalent) antigens with a population of cell-surface IgE antibodies when microfluidics is used to control the antigen concentration in solution (i.e., to expose cells to square-wave antigen pulses) | Harmon et al., 2017 |
| `kinetic_proofreading_hlavacek2001.bngl` | A model for kinetic proofreading in IgE receptor signaling | Hlavacek et al., 2001; 2002 |
| `lac_operon_dreisigmeyer2008.bngl` | A model for lactose-triggered induction of the lac operon | Dreisigmeyer et al., 2008 |
| `steric_effects_hlavacek1999.bngl` | A model for multivalent ligand-monovalent receptor interaction with steric effects | Hlavacek et al., 1999 |
| `blbr_dembo1978.bngl` | A model for bivalent ligand-bivalent cell-surface receptor interaction | Dembo & Goldstein, 1978; Perelson & DeLisi, 1980 |
| `tlbr_yang2008.bngl` | A model for trivalent ligand-bivalent cell-surface receptor interaction | Goldstein & Perelson, 1984; Yang et al., 2008 |
| `tlbr_solution_macken1982.bngl` | A model for trivalent ligand-bivalent receptor interaction in solution | Macken & Perelson, 1982 |

## Setup

Requires Python 3.12+.

```sh
uv sync
```
