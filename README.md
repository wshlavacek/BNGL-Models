# BNGL-Models

A collection curated models (in `models/`) drawn from the literature. Models are formulated using the BioNetGen language (BNGL). Provided BNGL files include annotation and simulation protocols. Model contributions include provenance and reference simulation results. Simulation results generated from BNGL files are independently verified.

## Models

| File(s) | Description of Primary Model | Reference(s) |
|-------|-------------------|---------------|
| `blbr_rings_posner1995/`<br>`blbr_rings_posner1995.bngl`<br>`blbr_rings_posner1995_no_rings.bngl`<br>`blbr_rings_posner1995_alt_formulation.bngl` | A model for interaction of a population of soluble bivalent ligands with a population of bivalent cell-surface receptors when cyclic receptor aggregates (e.g., a cyclic receptor dimer with 2:2 ligand:receptor stoichiometry) can form | Dembo & Goldstein, 1978; Posner et al., 1995 |
| `blbr_cooperativity_posner2004/`<br>`blbr_cooperativity_posner2004.bngl`<br>`blbr_cooperativity_posner2004_rings.bngl` | A model for interaction of a population of soluble bivalent ligands (e.g., anti-IgE IgG antibodies) with a population of bivalent cell-surface receptor (e.g., cell-surface IgE antibodies) when ligand-receptor binding is (negatively) cooperative | Wofsy & Goldstein, 1987; Posner et al., 2004 |
| `genetic_switch_gardner2000/`<br>`genetic_switch_gardner2000.bngl`<br>`genetic_switch_gardner2000_iptg.bngl` | A model for a genetic toggle switch | Gardner et al., 2000 |
| `blbr_heterogeneity_goldstein1980/`<br>`blbr_heterogeneity_goldstein1980.bngl` | A model for interaction of a population of soluble bivalent ligands with multiple populations of bivalent cell-surface receptors (i.e., receptors with distinct ligand-binding properties) | Goldstein & Wofsy, 1980 |
| `antigen_pulses_harmon2017/`<br>`antigen_pulses_harmon2017.bngl`<br>`antigen_pulses_harmon2017_simplified.bngl` | A model for interaction of a population of soluble antigens with a population of cell-surface IgE antibodies when microfluidics is used to control the antigen concentration in solution (i.e., to expose cells to square-wave antigen pulses) | Harmon et al., 2017 |
| `kinetic_proofreading_hlavacek2001/`<br>`kinetic_proofreading_hlavacek2001.bngl`<br>`kinetic_proofreading_hlavacek2001_messenger.bngl` | A model for kinetic proofreading in IgE receptor signaling | Hlavacek et al., 2001; 2002 |
| `mapk_erk_system_regulation_kocieniewski2026/`<br>`mapk_erk_system_regulation_kocieniewski2026.bngl` | A rule-based model of MAPK/ERK signaling regulation with EGFR, SOS, RAS, RasGAP, ARAF, BRAF, CRAF, MEK1, MEK2, ERK, and 14-3-3-mediated RAF regulation | Kocieniewski and Lipniacki, 2026 |
| `mek_isoform_erk_feedback_miller2026/`<br>`mek_isoform_erk_feedback_miller2026.bngl`<br>`mek_isoform_erk_feedback_miller2026_kocieniewski2013.bngl` | A wild-type MEK1/MEK2 ERK cascade model with ERK-mediated feedback to SOS1 and MEK1 Thr292; the primary file uses the Miller et al. PyBioNetFit maximum-likelihood parameterization, and the variant uses the original 2013 parameterization | Kocieniewski and Lipniacki, 2013; Miller et al., 2026 |
| `early_fceri_signaling_faeder2003/`<br>`early_fceri_signaling_faeder2003.bngl` | A detailed rule-based model for early FceRI signaling after IgE-dimer-induced receptor aggregation, including Lyn-mediated receptor phosphorylation, Syk recruitment, and Syk activation-loop phosphorylation | Faeder et al., 2003 |
| `lac_operon_dreisigmeyer2008/`<br>`lac_operon_dreisigmeyer2008.bngl` | A model for lactose-triggered induction of the lac operon | Dreisigmeyer et al., 2008 |
| `steric_effects_hlavacek1999/`<br>`steric_effects_hlavacek1999.bngl`<br>`steric_effects_hlavacek1999_nosteric.bngl`<br>`steric_effects_hlavacek1999_combinatoric.bngl` | A model for multivalent ligand-monovalent receptor interaction with steric effects | Hlavacek et al., 1999 |
| `blbr_dembo1978/`<br>`blbr_dembo1978.bngl`<br>`blbr_dembo1978_monovalent_inhibitor.bngl`<br>`blbr_dembo1978_with_rings.bngl` | A model for bivalent ligand-bivalent cell-surface receptor interaction | Dembo & Goldstein, 1978; Perelson & DeLisi, 1980 |
| `tlbr_yang2008/`<br>`tlbr_yang2008.bngl` | A model for trivalent ligand-bivalent cell-surface receptor interaction | Goldstein & Perelson, 1984; Yang et al., 2008 |
| `tlbr_solution_macken1982/`<br>`tlbr_solution_macken1982.bngl` | A model for trivalent ligand-bivalent receptor interaction in solution | Macken & Perelson, 1982 |

## Setup

Requires Python 3.12+.

```sh
uv sync
```
