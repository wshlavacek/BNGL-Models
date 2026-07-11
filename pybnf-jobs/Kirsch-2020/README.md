# Kirsch-2020 — ATF2 phosphoswitch, JNK & p38 co-regulation (PyBNF fitting jobs)

PyBNF edition-2 parameter-fitting jobs derived from one paper:

> Kirsch K, Zeke A, Tőke O, Sok P, Sethi A, Sebő A, Kumar GS, Egri P, Póti ÁL,
> Gooley P, Peti W, Bento I, Reményi A, Alexa A.
> **"Co-regulation of the transcription controlling ATF2 phosphoswitch by JNK and p38."**
> *Nat Commun* 2020; **11**(1):5769.
> PMCID: [PMC7666158](https://pmc.ncbi.nlm.nih.gov/articles/PMC7666158/) ·
> DOI: [10.1038/s41467-020-19582-3](https://doi.org/10.1038/s41467-020-19582-3)

Built with the `curate-pybnf-job` skill. Each job below is a **self-contained folder** — its
own model, conf, data/constraints, reproduction figure, and README with the exact
adaptations from the published model, verification results, and a ready-to-paste
`_manifest.py` snippet.

## The shared model

Both jobs reconstruct the **same rule-based model**, which the authors provide as
**Supplementary Software 1** (`Bionetgen_JNK-p38-ATF2_model.bngl`; its actions block
generates the paper's Fig. 7b). JNK and p38 dock the ATF2 transactivation domain at
distinct sites — JNK at the D-site, p38 bipartitely at the D-site **and** the F-site
(SPFENEF / S90 region) — and distributively phosphorylate the T69/T71 "phosphoswitch". JNK
additionally phosphorylates the F-site **S90**, and **S90-P sterically blocks p38's F-site
recruitment** — the paper's central mechanism. Deterministic ODE, **48 species, 152
reactions**. The edition-2 reconstruction reproduces the authors' trajectories **exactly**
(max relative difference 0.0).

## The jobs

| slug | fits | flavor | data source | status |
|---|---|---|---|---|
| [`fig7b_timecourse`](fig7b_timecourse/) | WT pp-ATF2(T69/T71) anisomycin time course (8 cell params) | quantitative, **native-only** (`normalization=init`) | Source Data, Fig. 4a (with SD → `chi_sq`) | ✅ tier-1 + fit **chi_sq ≈ 64.9** (8.4 % median rel-err) |
| [`phosphoswitch_bpsl`](phosphoswitch_bpsl/) | S90 phosphoswitch → p38 recruitment orderings (4 cell params) | **BPSL** constraints, **native-only** | Suppl. Table 2 mutants; Figs. 3c/4b binding | ✅ tier-1 + **`check` 6/6 satisfied** |

`fig7b_timecourse` is the quantitative reconstruction of the authors' headline fit;
`phosphoswitch_bpsl` is its constraint-bearing sibling, expressing the paper's central
*qualitative* claim (JNK's S90 phosphorylation diminishes p38 binding: **S90N > JNK
inhibitor > WT > MUT4** for p38:ATF2 recruitment) as PyBNF BPSL `.prop` constraints.

Both jobs are **native-only** (not PEtab-v2-exportable): `fig7b_timecourse` because the
data are treated/untreated fold changes (`normalization`), `phosphoswitch_bpsl` because
BPSL constraints have no PEtab representation. Each is verified accordingly (bounded fit +
figure reproduction; and `job_type = check` satisfaction for the BPSL job) rather than with
a PEtab round-trip.

## Source materials

All from the article's supplementary files (Springer CDN); no public GitHub / BioModels
deposition exists.

| file | content |
|---|---|
| [Supplementary Software 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19582-3/MediaObjects/41467_2020_19582_MOESM4_ESM.zip) | the authors' BioNetGen `.bngl` (WT parameters) |
| [Supplementary Information PDF](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19582-3/MediaObjects/41467_2020_19582_MOESM1_ESM.pdf) | **Supplementary Table 2** (all fitted parameters + data provenance); Fig. S7 in-vitro kinetics |
| [Source Data workbook](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19582-3/MediaObjects/41467_2020_19582_MOESM6_ESM.xlsx) | sheet `Figure_4` = the pp-ATF2(T69/T71) time courses (5 conditions, with SD) |

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Kirsch-2020/<slug>
pybnf -c <slug>.conf
```

See each slug's `README.md` for its full write-up.
