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

All three jobs reconstruct the **same rule-based model**, which the authors provide as
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
| [`p38atf2_binding`](p38atf2_binding/) | WT p38-ATF2 NanoBit **binding** (`p38ATF2all`; free `keq6,kstim6,dp2,dp4`) | quantitative, **native-only** (`normalization=init`, `chi_sq`) | Source Data `Figure_4` (with SD) | ✅ **validated 91/100** — Gate 3a 5.8 %, Gate 3b ≤1.3× |
| [`ppatf2_phospho`](ppatf2_phospho/) | WT pp-ATF2(T69/T71) **phosphorylation** (`pT69pT71`; free `dp3`; CTR + JNK-IN-8) | quantitative, absolute µM (`sos`), **PEtab-exportable** | **digitized** Fig. 7b pp-ATF2 panel | ✅ **validated** — Gate 3a 6.1 %, Gate 3b `dp3` 1.06×, PEtab round-trip |
| [`phosphoswitch_bpsl`](phosphoswitch_bpsl/) | S90 phosphoswitch → p38 recruitment orderings (4 cell params) | **BPSL** constraints, **native-only** | Suppl. Table 2 mutants; Figs. 3c/4b binding | ✅ tier-1 + **`check` 6/6 satisfied** (build-verified; not yet primary-source-audited) |

`p38atf2_binding` and `ppatf2_phospho` are the two quantitative pieces of the authors'
Fig. 7b parameter determination (the paper fit its 8 cell params by **decomposing** across
datasets — see below). `phosphoswitch_bpsl` is their constraint-bearing sibling, expressing
the paper's central *qualitative* claim (JNK's S90 phosphorylation diminishes p38 binding:
**S90N > JNK inhibitor > WT > MUT4** for p38:ATF2 recruitment) as PyBNF BPSL `.prop`
constraints.

> **Validation note (`validate-pybnf-job`, 2026-07-12).** The two quantitative slugs
> replace an earlier single `fig7b_timecourse` slug, which fit the p38-ATF2 **binding**
> data to the ATF2-**phosphorylation** observable `pT69pT71` by mistake. The audit
> (`p38atf2_binding/VALIDATION.md`) re-scoped it: `p38atf2_binding` fits the binding curve
> to `p38ATF2all` (the correct observable), and `ppatf2_phospho` fits the real
> pp-ATF2(T69/T71) curve. The shared model is **byte-identical** to the authors' BNGL
> (48/48 species, 152/152 reactions).

Per Supplementary Table 2, the 8 cell params come from different data: pp-JNK WB →
`keq7,kstim7,dp1`; p38-ATF2 NanoBit → `keq6,kstim6,dp2,dp4` (`p38atf2_binding`); pp-ATF2 WB
→ `dp3` (`ppatf2_phospho`). `p38atf2_binding` is **native-only** (treated/untreated fold
changes → `normalization`) and `phosphoswitch_bpsl` is native-only (BPSL), so both are verified
by a bounded fit + figure reproduction (`job_type = check` for BPSL) rather than a PEtab
round-trip. `ppatf2_phospho` is absolute-µM `sos` and **PEtab-v2-exportable**: it round-trips
through PEtab v2 (the two-phase `preequilibrate → stimulate` protocol renders as multi-period
experiments) and ships a committed [`petab/`](ppatf2_phospho/petab/) bundle, regenerated from the
conf by `make_petab.py`.

## Source materials

All from the article's supplementary files (Springer CDN); no public GitHub / BioModels
deposition exists.

| file | content |
|---|---|
| [Supplementary Software 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19582-3/MediaObjects/41467_2020_19582_MOESM4_ESM.zip) | the authors' BioNetGen `.bngl` (WT parameters) |
| [Supplementary Information PDF](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19582-3/MediaObjects/41467_2020_19582_MOESM1_ESM.pdf) | **Supplementary Table 2** (all fitted parameters + data provenance); Fig. S7 in-vitro kinetics |
| [Source Data workbook](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19582-3/MediaObjects/41467_2020_19582_MOESM6_ESM.xlsx) | sheet `Figure_4` = the Fig. 4a **p38-ATF2 NanoBit binding** time courses (treated/untreated, 5 conditions, with SD). The pp-ATF2(T69/T71) *phosphorylation* curve is **not** here — it is published only in the Fig. 7b panel (digitized for `ppatf2_phospho`). |

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Kirsch-2020/<slug>
pybnf -c <slug>.conf
```

See each slug's `README.md` for its full write-up.
