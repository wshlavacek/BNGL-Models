# Korwek-2023 — innate immune response to IFN-β and poly(I:C) (PyBNF fitting jobs)

PyBNF edition-2 parameter-fitting jobs derived from one paper:

> Korwek Z, Czerkies M, Jaruszewicz-Błońska J, Prus W, Kosiuk I, Kochańczyk M, Lipniacki T.
> **"Nonself RNA rewires IFN-β signaling: A mathematical model of the innate immune response."**
> *Sci Signal* 2023; **16**(815):eabq1173.
> DOI: [10.1126/scisignal.abq1173](https://doi.org/10.1126/scisignal.abq1173)

Built with the `curate-pybnf-job` skill. Each job below is a **self-contained folder** — its own
model, conf, data, reproduction figure, and README with the exact adaptations from the published
model, verification results, and a ready-to-paste `_manifest.py` snippet.

## The shared model

Korwek et al. provide their model as **data file S2** (`innate_immunity.bngl`). It couples five
regulatory modules — poly(I:C) sensing (RIG-I:MAVS, PKR–eIF2α, OAS3–RNase L), NF-κB (TAK1, IKK,
IκBα, A20), IRF3 (TBK1, IRF3), IFN-β (IFNB1 transcription, IFNAR) and STAT1/2 — into a
**53-species, 96-reaction** deterministic ODE network. The model is nondimensionalized in amount
(resting pools are 1) and its 38 independent parameters were fit with PyBioNetFit to 2915 data
points spanning several stimulation protocols.

The reconstruction each job uses is a fitting-ready copy of
[`models/innate_immune_response_korwek2023`](../../models/innate_immune_response_korwek2023),
whose primary file reproduces the authors' data file S2 **bit-for-bit** (maximum absolute
difference 0 over 43 observables × 4321 time points). That is what makes the published
parameterization a strong oracle here: "published" means the authors' actual fit, not a
transcription of it.

One parameter is added relative to data file S2: **`h_A20_gene`**, a gene-presence switch on
TNFAIP3 transcription in the idiom of the authors' own `h_Mavs` / `h_Pkr_gene` /
`h_Rnasel_gene`, so that the A20 KO cell line the paper simulates can be expressed as a
`condition:`. It is 1, and therefore inert, in the wild-type arm.

## The jobs

| slug | fits | flavor | data source | status |
|---|---|---|---|---|
| [`nfkb_tnfa`](nfkb_tnfa/) | the 15 NF-κB-module rate constants marked FITTED in table S1, against the TNF-α protocol; WT + A20 KO | quantitative, `chi_sq`, **PEtab-exportable** | **digitized** figs. S12A and S12B (111 points, 4 blots × 2 replicates) | ✅ tier-1 + PEtab round-trip + fit |

The paper states that TNF-α "is used to calibrate kinetic rate constants in the NF-κB pathway
module that follows Lipniacki et al. (2004) and Tay et al. (2010)", so `nfkb_tnfa` reproduces a
documented calibration step rather than inventing a fitting problem. TNF-α reaches the NF-κB
module alone — the poly(I:C), IRF3, IFN-β and STAT1/2 modules stay quiescent under it — so the
protocol isolates that module inside the complete five-module model.

## Data provenance

Korwek et al. tabulate none of their data: the quantified Western blots and the fitted model's
trajectories exist only as markers in the supplementary figures. Both were recovered by
[`models/innate_immune_response_korwek2023/digitize_korwek2023.py`](../../models/innate_immune_response_korwek2023/digitize_korwek2023.py),
which pulls each figure's bitmap out of the PDF at native resolution, calibrates every log panel
on its own axis tick marks (≈53 px per decade, so ≈4% per pixel), reads the model's open circles
and each blot's filled dots by colour, and flags every marker that sits on the axis floor or is
clipped at the top. The committed CSVs carry those flags; `nfkb_tnfa/make_exp.py` turns them into
`.exp` tables, writing a flagged marker as `NaN` because its value is not readable.

## Source files

- Paper and supplement: `dev/papers/Korwek-2023/` (primary sources, not redistributed)
- Authors' model: data file S2, `Korwek_et_al__Software/innate_immunity.bngl`
- Curated model collection: [`models/innate_immune_response_korwek2023`](../../models/innate_immune_response_korwek2023)
