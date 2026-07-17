# ppatf2_phospho — pp-ATF2(T69/T71) phosphorylation after anisomycin (PyBNF edition-2 job)

A PyBNF edition-2 parameter-fitting job derived from:

> Kirsch K, Zeke A, Tőke O, Sok P, Sethi A, Sebő A, Kumar GS, Egri P, Póti ÁL,
> Gooley P, Peti W, Bento I, Reményi A, Alexa A.
> **"Co-regulation of the transcription controlling ATF2 phosphoswitch by JNK and p38."**
> *Nat Commun* 2020; **11**(1):5769.
> PMCID: [PMC7666158](https://pmc.ncbi.nlm.nih.gov/articles/PMC7666158/) ·
> DOI: [10.1038/s41467-020-19582-3](https://doi.org/10.1038/s41467-020-19582-3)

> **Created by the `validate-pybnf-job` audit (2026-07-12)** as the sibling of
> [`p38atf2_binding`](../p38atf2_binding/). The former single `fig7b_timecourse` slug fit
> the p38-ATF2 **binding** data to the phosphorylation observable by mistake; this slug is
> the *actual* pp-ATF2(T69/T71) phosphorylation fit the paper reports (Fig. 7b pp-ATF2
> panel, absolute µM). See [`VALIDATION.md`](../p38atf2_binding/VALIDATION.md) for the
> full audit.

## The model

Identical to the sibling slug: the authors' JNK/p38/ATF2 TAD phosphoswitch model
(Supplementary Software 1), deterministic ODE, **48 species, 152 reactions**, generated
network byte-identical to the authors' file. JNK (D-site) and p38 (bipartite D+F-site)
distributively phosphorylate T69/T71; JNK phosphorylates S90, whose phosphorylation evicts
p38 from the FRS.

## What is fit

WT **pp-ATF2(pT69/pT71)** phosphorylation level (**absolute µM**) after **anisomycin**
stimulation, **0–40 min**, in **two conditions** (Fig. 7b "pp-ATF2 (pT69/pT71)" panel):

- **CTR** (control) — model observable `pT69pT71` at nominal params.
- **JNK-IN-8** (JNK inhibitor) — `pT69pT71` with **`k1 = k2 = 0`** (JNK catalytically dead,
  so only p38 phosphorylates T69/T71). A second experiment sharing the single free `dp3`.

`pT69pT71` rises **monotonically** to a plateau (constant `kstim` → no adaptation). This is
the phosphorylation counterpart to the peak-decay p38-ATF2 **binding** curve fit by
[`p38atf2_binding`](../p38atf2_binding/).

**Experimental design.** Pre-equilibrate at basal MAPK activity (`Stim = 0`) → apply
anisomycin (`Stim = 1`) → measure. The JNK-IN-8 experiment carries the inhibitor through
both phases.

## Data source & provenance — digitized (no Source Data table)

The pp-ATF2(T69/T71) **phosphorylation** time course is **not** in the paper's Source Data
workbook (which has sheets `Figure_1..6` + `S*`, but no `Figure_7`); in Fig. 4 it appears
only as blot *images*. It is published quantitatively only as the **Fig. 7b pp-ATF2
(pT69/pT71) panel** (absolute µM). `digitize.py` records the digitization:

- p.10 of the article PDF rendered at 400 dpi (`pdftoppm`); pp-ATF2 panel cropped.
- Linear axes; points at t = 0, 10, 20, 40 min; marker centers read to ~±0.03 µM.
- No error bars on this panel → objective **`sos`** (unweighted least squares).

| condition | t = 0 | 10 | 20 | 40 min |
|---|---|---|---|---|
| CTR (`ppatf2_phospho_ctr.exp`) | 0.21 | 0.55 | 0.72 | 0.80 µM |
| JNK-IN-8 (`ppatf2_phospho_jnki.exp`) | 0.085 | 0.24 | 0.42 | 0.42 µM |

## Free parameter (1) — `dp3`, fit to the pp-ATF2 WB

Per Supplementary Table 2, only **`dp3`** (the ATF2 pT69/pT71 phosphatase) was fit to the
pp-ATF2 signal ("data = Fig7b - pp-ATF2, WB"). Every other constant is held at its published
value (pp-JNK arm from pp-JNK WB, p38 arm from the NanoBit, `kon*/koff*/k1-k4` in vitro).
`dp3` sets the plateau height of the monotonic rise.

| id | published | PyBNF best-fit | ratio | role |
|---|---|---|---|---|
| `dp3` | 9.54e-3 | 1.008e-2 | 1.06× | phosphatase on ATF2 pT69 / pT71 |

## Verification (all pass)

- **Tier-1** (`check_conf.py`): edition 2, `job_type=de`, 2 experiments bound, 1 free param
  binds by id, no `__FREE`. **PASS.**
- **PEtab round-trip** (`petab_roundtrip.py --job-type de --method ode`): export → PEtab v2 lint
  clean → re-import (the two-phase protocol renders as multi-period experiments). **PASS.**
- **Gate 3a — reproduce at published params** (BNG2.pl): CTR to **median 6.1 %** rel. err,
  JNK-IN-8 to **median 12.5 %** (the lone 10-min JNK-IN-8 point sits below the authors' own
  calc curve too — see the paper's panel). Both match the authors' `calc-CTR` / `calc-JNK-IN-8`
  lines.
- **Gate 3b — recover by fitting** (`de`, ~8 s): `dp3` recovered **1.008e-2 vs 9.54e-3 =
  1.06×** of published; best-fit improves JNK-IN-8 to median 7.3 %.
- **Reproduction figure** `ppatf2_phospho_reproduction.png` (both conditions, published &
  best-fit vs digitized data).

## ✅ PEtab-v2-exportable · committed `petab/` bundle

The fit uses `sos` with no normalization (absolute µM), so it stays inside the PEtab v2
exportable subset. `scripts/petab_roundtrip.py --job-type de --method ode` exports → lints
clean → re-imports: the two-phase `preequilibrate → stimulate` protocol maps cleanly to PEtab
v2 **multi-period experiments** (each fit condition becomes a two-period experiment), so the
export faithfully renders CTR and the JNK-IN-8 `k1 = k2 = 0` inhibitor. (Contrast the two
native-only siblings: [`p38atf2_binding`](../p38atf2_binding/) uses `normalization`,
[`phosphoswitch_bpsl`](../phosphoswitch_bpsl/) uses BPSL `.prop` constraints.)

The runnable `ppatf2_phospho.conf` is the single source of truth; the exported problem is
committed under [`petab/`](petab/) as an illustrative byproduct (the same pattern as PyBNF's own
`examples/tutorial/12_petab_roundtrip/petab/`) for consumers who want PEtab v2 directly without
running PyBNF. It is **generated by `make_petab.py`, never hand-written**, so it can't drift from
the conf — re-run that script after any change. The bundle: 1 observable (`observableFormula` =
the bare `pT69pT71` function name), 1 estimated parameter (`dp3`) and 8 measurements, with
`conditions.tsv` + `experiments.tsv` capturing the two two-period experiments (4
experiment-period rows over 8 long-format condition-change rows: `Stim` 0→1 plus the `k1`/`k2`
knockouts).

| file | role |
|---|---|
| `make_petab.py` | regenerates the committed `petab/` bundle from `ppatf2_phospho.conf` (`export_job`); **run after any conf/model/exp change** |
| `petab/` | **generated** PEtab v2 bundle (do not hand-edit) — `problem.yaml` + `parameters.tsv` / `observables.tsv` / `measurements.tsv` + `conditions.tsv` / `experiments.tsv` + a copy of `ppatf2_phospho.bngl` |

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Kirsch-2020/ppatf2_phospho
pybnf -c ppatf2_phospho.conf
python make_petab.py                 # regenerate the committed petab/ bundle from the conf
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='ppatf2_phospho', conf='ppatf2_phospho.conf', simulator='ode',
    observables=('pT69pT71',),
    system='JNK/p38/ATF2 TAD phosphoswitch, WT pp-ATF2(T69/T71) phosphorylation anisomycin '
           'time course, CTR + JNK-IN-8 (Kirsch 2020, PMC7666158); ODE, pre-equilibrate -> '
           'stimulate, 40 min; digitized Fig. 7b, sos, absolute uM; PEtab-exportable'),
```
