# p38atf2_binding — p38-ATF2 TAD binding after anisomycin (PyBNF edition-2 job)

A PyBNF edition-2 parameter-fitting job derived from:

> Kirsch K, Zeke A, Tőke O, Sok P, Sethi A, Sebő A, Kumar GS, Egri P, Póti ÁL,
> Gooley P, Peti W, Bento I, Reményi A, Alexa A.
> **"Co-regulation of the transcription controlling ATF2 phosphoswitch by JNK and p38."**
> *Nat Commun* 2020; **11**(1):5769.
> PMCID: [PMC7666158](https://pmc.ncbi.nlm.nih.gov/articles/PMC7666158/) ·
> DOI: [10.1038/s41467-020-19582-3](https://doi.org/10.1038/s41467-020-19582-3)

> **Re-scoped by the `validate-pybnf-job` audit (2026-07-12).** This slug was formerly
> `fig7b_timecourse`, which fit its data to the **wrong observable**: it mapped the
> **p38-ATF2 NanoBit binding** curve to the ATF2-phosphorylation observable `pT69pT71`.
> The audit re-scoped it to the correct observable **`p38ATF2all`** (the p38:ATF2 complex),
> and split off the true pp-ATF2(T69/T71) phosphorylation fit into the sibling slug
> [`ppatf2_phospho`](../ppatf2_phospho/). See [`VALIDATION.md`](VALIDATION.md) for the full
> gate-by-gate evidence.

## The model

JNK and p38 dock the ATF2 transactivation domain at distinct sites — JNK at the D-site;
p38 bipartitely at the D-site **and** the F-site (SPFENEF / S90 region) — and
distributively phosphorylate the T69/T71 "phosphoswitch". JNK additionally phosphorylates
the F-site **S90**, and **S90-P sterically blocks p38's F-site (FRS) binding** — the
paper's central mechanism. MKK6/MKK7 activate p38/JNK and phosphatases reverse everything.
Deterministic ODE; **48 species, 152 reactions**. Faithful reconstruction of the authors'
Supplementary Software 1 (`Bionetgen_JNK-p38-ATF2_model.bngl`) — the generated network is
**byte-identical** to theirs (`model_diff.py`: 48/48 species, 152/152 reactions).

## What is fit

WT **p38-ATF2 TAD binding** time course after **anisomycin** stimulation of HEK293T cells,
**0–41 min** (19 points, mean ± SD, n=3), fit by the model observable **`p38ATF2all`**
(the p38:ATF2 complex — DRS+FRS bipartite plus DRS-only forms).

The signal **rises then falls** (peaks ~9 min): as pp-JNK accumulates it phosphorylates
S90 (`F~p`), and S90-P sterically evicts p38 from the FRS, so the bipartite complex is
transient. This peak-decay is the phosphoswitch, and it is exactly what the p38-ATF2
NanoBit assay measures. **The ATF2-phosphorylation observable `pT69pT71` is monotonic and
cannot produce this shape** — that is the sibling slug's job.

**Experimental design.** Pre-equilibrate at basal MAPK activity (`Stim = 0` → `keq6/keq7`),
then apply anisomycin (`Stim = 1` → `kstim6/kstim7`) and measure over 0–2466 s. PyBNF
synthesizes the two-phase protocol from `preequilibrate: basal` → `condition: stim`.

## Data source & provenance

The **Source Data** workbook
([MOESM6](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19582-3/MediaObjects/41467_2020_19582_MOESM6_ESM.xlsx)),
sheet `Figure_4` (Fig. 4a), WT column (mean rows 6–24, SD rows 29–47). That numeric block
is the Fig. 4a **left-panel NanoBit p38-ATF2 binding** curve ("treated/untreated"): WT
peaks ~1.6× at 9 min and decays; S90N rises to ~4×; MUT4 stays ~1 — matching that plot
exactly. The same data are overlaid with the model calc in the Fig. 7b "p38 + ATF2
interaction" panel. `extract_exp.py` transcribes it to `p38atf2_binding.exp` (t in seconds)
and reproduces the shipped file byte-for-byte. Model provided as
[Supplementary Software 1](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-19582-3/MediaObjects/41467_2020_19582_MOESM4_ESM.zip);
published parameters in **Supplementary Table 2**.

## Model observable → data mapping

`p38ATF2all` is a `Molecules` observable (plain name in the `.exp` header). The data are
**treated/untreated** ratios (anisomycin-stimulated ÷ unstimulated control), so the conf
compares them to `p38ATF2all` under `normalization = init` (each simulated series ÷ its own
basal/t=0 value = treated/untreated).

## Changes from the published model

- **No `begin actions` block** — the protocol is synthesized from the conf.
- **`Stim` gate** (0/1) selects basal (`keq*`) vs stimulated (`kstim*`) MAPK-activation
  rates through two `functions`, replacing the authors' `setParameter(kstim7, keq7)`
  equilibration idiom. Reproduces the authors' trajectories exactly (max rel. diff 0.0).
- Fitted constants are bare `id nominal` declarations at the published values; the conf's
  `*_var` free params bind to them **by name** (ADR-0034, no `__FREE`).

## Free parameters (4) — the p38-arm, fit to the NanoBit data

Per Supplementary Table 2, the p38-ATF2 NanoBit signal determined `keq6, kstim6, dp2, dp4`.
The pp-JNK-arm (`keq7, kstim7, dp1`, fit to pp-JNK WB), `dp3` (fit to pp-ATF2 WB), and the
in-vitro constants are **held at their published values** — refitting them here would break
the paper's dataset decomposition.

| id | published | PyBNF best-fit | ratio | role |
|---|---|---|---|---|
| `keq6`   | 1.74e-5 | 1.98e-5 | 1.14× | basal MKK6 → p38 activation |
| `kstim6` | 1.16e-4 | 1.54e-4 | 1.33× | stim. MKK6 → p38 activation |
| `dp2`    | 1.76e-3 | 2.29e-3 | 1.30× | phosphatase on pp-p38 |
| `dp4`    | 4.50e-3 | 4.52e-3 | 1.00× | phosphatase on ATF2 pS90 (drives p38 eviction = decay) |

## ⚠️ Native-only (not PEtab-v2-exportable)

The data are treated/untreated ratios → `normalization = init`, which PyBNF's `export_job`
cannot express in PEtab v2 (`NotImplementedError`), like `igf1r.conf` and the
`Rukhlenko-2022/cstar_trka` example. Verified **without** the PEtab round-trip.

## Verification (all pass)

- **Tier-1** (`check_conf.py`): edition 2, `job_type=de`, data bound, 4 free params bind by
  id, no `__FREE`. **PASS.**
- **Gate 3a — reproduce at published params** (BNG2.pl): `p38ATF2all` ÷basal reproduces the
  Fig. 4a WT binding curve to **median 5.8 % rel. err** (max 15.3 %), capturing the
  peak-decay transient.
- **Gate 3b — recover by fitting** (`de`, ~19 s): all 4 params recovered within ~1.3× of
  published (table above); best-fit reproduces to **median 4.2 %** (obj 16.7 < published's).
- **Reproduction figure** `p38atf2_binding_reproduction.png` (published vs best-fit vs data).

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Kirsch-2020/p38atf2_binding
pybnf -c p38atf2_binding.conf
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='p38atf2_binding', conf='p38atf2_binding.conf', simulator='ode',
    observables=('p38ATF2all',),
    system='JNK/p38/ATF2 TAD phosphoswitch, WT p38-ATF2 NanoBit binding anisomycin time '
           'course (Kirsch 2020, PMC7666158); ODE, pre-equilibrate -> stimulate, 41 min; '
           'normalization=init -> NATIVE-ONLY (not PEtab-exportable)'),
    # native-only: assert export_job raises NotImplementedError instead of a PEtab-lint test.
```
