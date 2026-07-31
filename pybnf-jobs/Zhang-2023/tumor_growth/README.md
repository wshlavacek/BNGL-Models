# tumor_growth — angiogenesis-driven tumor growth, ODE (PyBNF edition-2 job)

A PyBNF edition-2 parameter-fitting job that fits the **six parameters the Methods name as
fitted** in the angiogenesis-driven tumor growth model of:

> Zhang Y, Popel AS, Bazzazi H. **"Combining Multikinase Tyrosine Kinase Inhibitors Targeting
> the Vascular Endothelial Growth Factor and Cluster of Differentiation 47 Signaling Pathways
> Is Predicted to Increase the Efficacy of Antiangiogenic Combination Therapies."**
> *ACS Pharmacol Transl Sci* 2023; **6**:710–726.
> DOI: [10.1021/acsptsci.3c00008](https://doi.org/10.1021/acsptsci.3c00008)
> — model = Supporting Information **File S4** + the two equations of the *Computational Model
> of Angiogenesis-Driven Tumor Growth* section; fit target = **Fig. 4D**; reported values =
> **Table S1**.

to the mouse xenograft renal cell carcinoma growth curves of:

> Bridgeman VL, Wan E, Foo S, Nathan MR, Welti JC, Frentzas S, Vermeulen PB, Preece N,
> Springer CJ, Powles T, Nathan PD, Larkin J, Gore M, Vasudev NS, Reynolds AR. **"Preclinical
> evidence that trametinib enhances the response to antiangiogenic tyrosine kinase inhibitors
> in renal cell carcinoma."** *Mol Cancer Ther* 2016; **15**:172–183.
> DOI: [10.1158/1535-7163.MCT-15-0170](https://doi.org/10.1158/1535-7163.MCT-15-0170)

Built with the `curate-pybnf-job` skill. The **library-model** sibling — the same model with an
actions block, for reference simulation rather than fitting — is
[`models/endothelial_vegfr2_and_cd47_signaling_zhang2023/…_tumor_growth.bngl`](../../../models/endothelial_vegfr2_and_cd47_signaling_zhang2023/),
alongside the 870-species endothelial signaling model whose output drives it.

> ✅ **The fit reproduces all four arms of Fig. 4D that no published parameter set does:
> SSE 93.0 → 5.39 over 34 points (RMSE 1.65 → 0.40 fold).** And it lands on Table S1's gate:
> `kTD` **30.5** vs. the reported 29.28, `EC50TD` **0.573** vs. 0.584, `w_OR` **0.34** vs. 0.3,
> with `kg` **0.139** vs. File S4's 0.146.

## The model

Two states and two rules.

- `I() -> I() + TD()` with `d_TD()`: a delayed angiogenic driver **TD** relaxes at rate `r5`
  toward a Hill-type activation of ribosomal protein S6 by an ERK-weighted combination of the
  normalized endothelial signals, `w_OR·Erk + (1−w_OR)·Akt`.
- `TD() -> TD() + Cells()` with `g_Cells()`: tumor size **Cells** (a fold change, seeded at 1)
  grows exponentially at `kg`, softly capped at the linear rate `klinear` by a smooth minimum
  of exponent `psi = 20`, gated by `TD^kTD/(TD^kTD + EC50TD^kTD)`, and shrunk by first-order
  kill at `kkill`.

Both rules are in the authors' catalytic form, so BioNetGen multiplies each law by its
reactant — **the growth law carries an extra factor of `TD`** that the Methods equation does
not show. That factor is explicit in the authors' own SBML export (File S3: `rateLaw2 * S2`)
and is what makes the untreated arm come out right.

Network: **3 species, 2 reactions**, milliseconds per simulation. Not heavy.

## What is fit

Four time courses — the four arms of one xenograft study, differing only in the two normalized
endothelial signals that drive the growth model:

| experiment | condition | Erk | Akt | data | days |
|---|---|---|---|---|---|
| `control` | (nominals) | 1.000 | 1.000 | `control.exp` | 5–23 |
| `sunitinib` | `sunitinib` | 0.736 | 0.793 | `sunitinib.exp` | 5–41 |
| `trametinib` | `trametinib` | 0.580 | 1.000 | `trametinib.exp` | 8–28 |
| `combo` | `combo` | 0.355 | 0.793 | `combo.exp` | 0–41 |

`Erk` and `Akt` are the peak endothelial pERK1/2 and ppAKT of the companion signaling model,
normalized to untreated, digitized from the simulation bars of Fig. 4B and Fig. 4C. The fit
observable `Obs_Cells` is a `Molecules` observable (no parens); no per-point uncertainty is
reported, so the objective is plain **`sos`** (PyBNF reports it as half the residual sum of
squares).

**Data provenance.** All 34 points are digitized from the open circles of Fig. 4D: page 7 of
the article rendered at 1200 dpi, axis calibration from the plotted ticks (30.90 px/day,
126.0 px/fold), circle centers by annulus template matching. The extraction is documented and
reproducible in the model folder's `verify_zhang2023.ipynb` and archived as
`reference/zhang2023_fig4D_*_digitized.csv` there. The combination arm has no day-20 point —
its circle is crossed by the day-20 axis tick and the template match falls below threshold.

## Why this job exists

**No published parameter set reproduces Fig. 4D.** The two sources disagree, and each fails
in a different way:

| id | File S4 | Table S1 | this fit |
|---|---|---|---|
| `w_OR` | 0.950956 | 0.3 | **0.3423** |
| `kTD` | 2.999985 | 29.2812 | **30.505** |
| `EC50TD` | 0.000278 | 0.5840 | **0.5733** |
| `kg` | 0.146 | 0.2993 | **0.1388** |
| `klinear` | 0.3 | 1.9424 | **8.944** |
| `kkill` | 0 | *not reported* | **0.02641** |

- **File S4** gets the untreated arm right (RMSE 0.29 fold) and every treated arm wrong. Its
  gate is saturated: with `kTD ≈ 3` and `EC50TD ≈ 2.8e-4`, `TD^kTD/(TD^kTD + EC50TD^kTD) ≈ 1`
  in every arm, so the per-capita growth rate is `TD·(kg − kkill)` and has the **same sign in
  every arm**. No `kkill` can let the untreated tumor grow while the combination tumor shrinks.
  (The File S4 values of `kTD` and `EC50TD` appear verbatim in the same file as the
  commented-out `kg`/`taug` of a discarded logistic growth law, which is where they look to
  have come from.)
- **Table S1**'s sharp gate does separate the arms, but Table S1 reports no `kkill` and File S4
  sets it to 0 — and a growth law with `kkill = 0` cannot shrink a tumor. At `kkill = 0` its
  growth rates are far too fast: SSE 2622.

The fit resolves both: it recovers Table S1's gate (`kTD`, `EC50TD`, `w_OR`) together with a
`kkill` of 0.026 /day and File S4's `kg`.

`klinear` runs up to 8.9 because it is **weakly identified** — with `kg·Cells` never reaching
`klinear` over the plotted range, the smooth minimum never binds, so the linear cap is
unconstrained above about 1.3. Any value above that fits equally well.

## Free parameters (6)

| id | box | prior | fit | role |
|---|---|---|---|---|
| `w_OR` | 0 – 1 | uniform | **0.3423** | weight of the ERK arm vs. the AKT arm |
| `kTD` | 1 – 100 | loguniform | **30.505** | Hill coefficient of the TD growth gate |
| `EC50TD` | 0.01 – 1 | loguniform | **0.5733** | half-maximal TD for the growth gate |
| `kg` | 0.01 – 1 | loguniform | **0.1388** | exponential tumor growth rate (/day) |
| `klinear` | 0.05 – 10 | loguniform | **8.944** | linear tumor growth rate (/day) — weakly identified |
| `kkill` | 0.001 – 1 | loguniform | **0.02641** | first-order tumor kill rate (/day) |

`S6b`, `S6t`, `tau6`, `k6`, `r5` and `psi` stay pinned at their File S4 nominals, as in the
paper. Boxes are wide because the two published sets differ by three orders of magnitude on
`EC50TD`; `de` samples the whole box rather than starting from its center, so a wide box costs
iterations, not correctness. `klinear` is floored at 0.05 to keep `(kg·Cells/klinear)^psi`
inside double precision.

## ✅ PEtab.v2-compliant · not heavy

`sos` over four plain `.exp` time courses with `condition:` perturbations — entirely inside the
PEtab-exportable subset. The round-trip (export → `petab.v2` lint → import) is clean.

## Verification (see [`VALIDATION.md`](VALIDATION.md))

- **Tier-1** (`scripts/check_conf.py`): edition 2, `job_type=de` resolves, data bound, 6 free
  params bind by id, no `__FREE`. **PASS.**
- **PEtab.v2 round-trip** (`scripts/petab_roundtrip.py`): export → lint clean → import.
  **PASS.**
- **Real bngsim fit:** `pybnf -c tumor_growth.conf` converges to a **finite** objective
  **2.696** (`sos`, i.e. SSE 5.392 over 34 points) in about two minutes.
- **Reproduction** (`tumor_growth_reproduction.png`): SSE **93.0 → 5.39**, RMSE
  **1.65 → 0.40 fold**, all four arms; parameters recovered near Table S1's reported values.
- **Independent check:** the best-fit `.bngl` re-run through BNG2.pl agrees with an independent
  SciPy integration of the paper's two equations to **1.4e-6**.

## Run

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder with BNG2.pl
cd pybnf-jobs/Zhang-2023/tumor_growth
pybnf -c tumor_growth.conf        # differential evolution, ~2 min
python make_reproduction.py       # figure + metrics: File S4 vs Table S1 vs this fit
```

## `_manifest.py` entry (if promoted to the PyBNF real-world corpus)

```python
RealWorldExample(
    folder='tumor_growth', conf='tumor_growth.conf', simulator='ode',
    observables=('Obs_Cells',),
    system='Angiogenesis-driven tumor growth driven by normalized endothelial pERK1/2 and '
           'ppAKT (Zhang 2023, DOI 10.1021/acsptsci.3c00008, Fig. 4D); 2-state ODE with a '
           'delayed angiogenic driver and a Hill-gated exponential/linear growth law; four '
           'xenograft arms as conditions on the two signal inputs; sos on fold-change tumor '
           'volume digitized from Fig. 4D (Bridgeman 2016). PEtab-exportable. GLOBAL (de) fit '
           'of the six parameters the paper reports as fitted; no published set reproduces '
           'Fig. 4D and this one does (SSE 93.0 -> 5.39). See VALIDATION.md.'),
    stochastic=False,
    recover={'kTD': 29.2812, 'EC50TD': 0.5840}, tol=0.10,
)
```
