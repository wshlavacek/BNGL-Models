# VALIDATION — Zhang-2023/tumor_growth

Primary-source validation of the PyBNF job `pybnf-jobs/Zhang-2023/tumor_growth/`.
Confidence is **earned from the gate evidence below**.

> **Confidence: 75 / 100** — the model structure is the authors' own deposited File S4, the
> experimental design comes from the paper's own figures, all six free parameters are exactly
> the six the Methods name as fitted, and the fit reproduces all four arms of Fig. 4D
> (SSE **93.0 → 5.39**, RMSE **1.65 → 0.40 fold**) while landing on Table S1's reported gate
> (`kTD` 30.5 vs. 29.28, `EC50TD` 0.573 vs. 0.584, `w_OR` 0.34 vs. 0.3). Deductions: **both the
> fit data and the two model inputs are digitized from figures**, not obtained from source
> tables; the paper's own three-arm calibration protocol **does not** work on these data and
> had to be replaced by a four-arm fit (documented below); and `klinear` is not identifiable
> from these data.

Primary sources:
- Paper: Zhang Y, Popel AS, Bazzazi H. "Combining Multikinase Tyrosine Kinase Inhibitors
  Targeting the Vascular Endothelial Growth Factor and Cluster of Differentiation 47 Signaling
  Pathways Is Predicted to Increase the Efficacy of Antiangiogenic Combination Therapies."
  *ACS Pharmacol Transl Sci* 2023; 6:710–726. DOI 10.1021/acsptsci.3c00008.
  `dev/papers/Zhang2023/ACSPharmTranslSci.pdf`.
- Authors' deposited model: Supporting Information **File S4** (`S4.bngl`) and **File S3**
  (`S3.xml`, the SBML export of the same model), in `pt3c00008_si_002.zip`.
- Reported parameter values: Supporting Information **Table S1** (`pt3c00008_si_001.pdf`).
- Fit data: Bridgeman VL et al. *Mol Cancer Ther* 2016; 15:172–183,
  DOI 10.1158/1535-7163.MCT-15-0170, as replotted in Fig. 4D of the primary paper.

"The paper's result" for this job = **the four tumor growth curves of Fig. 4D**.

---

## Gate 0 — Materials inventory

| needed | present? | path / note |
|---|---|---|
| paper PDF | ✅ | `dev/papers/Zhang2023/ACSPharmTranslSci.pdf` |
| authors' model file | ✅ | `dev/papers/Zhang2023/pt3c00008_si_002/S4.bngl` (BNGL) and `S3.xml` (SBML) |
| authors' parameter table | ✅ | `pt3c00008_si_001.pdf`, Table S1 (11 of the 12 growth parameters) |
| authors' fit result | ⚠️ | **inconsistent** — File S4 and Table S1 disagree on 5 parameters, and neither reports `kkill` |
| fit data as a table | ❌ | only plotted, in Fig. 4D — digitized (Gate 1) |
| model inputs (Erk, Akt) as a table | ❌ | only plotted, in Fig. 4B/4C — digitized (Gate 1) |

**Verdict:** PARTIAL — the model structure is unambiguous, the parameterization is not, and
both the data and the two per-arm inputs had to be digitized.

## Gate 1 — Data provenance

All numbers are digitized from the article's own raster figures. The extraction code is
committed and re-runnable in the library-model folder's
[`verify_zhang2023.ipynb`](../../../models/endothelial_vegfr2_and_cd47_signaling_zhang2023/verify_zhang2023.ipynb),
and the extracted values are archived beside it as `reference/zhang2023_fig4*_digitized.csv`.

| `.exp` | source | method | units | n | verdict |
|---|---|---|---|---|---|
| `control.exp` | open circles, Fig. 4D | 1200-dpi render; ticks give 30.90 px/day and 126.0 px/fold; centers by 20–30 px annulus template match on the red mask | fold change of tumor volume | 6 | PASS |
| `sunitinib.exp` | idem, black mask | idem | idem | 11 | PASS |
| `trametinib.exp` | idem, blue mask | idem | idem | 6 | PASS |
| `combo.exp` | idem, magenta mask | idem | idem | 11 | PASS (see note) |

**Calibration check.** The recovered marker abscissae fall on integer days to within 0.06 d
(5.00, 8.00, 12.99, 15.97, 19.98, 22.96, 28.01, 30.99, 34.00, 37.01, 40.99), which validates
the x calibration independently of the y calibration.

**Known gap.** The combination arm has no **day-20** point: its circle is crossed by the day-20
axis tick, and the template match falls below threshold. 34 of a possible 35 points are used.

**Accuracy.** Marker radius is 30 px = 0.24 fold, so a center is good to roughly ±0.05 fold.
This is smaller than the cohort-to-cohort scatter of a xenograft experiment, which the figure
does not show at all (no error bars).

**Model inputs.** `Erk` and `Akt` per arm are the bar heights of Fig. 4B (simulation bars) and
Fig. 4C, each divided by the untreated bar, which is 1 by construction:

| arm | Erk (Fig. 4B sim) | Akt (Fig. 4C sim) |
|---|---|---|
| control | 1.000 | 1.000 |
| sunitinib | 0.736 | 0.793 |
| trametinib | 0.580 | 1.000 |
| combo | 0.355 | 0.793 |

Two internal consistency checks pass: Fig. 4C gives trametinib **exactly** the control value
(a MEK inhibitor does not touch AKT in this model) and gives the combination arm **exactly**
the sunitinib value.

## Gate 2 — Model provenance

`tumor_growth.bngl` is the authors' File S4 with:

1. the **simulation actions removed** (`generate_network`, `writeMexfile`, `writeSBML`) — PyBNF
   synthesizes the run from the conf, and the network is trivially finite (3 species,
   2 reactions), so no `max_stoich`-style directive needs to be retained;
2. `Erk` and `Akt` raised from File S4's placeholder **0** (no signal, hence no growth at all)
   to the untreated **1**, with the treated values applied as `condition:` perturbations;
3. dead code dropped — File S4 carries a commented-out alternative growth law and a
   commented-out alternative S6 law, neither of which is used;
4. house-style formatting, units comments and annotation.

**Rule semantics preserved.** The growth rule is the authors' `TD() -> Cells + TD()`, so
BioNetGen multiplies the growth law by the reactant `TD` in addition to the law's own factors.
This is *not* visible in the Methods equation, which is written as `dCells/dt = kg·Cells·(…)`,
but it **is** explicit in the authors' own SBML export (File S3: `<kineticLaw> rateLaw2 * S2`,
`S2 = TD`), and it is required for the untreated arm to match Fig. 4D. It is preserved here and
documented in the model file.

**Independent check.** The best-fit `.bngl` re-run through BNG2.pl agrees with an independent
SciPy integration of the paper's two equations (written from the Methods text, with the `TD`
factor) to **max |Δ| = 1.4e-6** on the untreated arm.

**Verdict:** PASS with the parameterization caveat of Gate 3.

## Gate 3 — Why the published parameter values are not used as-is

The paper reports the growth model's parameters twice, and the two reports disagree:

| id | File S4 | Table S1 |
|---|---|---|
| `w_OR` | 0.950956148 | 0.3 |
| `kTD` | 2.999985128 | 29.2812 |
| `EC50TD` | 0.000278233 | 0.5840 |
| `kg` | 0.146 | 0.2993 |
| `klinear` | 0.3 | 1.9424 |
| `kkill` | 0 | *absent* |
| `tau6`, `k6`, `r5`, `psi`, `S6t` | agree | agree |

(Table S1 also lists `S6b = 1`, which with `S6t = 1` would make S6 identically 1 and the whole
signal-dependence vanish; File S4's `S6b = 0` is clearly the intended value.)

Measured against the digitized Fig. 4D data (34 points, all four arms):

| parameter set | SSE | RMSE (fold) | which arms agree |
|---|---|---|---|
| File S4 | **93.04** | 1.654 | untreated only (RMSE 0.29) |
| Table S1 with `kkill = 0` | **2622.55** | 8.783 | none |
| this fit | **5.39** | 0.398 | all four |

The two failures have different causes, and both are structural:

- **File S4's gate is saturated.** With `kTD ≈ 3` and `EC50TD ≈ 2.8e-4`, the Hill factor
  `TD^kTD/(TD^kTD + EC50TD^kTD)` is ≈ 1 for every arm's steady-state `TD` (0.75 untreated down
  to 0.13 for the combination). The growth law then reduces to `TD·Cells·(kg − kkill)` in the
  exponential regime, whose **sign is arm-independent**: no single `kkill` can grow the
  untreated tumor and shrink the combination tumor, which is exactly what Fig. 4D shows.
  Circumstantially, `2.999985128` and `0.000278233` appear verbatim in the same File S4 as the
  commented-out `kg` and `taug` of a discarded logistic growth law, so they look like stale
  copy-paste rather than the fitted gate.
- **Table S1 has no `kkill`.** Its `kTD = 29.28`, `EC50TD = 0.584` gate does separate the arms
  (`TD` = 0.75 untreated vs. 0.47 for the combination, straddling `EC50TD`), but with the
  `kkill = 0` of File S4 nothing can shrink, and its faster `kg`/`klinear` overshoot every
  growing arm by an order of magnitude.

**The fit resolves both**, and does so by moving toward Table S1 rather than away from it:

| id | Table S1 | this fit | Δ |
|---|---|---|---|
| `kTD` | 29.2812 | 30.505 | +4 % |
| `EC50TD` | 0.5840 | 0.5733 | −2 % |
| `w_OR` | 0.3 | 0.3423 | +14 % |
| `kg` | 0.2993 (File S4: 0.146) | 0.1388 | File S4 −5 % |
| `klinear` | 1.9424 | 8.944 | not identifiable — see below |
| `kkill` | — | 0.02641 | first value for this parameter |

**`klinear` is not identifiable from these data.** The smooth minimum
`kg·Cells/(1+(kg·Cells/klinear)^psi)^(1/psi)` only binds once `kg·Cells` reaches `klinear`.
At the fitted `kg = 0.1388` and a maximum observed `Cells` of 7.24, `kg·Cells ≤ 1.005`, so any
`klinear` above about 1.3 gives an identical trajectory. The fitted 8.9 should be read as
"the linear cap is never reached", not as an estimate.

**Verdict:** PASS — the discrepancy is a property of the published values, is quantified here,
and the recalibration recovers the reported gate.

## Gate 4 — The paper's own calibration protocol does not work on these data

The Methods state that the six parameters were fitted "to the experimentally measured tumor
growth curves … without treatment, with sunitinib … and with trametinib", and that "the tumor
growth curve treated with the combination of sunitinib and trametinib is used for validation."

That protocol was run first, as a **three-arm calibration with the combination arm held out**
(same conf, `combo` removed from the `experiment:` list). Result:

| | three-arm fit | four-arm fit (this job) |
|---|---|---|
| objective on its own arms | `sos` **1.906** (SSE 3.81, n = 23) | `sos` 2.696 (SSE 5.39, n = 34) |
| control / sunitinib / trametinib SSE | 1.88 / 0.28 / 1.66 | 1.90 / 0.23 / 1.57 |
| **held-out combination arm** | SSE **16.67**, RMSE 1.23 fold | SSE 1.69, RMSE 0.39 fold |
| combination arm at day 41 | predicted **2.15**-fold, observed **0.35** | predicted 0.62, observed 0.35 |
| `kTD` / `EC50TD` | 4.45 / 0.739 | 30.5 / 0.573 |

The three-arm calibration fits its own arms slightly *better* and then predicts the held-out
arm **growing** where the data shrink. The reason is identifiability, not optimization: the
three calibration arms all grow, so all three sit on the same side of the `TD` gate and none of
them constrains where the gate is or how strong `kkill` is. Only the shrinking combination arm
does — and indeed the three-arm fit lands at `kTD = 4.45`, i.e. back near File S4's flat gate,
while the four-arm fit finds Table S1's sharp one.

This job therefore **fits all four arms**. The held-out run is reported here rather than
shipped; it is reproduced by deleting the `combo` experiment line from `tumor_growth.conf`.

**Verdict:** DOCUMENTED DEVIATION — the deviation from the paper's stated protocol is
deliberate, and the evidence for it is above.

## Gate 5 — Job mechanics

| check | command | result |
|---|---|---|
| tier-1 parse / well-formed | `scripts/check_conf.py tumor_growth.conf` | **PASS** — edition 2, `job_type=de`, data bound, 6 free params bind by id, no `__FREE` |
| PEtab.v2 round-trip | `scripts/petab_roundtrip.py tumor_growth.conf --job-type de` | **PASS** — export → `petab.v2` lint clean → import |
| real bngsim fit | `pybnf -c tumor_growth.conf` | **PASS** — finite objective **2.696**, converged at iteration 100 of 300, ~2 min, population 60 |
| heavy? | 3 species / 2 reactions, ms per simulation | **no** |

`sos` in PyBNF is **half** the residual sum of squares, so the printed objective 2.696
corresponds to SSE 5.392 over the 34 points.

## Gate 6 — Reproduction of the paper's result

`make_reproduction.py` simulates all four arms through **BNG2.pl** (not bngsim, so the figure
reproduces without the fitting toolchain) at each of the three parameter sets and writes
`tumor_growth_reproduction.png`.

| parameter set | control | sunitinib | trametinib | combo | ALL (SSE / RMSE fold) |
|---|---|---|---|---|---|
| File S4 | 0.51 | 60.19 | 13.69 | 18.65 | **93.04** / 1.654 |
| Table S1, `kkill = 0` | 916.38 | 553.35 | 1148.61 | 4.21 | **2622.55** / 8.783 |
| this fit | 1.90 | 0.23 | 1.57 | 1.69 | **5.39** / 0.398 |

Per-arm median relative error at the fit: control 0.146, sunitinib 0.068, trametinib 0.123,
combination 1.133. The combination arm's *relative* error is large because its data sit near
0.3 fold; in absolute terms its RMSE is 0.39 fold, comparable to the other arms and to the
±0.05 fold digitization uncertainty plus the unreported cohort scatter.

**Tolerance.** The stated bar is that the fit brings every arm inside ~0.6 fold RMSE, which is
about 2.5 marker radii on the source figure. All four arms clear it; no published set does for
more than one arm.

**Verdict:** PASS.

## Summary

| gate | verdict |
|---|---|
| 0 materials | PARTIAL — model unambiguous, parameterization inconsistent, data only plotted |
| 1 data provenance | PASS — digitized, calibration cross-checked, 34/35 points, method committed |
| 2 model provenance | PASS — File S4 structure preserved, including the `TD` reactant factor |
| 3 published values | PASS — both published sets quantified and shown not to reproduce Fig. 4D |
| 4 calibration protocol | DOCUMENTED DEVIATION — four-arm fit replaces the paper's three-arm holdout |
| 5 job mechanics | PASS — tier-1, PEtab round-trip, real fit, not heavy |
| 6 reproduction | PASS — SSE 93.0 → 5.39 over all four arms |

**Confidence: 75 / 100.**
