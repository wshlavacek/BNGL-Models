# VALIDATION — Suofu-2017/mito_camp

Primary-source validation of the PyBNF job `pybnf-jobs/Suofu-2017/mito_camp/`, per the
`validate-pybnf-job` skill. Confidence is **earned from the gate evidence below**.

> **Confidence: 65 / 100** — the reaction scheme, geometry and initial quantities are the
> paper's own (Supplementary Data Tables 1 and 2), the reconstruction is checked against a
> hand-written integration of that table to **3.8e-8**, the job mechanics all pass
> including a clean PEtab.v2 round-trip, and the fit **beats the authors' own plotted fit
> on the authors' own data** (SSE 4405.6 → 3684.6 over the same 77 points). Deductions:
> the authors' model file is **not deposited**, so there is no Gate-2 diff — only a prose
> specification plus two tables; **all 77 fit points are digitized from raster figures**;
> and Gate 3b (recovering the published parameters) is **impossible in principle**, because
> the paper reports its fit only as unit-less log10 histograms. Five of the fourteen free
> parameters are unidentified by these data.

Primary sources (in the untracked `dev/papers/Soufu2017/`; not redistributed — note the
folder name transposes the first author's name, which is Suofu):
- Model + data paper: Suofu Y, Li W, Jean-Alphonse FG, Jia J, Khattar NK, Li J, *et al.*
  "Dual role of mitochondria in producing melatonin and driving GPCR signaling to block
  cytochrome c release." *Proc Natl Acad Sci USA* 2017; 114:E7997–E8006.
  DOI 10.1073/pnas.1705768114, PMCID PMC5617291. `PNAS.pdf`.
- Supplement: `pnas.1705768114.sapp.pdf` — SI Appendix *Model constructions and analysis*,
  Supplementary Data Table 1 (reactions), Table 2 (initial quantities), Fig. S7
  (posteriors).
- Authors' model/job files: **none deposited**.

"The paper's result" for this job = the four simulated curves of **Fig. 4J**, fitted to the
measured time courses of **Fig. 4H** and **Fig. 4I**. There is no parameter table.

---

## Gate 0 — Materials inventory

| needed | present? | path / note |
|---|---|---|
| model paper PDF | ✅ | `dev/papers/Soufu2017/PNAS.pdf` |
| SI / supplementary tables | ✅ | `pnas.1705768114.sapp.pdf` — Tables 1 and 2, Fig. S7, the modelling section |
| authors' model file (`.bngl`/SBML) | ❌ | not deposited; the model is reconstructed from the prose spec + Tables 1 and 2 |
| authors' fitted parameter values | ❌ | reported only as twelve log10 posterior histograms (Fig. S7), with no unit convention |
| fit data as a table | ❌ | only plotted, in Fig. 4H / Fig. 4I — digitized (Gate 1) |
| authors' simulated curves as a table | ❌ | only plotted, in Fig. 4J — digitized, used as a comparator (Gate 5) |

**Verdict:** PARTIAL — the model *specification* is complete and unambiguous, but nothing
executable and no numeric parameterization was published.

## Gate 1 — Data provenance

All 77 points are digitized from the article's own raster figures. The extraction code is
committed and re-runnable in the library-model folder's
[`verify_suofu2017.ipynb`](../../../models/mitochondrial_mt1_camp_signaling_suofu2017/verify_suofu2017.ipynb);
re-running it against the PDF reproduces the committed CSV **exactly** (max difference 0).
The extracted values are archived beside it as
`reference/suofu2017_fig4HI_camp_timecourses_digitized.csv`.

| `.exp` column | source | method | units | n | verdict |
|---|---|---|---|---|---|
| `melatonin.exp` `camp_CY_pct()` | Fig. 4H, blue markers (PM sensor) | 1200-dpi render; axis calibration from the plotted ticks (61.2 px/min, 9.12 px per percentage point); centres by matched filter with a 16 px disk + greedy peak picking at 22 px exclusion | % of pre-agonist max | 29 | PASS |
| `melatonin.exp` `camp_ML_pct()` | Fig. 4H, black markers (OMM sensor) | idem | idem | 19 | PASS |
| `damgo.exp` `camp_CY_pct()` | Fig. 4I, blue markers | idem | idem | 19 | PASS |
| `damgo.exp` `camp_ML_pct()` | Fig. 4I, black markers | idem | idem | 10 | PASS |

**Why not connected components.** The markers are filled disks ~36 px across, drawn on top
of the mono-exponential line Prism fits through them and, in the plateau, spaced ~24 px
apart — closer than their own diameter. Labelling the colour mask merges up to a dozen of
them into one blob (a first attempt recovered 2 of the 19 black markers in Fig. 4H).
Matched filtering plus greedy non-maximum suppression separates them.

**Visual check.** The notebook plots the picked centres over the source panels. Every pick
sits on a marker — no false positives from the fitted line, the error bars (which are a
paler colour and fall outside both masks) or the axis furniture.

**Known gap.** Coverage is not provably complete: two or three of the most crowded blue
markers in the Fig. 4H plateau may be absorbed by a neighbour's exclusion disk. With 29
points on that curve and a plateau that is flat, the effect on the fit is negligible.

**Accuracy.** Marker radius is 18 px = 2.0 percentage points, so a centre is good to about
±0.5 pp. That is far smaller than the SEM the figure plots (bars spanning 10–30 pp in the
plateau), which is the dominant uncertainty and is *not* carried into the `.exp` — the bars
overlap too densely to attribute per point. Consequence: the objective is `sos`, not
`chi_sq`, and it weights the noisy plateau as heavily as the well-determined onset.

**Normalization.** The figures plot "cAMP (% of max)", i.e. each averaged trace divided by
its own pre-agonist maximum. The observables `camp_CY_pct()` / `camp_ML_pct()` are cAMP as
a percentage of the **Supplementary Data Table 2** basal quantity, and the experiments are
pre-equilibrated receptor-silent, so matching the data also forces the model's basal state
onto Table 2's level. This is stricter than a free scale factor and keeps the job
PEtab-exportable; see Gate 5 for how well the model manages it.

## Gate 2 — Model fidelity

Reference compared against: the SI Appendix prose specification + Supplementary Data
Table 1 + Table 2. **No author file exists to diff against**, so `model_diff.py` cannot be
run; the substitute is an independent re-implementation.

| aspect | paper | our `.bngl` | verdict |
|---|---|---|---|
| compartments | 1 pL spherical cell, 100 mitochondria of radius 0.3 µm; cytosol, mitochondrial lumen, PM, OMM | `EC`/`PM`/`CY`/`OMM`/`ML`, V_cy 1000 µm³, V_ml 11.3097 µm³, A_pm 483.6 µm², A_omm 113.1 µm² | match (`EC` added: BioNetGen needs a 3D root for `PM`) |
| reactions | 8 types, Supplementary Data Table 1 | 20 rules = the 8 types written per compartment pair, + nucleotide permeation | equivalent |
| rate constants | 12, named but not valued | same 12 names | names match; values refit (Gate 3) |
| initial quantities | Supplementary Data Table 2 | same, with luminal = cytosolic × `phi` (which is what Table 2's rounded values are) | match |
| species | "30 molecular species … 30 ODEs" | **26** | see below |
| observables ↔ measured | cytosolic cAMP ↔ PM sensor; luminal cAMP ↔ OMM sensor | `camp_CY_pct()`, `camp_ML_pct()` on free cAMP | match |
| network cap | none needed | none; the network is finite | match |

**Four deliberate deviations, all documented in the model's `#@note:`:**

1. **26 species, not 30.** The four not represented are µ-OR and its Gi complex *on the
   OMM* — which the paper's model carries at zero copies, and which the
   compartment-explicit rules never create — and the two agonists, which are applied as a
   switch on receptor activity rather than as bound ligand. The paper does the same thing
   ("plus 10⁶ *active* copies of µOR on the PM"), so the agonists are bookkeeping there
   too.
2. **ATP, cAMP and AMP are three states of one molecule type `Nuc`.** This conserves the
   adenine-nucleotide pool by construction and lets one rule transport all three across the
   OMM. The generated network is equivalent to three separate types.
3. **Nucleotide permeation obeys detailed balance across the volume ratio.** The paper
   names one constant, `kdiff`. It is implemented as the efflux (ML → CY) constant, with
   influx = `kdiff × phi`, so the compartments equilibrate at equal *concentration*. Equal
   constants both ways would equilibrate at equal *count*, i.e. 88× the cytosolic
   concentration inside the mitochondrion, contradicting the paper's stated design that
   "protein concentrations in the ML were identical to those in the cytosol". The
   equilibration protocol confirms the intended behaviour: `cAMP_ML / (phi · cAMP_CY)` =
   1.000000 at basal.
4. **AC.ATP.Gi\* is never formed.** Rule 1 requires AC's Gi site free and rule 8 requires
   its substrate site free, so Gi\* binding is competitive with ATP. Table 1 lists AC.ATP
   and AC.Gi\* as the only AC complexes, and inhibition is the point of rule 8.

**Independent check.** The library-model sibling, run through BNG2.pl at these parameters,
agrees with a **hand-written** 26-state SciPy integration of Supplementary Data Table 1 —
written from the table, not parsed from the generated `.net` — to **max relative error
3.8e-8** across all three protocols and nine observables. This is the check that
BioNetGen's compartmental volume scaling (divide a surface–volume bimolecular rate by the
*volume* compartment's size) is doing what the paper's construction assumes.

**Verdict:** PASS (equivalent) — every deviation is named, and the reconstruction is
verified against an independent implementation of the published reaction table.

## Gate 3 — The published parameter values cannot be used, or recovered

**Gate 3a (reproduce at the paper's parameters) is not executable, and Gate 3b (recover
them by fitting) is not decidable.** Suofu et al. fit by replica-exchange MCMC (eight
chains, 10,000 swaps) and report the result *only* as twelve log10 posterior histograms,
Supplementary Data Fig. S7. There is no table, and no statement of the unit convention
behind the axes — with molecule counts and BioNetGen compartmental scaling, a second-order
constant differs by whatever the compartment size is in whatever unit they chose, and the
paper does not say. So the histograms cannot be read back into a runnable model.

They are also not a point estimate:

| Fig. S7 panel | prior box (red outline) | posterior |
|---|---|---|
| `log_k1f` | −11.5 … −4 | narrow spike at ≈ −11.9, **against the lower edge** |
| `log_k2` (panel 2) | −3 … 0 | ≈ −1.5 |
| `log_k2` (panel 3 — same label; presumably `log_k1r`) | −3 … 0 | ≈ −1.5 |
| `log_k3f` | −8 … −4 | narrow spike at ≈ −7.6, **against the lower edge** |
| `log_k3r` | −8 … 1 | bimodal, ≈ −2.7 and ≈ +0.5 |
| `log_k4` | −2 … 6 | bimodal, ≈ 2 and ≈ 4.5 |
| `log_k5` | −7 … 1 | bimodal, ≈ −6 and ≈ −1 |
| `log_k6` | −8 … 0 | bimodal, ≈ −7 and ≈ −6 |
| `log_k7` | −2 … 6 | ≈ −1.7, against the lower edge, long right tail |
| `log_k8f` | −2 … 6 | ≈ 4.3 |
| `log_k8r` | −2 … 6 | bimodal, ≈ 3 and ≈ 4.7 |
| `log_kdiff` | −3 … 3 | narrow spike at ≈ −2.6, **against the lower edge** |

(Modes read off the rendered figure; two axes carry the same label.) Three parameters are
pinned against a prior bound and four span a whole prior box.

**A suggestive but inconclusive comparison.** If the published axes are read as
*count-based* constants (rate = k·N₁·N₂, no volume division), our fitted second-order
constants map as log10(k/V_cy). On that reading four of the twelve land within about a
decade of the published mode — `k3f` −7.35 vs ≈ −7.6, `kdiff` −2.02 vs ≈ −2.6, `k2` −0.60
vs ≈ −1.5, `k1r` −2.52 vs ≈ −1.5 — and the rest do not, `k8f` by nine decades. Four out of
twelve is not evidence of a shared convention; it is recorded so a reader with the authors'
code can check.

**Verdict:** NOT DECIDABLE — recorded, not scored. This is the single largest reason the
confidence here is 65 and not 80.

## Gate 3c — A contradiction in the paper, resolved by the fit

The paper states the OMM receptor number twice, and the two statements differ by
`A_omm/(A_pm·phi)` = **20.68**:

| source | MT1 on the OMM | reading |
|---|---|---|
| Supplementary Data Table 2 | 1.13e4 copies | = 1e6 × `phi`: the same **concentration** as a PM carrying 1e6 copies |
| SI Appendix text | "2068 MT1 per square µm of OMM" = 2.34e5 copies | the same **areal density** as that PM (1e6 / A_pm = 2068 /µm², which is where the number comes from) |

Which is operative decides how much more efficacious mitochondrial MT1 is than
plasma-membrane µ-OR — the paper's central quantitative claim. The job boxes the two
readings as `MT1_OMM_amp` ∈ [1, 20.68] and fits it. Result: **1.130**, with all 1268
archived sets within 1% of the best objective lying in [1.00, 1.27]. The data pick
**Table 2**, and exclude the SI text's reading by more than an order of magnitude.

Read as densities, the fit puts 8.2e4 MT1 on the PM (170 /µm²) and 1.28e4 on the OMM
(113 /µm²) — a concentration ratio of 0.073 but comparable areal densities, which is the SI
text's qualitative picture at ~20× lower absolute density.

**Verdict:** PASS — the ambiguity is quantified and settled by the data.

## Gate 4 — Job mechanics

| check | command | result |
|---|---|---|
| tier-1 parse / well-formed | `scripts/check_conf.py mito_camp.conf` | **PASS** — edition 2, `job_type=de`, data bound, 14 free params bind by id, no `__FREE` |
| PEtab.v2 round-trip | `scripts/petab_roundtrip.py mito_camp.conf --job-type de` | **PASS** — export → `petab.v2` lint clean → import |
| committed PEtab bundle | `make_petab.py` | **PASS** — `petab/`: 2 observables, 14 estimated parameters, 77 measurements, 6 condition rows, 2 experiments; lints clean |
| real bngsim fit | `pybnf -c mito_camp.conf` | **PASS** — finite objective **1842.28**, 300 of 300 DE generations × 130 individuals ≈ 39,000 parameter sets, ~28 min on 16 cores, **0 failed simulations** |
| heavy? | 26 species / 30 reactions, ms per simulation | **no** |

`sos` in PyBNF is **half** the residual sum of squares, so 1842.28 corresponds to
SSE 3684.56 over 77 points (RMSE 6.92 percentage points).

**Two simulation paths agree.** `make_reproduction.py` drives BNG2.pl and reaches the basal
state by a 10⁷ s unstimulated integration, where PyBNF uses bngsim's steady-state solve for
`preequilibrate:`. Scored on the same points they give SSE 3684.6 and 2 × 1842.28 =
3684.56 — agreement to five digits, which validates the `preequilibrate:` path against an
ordinary integration.

## Gate 5 — Reproduction of the paper's result

`make_reproduction.py` runs both stimulations through BNG2.pl at the fitted parameters and
writes `mito_camp_reproduction.png`. Two comparisons.

**(a) Against the measurements, scored on the same 77 points as the authors' own curves:**

| curve | n | RMSE, this fit (pp) | RMSE, Suofu Fig. 4J (pp) |
|---|---|---|---|
| melatonin, PM sensor (cytosol) | 29 | 6.95 | 5.70 |
| melatonin, OMM sensor (lumen) | 19 | 7.10 | 6.27 |
| DAMGO, PM sensor (cytosol) | 19 | **6.42** | 10.14 |
| DAMGO, OMM sensor (lumen) | 10 | **7.37** | 8.72 |
| **ALL** | **77** | **6.92** (SSE 3684.6) | 7.56 (SSE 4405.6) |

The two fits trade off: the authors' does better on the melatonin curves, this one does
better on both DAMGO curves, and the total favours this one by 16% in SSE.

**(b) Against the authors' plotted trajectories, Fig. 4J**, in fraction-of-basal units:

| curve | max \|Δ\| | median \|Δ\| | this fit @600 s | Fig. 4J @600 s |
|---|---|---|---|---|
| melatonin, PM | 0.304 | 0.234 | 0.536 | 0.441 |
| melatonin, OMM | 0.161 | 0.129 | 0.280 | 0.177 |
| DAMGO, PM | 0.141 | 0.066 | 0.009 | 0.034 |
| DAMGO, OMM | 0.216 | 0.175 | 0.496 | 0.299 |

Digitization of Fig. 4J is worth about ±0.02 (line half-width ~8 px against 912 px per unit
fraction, plus the grey credible band overdrawing the black curves early), so these
differences are real, not extraction noise. They are two fits of a sloppy model to the same
data, and neither targets the other. **The qualitative result the paper draws is
reproduced**: under melatonin the lumen loses more cAMP than the cytosol, and under DAMGO
the order reverses.

The largest single disagreement, the DAMGO luminal plateau (0.50 here vs. 0.30 plotted), is
mechanistic rather than accidental. With no receptor on the OMM, luminal cAMP under DAMGO
falls only as fast as the emptying cytosol drains it, so the same `kdiff` that lets the
lumen reach 0.18 under melatonin — where the cytosol stays near 0.44 and props it up — also
keeps it near 0.5 under DAMGO. A search over 39,000 parameter sets did not find a set that
reaches both plotted plateaus, and the measurements side with this fit on that curve
(RMSE 7.4 vs 8.7 pp).

**Where the model strains: the basal level.** The receptor-silent steady state comes out at
**7.61e6** free cytosolic cAMP against Supplementary Data Table 2's 1.0e7 — 24% low — so
every simulated curve starts near 76% rather than 100%. Since Table 2's cAMP is the
observable's fixed reference, this cost is paid in the objective rather than hidden by a
scale factor, and it is visible as the offset at t = 0 in both reproduction figures. The
digitized data's own earliest points are at 83–94%, so part of the gap is the experimental
dead time between switching the perfusion and the first frame; the rest is the model's.

**Tolerance.** The bar set here is that the fit must score no worse on the 77 measured
points than the authors' own published curves do. It clears it (6.92 vs 7.56 pp RMSE),
which is the strongest statement available given that the published parameter values do not
exist in usable form.

**Verdict:** PASS.

## Gate 6 — Identifiability

Range covered by the 1268 archived parameter sets within 1% of the best objective:

| well determined (< ½ decade) | `k3f` 0.01 · `k2` 0.04 · `kdiff` 0.05 · `MT1_OMM_amp` 0.10 · `ratio_MT1` 0.18 · `k6` 0.21 · `k1f` 0.41 · `k8f` 0.43 |
|---|
| **not determined (> 1 decade)** | `k1r` 3.3 · `k3r` 2.6 · `k5` 1.2 · `k7` 1.2 · `k8r` 1.2 |

The undetermined five are all steps whose rate only has to be *fast enough* — the reverse
legs of the two enzyme-substrate bindings, AMP recycling, G-protein release and AC.Gi\*
dissociation — so the objective is flat above a threshold and the fitted values should be
read as lower bounds, not estimates. The authors' own posteriors are broad in overlapping
places (`k3r`, `k4`, `k5` bimodal across a whole prior box), so this is a property of the
model and the data, not of the search.

**Verdict:** DOCUMENTED — nine of fourteen parameters are pinned; five are not, and are
labelled as such in the job README.

---

## Divergence & corrections

- Scope vs. paper: matches. Both simulations of Fig. 4J are fitted jointly, as the paper
  did ("model parameters were simultaneously fit to data from both experiments").
- One free parameter is added beyond the paper's twelve rate constants + MT1 ratio:
  `MT1_OMM_amp`, which resolves a contradiction *in the paper* rather than relaxing it
  (Gate 3c). It is boxed between the two published readings, so it cannot escape them.
- Corrections applied to job files: none after the first successful run; the only mid-course
  change was adding `MT1_OMM_amp` after Gate 3c was discovered, and the fit was restarted
  from scratch.

## Bottom line

Solid: the reaction scheme and its BioNetGen implementation (verified to 3.8e-8 against a
hand-written integration of the published table), the job mechanics, the digitization, and
a fit that outperforms the authors' own plotted curves on the authors' own data while
resolving a 20.7-fold ambiguity in their reported receptor density. Residual risk: the
model was reconstructed from prose because nothing was deposited, so an undetected
structural difference from the authors' file cannot be excluded; the fit is not unique
(five parameters unidentified); and the basal cAMP level sits 24% below the tabulated
value. The single most valuable next step would be a `job_type = pt` run reproducing the
authors' replica-exchange posteriors rather than a point estimate — identifiability, not
the point estimate, is what this model's parameterization is really about.

| gate | verdict |
|---|---|
| 0 materials | PARTIAL — complete specification, nothing executable, no numeric parameters |
| 1 data provenance | PASS — 77 points digitized, method committed and re-runnable, visually checked |
| 2 model fidelity | PASS (equivalent) — 4 named deviations; independent SciPy check to 3.8e-8 |
| 3 published values | NOT DECIDABLE — reported only as unit-less log10 histograms |
| 3c OMM density contradiction | PASS — 20.68-fold ambiguity resolved in favour of Table 2 |
| 4 job mechanics | PASS — tier-1, PEtab round-trip, 39,000-set fit, 0 failures, not heavy |
| 5 reproduction | PASS — SSE 4405.6 → 3684.6 vs. the authors' own curves on the same points |
| 6 identifiability | DOCUMENTED — 9 of 14 pinned, 5 unidentified |

**Confidence: 65 / 100.**
