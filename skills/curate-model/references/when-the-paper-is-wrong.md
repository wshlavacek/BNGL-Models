# When the paper is wrong

Curation is reproduction, so the interesting cases are the ones that will not reproduce. A
published number that violates the paper's own conservation law, a distributed model file that
disagrees with its own table, an unstated convention that changes the answer by a factor of two —
these are not obstacles to the curation. **They are the curation's most valuable output**, because
nobody else has checked, and the collection is often the first place the discrepancy is recorded.

This document fixes the doctrine, the triage, and the house wording. It is generalized from a
dozen worked cases in `models/`, indexed in §5.

## Contents
1. Doctrine
2. Triage — what kind of wrong is it?
3. The five responses
4. Wording
5. Worked cases

---

## 1. Doctrine

Four rules, in priority order.

1. **Reproduce what was published; do not improve it.** The deliverable is the published model,
   not the best model. If it disagrees with the paper's own figure, that disagreement is the
   result.
2. **Never silently repair anything.** A changed rate, a corrected rule, a chosen convention — each
   is recorded in the `.bngl` header, `metadata.yaml`, the notebook, and the README row. A reader
   must be able to recover exactly what the paper said from what you shipped.
3. **Prefer the smallest repair that closes the check.** Mu 2010's X5P pool balance does not close
   at the precision the fluxes are printed to; the primary file carries a two-value correction and
   the notebook shows it is "the smallest repair available". Establishing minimality is part of
   justifying the repair.
4. **Do not resolve an ambiguity by fiat — ship both readings.** If the paper genuinely permits two
   interpretations and they differ observably, that is a two-file collection, not a judgement call
   (§3, "bracket it").

## 2. Triage — what kind of wrong is it?

Work down this list. The response follows from the category, not from taste.

1. **Is it *you*?** Assume so until the deterministic check passes. An independent implementation
   agreeing to ~1e-6 is what earns you the right to call a residual disagreement the paper's.
2. **Does it violate something the paper itself asserts?** A conservation law, a detailed-balance
   loop, a stated steady state, a mass balance. This is the strongest possible evidence, because
   the paper is refuting itself — no external judgement is involved. *(Samoilov's N; Blinov's
   detailed-balance loops; Mu's X5P balance.)*
3. **Do two published artifacts disagree?** Paper text vs. deposited model, table vs. distributed
   file, main text vs. SI, Eq. 15.16 vs. Eqs. 15.27–15.30. Then the question is only *which one to
   follow*, and you must say why. *(Michalski; Blinov; Zhang; Mu.)*
4. **Is something simply unstated?** An initial condition, a propensity convention, a fitting grid.
   Not an error — an underspecification. It needs evidence or brackets, not a correction.
   *(McKane's ICs; Samoilov's homodimer counting; Rohrs's Prism sample times.)*
5. **Is the published result unreachable from the published values?** Then say exactly that, and
   demonstrate it across the plausible parameter sets rather than asserting it. *(Zhang's treated
   arms.)*

## 3. The five responses

### Correct it — but only on two independent lines of evidence

Change a distributed value **only** when two things independently point the same way. Precedent:
`combinatorial_egfr_signaling_blinov2006` sets `km21 = 0.1 s⁻¹` (Table 1's `k_-21`) against the
`0.01` in the distributed BNGL, because (a) the correction restores detailed-balance loops `loop2`
and `loop4` to unity, and (b) the digitized panel-C steady state (~64 nM) independently matches the
corrected value, not the uncorrected one (~77 nM). One line of evidence would have been a guess;
two make it a finding. The header, metadata and README all name the correction.

Where the correction defines a distinct published model rather than a repair to the primary one,
carry it in the **variant** — `scaffolded_mapk_cascade_kocieniewski2012_scaff11.bngl` corrects the
Appendix A rule-7 typo (an unphosphorylated left-hand side that abolishes the prozone effect), and
the corrected rule reproduces the paper's own Fig. 2.

### Document it and preserve the model

When the discrepancy is in a *reported number* rather than the mechanism, keep the mechanism and
record the number. `noise_induced_bistable_futile_cycle_samoilov2005` reproduces six of the seven
values of the Supporting Text fixed point to ≤ 5e-5; the seventh, `N = 2.213`, violates the paper's
own conserved pool `E+ + C+ + N = 30`, which `N = 2.132` satisfies exactly. That is reported as a
transcription error in the paper. Nothing in the model changes.

Likewise a **uniform transformation** between the figures and the model: `mckane2005` finds the
published figures drawn on a time base about twice as slow as the model's, shows the same factor in
*both* Fig. 1 and Fig. 2, confirms every verified quantity is invariant under it, and preserves the
supplied propensities. Two tests make a transformation a finding rather than a fudge — it must
appear in **more than one** figure, and the quantities you verify must be invariant under it.

### Bracket it — ship both readings

When the paper leaves a convention unstated and the readings differ observably, ship a file for
each and let the evidence assign them. `samoilov2005` reads the homodimeric driver step
`N + N -> E+ + N` as `k21·N·(N−1)` in the primary file — the reading that matches Eqs. S5 and the
reported fixed point — and as half that in `_unordered_pair.bngl`, the reading whose switching
statistics and stationary modes match the digitized Fig. 3A. Neither is "the" answer; the pair is,
and the README says which evidence pins which.

### Surface it and leave the published rules alone

Where the published rules have a semantic consequence the authors may not have intended, state it
and do not rewrite them. `shp2_regulation_and_function_barua2007`: occlusion suppresses
C-SH2/PTP avidity and so is "more restrictive than a true N-SH2 deletion. We surface this
discrepancy rather than alter the published rules."

Use the same response for an underspecified *analysis* choice. `car_cd3zeta_phosphorylation_rohrs2018`:
"The paper does not say which sample times were entered into Prism, and that choice matters" —
refitting the same curve on a log grid, a linear grid, or the reported points moves t½. The
notebook quantifies the spread instead of picking one and calling it agreement, and marks the Fig. 3
sequential-order summaries **not reproduced** from Data S1.

### Declare it unreachable — and prove that

The strongest claim, so it needs the most work: show the target cannot be reached from *any* of the
published parameter sets, not just the one you tried. `endothelial_vegfr2_and_cd47_signaling_zhang2023`:
"The three treated arms of Fig. 4D cannot be reproduced with either published parameter set, and
this is a property of the published values rather than of the curation" — established across the
File S4 and Table S1 values, which themselves disagree.

## 4. Wording

The house style is flat and falsifiable. Name the artifact, name the quantity, give the number,
say what you did.

- **Attribute precisely.** "the reported 2.213 violates the paper's own conserved pool
  `E+ + C+ + N = 30`" — not "the paper is wrong here".
- **Quantify the disagreement.** A discrepancy with a number attached is a result; without one it
  is a complaint. `insulin_signalling_and_oxidative_stress_smith2013` records "the one genuine
  discrepancy" as an overshoot "by up to about 0.03 of full scale between sixteen and thirty-five
  minutes, while matching its shape and its late decay."
- **Say when something is not reproduced.** "not reproduced (documented)" is an acceptable and
  honest outcome. Burying it is not.
- **Do not speculate about cause.** Record what disagrees and by how much. Whether it was a typo, a
  different code path, or an undocumented step is usually unknowable from outside — and an
  attribution you cannot support is worse than none, because it reads as a diagnosis. The canonical
  case is `lambda_switch_arkin1998`, whose full circuit saturates to ~100% lysogeny at high MOI
  against Arkin's Fig. 6a "Full" plateau of ~82%. An earlier draft attributed this to the
  mean-field promoter treatment; that attribution was **removed** because there was no evidence for
  it and Arkin's original code — the only thing that could diagnose it — cannot be found. The
  shipped wording is now "high-MOI over-prediction (cause unidentified)".
- **Separate "unexplained" from "unreasonable".** The same model records that ~100% saturation is
  a *supported* behavior: Cortes et al. (2017) likewise report a lysogeny probability saturating at
  1, in agreement with experiment, so the lower published plateau may be the outlier. Noting that a
  discrepancy is defensible is not the same as explaining it, and both belong in the record.
- **Do not soften a real match either.** If it reproduces, say so plainly.

Every discrepancy lands in four places: the `.bngl` `#@note`, `metadata.yaml`, the verification
notebook (with the number), and the README Models-table row. The README row is what a reader sees
first and is where the one-sentence version belongs.

## 5. Worked cases

| model | category | what was found |
|---|---|---|
| `noise_induced_bistable_futile_cycle_samoilov2005` | violates own law; unstated convention | reported `N = 2.213` breaks the paper's conserved pool (`2.132` satisfies it); homodimer counting bracketed by two files |
| `combinatorial_egfr_signaling_blinov2006` | distributed file vs. table | `km21` `0.01` → `0.1`; confirmed by detailed balance *and* digitized panel C |
| `scaffolded_mapk_cascade_kocieniewski2012` | typo in published rules | Appendix A rule 7 had an unphosphorylated LHS abolishing the prozone effect; corrected in the Scaff-11 variant |
| `demographic_noise_predator_prey_cycles_mckane2005` | uniform transformation; unstated IC | figures on a ~2× slower time base than the model, in both Fig. 1 and Fig. 2; Fig. 1 ICs never stated, 800/800 confirmed from the transient dip depth to ~±25 |
| `endothelial_vegfr2_and_cd47_signaling_zhang2023` | unreachable from published values | Fig. 4D treated arms reproduce from neither File S4 nor Table S1, which disagree with each other |
| `camkii_holoenzyme_activation_michalski2012` | text vs. deposited model | whether `Dpp` as well as `Dpu` drives autophosphorylation at `r3`; the difference is quantified rather than decided |
| `car_cd3zeta_phosphorylation_rohrs2018` | underspecified analysis | Prism sample times unstated and outcome-relevant; Fig. 3 sequential-order summaries not reproduced from Data S1 |
| `shp2_regulation_and_function_barua2007` | semantic consequence of published rules | occlusion is more restrictive than a true N-SH2 deletion; surfaced, rules untouched |
| `insulin_signalling_and_oxidative_stress_smith2013` | genuine model-data mismatch | Fig. 2C IRS1-Yp overshoot quantified and reported |
| `calvin_cycle_isotopomer_labeling_mu2010` | internal inconsistency; imprecise printing | Eq. 15.16 vs. Eqs. 15.27–15.30 on the `S → M₂` carbon transposition; X5P balance does not close at printed precision, repaired minimally |
| `lambda_switch_arkin1998` | unexplained standing discrepancy | ~100% lysogeny at high MOI vs. the ~82% Fig. 6a plateau; an unsupported mechanistic attribution was removed, leaving "cause unidentified", and the saturation independently noted as defensible against Cortes et al. (2017) |
