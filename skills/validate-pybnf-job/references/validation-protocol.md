# Validation protocol — the five gates, pass criteria, and the confidence rubric

The spine of `validate-pybnf-job`. Each gate has a definition, a **pass criterion**, the
**evidence** to record in `VALIDATION.md`, and how it feeds the confidence score. Run the gates in
order; a hard failure early (missing sources, wrong data, wrong model) makes later gates moot.

Every gate resolves to one of: **PASS** · **PASS (equiv)** — functionally equivalent, not identical
· **PARTIAL** — mostly grounded, a documented gap · **FAIL** — contradicts the source · **BLOCKED**
— cannot assess (source or run missing). Record the verdict *and its evidence* (file:line, table
cell, figure panel, metric, command).

---

## Gate 0 — Materials inventory

**Definition.** Confirm the parking garage (`dev/papers/<FirstAuthor>-<Year>/`) holds the primary
sources this job needs: the **model paper**, the **data paper** (may differ — e.g. tlbr: model
Monine, data Posner), the **SI/Source-Data**, and ideally the **authors' own model/job files**.

**Pass criterion.** Every source the job depends on is physically present and readable.

**Evidence.** A checklist: each needed item → present? path? Note which are the authors' original
files (enables a direct Gate-2 diff) vs. only prose specs.

**Score effect.** Any missing *required* source → the gates that depend on it are **BLOCKED**, and
the overall confidence is **capped low** (see rubric). Do not paper over a missing source.

---

## Gate 1 — Data provenance

**Definition.** Every column/value in each `.exp` traces to a primary source.

**Pass criterion.**
- **Tabular source:** the `.exp` values match a paper/SI table (transcription) within rounding.
- **Figure source:** re-digitizing the cited panel from the PDF reproduces the `.exp` within a
  stated tolerance (see `data-provenance.md`). Record panel, axis calibration, any normalization /
  unit / scale factor, and the digitization tolerance.
- Independent variable, observable mapping, units, and any `_SD` all match the source.

**Evidence.** A provenance table: `.exp` column → source (table N / Fig Np) → method (transcribed /
digitized) → diff metric. Use `scripts/compare_data.py <redigitized.csv> <shipped.exp>` for the
figure case (reports max/median abs & rel difference, aligned on the independent variable).

**Common failure modes.** Data actually came from a different paper than cited; the shipped `.exp`
is a *model output* (synthetic) not measured data — legitimate only if the job is explicitly a
synthetic-recovery job; normalization applied to the data but undocumented; digitization off due to
log axes or a legend scale factor.

**Score effect.** FAIL (data don't match any primary source) **caps confidence at ~30** — the job
fits the wrong numbers. BLOCKED (figure un-digitizable, no table) caps mid.

---

## Gate 2 — Model fidelity

**Definition.** The `.bngl` is **identical** or **functionally equivalent** to the published model.

**Pass criterion.**
- *Identical* — same molecule types, seed species, observables, reaction rules, rate-law
  expressions, initial conditions, and units as the author's file (modulo cosmetic edits:
  whitespace, comments, `__FREE`→bind-by-id, stripped simulation actions, a retained
  `generate_network` cap).
- *Functionally equivalent (PASS (equiv))* — different syntax but the **same generated reaction
  network** (species set + reaction topology) and the same observable definitions; rate-constant
  *values* may differ (they are being fit) but rate-law *structure* and parameter *roles* match.

**Evidence.**
- If the author's `.bngl` is in the parking garage: `scripts/model_diff.py <our.bngl>
  <author.bngl>` — compares generated networks (species + reaction topology, tolerant of rate
  values) for network-generating models, or rules-level for network-free. Report species-only-in-A,
  reactions-only-in-A, etc.
- If only a paper/SI spec: a hand diff — a table of (species, rules, rate laws, ICs, units) paper
  vs. `.bngl`, each marked match / equivalent / differs. Confirm the model builds through BNG2.pl
  and that observables map to the paper's measured quantities.
- Note explicitly whether the **network cap** (`max_stoich`/`max_agg`/`max_iter`) and **units**
  match the paper; both are frequent silent-divergence points.

**Score effect.** A real structural difference (missing/extra rule, wrong rate law, wrong ICs,
wrong units) is a **FAIL** and caps confidence at ~35 until corrected. PASS (equiv) is full credit.

---

## Gate 3a — Reproduce the paper's result at the paper's parameters

**Definition.** With the model set to the **paper's reported best-fit** parameters, it reproduces
the paper's figure/data.

**Pass criterion.** Overlay of model-at-published-params vs. the (Gate-1) data meets a stated
quantitative bar (median/max relative error, peak amplitude/timing, or nearest-curve distance),
justified from the data precision.

**Evidence.** Reuse the job's `make_reproduction.py` **pointed at the paper's numbers** (not a
RuleHub/corpus re-fit, not a self-produced fit). Record the params used, their source (paper table
/ figure), the metric, and the figure. For a synthetic-recovery job, "published params" = the
ground truth, and 3a shows the ground-truth model regenerates the synthetic data (sanity-check the
t=0 values match the data's t=0 — a classic reproduction bug).

**Score effect.** Strong contributor. FAIL (model at published params does NOT match the paper's
own figure) is a serious red flag about Gate 1 or Gate 2 and should send you back.

---

## Gate 3b — Recover the paper's parameters by fitting

**Definition.** A PyBNF run recovers the paper's reported estimates within a stated tolerance.

**Pass criterion.** Fitted parameters land within tolerance of the published values (loose is fine
if justified — many rate-constant fits are sloppy/non-identifiable; then the *identifiable*
combinations or the reproduced curve is the real check). `scripts/compare_params.py
<published.json> <fit_sorted_params.txt>` tabulates published vs. recovered (log10 ratio,
within-tol flag) and summarizes.

**Evidence.** The pybnf run command + budget, the recovered params, the comparison table, and an
identifiability note (multiple basins? bounds hit?). If the job is **heavy** (cluster-scale
network / SSA), a partial or seeded run that shows the objective and the estimates *trending*
toward the published values is acceptable — say so and mark it PARTIAL, not PASS.

**Score effect.** PASS = strong. PARTIAL (trending, or only identifiable combos recovered, or heavy
→ not fully run) = partial credit. A job whose scope *cannot* reach the published params (a reduced
parameterization) FAILS 3b as "the paper's fit" — trigger the divergence policy.

---

## Divergence policy (restated)

Default: **make the job match what the paper reported.** Multiple published models/fits → multiple
slugs. Never substitute a simpler problem or a re-fit for the paper's own result without saying so.
If a reduced scope is kept (only when re-scoping is infeasible), relabel it explicitly as a reduced
variant and cap its confidence accordingly (it is not "the paper's fit").

---

## Confidence rubric (0–100, earned)

Start from the gate verdicts; the score must be reconstructable from them.

| Situation | Confidence band |
|---|---|
| Gates 0–2 PASS/equiv, 3a PASS, 3b PASS (params recovered within stated tol) | **90–100** |
| Gates 0–2 PASS/equiv, 3a PASS, 3b PARTIAL (heavy→trending, or only identifiable combos, or sloppy but curve reproduced) | **75–89** |
| Data + model grounded, 3a PASS, but 3b not run or scope-limited (documented) | **60–74** |
| A material gap: one source BLOCKED, or model only PASS (equiv) via prose (no author file), 3a ok | **45–59** |
| Data or model **FAIL/uncorrected**, or job scope diverges from the paper and is not re-scoped | **≤35** |
| Required sources missing (Gate 0 BLOCKED) — cannot validate | **≤25, explicitly "unassessed"** |

Adjust within a band from the strength of the evidence (author file vs. prose; digitization
tolerance; identifiability). Always write the one-line rationale naming the gates that set the band.
A score is only as good as the worst *grounding* gate (0/1/2): you cannot buy back bad data or a
wrong model with a pretty reproduction.
