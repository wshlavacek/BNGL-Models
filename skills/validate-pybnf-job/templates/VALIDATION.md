# VALIDATION — <FirstAuthor>-<Year>/<slug>

Primary-source validation of the PyBNF job `pybnf-jobs/<FirstAuthor>-<Year>/<slug>/`, per the
`validate-pybnf-job` skill. Confidence is **earned from the gate evidence below**.

> **Confidence: NN / 100** — <one-line rationale naming the gates that set the band>.

Primary sources (in the untracked `dev/papers/<FirstAuthor>-<Year>/`; not redistributed):
- Model paper: <cite + PMCID/DOI>
- Data paper (if different): <cite + PMCID/DOI>
- Supplement / author files used: <SI, Source Data, author `.bngl`/`.conf`/`.exp` — or "none, prose spec only">

"The paper's result" for this job = **<Figure/Table N>** + parameter table **<Table N>**.

---

## Gate 0 — Materials inventory

| needed | present? | path / note |
|---|---|---|
| model paper PDF | ✅/❌ | |
| data paper PDF | ✅/❌ | |
| SI / Source Data | ✅/❌ | |
| author model/job files | ✅/❌ | (enables a direct Gate-2 diff) |

**Verdict:** PASS / BLOCKED — <what's missing>.

## Gate 1 — Data provenance

| `.exp` | source (table/fig, series) | method | normalization/units | diff vs. shipped | verdict |
|---|---|---|---|---|---|
| `<name>.exp` | <Fig 5B, WT> | transcribed / digitized | <none / to-t0 / %max> | median X%, max Y% | PASS/… |

**Verdict:** … (`scripts/compare_data.py` output summarized).

## Gate 2 — Model fidelity

Reference compared against: <author `.bngl` in garage / paper SI spec>.

| aspect | paper | our `.bngl` | verdict |
|---|---|---|---|
| molecule types / seed species | | | match / equiv / differs |
| reaction rules + rate laws | | | |
| initial conditions / totals | | | |
| units | | | |
| network cap (max_stoich/agg/iter) | | | |
| observables ↔ measured quantities | | | |

`scripts/model_diff.py` result: <species-only-in-X, reactions-only-in-X, …>.

**Verdict:** PASS (identical) / PASS (equiv) / FAIL — <every deviation named>.

## Gate 3a — Reproduce at the paper's parameters

- Published params used: <values + source: Table N / Fig caption>.
- Reproduction: `make_reproduction.py` at those params → `<slug>_reproduction.png`.
- Metric: <median/max rel err; peak amp/timing; nearest-curve distance>.

**Verdict:** PASS/… — model at the paper's own params <does/does not> reproduce <Figure N>.

## Gate 3b — Recover the paper's parameters by fitting

- Run: `pybnf -c <slug>.conf` (<budget>; heavy? partial?).
- `scripts/compare_params.py <published.json> <fit_sorted_params.txt>`:

| param | published | recovered | log10 ratio | within tol? |
|---|---|---|---|---|
| | | | | |

- Identifiability: <sloppy? multiple basins? bounds hit?>.

**Verdict:** PASS / PARTIAL (trending / heavy / identifiable-combos-only) / FAIL (scope can't reach).

---

## Divergence & corrections

- Scope vs. paper: <matches / job was reduced — action taken>.
- Corrections applied to job files: <list, or "none">.
- Re-run after corrections: <gates re-run>.

## Bottom line

<2–3 sentences: what is solid, what is the residual risk, and the single most valuable next step to
raise the confidence.>
