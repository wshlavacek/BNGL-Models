---
name: validate-pybnf-job
description: Use when confirming that an existing PyBNF parameter-fitting job setup (one produced by curate-pybnf-job, e.g. a folder under pybnf-jobs/<FirstAuthor>-<Year>/<slug>/) is genuinely FAITHFUL TO ITS PRIMARY SOURCES — that its data trace to the paper's tables/figures, its BNGL model is identical or functionally equivalent to the published model, and a PyBNF run reproduces the estimates the paper reports. Grounds each job in the original PDF + supplemental materials (staged untracked in dev/papers/<FirstAuthor>-<Year>/), extracts/digitizes reported data when needed, diffs the model against the authors' own files, runs the fit, and emits a committed VALIDATION.md scorecard with an EARNED confidence score and a list of corrections. Trigger whenever the user wants to validate / audit / ground / fact-check a pybnf-jobs entry against the literature, verify a job matches the paper, re-scope a job that drifted from what was published, or assign a defensible confidence to a job — even if they don't say "validate" or "primary source".
---

# Validate PyBNF Job

The auditor sibling of `curate-pybnf-job`. That skill *builds* a job by mirroring/reconstructing;
this one *proves the job is faithful to the primary literature* and assigns an **earned** confidence
(not a guessed one). The deliverable is a committed **`VALIDATION.md`** scorecard in the job folder
plus any corrections the audit forces — including **re-scoping a job that was simplified away from
what the paper reported** (the default is always: match the literature).

Motivation: the `pybnf-jobs/` corpus was first built by mirroring PyBNF's `examples/real-world/` +
RuleHub, not from primary reconstruction. Running/exporting cleanly proves a job is *well-formed*;
it does **not** prove the data, model, and results match the paper. This skill closes that gap.

## Hard prerequisite: the primary sources must be in hand

You cannot validate against sources you do not have. **Before doing anything else, confirm the
parking garage holds the materials** for this job:

```
dev/papers/<FirstAuthor>-<Year>/     # git-IGNORED (whole dev/ tree); never commit or redistribute
   <paper>.pdf                        # the primary paper(s) — model paper AND data paper if different
   <supplement>.pdf / .xlsx / .txt    # SI, Source Data, author BNGL / BioNetFit / PyBNF job files
   NEEDED.md                          # the collection checklist for this job
```

If the required PDFs / supplements are **missing**, STOP: list exactly what is needed (use the
dir's `NEEDED.md`), ask the user to fetch them into the parking garage, and do not proceed by
guessing or by re-mirroring the corpus. A validation built on absent sources is worthless.
Prefer the **authors' own model/job files** when the paper ships them (SI archive, GitHub) — they
turn Gate 2 from a judgement call into a direct diff.

## Required reading

Before auditing, read (source of truth, not memory):
1. `references/validation-protocol.md` — the gates, their pass criteria, the divergence policy,
   and the confidence rubric. This is the spine of the skill.
2. `references/data-provenance.md` — how to trace an `.exp` to a primary table, and how to
   **re-digitize** a figure from the PDF faithfully (panel, axis calibration, scale factors).
3. `../curate-pybnf-job/references/` — `real-world-anatomy.md`, `edition2-conf-reference.md`,
   `petab-compliance.md`, `bpsl-constraints.md` (what the job's files are supposed to be).
4. `../curate-model/SKILL.md` §"Reported Simulation Data Verification" — the digitization +
   quantitative-agreement discipline this skill reuses.
5. The job's own `README.md` and the parking-garage `NEEDED.md`.

## The five gates (run in order; record every one in VALIDATION.md)

Full criteria in `references/validation-protocol.md`. In brief:

- **Gate 0 — Materials inventory.** Enumerate what's in the parking garage vs. what the job needs.
  Identify: the model paper, the data paper (may differ), the SI, and any author model/job files.
  A gate you can't ground because a source is missing is `BLOCKED`, not `PASS`.
- **Gate 1 — Data provenance.** Every `.exp` value traces to a primary **table** (transcribe and
  diff) or a **figure** (re-digitize from the PDF, `scripts`/`curate-model` method, and diff
  against the shipped `.exp`). State the extraction method and tolerance. `scripts/compare_data.py`
  diffs a re-digitized CSV against the shipped `.exp`.
- **Gate 2 — Model fidelity.** The `.bngl` is **identical or functionally equivalent** to the
  published model (species, molecule types, rules, rate laws, initial conditions, UNITS,
  conserved observables). Diff against the author's file if available; else against the paper/SI
  spec. `scripts/model_diff.py` compares two models by their generated networks (species +
  reaction topology, tolerant of the fitted rate *values*). Record **every** deviation and whether
  it is cosmetic, functionally-equivalent, or a real difference.
- **Gate 3a — Reproduce the paper's result at the paper's parameters.** Set the model to the
  **published best-fit** values and confirm it reproduces the paper's figure/data (overlay +
  quantitative metric). Reuse the job's `make_reproduction.py`, but pointed at the **paper's**
  numbers — never a corpus re-fit or a self-produced fit standing in for "the paper".
- **Gate 3b — Recover the paper's parameters by fitting.** Run PyBNF and confirm the fit lands the
  paper's reported estimates within a stated tolerance. `scripts/compare_params.py` tabulates
  published vs. recovered (log-ratio, within-tol). If the job is heavy, a partial/seeded run that
  shows the estimates *trending* to the reported values is acceptable — say so.

**Gate 3a and 3b routinely diverge, and the divergence is the finding.** Define "the paper's
result" per job *before* running (a specific figure, a specific parameter table). For synthetic-data
jobs (e.g. Gupta fceri) "the paper's result" is **ground-truth recovery**, which Gate 3b must
actually demonstrate, not assert.

## Divergence policy — match the literature

If the audit finds the job was **scoped differently from the paper** — a simplified
parameterization, a subset of the reported fits, a re-fit substituted for the paper's own values —
the default is to **correct the job to match what was published**, not to keep the simplification.

- A paper with **multiple models / multiple fits** yields **multiple slugs**; validate (and, if
  missing, build) each one that is tractable. Do not collapse them into one easier problem.
- Worked example: `Erickson-2019/igf1r` was a reduced 3-parameter (`K1`/`K2`/`K1prime`), F5B-only
  distillation; the paper (its SI ships the authors' own BioNetFit files) reports a fuller
  7-rate-constant fit to three datasets. It was **re-scoped** to the published model/fit — the
  authors' model verbatim, three datasets, detailed balance — and required **legacy (edition-1)**
  mode because edition-2 `experiment:` directives cannot express the F5D preincubate/wash/reset
  protocol. (If re-scoping is genuinely infeasible, the fallback is to keep the reduced form but
  **relabel it explicitly** as a reduced variant, not "the paper's fit", and score Gate 3b down.)

When a correction is required, apply it to the job files (model/conf/exp/README + reproduction),
re-run the affected gates, and note the change in both `VALIDATION.md` and the job `README.md`.

## Confidence score (earned, not guessed)

`VALIDATION.md` ends with a 0–100 confidence and a one-line rationale, computed from the gate
outcomes per the rubric in `references/validation-protocol.md` (each gate contributes; a `BLOCKED`
or `FAIL` on data/model caps the score hard). The number must be *defensible from the recorded gate
evidence* — if a reader disagrees, they should be pointing at a specific gate, not a vibe.

## Environment

Same as `curate-pybnf-job` (its scripts + BNG2.pl); this skill's scripts also use only `numpy` +
BNG2.pl (via `BNGPATH`) unless noted:

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"     # folder containing BNG2.pl
PY=~/Code/PyBNF/.venv/bin/python
VSK=~/Code/BNGL-Models/skills/validate-pybnf-job
JOB=~/Code/BNGL-Models/pybnf-jobs/<FirstAuthor>-<Year>/<slug>
PAPERS=~/Code/BNGL-Models/dev/papers/<FirstAuthor>-<Year>   # parking garage (untracked)

$PY $VSK/scripts/compare_data.py   <redigitized.csv> <shipped.exp>        # Gate 1
$PY $VSK/scripts/model_diff.py     <our.bngl> <author_reference.bngl>     # Gate 2
$PY $VSK/scripts/compare_params.py <published.json> <fit_sorted_params.txt>  # Gate 3b
```

Run `pybnf` and the job's `make_reproduction.py` from inside the slug folder (relative paths).

## Workflow

1. **Gate 0.** Inventory `dev/papers/<Author-Year>/` against `NEEDED.md`. Missing sources → report
   and stop.
2. **Read the primary sources.** Read the paper (and SI) for: the model spec, the exact
   data/figure the job fits, the reported parameters, and any qualitative claims (for BPSL jobs).
3. **Gate 1 (data).** Trace/re-digitize; diff vs. the shipped `.exp`.
4. **Gate 2 (model).** Diff the `.bngl` vs. the author's file or the paper spec.
5. **Decide "the paper's result"** for this job (figure + parameter table), and whether the job's
   scope matches it. Apply the divergence policy.
6. **Gate 3a / 3b.** Reproduce at published params; fit and compare recovered vs. published.
7. **Apply corrections** the audit forces; re-run affected gates.
8. **Write `VALIDATION.md`** (from `templates/VALIDATION.md`): per-gate PASS/FAIL/BLOCKED with
   evidence, the model diff, the data-provenance table, 3a/3b metrics, the earned confidence + one
   corrections list. Update the job `README.md` if anything changed. Commit only derived artifacts
   (never the PDFs).
9. **Report.** Summarize each gate, the confidence, and what (if anything) was corrected.

## Completion criteria

Not complete until, for the job (or each slug):
- Gate 0 materials are present (or the job is explicitly `BLOCKED` with the missing list);
- Gate 1 data provenance is established with a quantitative `.exp`-vs-source diff;
- Gate 2 records the model as identical / functionally-equivalent / different, with every deviation
  named;
- "the paper's result" is defined, and Gate 3a (reproduce at published params) + Gate 3b (recover
  by fitting) are run, or their infeasibility is documented with the reason;
- the divergence policy has been applied (job re-scoped to the literature, or the reduced scope
  explicitly relabeled);
- `VALIDATION.md` exists with per-gate evidence and an earned confidence score, and any forced
  corrections are applied and reflected in the job `README.md`.

If a gate cannot be closed, state exactly which one, why, and what source or run would close it.
