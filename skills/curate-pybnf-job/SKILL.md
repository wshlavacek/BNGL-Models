---
name: curate-pybnf-job
description: Use when turning a published systems-biology / rule-based modeling paper (a PubMed Central article, PMCID/DOI, or a dev/papers/ folder) into a PyBNF edition-2 (new-era), PEtab.v2-compliant parameter-fitting job setup for lanl/pybnf's examples/real-world/ collection. Reconstructs the BNGL model from the paper, extracts the fitting data into .exp files (and qualitative properties into BPSL .prop/.con constraint files), authors the annotated .conf, verifies it end-to-end (tier-1 parse, PEtab export/lint/import round-trip for quantitative jobs or a job_type=check satisfaction pass for constraint-bearing ones, a bounded bngsim fit, and reproduction of the paper's reported fit), and registers it in _manifest.py + README + tests. Trigger whenever the user wants to add a real-world PyBNF fitting example from a paper, expand examples/real-world/, build a PEtab v2 parameterization problem from a model + data, fit qualitative/BPSL constraints, or "make a job setup" from a paper — even if they don't say "PEtab", "edition 2", or "constraints".
---

# Curate PyBNF Job

Use this skill to turn a published model + data into a **PyBNF edition-2 (new-era),
PEtab.v2-compliant parameter-fitting job setup**, and add it as a self-contained
example under `~/Code/PyBNF/examples/real-world/`. It is the fitting-job sibling of
`curate-model`: same "point at a paper, reconstruct faithfully, verify quantitatively"
discipline, but the deliverable is a runnable PyBNF fit rather than a static reference
simulation.

Prefer an explicit source in the request — a PMCID/URL, a DOI, or a local paper folder:

```text
Use curate-pybnf-job for PMC5334499  (or: for dev/papers/<Folder>)
```

## Inputs & where things go

- **Input:** a paper in PubMed Central (PMCID or URL), a DOI, or a local paper folder
  (e.g. `dev/papers/<Folder>` with a PDF + any author-provided model/data files).
- **Output:** job folders are **filed by source paper under a `<FirstAuthor>-<Year>/`
  directory** — one such directory per paper, holding one or more job **slugs**, because a
  single paper often yields several fitting problems (e.g.
  `Rukhlenko-2022/{cstar_trka, cstar_trkb, cstar_skmel133, cstar_skmel133_bpsl}/`). Unless
  the user names another location, create it under the PyBNF repo's real-world corpus:
  `~/Code/PyBNF/examples/real-world/<FirstAuthor>-<Year>/<slug>/`, plus edits to that repo's
  `_manifest.py`, `README.md`, and test suite. Do **not** commit or push unless asked; leave
  a clean, PR-ready set of changes.
- **Paper-level landing README:** when the grouping directory holds more than one slug, add
  a `<FirstAuthor>-<Year>/README.md` that ties them together — the full paper citation, the
  shared model (what was reconstructed, from which supplementary file), a one-row-per-slug
  table (what each fits · flavor · data source · verification status), and the source-file
  links. It complements the corpus-level coverage matrix, and *is* the entry point when the
  output is a self-contained per-folder set (the pattern below).

## Required reading

Before authoring anything, read (they are the source of truth, not your memory):

1. `references/edition2-conf-reference.md` — the full edition-2 `.conf` surface
   (`job_type`, `objective`, `experiment:`/`condition:`, `*_var` free params). What is
   **[E2]** vs rejected **[LEGACY]**.
2. `references/real-world-anatomy.md` — exactly what a real-world folder,
   `_manifest.py` entry, and the two test tiers require.
3. `references/petab-compliance.md` — what "PEtab.v2-compliant" means (a *verified
   round-trip property* of the native conf) and the **PEtab-exportable subset** you
   must stay inside.
4. `references/bpsl-constraints.md` — BPSL `.prop`/`.con` constraint files: the
   grammar, how to attach them (data fusion or constraint-only), `job_type = check`,
   and the fact that a constraint-bearing job is **native-only (not PEtab-exportable)**.
   Read this whenever the paper reports *qualitative* facts (orderings, thresholds,
   monotonic dose-response, "peaks before", bistability) — BPSL is PyBNF's signature
   capability and a strong reason to add an example.
5. `skills/bngl/skill.md` and `skills/bngl/templates/model_skeleton.bngl` — house style
   for the BNGL model you reconstruct. Read `skills/nfsim/SKILL.md` if the model is
   network-free (crosslinking, aggregation, cyclic complexes).
6. Two or three existing examples for the idioms to imitate:
   `examples/real-world/receptor/` (simplest: ODE + pre-equilibration + `sos`),
   `examples/real-world/egfr_ode/` (time course + dose-response scan + `chi_sq`), and,
   for a stochastic paper, `examples/real-world/tlbr/` (NFsim scan).

## Environment

The verification scripts and `pybnf` run in the **PyBNF** environment. Use its
interpreter and set `BNGPATH` (BNG2.pl is needed even to parse a model). The scripts
`chdir` into the conf's folder themselves, so they may be invoked from anywhere with
absolute paths:

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder containing BNG2.pl
PY=~/Code/PyBNF/.venv/bin/python
SKILL=~/Code/BNGL-Models/skills/curate-pybnf-job     # this skill's directory
CONF=~/Code/PyBNF/examples/real-world/<FirstAuthor>-<Year>/<slug>/<slug>.conf
```

If BNGPATH/BNG2.pl or the PyBNF venv can't be found, ask the user rather than guessing.

## Workflow

1. **Read the paper and extract the fitting problem.** Identify: the biochemical model
   (species, rules, parameters, initial conditions); the **data to fit** and which
   figure/table it lives in; the **experimental design** (time course? dose-response?
   pre-equilibration/washout?); the **observables** (what each measured quantity maps to
   in the model); the **free parameters** and their published values / plausible search
   ranges; whether the dynamics are deterministic (ODE), stochastic (SSA), or
   network-free (NFsim); and any **qualitative properties** the paper asserts that a
   time-series doesn't capture (orderings, thresholds, steady-state levels, monotonic
   dose-response, "peaks before", bistability, oscillation) — these become BPSL `.prop`
   constraints. Write this down before building.

   Then **decide the example's flavor**, because it sets the verification path
   (`references/bpsl-constraints.md`):
   - **quantitative** (`.exp` only) → PEtab-exportable → the full PEtab round-trip applies;
   - **data fusion** (`.exp` + `.prop`) or **constraint-only** (`.prop` only) → uses BPSL →
     **native-only, not PEtab-exportable** → verify with `job_type = check` instead.
   Prefer capturing the paper's qualitative claims as constraints when they carry real
   information — a BPSL example is often a *more* compelling addition than a plain fit.

2. **Name and create the folder.** File the job under its source paper: a
   `<FirstAuthor>-<Year>/` directory (e.g. `Rukhlenko-2022/`) holding a short, descriptive
   slug like the corpus (`receptor`, `egfr_ode`, `tlbr`) — e.g. `Rukhlenko-2022/cstar_trka/`.
   One author-year directory per paper; for each additional fitting problem the same paper
   yields (a second cell line, a dose-response, a BPSL landscape job, …) add a new slug
   *beside* the others in the same directory. Append `b`/`c` to the year (`Smith-2020b`) only
   to disambiguate two papers sharing a first author and year. Create the folder at the
   output location from "Inputs & where things go". Once the directory holds more than one
   slug, add (or update) the paper-level `<FirstAuthor>-<Year>/README.md` landing page
   described in "Inputs & where things go".

3. **Reconstruct the BNGL model from scratch**, fitting-ready and on the edition-2
   surface (differences from a curate-model reference model):
   - **No `begin actions` block** — the simulation is synthesized from the conf.
   - Fitted rate constants are bare `id nominal` declarations; the free-parameter names
     in the conf bind to these ids **by name** (no `__FREE` alias, ADR-0034).
   - **Observable/function names are the contract with the `.exp` header *and* any
     `.prop` constraint** — a `Molecules`/`Species` observable is a plain name; a
     `functions` entry appears **with parentheses** in the exp header. Every name a
     `.prop` references must exist as a model observable.
   - For pre-equilibration, add a boolean gate parameter (e.g. `Ligand_isPresent 0`)
     toggled by two `condition:` states.
   Follow `skills/bngl/skill.md` house style throughout. Confirm the model builds
   through BNG2.pl and behaves as the paper describes (an independent check of the
   network/dynamics, in the spirit of curate-model's model-specification verification).

4. **Extract the data into `.exp` file(s).** First column is the independent variable
   (`time` for a time course; the **model-parameter name** for a dose/scan column);
   remaining columns are observables named exactly as in the model. Add `<obs>_SD`
   columns only if the paper reports per-point uncertainty (they switch the objective to
   `chi_sq`); use `NaN` for missing points. Prefer source tabular/supplementary data;
   if only plotted curves exist, **digitize** the panel (record figure/panel, extraction
   method, axis calibration, any legend scale factors) exactly as `curate-model`
   prescribes — the digitized fit-target data is what "reproducing the paper" is measured
   against, so keep it faithful and documented.

   **If the example uses BPSL** (data-fusion or constraint-only), also author the
   `.prop` file(s): one qualitative statement per line in the BPSL grammar
   (`<obs> <op> <obs|const> <always|once|at …|between …> weight <w>`), using dotted
   `suffix.obs` to compare across experiments/mutants. Translate each qualitative claim
   in the paper into one line; weight stronger claims higher. Full grammar and real
   examples: `references/bpsl-constraints.md`.

5. **Author `<name>.conf`** from `templates/job_setup.conf`, on the edition-2 surface:
   `edition = 2`, `bngl_backend = bngsim`, `model:`, the `condition:`/`experiment:`
   lines that encode the design, a `job_type` (a global metaheuristic — `de`/`ss`/`pso`
   — for a first paper fit), the `objective` (`sos` without `_SD`, `chi_sq` with), a
   search budget, and one `*_var` free parameter per fitted id, bracketing the published
   value. Open the conf with a banner comment citing the paper (PMCID/DOI) and the figure
   each `.exp` came from.
   - **Quantitative example:** bind only `.exp` (`experiment: <name>, data: <name>.exp`)
     and **stay inside the PEtab-exportable subset** — no `normalization`, `cumulative`,
     `neg_bin`/`lognormal` noise (`references/petab-compliance.md`).
   - **BPSL example:** attach the constraints via the experiment's `data:` — data fusion
     (`experiment: <name>, data: <name>.exp, <name>.prop`) or constraint-only
     (`experiment: <name>, t_end: <T>, n_steps: <N>, data: <name>.prop`, where `t_end` is
     required). Optionally set `constraint_scale` to balance the qualitative penalty
     against the quantitative objective. This job is native-only — that's expected, not a
     defect (`references/bpsl-constraints.md`).

6. **Verify end-to-end** (the completion bar), with BNGPATH set. The `scripts/*`
   (`check_conf.py`, `petab_roundtrip.py`) `chdir` into the conf's folder themselves, so run
   them from anywhere with an absolute `$CONF`. Run `pybnf` itself **from inside the job
   folder** (`cd <FirstAuthor>-<Year>/<slug>`), since its relative paths — `model:`, `data:`,
   `output_dir` — resolve against the working directory:
   - **Tier-1 (parses & well-formed):** `$PY $SKILL/scripts/check_conf.py $CONF`
     — edition 2, `job_type` resolves, data and/or constraints bound, free params bind by
     id. It reports whether the job carries BPSL constraints (i.e. is native-only).
   - **PEtab.v2 compliance — *quantitative examples only* (the round-trip):**
     `$PY $SKILL/scripts/petab_roundtrip.py $CONF --job-type <jt>`
     — export → `petab.v2` lint clean → import. If it reports a *non-exportable feature*
     (e.g. `normalization`), rework the conf (step 5). If it reports *BPSL constraint
     data*, that is expected for a BPSL example — skip this check and use the `check`
     verification below instead.
   - **BPSL satisfaction — *constraint-bearing examples only* (`job_type = check`):** run
     a `check` job at the fitted (or published) parameters and confirm it prints
     `Satisfied M out of M constraints` (or document which the paper itself does not
     require). This is the native-only analog of the PEtab round-trip.
   - **A real bngsim fit reaches a finite objective:** run a short bounded fit
     (`pybnf -c <name>.conf` with a small `max_iterations`/`population_size`, or the
     recovery-tier path in `test_real_world_examples.py`) and confirm the
     simulate→score→propose loop yields a finite score. If the network is cluster-scale
     (minutes to generate, or NFsim on ~10³ molecules), don't force a full fit — mark the
     example `heavy=True` and document it.
   - **Reproduce the paper's reported fit:** compare the fitted model against the
     digitized/extracted target data quantitatively (curate-model metrics: max/median
     relative error, peak amplitude/timing, or nearest-curve distance), and — when the
     paper reports point estimates — check the recovered parameters land in the right
     ballpark. State the tolerance and justify it from the data precision. Keep this
     verification (a short notebook/script + a PNG) in your scratchpad workspace and
     summarize the result; it is evidence for the PR, not a committed real-world file.

7. **Register the example (PR-ready).**
   - Add the `_manifest.py` entry (`templates/manifest_entry.py`): `folder`, `conf`,
     `simulator`, `observables` (the `.exp` columns; functions without parens),
     `system` (biology + paper mapping; for a BPSL example, note the constraints), and
     `stochastic` (True iff ssa/nf), `heavy` if cluster-scale, optional
     `recover={param: published_value}` + `tol`.
   - Update `README.md`: add the coverage-matrix row (status ✅ if it runs in the
     recovery tier, 🔶 if too heavy for routine CI); add a Known-limitations bullet if
     heavy. For a BPSL example, note in the row that it is native-only (not
     PEtab-exportable).
   - Add a **test assertion** for the new example (the gap the research flagged):
     - *Quantitative example:* a PEtab export/lint check — extend
       `tests/test_real_world_examples.py` (or a sibling) with a parametrized test that
       runs `pybnf.petab.export_job` + `petab.v2` `lint_problem`, mirroring
       `scripts/petab_roundtrip.py`. Gate it on the `petab` extra being present.
     - *BPSL example:* **do not** assert PEtab (it will correctly refuse). Assert instead
       that `export_job` raises `NotImplementedError` (a guard that the constraint stays
       native-only), and/or that a `check` run reports the expected satisfaction.

8. **Report.** Summarize the new folder, the verification results (tier-1; PEtab
   round-trip *or* BPSL `check` satisfaction; fit objective; paper-reproduction metric),
   and the exact `_manifest.py` / `README.md` / test edits — everything needed to open
   the PR to lanl/pybnf.

## Two guardrails that bite

- **PEtab-exportable subset (quantitative examples).** `export_job` fails loud on
  features PEtab v2 can't express. The frequent trap: `normalization = init/peak/…` is
  **not** exportable (the shipped `igf1r.conf` uses it and is non-compliant). If the
  paper's data is relative-to-a-reference, encode the reference in the model/observable
  instead. Full list and substitutions: `references/petab-compliance.md`.
- **BPSL constraints are native-only.** Any `.prop`/`.con` on an experiment's `data:`
  makes the job non-exportable — `export_job` raises `NotImplementedError`, and there is
  no round-trip. This is a *property of the flavor you chose*, not a bug: verify a BPSL
  example with `job_type = check`, not `petab_roundtrip.py`, and register it as
  native-only. Full grammar and rationale: `references/bpsl-constraints.md`.

## Deliverables

```text
~/Code/PyBNF/examples/real-world/<FirstAuthor>-<Year>/       # e.g. Rukhlenko-2022/
├── README.md            # paper-level landing page (when >1 slug): citation, shared model,
│                        #   one-row-per-slug table (fits · flavor · data · status), sources
└── <slug>/                                                  # e.g. cstar_trka/
    ├── <slug>.bngl      # edition-2, fitting-ready, no actions block
    ├── <slug>.conf      # edition-2 job setup, banner-commented
    ├── <data>.exp       # column headers == model observable names (quantitative / fusion)
    └── <slug>.prop      # optional: BPSL constraints (data-fusion or constraint-only example)
# plus, in the PyBNF repo:
#   _manifest.py         # one RealWorldExample(...) entry
#   README.md            # coverage-matrix row (+ limitations if heavy; note native-only if BPSL)
#   tests/…              # a PEtab export/lint assertion (quantitative) OR an export-refused /
#                        # check-satisfaction assertion (BPSL)
# plus, in your scratchpad (evidence, not committed to PyBNF):
#   verification notebook/script + PNG comparing the fit to the paper
```

## Completion criteria

Not complete until:
- `scripts/check_conf.py` passes (edition 2, job_type resolves, data and/or constraints
  bound, free params bind by id with no `__FREE`);
- **either** `scripts/petab_roundtrip.py` passes (quantitative example: export → lint
  clean → import) **or** a `job_type = check` run reports the intended constraint
  satisfaction (BPSL example) and the example is registered native-only;
- a real bngsim fit reaches a finite objective, **or** the example is justifiably marked
  `heavy=True` with the reason documented;
- the fit is compared quantitatively to the paper's reported data/parameters, with the
  metric and tolerance stated (or it is explicitly documented why the reported data
  could not be extracted);
- the `_manifest.py` entry, `README.md` coverage row, and the matching test assertion
  (PEtab-lint for quantitative, export-refused/`check` for BPSL) are added and consistent
  with the folder;
- the conf header and manifest `system` field cite the source paper (PMCID/DOI) and the
  figure/table the data (and any qualitative properties) came from.

If a required artifact cannot be produced, state exactly which one and why.
