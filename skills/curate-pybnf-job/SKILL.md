---
name: curate-pybnf-job
description: Use when turning a published systems-biology / rule-based modeling paper (a PubMed Central article, PMCID/DOI, or a dev/papers/ folder) into a PyBNF edition-2 (new-era), PEtab.v2-compliant parameter-fitting job setup for lanl/pybnf's examples/real-world/ collection. Reconstructs the BNGL model from the paper, extracts the fitting data into .exp files, authors the annotated .conf, verifies it end-to-end (tier-1 parse, PEtab export/lint/import round-trip, a bounded bngsim fit, and reproduction of the paper's reported fit), and registers it in _manifest.py + README + tests. Trigger whenever the user wants to add a real-world PyBNF fitting example from a paper, expand examples/real-world/, build a PEtab v2 parameterization problem from a model + data, or "make a job setup" from a paper — even if they don't say "PEtab" or "edition 2".
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
- **Output:** a new folder `~/Code/PyBNF/examples/real-world/<name>/` (in the *PyBNF*
  repo, a separate working directory), plus edits to that repo's `_manifest.py`,
  `README.md`, and test suite. Do **not** commit or push unless asked; leave a clean,
  PR-ready set of changes.

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
4. `skills/bngl/skill.md` and `skills/bngl/templates/model_skeleton.bngl` — house style
   for the BNGL model you reconstruct. Read `skills/nfsim/SKILL.md` if the model is
   network-free (crosslinking, aggregation, cyclic complexes).
5. Two or three existing examples for the idioms to imitate:
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
CONF=~/Code/PyBNF/examples/real-world/<name>/<name>.conf
```

If BNGPATH/BNG2.pl or the PyBNF venv can't be found, ask the user rather than guessing.

## Workflow

1. **Read the paper and extract the fitting problem.** Identify: the biochemical model
   (species, rules, parameters, initial conditions); the **data to fit** and which
   figure/table it lives in; the **experimental design** (time course? dose-response?
   pre-equilibration/washout?); the **observables** (what each measured quantity maps to
   in the model); the **free parameters** and their published values / plausible search
   ranges; and whether the dynamics are deterministic (ODE), stochastic (SSA), or
   network-free (NFsim). Write this down before building.

2. **Name and create the folder.** Choose a short, descriptive slug like the corpus
   (`receptor`, `egfr_ode`, `tlbr`); disambiguate with `_<firstauthor><year>` only if
   needed. Create `~/Code/PyBNF/examples/real-world/<name>/`.

3. **Reconstruct the BNGL model from scratch**, fitting-ready and on the edition-2
   surface (differences from a curate-model reference model):
   - **No `begin actions` block** — the simulation is synthesized from the conf.
   - Fitted rate constants are bare `id nominal` declarations; the free-parameter names
     in the conf bind to these ids **by name** (no `__FREE` alias, ADR-0034).
   - **Observable/function names are the contract with the `.exp` header** — a
     `Molecules`/`Species` observable is a plain name; a `functions` entry appears
     **with parentheses** in the exp header.
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

5. **Author `<name>.conf`** from `templates/job_setup.conf`, on the edition-2 surface:
   `edition = 2`, `bngl_backend = bngsim`, `model:`, the `condition:`/`experiment:`
   lines that encode the design, a `job_type` (a global metaheuristic — `de`/`ss`/`pso`
   — for a first paper fit), the `objective` (`sos` without `_SD`, `chi_sq` with), a
   search budget, and one `*_var` free parameter per fitted id, bracketing the published
   value. **Stay inside the PEtab-exportable subset** (no `normalization`, `cumulative`,
   `neg_bin`/`lognormal` noise, or `.prop` constraints — see `petab-compliance.md`).
   Open the conf with a banner comment citing the paper (PMCID/DOI) and the figure each
   `.exp` came from.

6. **Verify end-to-end** (the completion bar). Run from the PyBNF repo root with BNGPATH
   set:
   - **Tier-1 (parses & well-formed):** `$PY $SKILL/scripts/check_conf.py $CONF`
     — edition 2, `job_type` resolves, data bound, free params bind by id.
   - **PEtab.v2 compliance (the round-trip):**
     `$PY $SKILL/scripts/petab_roundtrip.py $CONF --job-type <jt>`
     — export → `petab.v2` lint clean → import. If it reports a non-exportable feature,
     rework the conf (step 5) rather than committing a non-compliant example.
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
     `system` (biology + paper mapping), `stochastic` (True iff ssa/nf), `heavy` if
     cluster-scale, and optionally `recover={param: published_value}` + `tol`.
   - Update `README.md`: add the coverage-matrix row (status ✅ if it runs in the
     recovery tier, 🔶 if too heavy for routine CI); add a Known-limitations bullet if
     heavy.
   - Add a **PEtab export/lint test assertion** for the new example (the gap the
     research flagged — the manifest doesn't yet assert PEtab). Extend
     `tests/test_real_world_examples.py` (or a sibling test) with a parametrized check
     that runs `pybnf.petab.export_job` + `petab.v2` `lint_problem` for the example,
     mirroring `scripts/petab_roundtrip.py`. Gate it on the `petab` extra being present.

8. **Report.** Summarize the new folder, the verification results (tier-1, PEtab
   round-trip, fit objective, paper-reproduction metric), and the exact `_manifest.py` /
   `README.md` / test edits — everything needed to open the PR to lanl/pybnf.

## PEtab-exportable subset (the one guardrail that bites)

`export_job` fails loud on features PEtab v2 can't express. The frequent trap:
`normalization = init/peak/…` is **not** exportable (the shipped `igf1r.conf` uses it
and is non-compliant). If the paper's data is relative-to-a-reference, encode the
reference in the model/observable instead. Full list and substitutions:
`references/petab-compliance.md`.

## Deliverables

```text
~/Code/PyBNF/examples/real-world/<name>/
├── <name>.bngl          # edition-2, fitting-ready, no actions block
├── <name>.conf          # edition-2 job setup, PEtab-exportable, banner-commented
└── <data>.exp           # one or more; column headers == model observable names
# plus, in the PyBNF repo:
#   _manifest.py         # one RealWorldExample(...) entry
#   README.md            # coverage-matrix row (+ limitations if heavy)
#   tests/…              # a PEtab export/lint assertion for the example
# plus, in your scratchpad (evidence, not committed to PyBNF):
#   verification notebook/script + PNG comparing the fit to the paper
```

## Completion criteria

Not complete until:
- `scripts/check_conf.py` passes (edition 2, job_type resolves, data bound, free params
  bind by id, no `__FREE`);
- `scripts/petab_roundtrip.py` passes (export → lint clean → import) — the example is
  inside the PEtab-exportable subset;
- a real bngsim fit reaches a finite objective, **or** the example is justifiably marked
  `heavy=True` with the reason documented;
- the fit is compared quantitatively to the paper's reported data/parameters, with the
  metric and tolerance stated (or it is explicitly documented why the reported data
  could not be extracted);
- the `_manifest.py` entry, `README.md` coverage row, and the PEtab test assertion are
  added and consistent with the folder;
- the conf header and manifest `system` field cite the source paper (PMCID/DOI) and the
  figure/table the data came from.

If a required artifact cannot be produced, state exactly which one and why.
