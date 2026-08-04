---
name: curate-model
description: Use when converting a source paper folder such as dev/papers/Faeder2003 into a complete curated BNGL model collection under models/, including house-style BNGL files, reference simulation data, independent verification against the model specification and reported simulation data, metadata.yaml, and README Models-table update.
---

# Curate Model

Use this skill to add a published model, or a collection of closely related model
variants, to the BNGL-Models library from a source folder under `dev/papers/`.

Prefer an explicit source path in the user request:

```text
Use curate-model for dev/papers/Faeder2003.
```

If the user supplies only a folder name such as `Faeder2003`, infer
`dev/papers/Faeder2003`.

## Required Reading

Before creating or editing model artifacts:

1. Read `skills/bngl/skill.md`. It is authoritative for **content** — `.bngl`
   structure and formatting, naming, units, annotation, the `metadata.yaml`
   schema, the overlay-plot style (§2), and the lint rules (§9). **This skill is
   authoritative for the workflow**: which artifacts exist, the two verification
   levels, the verification-artifact shape, and the completion criteria. Where
   the two ever disagree on process, this skill wins.
2. Read `skills/bngl/skill.md` §5 for the `metadata.yaml` schema, the source-tag
   vocabulary, and the file-role vocabulary, and §9 for the lint rules a finished
   folder must pass.
3. Read `skills/bngl/templates/model_skeleton.bngl` before writing a new BNGL
   file.
4. Read `skills/nfsim/SKILL.md` when NFsim behavior, crosslinking, molecularity,
   species observables, cyclic complexes, or network-free simulation may affect
   correctness.
5. Read `references/stochastic-verification.md` when the model's defining result is
   noise — a stationary distribution, switching statistics, a spectrum, bursts — or
   whenever the simulation protocol is SSA or NFsim. A stochastic model is not
   verified by comparing trajectories, and that reference fixes what to compare
   instead, the independence rule, and the tolerances in current use.
6. Read `references/when-the-paper-is-wrong.md` as soon as anything fails to
   reproduce. A published result that will not reproduce is a finding, not a
   blocker; that reference gives the triage, the five permitted responses, and the
   house wording. Read it *before* changing any published value.
7. Review existing `models/` folders for naming, metadata, verification notebook,
   PNG, reference-data, and README conventions.
8. Read all relevant files in the requested `dev/papers/<folder>/`, including the
   PDF and any existing source BNGL or model files.

## Workflow

1. Identify the source paper, existing source model files, model scope, and the
   published figure or figures to reproduce. If no figure is specified, choose a
   figure that can be reproduced faithfully and explain the choice.
2. Create `models/<canonical_model_name>/`, where the folder name follows the
   BNGL skill naming rule: `<descriptive_slug>_<firstauthor><year>`.
3. Create `reference/` inside the new model folder.
4. Write the primary BNGL file as `<canonical_model_name>.bngl`, following
   `skills/bngl/skill.md`. Add complete variant BNGL files when needed, using
   `<canonical_model_name>_<variant>.bngl`. Two content rules drift most often and
   are worth checking explicitly before moving on: every parameter carries an
   inline unit comment (`skills/bngl/skill.md` §1.3.1, `# /s`, `# M`,
   `# dimensionless`), and the header **declares the simulation intent** — whether
   the model is population-based (molecule counts) or concentration-based with
   conversion through `NA` and a volume (§1.3.4). About half the collection is
   missing the second one; do not add to that.
5. Run the BNGL file or files with BioNetGen. Copy generated reference outputs
   (`.gdat`, `.cdat`, `.net`, `.scan`, `.xml`, `.species`, and scan output
   directories) into `reference/`.
6. Extract reported simulation data for the reproduced figure or figures.
   Prefer source tabular data when available. If the paper only reports curves
   in figures, digitize the relevant panel(s) from the source PDF or image,
   calibrate the plotted axes, and save the digitized data in `reference/`.
7. Create the verification artifact — normally `verify_<author><year>.ipynb`,
   or a driver script for a campaign a notebook cannot hold (see "Verification
   Artifact Shape"). Whichever shape, it must run BioNetGen, independently
   implement the paper's equations or expected dynamics in Python/SciPy, compare
   BioNetGen output against that independent implementation quantitatively,
   compare the curated simulation against the reported simulation data
   quantitatively, and save `verify_<author><year>.png` — which must make that
   agreement *readable*, per "The Verification Figure": both comparisons shown,
   residuals or parity as well as overlays, the metrics printed on the figure, and
   the verdict stated.
8. Create `metadata.yaml` following `skills/bngl/skill.md` section 5. List every
   deliverable, including all reference files. A scan output directory may be
   listed as a single `reference/<prefix>_scan/` entry rather than one entry per
   file inside it; list everything else individually. Declare each file's intended
   annotation depth with `documentation_target:`; there is no numeric `rating:`
   field — it was removed 2026-08-04, having never been computed.
9. After all curation artifacts are complete and verified, update the Models table
   in `README.md` with the new collection. Include the folder and BNGL file names,
   a concise description of the primary model, and source reference(s).
10. **Link the fitting-job sibling, if there is one.** Check `pybnf-jobs/` for a
    `<FirstAuthor>-<Year>/` directory from the same paper. When one exists, name it
    in this model's README row and in `metadata.yaml`, and add the reverse link to
    the job's README — the two halves of one paper should each say the other
    exists. Use the corpus's single canonical `<firstauthor><year>` key for the
    paper across `models/`, `pybnf-jobs/`, and `dev/papers/`; check for an existing
    spelling before coining one.

## Reported Simulation Data Verification

Curated models must be verified at two levels:

1. **Model-specification verification:** BioNetGen output from the curated BNGL
   file agrees quantitatively with an independent implementation of the same
   generated network, equations, or expected dynamics.
2. **Reported-data verification:** BioNetGen output from the curated BNGL file is
   quantitatively consistent with the simulation data reported in the source
   paper, supplementary information, or source model documentation.

When reported data are available only as plotted curves:

- render or otherwise obtain the original figure at sufficient resolution;
- record the figure/panel, rendering or extraction method, axis calibration, and
  any plotted scale factors (for example, "10x" legend entries);
- save digitized data under `reference/` using a descriptive filename such as
  `reference/<author><year>_<figure>_<condition>_digitized.csv`;
- compare the curated simulation against the digitized data with explicit,
  reproducible metrics in the verification notebook;
- choose metrics appropriate to the figure type, such as max/median relative
  error for tabular data, peak amplitude and timing errors for trajectories, or
  normalized nearest-curve distances for thick plotted curves;
- state the tolerance and justify it based on the precision of the reported data,
  the figure resolution, line thickness, axis scaling, and any manual or
  automated digitization uncertainty.

If reported simulation data cannot be digitized or otherwise extracted, document
the reason in the notebook and metadata. Do not treat visual qualitative
agreement as sufficient when a figure or table can be quantified.

## Reproducibility Of Committed Reference Output

`reference/` is a regression baseline, so "the fresh run matches the reference" needs a
definition. Re-running the final `.bngl` must reproduce the committed output under these
rules — the verification artifact should assert them, not eyeball them:

- **`.net` and `.xml`:** exact structural match, ignoring comment lines and whitespace.
- **`.gdat` and `.scan`:** for each value, either **relative error ≤ 1e-3**, **or**
  **absolute error ≤ `atol`**, where `atol = 1e-6 × max(|column|)`. The absolute-error
  fallback handles near-zero trajectories, where relative error is undefined or inflated.
- A **seeded stochastic** protocol is held to the same bar: same seed, same numbers. If a
  fresh run does not reproduce a committed seeded trajectory, the seed is not pinned.

Commit reference output for every uncommented simulation protocol — `.gdat`/`.scan` from
`simulate`/`parameter_scan`, the `.net` for any generate-first protocol, the `.xml` for any
network-free one, and the per-point files in scan subdirectories. Reference data for
commented-out protocols is welcome but not required.

## The Verification Figure

`verify_<author><year>.png` exists to answer exactly one question: **is the BioNetGen
output correct?** A reader must be able to answer it *from the figure alone*, without
opening the notebook. A figure that shows what the model does, rather than that the
model is right, has failed — however many panels it has.

### Required content

1. **Both comparisons must be visible and labeled as such.**
   - **BioNetGen vs. the independent implementation** (verification level 1)
   - **BioNetGen vs. the reported data** (verification level 2)

   If the notebook performs a level, the figure must show it. Computing a comparison
   and then leaving it out of the PNG is the most common defect in the collection —
   `innate_immune_response_korwek2023` ships `independent_korwek2023.py` and a
   sixteen-panel figure in which the independent implementation appears nowhere.
   Where a level genuinely cannot be performed, say so *in the figure*, not only in
   the notebook.

2. **Show the disagreement, not just the overlay.** Two curves drawn on top of each
   other look identical whether the error is 1e-6 or 5%, so an overlay alone proves
   nothing. Every comparison needs one of:
   - a **residual panel** — difference or relative difference against the reference,
     with the tolerance drawn as a horizontal band;
   - a **parity plot** — model against reference with the identity line, when the
     comparison is over many points rather than a time course;
   - an **inset** on the overlay, when space is tight.

3. **Put the numbers on the figure.** Every metric the notebook asserts belongs in the
   PNG — as a per-panel annotation (`max rel err 1.3e-6`) or a summary block. A reader
   should never have to open the notebook to learn how well it agreed.

4. **State the verdict.** Name the tolerance and whether it was met. One line is
   enough: `BNG vs SciPy: max rel err 1.3e-6 (tol 1e-4) PASS · BNG vs Fig 2A: median
   rel err 8.8% (tol 15%) PASS`.

### Layout

The default that satisfies the above:

```text
row 1   overlay: BNG + independent          overlay: BNG + reported data
row 2   residual vs independent (+ tol band) residual vs reported (+ tol band)
footer  metrics and verdict
```

Behavior panels — dose-response shapes, mechanism illustrations, anything that shows
what the model *does* — are optional, go last, and must not crowd out the four
verification panels. If the figure is getting wide, drop a behavior panel, never a
residual.

### When the paper has many panels

Reproducing a sixteen-panel published figure panel-for-panel dilutes the verification
into sixteen eyeball judgements. Keep the per-panel overlays if they are informative,
but **add one aggregate agreement panel** that collapses every point into a single
readable claim: a parity plot of all model values against all reference values with
the identity line, coloured by panel, annotated with the pooled metric. That one panel
is what makes correctness a glance instead of a survey.

### Worked examples

`mapk_adaptive_resistance_frohlich2023` is the closest current example of required
content 1 — one axes carrying curated BioNetGen, the independent SciPy integration,
the authors' own source-model trajectory, and the Fig. 2A measurements, each labeled.
It still lacks a residual view and prints no metric, so the reader cannot tell that
the SciPy agreement is 1.3e-6 rather than merely close; add rows 2 and the footer and
it would be exemplary.

## Verification Artifact Shape

The verification is a *contract*, not a file format: run BioNetGen, check it against an
independent implementation, check it against the reported data, emit
`verify_<author><year>.png`. Two shapes satisfy it.

**Notebook (default).** `verify_<author><year>.ipynb`, committed **with its outputs**
so a reader sees the numbers without rerunning. Use it whenever the campaign fits
in a notebook run.

Either shape emits `verify_<author><year>.png`, and the figure follows the overlay
convention in `skills/bngl/skill.md` §2: BioNetGen output as solid lines, the
independent solution as open markers subsampled every ~15th point, and a legend or
subtitle stating the convention. Perfect agreement must be *visible* as markers
sitting on lines, not as one curve hiding another.

**Driver script.** `run_<author><year>.py` plus the committed PNG, and no notebook.
Use this when the campaign cannot reasonably live in a notebook — a network-free
(NFsim/RuleMonkey) model, a large replicate ensemble, or any run measured in hours
rather than seconds. `p53_nhej_dolan2015` (50 loci, 2400-minute NFsim single-cell
run) and `tcr_signaling_chylek2014` are the two current examples. The script must
still perform and report every comparison a notebook would, and the folder must
still ship the PNG. Record in `metadata.yaml` and the README row *why* the notebook
shape was not used.

Do not stretch a notebook around a multi-hour job to satisfy the letter of the
default; do not reach for a driver script to avoid writing the comparisons.

**More than one verification file is fine.** Where a folder verifies genuinely
separate things, add a descriptive suffix rather than forcing one file:
`verify_<author><year>_<aspect>.ipynb` (e.g. `verify_michalski2012_six_state.ipynb`,
`verify_dembo1978_monovalent.ipynb`). Exactly one of them carries the plain
`verify_<author><year>` name and is the entry point.

### Helper scripts

Verification work factored out of the notebook is committed beside it, listed in
`metadata.yaml` with `role: verification`. Use these prefixes — they are the
conventions already in the collection, and they tell a reader what a file is
without opening it:

| prefix | purpose |
|---|---|
| `independent_<author><year>.py` | the level-1 independent implementation, when it is too long to sit inline in the notebook |
| `digitize_<author><year>.py` | recovers plotted curves from the source PDF; must record panel, extraction method, axis calibration, and any legend scale factors |
| `extract_<author><year>.py` | pulls tabular data out of supplementary material |
| `run_<author><year>.py` | drives the simulation campaign (the driver-script shape above, or a long run a notebook then reads) |
| `generator/build_<author><year>.py` | generates the `.bngl` when the model is written programmatically; put the generator and its modules in a `generator/` subdirectory (`lambda_switch_arkin1998`, `amyloid_beta_competing_aggregation_pathways_rana2020`) |

A committed digitizer or generator is what makes a derived artifact reproducible
rather than merely present. If a `reference/` CSV was digitized, the script that
produced it belongs in the folder.

## Point Of Contact

`metadata.yaml` requires a model-specific `point_of_contact`. Do not infer or
default this field from existing model folders, paper authors, commit history, or
the current user.

If `point_of_contact` is not supplied in the triggering prompt, continue paper
reading, model construction, simulation, and verification work when possible.
Before finalizing `metadata.yaml`, ask the user for:

- `name`
- `email`
- optional `orcid`
- optional `github`

Do not mark the model curation complete until the user supplies the required
`point_of_contact.name` and `point_of_contact.email`.

## Deliverables

A complete curated model folder must contain, at minimum:

```text
models/<canonical_model_name>/
├── <canonical_model_name>.bngl
├── metadata.yaml
├── verify_<author><year>.ipynb      # or run_<author><year>.py — see
│                                    #   "Verification Artifact Shape"
├── verify_<author><year>.png        # required in BOTH shapes
└── reference/
    ├── <BioNetGen-generated reference outputs>
    └── <reported or digitized simulation data used for verification>
```

Optional, and common:

```text
├── <canonical_model_name>_<variant>.bngl   # complete variant / related models
├── verify_<author><year>_<aspect>.ipynb    # additional verification files
├── independent_<author><year>.py           # role: verification
├── digitize_<author><year>.py              # role: verification
├── extract_<author><year>.py               # role: verification
├── run_<author><year>.py                   # role: verification
└── generator/                              # programmatically generated models
    └── build_<author><year>.py
```

Additional complete BNGL variant or related files may be added when they are
needed to reproduce the paper or represent closely related published protocols.

## Completion Criteria

The task is not complete until:

- all deliverables exist;
- the active BNGL simulation protocol runs, unless a missing dependency prevents
  execution;
- `reference/` contains committed reference outputs from the final BNGL file or
  files;
- the verification artifact — notebook or driver script — regenerates the PNG and
  reports quantitative agreement between BioNetGen and the independent
  implementation;
- **the PNG itself shows both comparisons, shows the disagreement (residual or
  parity, not overlay alone), prints the metrics, and states the verdict** — a
  reader must be able to judge correctness without opening the notebook;
- it reports quantitative agreement between BioNetGen output and the reported
  simulation data for the reproduced figure, or explicitly documents why reported
  simulation data could not be extracted;
- if the driver-script shape was used, the reason is recorded in `metadata.yaml`
  and the README row;
- a committed notebook carries its executed outputs;
- any reported or digitized simulation data used for comparison are committed in
  `reference/` and listed in `metadata.yaml`, together with the `digitize_`/
  `extract_` script that produced them when they were derived;
- `metadata.yaml` is complete and includes a user-supplied point of contact;
- `README.md` includes the completed model collection in the Models table.

If a required artifact cannot be produced, state exactly which artifact is
missing and why.
