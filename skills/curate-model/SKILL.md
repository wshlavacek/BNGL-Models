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

1. Read `skills/bngl/skill.md`; treat it as the authoritative workflow and
   house-style contract.
2. Read `skills/bngl/rating.md` for trust, annotation, formatting, and reference
   data expectations.
3. Read `skills/bngl/templates/model_skeleton.bngl` before writing a new BNGL
   file.
4. Read `skills/nfsim/SKILL.md` when NFsim behavior, crosslinking, molecularity,
   species observables, cyclic complexes, or network-free simulation may affect
   correctness.
5. Review existing `models/` folders for naming, metadata, verification notebook,
   PNG, reference-data, and README conventions.
6. Read all relevant files in the requested `dev/papers/<folder>/`, including the
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
   `<canonical_model_name>_<variant>.bngl`.
5. Run the BNGL file or files with BioNetGen. Copy generated reference outputs
   (`.gdat`, `.cdat`, `.net`, `.scan`, `.xml`, `.species`, and scan output
   directories) into `reference/`.
6. Extract reported simulation data for the reproduced figure or figures.
   Prefer source tabular data when available. If the paper only reports curves
   in figures, digitize the relevant panel(s) from the source PDF or image,
   calibrate the plotted axes, and save the digitized data in `reference/`.
7. Create `verify_<author><year>.ipynb`. The notebook must run BioNetGen,
   independently implement the paper's equations or expected dynamics in
   Python/SciPy, compare BioNetGen output against the independent implementation
   quantitatively, compare the curated simulation against the reported
   simulation data quantitatively, and save `verify_<author><year>.png`.
8. Create `metadata.yaml` following `skills/bngl/skill.md` section 5. List every
   deliverable, including all reference files.
9. After all curation artifacts are complete and verified, update the Models table
   in `README.md` with the new collection. Include the folder and BNGL file names,
   a concise description of the primary model, and source reference(s).

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
├── verify_<author><year>.ipynb
├── verify_<author><year>.png
└── reference/
    ├── <BioNetGen-generated reference outputs>
    └── <reported or digitized simulation data used for verification>
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
- the notebook regenerates the PNG and reports quantitative agreement between
  BioNetGen and the independent implementation;
- the notebook reports quantitative agreement between BioNetGen output and the
  reported simulation data for the reproduced figure or explicitly documents why
  reported simulation data could not be extracted;
- any reported or digitized simulation data used for comparison are committed in
  `reference/` and listed in `metadata.yaml`;
- `metadata.yaml` is complete and includes a user-supplied point of contact;
- `README.md` includes the completed model collection in the Models table.

If a required artifact cannot be produced, state exactly which artifact is
missing and why.
