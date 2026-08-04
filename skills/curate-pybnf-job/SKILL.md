---
name: curate-pybnf-job
description: Use when turning a published systems-biology / rule-based modeling paper (a PubMed Central article, PMCID/DOI, a dev/papers/ folder, or a PEtab problem) into a PyBNF edition-2 parameter-fitting job under pybnf-jobs/<FirstAuthor>-<Year>/<slug>/. Picks the archetype (hand-built BNGL, PEtab/SBML-imported, or BPSL constraint-bearing), reconstructs the model, extracts the fit data into .exp files (and qualitative claims into BPSL .prop constraints), authors the annotated .conf, and clears the acceptance bar — a declared reference objective J* plus a measured optimality gap OG = J_paper - J*, recorded in jstar.txt / nominal_check.json — then ships the slug README, make_reproduction.py and reproduction PNG. Trigger whenever the user wants to add a PyBNF fitting job from a paper, expand pybnf-jobs/, build a PEtab parameterization problem from a model + data, fit qualitative/BPSL constraints, score a job against a benchmark optimum, or "make a job setup" from a paper — even if they don't say "PEtab", "edition 2", "optimality gap", or "constraints".
---

# Curate PyBNF Job

Turn a published model + data into a **PyBNF edition-2 parameter-fitting job** and add it to
the corpus at `pybnf-jobs/`. It is the fitting-job sibling of `curate-model`: same "point at a
paper, reconstruct faithfully, verify quantitatively" discipline, but the deliverable is a
runnable fit whose quality is a **measured number**, not a claim.

Prefer an explicit source in the request — a PMCID/URL, a DOI, or a local paper folder:

```text
Use curate-pybnf-job for PMC5334499  (or: for dev/papers/<Folder>)
```

## Inputs & where things go

- **Input:** a paper in PubMed Central (PMCID or URL), a DOI, a local paper folder
  (`dev/papers/<FirstAuthor>-<Year>/` with a PDF + any author-provided model/data files), or an
  upstream PEtab problem.
- **Output — this repo is the corpus.** Jobs are filed by source paper under
  `pybnf-jobs/<FirstAuthor>-<Year>/<slug>/`, one directory per paper holding one or more job
  **slugs**, because a single paper often yields several fitting problems (e.g.
  `Rukhlenko-2022/{cstar_trka, cstar_trkb, cstar_skmel133, cstar_skmel133_bmra}/`). Each slug is
  **self-contained and self-documenting**: it carries its own `README.md`, its reproduction
  script and PNG, and its scoring provenance.
- **Upstreaming to PyBNF is a separate, optional, final step** (step 9). `~/Code/PyBNF/examples/
  real-world/` is a small teaching subset of this corpus, not its home. Do not write a new job
  there first, and do not treat `_manifest.py` as the registry — the corpus registry is the
  paper-level README plus this repo's root `README.md`.
- **Paper-level landing README:** every `<FirstAuthor>-<Year>/` directory gets a `README.md` —
  the full citation, the shared model (what was reconstructed, from which supplementary file),
  a one-row-per-slug table (what each fits · archetype · data source · **J\*, OG, status**), and
  the source links. When the paper also has a curated entry in `models/`, link it, and link back
  from that model's README row.
- **Never commit the PDFs.** `dev/` is git-ignored; it is the parking garage, not a deliverable.

## Pick the archetype first — it sets everything downstream

| | **A · hand-built BNGL** | **B · PEtab/SBML-imported** | **C · BPSL constraint-bearing** |
|---|---|---|---|
| when | the paper's model is rule-based, or you reconstruct it | an upstream PEtab problem exists (e.g. a benchmark collection) | the paper asserts *qualitative* facts a time series can't carry |
| model file | `<slug>.bngl` you author | `model_<slug>.xml` copied verbatim from upstream | `<slug>.bngl` (+ per-mutant variants) |
| conf | `bngl_backend = bngsim` | `sbml_backend = bngsim` | `bngl_backend = bngsim` |
| data | `.exp` you extract/digitize | `.exp` + `_measparams.tsv` emitted by the importer | `.prop` (± `.exp` for data fusion) |
| PEtab export | in-subset ⇒ round-trips | it *came from* PEtab | **refused** — native-only, by design |
| extra provenance | — | `upstream.json` pinning the source commit + per-file sha256 | — |
| acceptance | OG vs. J\* | OG vs. the benchmark's published J\* | OG **and** `job_type = check` satisfaction |

Archetypes B and C are not exotic: B is a third of the corpus (23 of 75 confs), C is the
capability that makes PyBNF worth benchmarking. A paper can yield slugs of more than one kind —
prefer that over collapsing them.

## Required reading

Before authoring anything, read (they are the source of truth, not your memory):

1. `references/og-acceptance.md` — **the acceptance bar.** What J\* is, the three provenance
   tiers it can come from, how OG is computed, the `OG < 1.92` threshold, the status vocabulary,
   and what to do when the objective has no likelihood. Read this *before* you choose an
   objective, because the choice determines whether the job can be scored at all.
2. `references/edition2-conf-reference.md` — the full edition-2 `.conf` surface
   (`job_type`, `objective`/`noise_model`, `experiment:`/`condition:`, `*_var` free params).
   What is **[E2]** vs. rejected **[LEGACY]**.
3. `references/real-world-anatomy.md` — exactly what a slug folder contains, file by file.
4. `references/petab-compliance.md` — what "PEtab.v2-compliant" means (a *verified round-trip
   property* of the native conf) and the **PEtab-exportable subset** you must stay inside for
   archetype A.
5. `references/bpsl-constraints.md` — BPSL `.prop`/`.con` files: the grammar, how to attach them,
   `job_type = check`, and the fact that a constraint-bearing job is **native-only**. Read this
   whenever the paper reports orderings, thresholds, monotonic dose-response, "peaks before",
   bistability, or oscillation.
6. `skills/bngl/skill.md` + `skills/bngl/templates/model_skeleton.bngl` — house style for any
   BNGL you author. Read `skills/nfsim/SKILL.md` if the model is network-free.
7. Two or three existing slugs for the idioms to imitate — pick by archetype:
   `Kozer-2013/egfr_ode/` (A: time course + dose-response scan, retained network cap),
   `Grein-2026-benchmark-subset-I/Bertozzi_PNAS2020/` (B: the complete OG provenance set),
   `Rukhlenko-2022/cstar_skmel133_bmra/` (C: data fusion + `check` conf).

`validate-pybnf-job` is the auditor sibling. It grounds a finished job in the primary literature
and writes `VALIDATION.md`. **Do not write `VALIDATION.md` from this skill** — curation produces
the job and its score; validation earns the confidence. Hand off when this skill's completion
criteria are met.

## Environment

The scripts and `pybnf` run in the **PyBNF** environment. `BNGPATH` is needed even to parse a
model. `scripts/*` `chdir` into the conf's folder themselves, so they take absolute paths from
anywhere:

```bash
export BNGPATH="$HOME/Simulations/BioNetGen-2.9.3"   # folder containing BNG2.pl
PY=~/Code/PyBNF/.venv/bin/python
SKILL=~/Code/BNGL-Models/skills/curate-pybnf-job
JOB=~/Code/BNGL-Models/pybnf-jobs/<FirstAuthor>-<Year>/<slug>
```

Run `pybnf` itself **from inside the slug folder** — `model:`, `data:` and `output_dir` resolve
against the working directory. If BNGPATH/BNG2.pl or the PyBNF venv can't be found, ask rather
than guessing.

## Workflow

1. **Read the paper and extract the fitting problem.** Identify: the model (species, rules,
   parameters, initial conditions); the **data to fit** and which figure/table it lives in; the
   **experimental design** (time course? dose-response? pre-equilibration/washout?); the
   **observables** (what each measured quantity maps to in the model); the **free parameters**
   and their published values / plausible ranges; whether the dynamics are ODE, SSA, or NFsim;
   any **qualitative properties** the paper asserts; and — critically — **what the paper reports
   about the fit itself**: an objective value, a parameter table, confidence intervals, a
   figure. Write this down before building. Then pick the archetype from the table above.

2. **Declare J\* and the reproduction target *before* you fit.** This is the discipline change:
   decide what "the paper's result" is while you still have the paper open, not after you have a
   number you like. Per `references/og-acceptance.md`, J\* comes from one of three tiers —
   **T1** a published/benchmark objective value, **T2** PyBNF's objective evaluated at the
   paper's published best-fit parameters, or **T3** the corpus's own best-known objective
   (a regression anchor, not an optimality claim). Record the tier; it is what the status badge
   means. A job whose J\* provenance you cannot state is a job you cannot score.

   **For the common T2 case, do not derive it by hand** — once the model and conf exist (steps
   4–6), `scripts/make_jstar.py` computes it:

   ```bash
   $PY $SKILL/scripts/make_jstar.py $JOB --write
   ```

   It evaluates PyBNF's objective at the model's shipped nominal parameters through
   `likelihood_information_criteria` — the same path that writes `information_criteria.txt` at
   the end of a real fit — so the number lands on the comparable `-log_likelihood` scale, and it
   writes both `jstar.txt` and `nominal_check.json`. It also gates on the *resolved* objective
   rather than the conf token, cross-checks its own scored-point count against `ic.n`, and
   refuses with a specific reason when T2 is structurally undefined (an estimated noise nuisance,
   a constraint-only BPSL job, a published point outside the conf's own search box). Treat a
   refusal as information about the job, not as a tooling failure — `references/og-acceptance.md`
   §3 lists the cases.

3. **Name and create the folder.** `pybnf-jobs/<FirstAuthor>-<Year>/<slug>/`, with a short
   descriptive slug (`egfr_ode`, `tlbr`, `cstar_trka`). One author-year directory per paper; each
   additional fitting problem the paper yields becomes a new slug *beside* the others. Append
   `b`/`c` to the year (`Smith-2020b`) only to disambiguate two papers sharing first author and
   year. **Use one canonical `<FirstAuthor>-<Year>` key per paper across `models/`,
   `pybnf-jobs/`, and `dev/papers/`** — check for an existing spelling before inventing one.
   Create or update the paper-level `README.md`.

4. **Build the model, fitting-ready and on the edition-2 surface.**

   *Archetype A — reconstruct the BNGL from scratch:*
   - **Strip the *simulation* actions, but KEEP *network-definition* directives.** Remove
     `simulate`/`parameter_scan`/`setConcentration` (synthesized from the conf), but **retain a
     `generate_network({...,max_stoich=>{...}})` line** whenever the network is only finite under
     a `max_stoich`/`max_agg`/`max_iter` cap. That cap is part of the *model*, not the
     experiment: strip it and pybnf falls back to a bare `generate_network({overwrite=>1})`
     (`pset.py:638`) → **unbounded network → silent hang**, or a wrong finite network. pybnf
     captures a retained line (`pset.py:617`) and the job stays PEtab-exportable. Network-free
     (NFsim) models keep no actions block.
   - Fitted rate constants are bare `id nominal` declarations; conf free-parameter names bind to
     these ids **by name** (no `__FREE` alias, ADR-0034).
   - **Observable/function names are the contract with the `.exp` header *and* any `.prop`** — a
     `Molecules`/`Species` observable is a plain name; a `functions` entry appears **with
     parentheses** in the exp header. Every name a `.prop` references must exist as an observable.
   - For pre-equilibration, add a boolean gate parameter (e.g. `Ligand_isPresent 0`) toggled by
     two `condition:` states.
   - Follow `skills/bngl/skill.md` house style. Confirm the model builds through BNG2.pl and
     behaves as the paper describes.

   *Archetype B — import, don't author:* run PyBNF's PEtab importer, copy the SBML **verbatim**,
   and write `upstream.json` pinning the source commit plus a per-file sha256 of LF-normalized
   content. Never edit an upstream file; everything you author sits beside it. State in the slug
   README which bytes are copied and which are yours.

5. **Extract the data into `.exp` file(s).** First column is the independent variable (`time`
   for a time course; the **model-parameter name** for a dose/scan column); remaining columns are
   observables named exactly as in the model. Add `<obs>_SD` columns only if the paper reports
   per-point uncertainty; use `NaN` for missing points. Prefer source tabular/supplementary data;
   if only plotted curves exist, **digitize** the panel — record figure/panel, extraction method,
   axis calibration, and any legend scale factors, exactly as `curate-model` prescribes, and
   commit the extraction script. The digitized target is what "reproducing the paper" is measured
   against, so keep it faithful and documented.

   **If the job uses BPSL** (archetype C), also author the `.prop` file(s): one qualitative
   statement per line in the BPSL grammar (`<obs> <op> <obs|const> <always|once|at …|between …>
   weight <w>`), using dotted `suffix.obs` to compare across experiments/mutants. Translate each
   qualitative claim in the paper into one line; weight stronger claims higher.

6. **Author `<slug>.conf`** from `templates/job_setup.conf`. Open it with a banner comment citing
   the paper (PMCID/DOI), the figure each `.exp` came from, and the declared J\* and its tier.
   Then:
   - **Objective — a modeling choice; check it is scoreable.** Under edition 2, `objective = sos`
     resolves to a Gaussian likelihood with σ ≡ 1, so it is scoreable as-is — use it when the
     paper's fit is unweighted least squares. Use `chi_sq` when the `.exp` carries `_SD`, and a
     fitted `sigma` (`noise_model = normal, sigma = fit sd_<obs>`) when the noise scale is
     unknown. What is *not* scoreable is a resolved objective with no per-point log-likelihood
     (`objective.py:73-96`) — in practice the **legacy edition-1** `objfunc = sos` path, which
     resolves to `SumOfSquaresObjective`. Confirm with `make_jstar.py`, which gates on the
     resolved object rather than the token. Stay inside the PEtab-exportable subset for archetype
     A — no `normalization`, no `cumulative`, no `neg_bin`/`lognormal`.
   - **`job_type` — choose by running candidates, not by assumption.** The corpus uses `de` (25),
     `gntr` (20), `ss` (10), `check` (6), `cmaes` (5), `lbfgs` (4), `am` (1). Defaults that have
     earned their place: **`gntr`** (Fisher/Gauss-Newton trust region, ADR-0068) is the gradient
     workhorse and the one method that handles an **estimated** noise scale — plain `trf` refuses
     it; **`cmaes`** with IPOP restarts where the gradient path refuses the problem or the
     landscape is strongly multimodal; **`de`/`ss`** for a first global sweep on a small model;
     **`am`** when the paper reports posteriors rather than point estimates; **`check`** for a
     constraint-satisfaction job. Gradient methods need `bngsim` and differentiable dynamics —
     a model with discrete events will be refused (`_require_differentiable_dynamics`).
   - **Budget is part of the result.** The collection default (20 starts × 500) is *not* tuned.
     If a slug needs 100 × 1000 to reach its basin, ship that budget and say so in the README.
   - **Free parameters:** one `*_var` per fitted id, bracketing the published value.
   - **BPSL:** attach constraints via the experiment's `data:` — fusion (`data: <n>.exp, <n>.prop`)
     or constraint-only (`data: <n>.prop` with a required `t_end:`). Ship a companion
     `<slug>_check.conf` with `job_type = check`.

7. **Verify, and score.** With BNGPATH set:
   - **Tier-1 (parses & well-formed):** `$PY $SKILL/scripts/check_conf.py $JOB/<slug>.conf` —
     edition 2, `job_type` resolves, data and/or constraints bound, free params bind by id.
   - **PEtab round-trip — archetype A only:**
     `$PY $SKILL/scripts/petab_roundtrip.py $JOB/<slug>.conf --job-type <jt>` (export → `petab.v2`
     lint clean → import). A reported *non-exportable feature* means rework the conf; a report of
     *BPSL constraint data* means you are archetype C and this check does not apply.
   - **BPSL satisfaction — archetype C only:** run the `check` conf at the fitted (or published)
     parameters and confirm `Satisfied M out of M constraints`, or document which the paper
     itself does not require.
   - **The acceptance bar — run the fit and measure OG:**
     ```bash
     cd $JOB && pybnf -c <slug>.conf
     $PY $SKILL/scripts/score.py $JOB output
     ```
     `score.py` reads `jstar.txt` and the run's `Results/information_criteria.txt`, computes
     `J_paper = -log_likelihood` and `OG = J_paper - J*`, and reports **solved** iff `OG < 1.92`
     (χ², α=0.05, 1 dof). Commit the provenance it depends on: `jstar.txt`,
     `best_fit_params.txt`, `information_criteria.txt`, and `nominal_check.json` recording the
     tier, the numbers, and a one-paragraph interpretation. Full rules, including what to do when
     OG is large and when a J\*-at-nominal check is the honest claim, are in
     `references/og-acceptance.md`.
   - **Reproduce the paper's figure:** write `make_reproduction.py` that simulates at the fitted
     (and, where the paper reports them, the published) parameters and overlays the target data,
     saving `<slug>_reproduction.png`. Report a quantitative metric — max/median relative error,
     peak amplitude/timing, or nearest-curve distance — and justify the tolerance from the data's
     precision. Both script and PNG are **committed** slug files, not scratchpad evidence.

8. **Write the slug `README.md`.** It is the entry point and must stand alone: what the job fits
   and why; the model and where it came from; the training data and its figure/table provenance;
   the free parameters; the archetype and its consequences (PEtab-exportable? native-only?);
   **a `**Run cost:**` line directly under the title** (see below); **the status line — J\*, its
   tier, OG, and the badge**; the exact run command; and a link to `VALIDATION.md` once
   `validate-pybnf-job` has been run. Update the paper-level README's slug table — including its
   **`run cost`** column — and the root `README.md` if the job pairs with a curated model.

9. **Report — and close the loop upstream.** Summarize the new slug, the verification results
   (tier-1; PEtab round-trip *or* BPSL `check`; **OG and status**; the paper-reproduction metric),
   and the README edits. **If curating this job exposed a bngsim or PyBNF defect, file it and cite
   the issue number in the slug README and `nominal_check.json`** — that is the corpus doing its
   job, and it is the single most valuable output after the job itself. Optionally, if the job is
   small, fast, and pedagogically clean, propose promoting it to `~/Code/PyBNF/examples/
   real-world/` with a `_manifest.py` entry (`templates/manifest_entry.py`) — a separate PR
   against a separate repo, never a substitute for landing it here.

## The acceptance bar, in one paragraph

A job is not done when it runs. It is done when it has a **declared reference objective J\*** with
a stated provenance tier, and a **measured optimality gap** `OG = J_paper − J*` where
`J_paper = −log_likelihood` from `Results/information_criteria.txt`. `OG < 1.92` is **solved**.
An OG measured at the model's nominal point rather than from a fit is **objective validated** —
it validates the model, the observables and the objective, and claims nothing about the optimizer;
mark it `†` and say so. A job that imports, simulates and scores but whose nominal point is not
the published optimum and which has not been driven to `OG < 1.92` is **setup only** — a
legitimate, useful deliverable, but it must not be described as reproducing the paper's fit. The
old bar ("a real bngsim fit reaches a finite objective") passes a broken model and is retired.

## Declare the scale — and that is where `heavy` comes from

Every slug states what it costs to run, twice: a `**Run cost: \`<class>\`** — <why>.` line
directly under the slug README's title, and a **`run cost`** column in the paper-level slug
table. Same four-value vocabulary as `models/` (`skills/bngl/skill.md` §5.6), so one word means
one thing across the repository. Give the *why* — the evaluation count and what each evaluation
costs — not just the word; that is what lets the next reader re-derive the class when the conf's
budget changes.

**Call the column `run cost`, not `scale`.** In `pybnf-jobs/` the word *scale* is already taken —
PEtab's observable measurement scale (`lin`/`log10`/`ln`), which the Grein subset tabulates per
slug and which is load-bearing for how its objective constants are restored. `models/` has no
such conflict and uses `scale:` there.

| scale | can I run this on a laptop? | budget for a full fit |
|---|---|---|
| `trivial` | yes, right now | ≲ 1 min · a `job_type = check` or `sim` conf, ~1 evaluation |
| `minutes` | yes, walk away | ≲ 1 core-hour |
| `hours` | technically, painfully | ≲ 48 core-hours — a workstation, or overnight |
| `cluster` | no — it needs a scheduler | > 48 core-hours, or many-node islands |

**`heavy` is not a field. `heavy ≡ scale ∈ {hours, cluster}`.** That is the whole definition. The
manifest's `heavy=True`, and the tier-2 exclusion it drives
(`references/real-world-anatomy.md`), are *consequences* of the class rather than a separate
boolean a curator can forget to flip — which is exactly how the flag went missing before.

**Assign it statically — never by running the fit to convergence to find out.**

```text
core-hours  ≈  evaluations × seconds-per-evaluation / 3600
evaluations  =  population_size × max_iterations   (× smoothing, for SSA/NF replicates)
```

- **`population_size` × `max_iterations` is the number that matters.** It spans four orders of
  magnitude across this corpus — 60 for a `Mallela-2024` lbfgs slug, 1,062,000 for
  `Lang-2024/v3_2_0`. Read it straight off the conf.
- **Seconds per evaluation** comes from the model, not the optimizer: inherit the `models/` entry's
  `scale:` when the paper has a curated sibling, otherwise use its reaction count and method by
  §5.6.1. Anchor: `Suofu-2017/mito_camp` records 39,000 evaluations in 28 minutes on 16 cores —
  **0.69 s/evaluation** for a 30-reaction ODE model. Network-free and SSA evaluations are one to
  three orders of magnitude worse, and `smoothing = N` multiplies every one of them by *N*.
- **`wall_time_sim` is a ceiling you declared, not a measurement.** Half this corpus sets it to
  3600 as a safety net on jobs that finish in seconds. Use it to break a tie, never to classify.
- **Free-parameter count is the honest tiebreaker at the top.** A 177-parameter problem is
  `cluster` even under a generous per-evaluation estimate, because the budget it *needs* to
  converge is not the budget the conf currently declares.

Sanity check your assignment against the two slugs the corpus already calls heavy in prose:
`Lang-2024/v3_2_0` (1.06M evaluations, 177 free parameters) and both `Miller-2026` slugs
(240k–300k evaluations; "the authors used an HPC cluster"). If your procedure does not put those
in `cluster`, it is miscalibrated.

Two costs, stated separately when they differ: a `cluster` fit whose **reproduction** script
replays committed best-fit parameters in seconds should say so, because that is the part a reader
without a cluster can actually run.

## Guardrails that bite

- **A missing `information_criteria.txt` does not mean "unscoreable".** PyBNF writes it at the end
  of a *fit*, so a job never run has none whatever its objective — in the Grein subset both
  `Smith_BMCSystBiol2013` (`sos`) and `Weber_BMC2015` (`chi_sq`) lack it because both are ⚪
  setup-only. Genuine unscoreability means the *resolved* objective carries no per-point
  log-likelihood, which under edition 2 essentially never happens: `objective = sos` resolves to
  `Gaussian(sigma=1)`. See `references/og-acceptance.md` §2.
- **PEtab-exportable subset (archetype A).** `export_job` fails loud on features PEtab v2 can't
  express. The frequent trap: `normalization = init/peak/…` is **not** exportable. If the paper's
  data is relative to a reference, encode the reference in the model/observable instead.
- **BPSL constraints are native-only.** Any `.prop`/`.con` on an experiment's `data:` makes the
  job non-exportable — `export_job` raises `NotImplementedError`. That is a property of the
  archetype you chose, not a defect: verify with `job_type = check`, register it native-only.
- **A stripped network cap is a silent hang, not an error.** See step 4.
- **Untuned budgets mislead in both directions.** A default-budget run that lands far from J\*
  says nothing about PyBNF's optimizers; a slug in the corpus converged *worse* than its own
  nominal point on 20 starts. Tune, or label the row honestly.

## Deliverables

```text
pybnf-jobs/<FirstAuthor>-<Year>/                 # e.g. Rukhlenko-2022/
├── README.md                    # paper landing page: citation, shared model, slug table
│                                #   (fits · archetype · data · J* · OG · status), source links
└── <slug>/                                      # e.g. cstar_trka/
    ├── README.md                # stands alone; carries the status line
    ├── <slug>.bngl              # archetype A/C — no simulation actions; KEEP a needed
    │                            #   generate_network(...max_stoich...) directive
    ├── model_<slug>.xml         # archetype B — SBML verbatim from upstream
    ├── upstream.json            # archetype B — pinned commit + per-file sha256 (LF-normalized)
    ├── <slug>.conf              # edition-2 job setup, banner-commented with paper + J* tier
    ├── <slug>_check.conf        # archetype C — job_type = check companion
    ├── <data>.exp               # headers == model observable names (A/B, and C data fusion)
    ├── <slug>.prop              # archetype C — BPSL constraints
    ├── make_<data>.py           # data extraction / digitization script, when data was derived
    ├── make_reproduction.py     # simulates at fitted/published params, overlays the target
    ├── <slug>_reproduction.png  # the overlay
    ├── jstar.txt                # the declared reference objective J*
    ├── best_fit_params.txt      # the fit that produced the reported OG
    ├── information_criteria.txt # k, n, log_likelihood, AIC/BIC/AICc from that fit
    └── nominal_check.json       # J* tier, J_paper, OG, k, n, optimizer, interpretation
# VALIDATION.md is written later, by validate-pybnf-job — not by this skill.
```

## Completion criteria

Not complete until:
- the archetype is stated, and the folder matches its column in the archetype table;
- `scripts/check_conf.py` passes (edition 2, job_type resolves, data and/or constraints bound,
  free params bind by id with no `__FREE`);
- **either** `scripts/petab_roundtrip.py` passes (A: export → lint clean → import) **or** a
  `job_type = check` run reports the intended satisfaction (C) and the slug is registered
  native-only;
- **`jstar.txt` exists** (for a T2 anchor, generated by `scripts/make_jstar.py`), **its
  provenance tier is stated, and `scripts/score.py` reports an OG**,
  with the status badge (solved / objective validated / setup only) recorded in
  `nominal_check.json`, the slug README, and the paper-level slug table — **or** it is documented
  precisely why no J\* can be declared for this paper;
- `make_reproduction.py` and `<slug>_reproduction.png` are committed, and the reproduction metric
  and its tolerance are stated;
- the slug `README.md` and the paper-level `README.md` exist, agree with the conf and the shipped
  artifacts, and use the corpus's one canonical `<FirstAuthor>-<Year>` key for this paper;
- **the scale is declared** in the slug README status line and the paper-level slug table, from
  `population_size × max_iterations` and the model's per-evaluation cost — and if it comes out
  `hours` or `cluster` (i.e. `heavy`), the README says what a reader without that hardware can
  still run;
- the conf banner cites the paper (PMCID/DOI) and the figure/table each `.exp` (and any
  qualitative property) came from;
- any bngsim/PyBNF defect this curation exposed is filed and cited.

If a required artifact cannot be produced, state exactly which one and why.
