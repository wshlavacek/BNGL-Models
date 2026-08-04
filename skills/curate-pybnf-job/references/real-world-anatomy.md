# Anatomy of a `pybnf-jobs/` slug folder

What a slug under `pybnf-jobs/<FirstAuthor>-<Year>/<slug>/` must contain to (a) run, (b) be
scored against its reference objective, and (c) stand alone for a reader. Grounded in the PyBNF
source; every schema claim cites `file:line`.

A slug is **self-contained and self-documenting**: `{model, conf, data}` + its own `README.md` +
the reproduction pair + the OG provenance set. It is not an entry in a central manifest — the
corpus registry is the paper-level `README.md` beside it.

> **On `~/Code/PyBNF/examples/real-world/`.** That directory is a small teaching subset of this
> corpus (8 papers, against 20 here), not its home. Its `_manifest.py` + two-test-tier machinery
> is documented in the **Appendix** and applies only when you *promote* a finished slug there —
> a separate PR against a separate repo. Do not author a new job there.

## Contents
1. File inventory
2. File formats (with the `receptor` template)
3. The slug `README.md`
4. Appendix: promoting a slug to PyBNF's `examples/real-world/`
   (the `_manifest.py` entry, the two test tiers, the README edits)

---

## 1. File inventory

Each slug lives in its own folder, **grouped by source paper under a `<FirstAuthor>-<Year>/`
directory** (e.g. `Rukhlenko-2022/cstar_trka/`; one paper's several jobs sit side by side there).
All paths inside the `.conf` resolve relative to the slug folder, so run `pybnf` from inside it.

**Core — every slug:**

| file | role | required |
|---|---|---|
| `README.md` | the entry point: what is fit, the model, the data provenance, free params, archetype, **the status line (J\*, tier, OG, badge)**, the run command | yes |
| `<slug>.conf` | the edition-2 job setup, banner-commented with the paper and the J\* tier | yes |
| `jstar.txt` | the declared reference objective (see `og-acceptance.md`) | yes |
| `nominal_check.json` | J\* tier + source, `J_paper`, `OG`, `k`, `n`, optimizer, status, and a prose `interpretation` | yes |
| `best_fit_params.txt` · `information_criteria.txt` | the fit that produced the reported OG, and its `k`/`n`/`log_likelihood` | yes, once a fit has been run |
| `make_reproduction.py` · `<slug>_reproduction.png` | simulate at fitted/published params, overlay the target data | yes |
| `make_<data>.py` | the extraction/digitization script, when the data was derived rather than transcribed | when data was derived |
| `VALIDATION.md` | written later by `validate-pybnf-job`, **not** by curation | no |

**Archetype A — hand-built BNGL:**

| file | role |
|---|---|
| `<slug>.bngl` | the model — **no *simulation* actions** (`simulate`/`parameter_scan`; synthesized from the conf). **RETAIN network-generation directives** (`generate_network` with `max_stoich`/`max_agg`/`max_iter`) when the network is only finite/correct under them — see §2 |
| `<data>.exp` | ≥1 data file; each `experiment:` binds one; column headers ARE model observable/function names |
| `*_ground_truth.bngl` | synthetic-data jobs only: the model at known-true params (documentation) |

**Archetype B — PEtab/SBML-imported:**

| file | role |
|---|---|
| `model_<slug>.xml` | the SBML, **copied verbatim** from upstream — never edited |
| `upstream.json` | the pinned upstream commit + per-file sha256 of LF-normalized content, so a re-copy cannot silently drift. `.gitattributes` stores `*.xml` as LF for exactly this reason |
| `<data>.exp` · `<data>_measparams.tsv` | PyBNF-format translations of the upstream measurement tables, emitted by the importer |

There is no `<slug>.bngl`: the model is the SBML, and the conf uses `sbml_backend = bngsim`.
State in the README which bytes are copied and which are yours.

**Archetype C — BPSL constraint-bearing:**

| file | role |
|---|---|
| `<slug>.prop` (`.con`) | BPSL constraint file(s), attached on an `experiment:`'s `data:`. Makes the job **native-only** (not PEtab-exportable). See `bpsl-constraints.md` |
| `<slug>_check.conf` | the `job_type = check` companion that reports constraint satisfaction |
| `<variant>.bngl` | per-mutant/per-condition model variants when constraints compare across them |

---

## 2. File formats (canonical template = `receptor/`)

### `<name>.bngl` (edition-2, fitting-ready)
- Standard `begin model … end model` with `parameters`, `molecule types`, `seed
  species`, `observables`, `functions`, `reaction rules`. Follow `skills/bngl/skill.md`
  house style.
- **Actions block: strip *simulation* actions, but KEEP *network-definition* directives.**
  The two kinds of action are not the same thing:
  - *Simulation/experiment actions* (`simulate`, `parameter_scan`, `setConcentration`,
    `t_end`, `method`) → **remove them**; they are synthesized from the conf's
    `experiment:`/`condition:` lines (`receptor.bngl:2-7`), PEtab-style.
  - *Network-generation directives* (`generate_network({...,max_stoich=>{...},max_agg=>...,
    max_iter=>...})`) → **KEEP them** when the model needs them. `max_stoich`/`max_agg`/
    `max_iter` are part of the *model's specification* (they are what make an
    aggregation/polymerization network **finite**), not the experiment design. Strip them
    and pybnf falls back to a bare `generate_network({overwrite=>1})` (`pset.py:638-639`) —
    **no cap → unbounded network → generation never terminates** (a silent hang), or, worse,
    a finite-but-*different* network → quietly wrong results. pybnf **captures** an existing
    `generate_network` line from the model and uses it in place of that default
    (`pset.py:617-619`; bngsim routes it to the `.net` backend,
    `bngsim_model/classification.py:291`), and the job stays PEtab-exportable, so keeping the
    directive is the supported fix. **Test:** if your model is network-generating and its
    rules do not *themselves* bound complex size, it needs a retained `generate_network` with
    the cap. Network-**free** (NFsim) models have no `generate_network` and correctly keep no
    actions block. Canonical example: `Kozer-2013-2014/egfr_ode` keeps only
    `generate_network({overwrite=>1,max_stoich=>{EGF=>4,EGFR=>4}})`.
- Fitted rate constants are bare `id nominal` declarations, e.g. `KD1 1.0`
  (`receptor.bngl:52`); the optimizer overrides them in place. The free-parameter names
  in the conf must exactly equal these ids (ADR-0034).
- **Observable/function names are the contract with the `.exp` header.** A
  `Molecules`/`Species` observable → plain name in the exp header (`receptor.bngl:140`
  `RLbonds`). A `functions` entry → name **with parentheses** in the exp header
  (`tlbr.bngl:78` `FL()` ↔ exp column `FL()`; the manifest lists it without parens).
- For pre-equilibration, a boolean gate parameter (e.g. `Ligand_isPresent 0`,
  `receptor.bngl:90`) toggled by the two conditions, used to gate a rate.

### `<name>.conf`
See `edition2-conf-reference.md`. The minimal ODE shape (`receptor.conf`):
`output_dir` → `edition = 2` → `bngl_backend = bngsim` → `model:` → `condition:`(s) →
`experiment:` → `job_type` → `objective` (`sos` if the exp has no `_SD`, else `chi_sq`)
→ budget → `*_var` free params bound by id.

### `<data>.exp` — exact format
Header line 1 begins with `#`; columns whitespace- or tab-separated; names match model
observable/function names, independent variable first.
- **First column = independent variable.** Time course ⇒ `time` (`receptor.exp:1`
  `# time RLbonds pR`). `parameter_scan` ⇒ the scanned **model parameter** name (the
  dose column), e.g. `LTconc`, `IGF1_cold_conc`.
- **Remaining columns = observables**, one per model observable/function.
- **`_SD` columns** (optional): `<obs>_SD` supplies per-point σ and switches the
  objective to `chi_sq`; absent ⇒ `sos`.
- **`NaN`** marks a missing point (value and its `_SD`), skipped in scoring.
- Numbers may be plain or scientific; values may be negative.

Two real headers:
```
# time RLbonds pR                                      # -> sos
#	IGF1_cold_conc	IGF1_hot_bound	IGF1_hot_bound_SD  # -> chi_sq (tab-separated)
```

---

## 3. The slug `README.md`

The slug README is the deliverable a reader meets first, and it must stand alone — assume they
have neither the paper nor this skill open. The corpus converged on this shape:

```markdown
# <slug> — <one-line biology>, <ODE|SSA|NFsim> (PyBNF edition-2 job)

<2-4 sentences: what this job fits and why it is worth having.>

## Status
J* = <value>  (tier <T1|T2|T3>, <source>) · OG = <value> · <✅ solved | 🟢 objective validated †
| 🔁 regression-anchored | ⚪ setup only> · <PEtab-exportable | native-only (BPSL)> · <heavy?>

## Reference
<full citation, PMCID/DOI, and which figure/table the data came from>

## The model
<what it is, where it came from (author file? reconstructed? imported?), size, key mechanisms>

## What is fit
<the experimental design: time course / dose-response / pre-equilibration; the observables and
 how they map to the paper's measured quantities>

## Free parameters (<k>)
<table or list: id, published value, search range, note>

## Optimizer
<job_type and WHY — which candidates were tried, what the budget is, what it costs>

## Verification
<tier-1, PEtab round-trip or BPSL check, the OG, the reproduction metric + tolerance;
 link to VALIDATION.md once validate-pybnf-job has run>

## Run
<the exact commands, from inside the slug folder>
```

Two rules that the corpus enforces and that reviews catch:

- **The README must agree with the conf and the shipped artifacts.** A README claiming `de` while
  the conf says `gntr`, or quoting a budget the conf does not carry, is a defect — the whole
  corpus was swept for exactly this in "Make every slug README agree with its own conf and
  artifacts".
- **Never let a declared number read as a measured one.** An OG evaluated at the nominal point
  carries a `†` and says so in words. See `og-acceptance.md` §4.

The paper-level `<FirstAuthor>-<Year>/README.md` above it carries the citation, the shared model,
and a one-row-per-slug table (what each fits · archetype · data source · J\* · OG · status).

---

## 4. Appendix: promoting a slug to PyBNF's `examples/real-world/`

Optional, and only for a slug that is small, fast, and pedagogically clean. This is a PR against
`lanl/PyBNF`, not part of landing the job here. Nothing below governs `pybnf-jobs/`.

### The `_manifest.py` entry

Each promoted example registers one frozen `RealWorldExample` in `EXAMPLES`
(`_manifest.py:35-94`):

| field | type | meaning | required |
|---|---|---|---|
| `folder` | str | subfolder name; also the pytest id and `example_by_folder` key | yes |
| `conf` | str | conf filename inside the folder | yes |
| `simulator` | `'ode'`/`'ssa'`/`'nf'` | the method the `experiment:` synthesizes | yes |
| `observables` | tuple | data-bound observable/function names (the `.exp` columns; functions without parens) | yes |
| `system` | str | one-line biology + **paper mapping** (cite the paper here) | yes |
| `stochastic` | bool (False) | `True` iff simulator is `ssa`/`nf`; **cross-checked in tier-1** | default |
| `heavy` | bool (False) | cluster-scale build/fit; **excludes from the executable tier-2 set** | default |
| `blocked` | str (`''`) | non-empty ⇒ can't complete through bngsim, with reason | default |
| `recover` | dict (`{}`) | optional `{param: truth}` for parameter-recovery assertion | default |
| `tol` | float (0.5) | relative tolerance for the `recover` check | default |

Minimal ODE, experimental-data entry:
```python
RealWorldExample(
    folder='<name>', conf='<name>.conf', simulator='ode',
    observables=('<obs1>', '<obs2>'),
    system='<biology> (<First-author Year>, PMCID); ODE, <protocol>'),
```
For SSA/NF add `stochastic=True`; add `heavy=True` if a single build/fit is
cluster-scale. Entries are grouped by simulator with section comments
(`_manifest.py:54,68,75`).

**`recover` for paper-derived (real-data) examples:** there's no synthetic truth, but
you *can* seed `recover` with the paper's reported best-fit values and a loose `tol`
to assert the fit lands in the right ballpark — or leave it `{}` and document the
comparison in the PR. Prefer populating it when the paper reports point estimates.

---

### What the two test tiers require (`tests/test_real_world_examples.py`)

**Tier 1 — backend-free, default CI** (`test_real_world_conf_is_wellformed`,
`:85-108`). Needs only `BNGPATH` set (to locate BNG2.pl for model parsing), no
simulation. Asserts: conf parses; `edition == 2`; `job_type` resolves in
`FIT_TYPE_REGISTRY`; `conf.exp_data` non-empty; ≥1 free param, none with `__FREE`;
resolved `model.stochastic` == manifest `stochastic`. **`scripts/check_conf.py`
reproduces this tier locally.**

NF-only backend-free guards (`:111-177`) apply if `simulator == 'nf'`: synthesis must be
network-free (no `resetConcentrations`, `generates_network` False); NF
pre-equilibration must carry `equil_t_end:`; `gml:`/`complex:` must ride into the
actions. A new NF example is covered automatically by the parametrized test.

**Tier 2 — opt-in `recovery`** (`test_real_world_runs_through_bngsim`, `:197-229`).
Runs only for **non-`heavy`** examples; needs BNG2.pl + bngsim. With a short bounded
fit (`max_iterations=2, population_size=6`) it asserts models build, and
`trajectory.best_score()` is finite (the whole simulate→score→propose loop ran). If
`recover` is set, each recovered param is within `tol` of truth.

To land in the ✅ executable tier: keep `heavy=False`, ensure the model builds fast
through BNG2.pl, and ensure a finite objective. If cluster-scale, set `heavy=True`;
it then stays backend-free-only (🔶).

---

### README + manifest + test update steps

Promoting an example touches the copied folder plus:

1. **`_manifest.py`** — one `RealWorldExample(...)` in the right simulator group. The
   test imports `EXAMPLES` dynamically, so both tiers pick it up automatically.
2. **`README.md`, two places:**
   - **Coverage matrix** (`README.md:43-51`): add a row
     `| [`<name>`](<name>/) | <paper mapping> | **<SIM>** | <features exercised> | <status> |`.
     Status: ✅ = validated end-to-end (runs in `recovery` tier), 🔶 = builds/runs but
     too heavy for routine CI.
   - Only if the example is cluster-scale, add a **Known limitations** bullet
     (`README.md:85-94`).
3. **Compliance test assertion (this skill's completion bar):**
   - *Quantitative example:* assert the new example round-trips through PEtab v2 (export →
     `petab.v2` lint clean → import). The real-world manifest does not yet assert this; add
     a parametrized test (or extend the existing module) that runs
     `pybnf.petab.export_job` + `lint_problem`, mirroring `scripts/petab_roundtrip.py`. See
     `petab-compliance.md`.
   - *BPSL example:* PEtab export is (correctly) refused, so assert instead that
     `export_job` raises `NotImplementedError` (a guard the constraint stays native-only)
     and/or that a `job_type = check` run reports the expected satisfaction. See
     `bpsl-constraints.md`.
