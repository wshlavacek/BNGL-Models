# Anatomy of a `examples/real-world/` job-setup folder

What a new example under `~/Code/PyBNF/examples/real-world/` must contain to (a) run
and (b) pass the two committed test tiers. Grounded in the PyBNF source (v1.6.0); every
schema claim cites `file:line`. The corpus is the 2019 PyBNF-paper case studies
re-expressed on the edition-2 surface; a new paper-derived example follows the same
mold.

A new example = **a folder of `{model, conf, exp}` + one `_manifest.py` entry + two
`README.md` edits** (+ a PEtab test assertion, per the skill's completion bar).

## Contents
1. File inventory
2. File formats (with the `receptor` template)
3. The `_manifest.py` entry (field-by-field)
4. What the two test tiers require
5. README + manifest + test update steps

---

## 1. File inventory

Each example lives in its own self-contained slug folder, **grouped by source paper under
a `<FirstAuthor>-<Year>/` directory** (e.g. `Rukhlenko-2022/cstar_trka/`; one paper's
several jobs sit side by side there — see the skill's "Name and create the folder" step).
All paths inside the `.conf` resolve relative to the slug folder (`README.md:98`; the test
`chdir`s in before parsing — `test_real_world_examples.py:57-60`), so the grouping directory
is transparent to the conf. Where a `_manifest.py`/README/test row names the folder, use the
paper-relative path `<FirstAuthor>-<Year>/<slug>`.

| file | role | required |
|---|---|---|
| `<name>.bngl` | the BioNetGen model — **no *simulation* actions** (`simulate`/`parameter_scan`; synthesized from the conf). **RETAIN network-generation directives** (`generate_network` with `max_stoich`/`max_agg`/`max_iter`) when the model's network is only finite/correct under them — see §2 | yes |
| `<name>.conf` | the edition-2 job setup (this *is* the test input) | yes |
| `<data>.exp` | ≥1 data file; each `experiment:` line binds one; column headers ARE model observable/function names | yes, unless the example is constraint-only |
| `<data>.prop` | BPSL constraint file(s) for a data-fusion or constraint-only example; attached on an `experiment:`'s `data:`. Makes the job native-only (not PEtab-exportable). See `bpsl-constraints.md` | only for BPSL examples |
| `*_ground_truth.bngl` | only for synthetic-data examples: the model at known-true params (documentation; not referenced by conf/tests) | no |

There is **no per-folder README**; documentation lives in `examples/real-world/README.md`
and metadata in `examples/real-world/_manifest.py`.

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

## 3. The `_manifest.py` entry

Each example registers one frozen `RealWorldExample` in `EXAMPLES`
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

## 4. What the two test tiers require (`tests/test_real_world_examples.py`)

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

## 5. README + manifest + test update steps

Adding an example touches the new folder plus:

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
