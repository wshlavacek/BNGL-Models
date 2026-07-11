# BPSL constraint files (`.prop` / `.con`) in edition-2 job setups

PyBNF's **Biological Property Specification Language (BPSL)** lets a fit use
*qualitative* facts a paper reports — orderings, thresholds, steady-state levels,
monotonic dose-response, "peaks before", bistability — not just quantitative
time-series. BPSL is the signature capability of the 2019 Mitra et al. paper the
`examples/real-world/` corpus derives from, so a constraint-bearing job is a
*first-class* reason to add a real-world example, not an edge case. Grounded in the
PyBNF source (v1.6.0); cites `file:line`.

`.prop` and `.con` are **synonyms** — one grammar, one loader
(`config.py:1566,2078`; `export.py:479`). The repo ships only `.prop`; prefer it.

## The one consequence to decide up front

A job that attaches **any** `.prop`/`.con` to an experiment is **native-PyBNF-only —
not PEtab-v2-exportable**. `export_job` raises `NotImplementedError` (`export.py:479-487`,
verified): BPSL has no core-PEtab representation, and there is **no round-trip**. So
before authoring, pick the example's flavor:

| flavor | data | PEtab round-trip | how you verify |
|---|---|---|---|
| **quantitative** | `.exp` only | **required** (`petab_roundtrip.py`) | tier-1 + PEtab + bounded fit + paper repro |
| **data fusion** | `.exp` + `.prop` | **N/A** (native-only) | tier-1 + bounded fit + `check` satisfaction + paper repro |
| **constraint-only** | `.prop` only | **N/A** (native-only) | tier-1 + `check` satisfaction (+ bounded fit if params are free) |

For a native-only example, don't run (or promise) `petab_roundtrip.py`; verify the
constraints with `job_type = check` instead (§3), and mark it non-exportable in the
manifest/README.

---

## 1. BPSL grammar (`constraint.py:240-284`)

One statement per line; blanks and `#` comments skipped (`constraint.py:90`):

```
<inequality> <enforcement> [<penalty-clause>] [# comment]
```

**Inequality** (`constraint.py:249-251`): `obs OP (obs | number)` or `number OP obs`
— at least one side is an observable. `OP` ∈ `< <= > >=` (no `==`/`!=`).
- **Observable token** (`constraint.py:242`): a bare name (`Activity`, `pERK`) resolves
  to *this experiment's own* simulation output; a **dotted** `suffix.Observable`
  (`mut.y2`, `KO.MEK_pRDS_KO`) references another experiment/suffix's output — how you
  compare across conditions/mutants (`constraint.py:446-485`).

**Enforcement** (temporal quantifier; `constraint.py:255-260`):

| form | meaning |
|---|---|
| `always` | holds at every output point (worst case) |
| `once` | holds at ≥1 output point (best case) |
| `at <crit> [everytime\|first] [before]` | holds where a condition is met |
| `between <crit>,<crit>` | holds on the interval between two conditions |
| `once between <crit>,<crit>` | holds at ≥1 point in the interval |

`<crit>` = `number` or `obs = number` (`constraint.py:253-254`): `at 300` (a time on the
independent variable) or `at y1=400` / `at FLAG_SPC=1` (a variable crossing). `everytime`
re-applies at every crossing; `before` uses the point just before it.

**Split-at** (`constraint.py:261-262`): compare two observables each at its own point —
`obs1 at <crit> OP obs2 at <crit>` (e.g. encode monotonicity: `pERK at time=300 > pERK at time=1800`).

**Penalty clause** — pick one family per line (mutually exclusive, `constraint.py:363-399`):
1. **Weighted / hinge (default):** `weight <num> [altpenalty <ineq>] [min <num>]`
   (`constraint.py:263-266`). No clause ⇒ `weight 1` (`constraint.py:156-159`). Penalty =
   `max(0, margin)·weight`. This is what every shipped `.prop` uses.
2. **Probit likelihood:** `confidence <num> [tolerance <num>]` or `pmin <num> pmax <num>
   [tolerance <num>]` (`constraint.py:267-271`).
3. **Logit / softplus** (Miller 2025): `logit scale <num> [pmin <num> pmax <num>]`
   (`constraint.py:272-276`); `ln(1 + e^{margin/scale})`.

### Real example lines (verbatim)

`constraint_raf/raf.prop` — threshold at a time:
```
Activity>13.74 at 0.000100372016637506	weight 0.03
Activity<16.74 at 3.1126883457236 weight 0.03
```
`constraint_advanced/wt.prop` — `once`, `at var=val`, `between`, cross-suffix, `altpenalty`, `min`:
```
y1<0 once weight 40
y2<77 at y1=400 weight 40
y2<y1 between y2=-25,y2=0 weight 100
y2>y1 at time=7 weight 1 altpenalty mut.y2<mut.y1
y2>26.9 at 0 weight 20 min 20
y1<=mut.y1 always weight 50
```
`Miller2025_MEK_Isoforms/.../WT.prop` — split-at monotonicity + cross-mutant ordering:
```
WT.MEK_pRDS_WT at time=300>WT.MEK_pRDS_WT at time=1800 weight 0.001
WT.MEK_pRDS_WT at time=1800<KO.MEK_pRDS_KO at time=1800 weight 0.001
```

Every observable named in a `.prop` must exist in the model (a bare name in the
experiment's own model; a dotted name in the referenced suffix's model).

---

## 2. Attaching a `.prop` in an edition-2 conf (`config.py:1334-1421,1556-1642`)

Constraints ride an experiment's `data:` list; files split by extension — `.exp` =
measurements, `.prop`/`.con` = constraints (`config.py:1566`). A **bare** observable in
the `.prop` binds to that experiment's own output (model + condition inherited); no
extra syntax.

**Data fusion — `.exp` + `.prop` on one experiment** (the Miller2025 / constraint_raf
pattern; the `.exp` sets the grid + quantitative objective, the `.prop` adds a penalty
on the same run):
```conf
experiment: meas, condition: withligand, data: meas.exp, meas.prop
```

**Constraint-only experiment — `.prop` with no `.exp`** (the constraints alone define
the job; the experiment must carry its own timing, `config.py:1575-1618`):
```conf
experiment: qual, t_end: 10, n_steps: 5, data: qual.prop
```
- `t_end:` is **required** (integration endpoint; error if absent).
- `t_start:` optional (default 0, must be `< t_end`); `n_steps:` optional output
  resolution (default step 1). `type:` must be `time_course` (scan is refused).

Both are edition-gated; the field list is in the parser error text
(`parse.py:865-874`).

---

## 3. `job_type = check` — model checking, no fitting (`algorithms/model_check.py`)

`check` evaluates a **single** parameter set (empty search) and **reports satisfaction**
— it prints `Objective value is …`, `Satisfied N out of M constraints`, and writes a
per-line itemized `<suffix>_constraint_eval.txt` (`model_check.py:66-104`). Use it to
confirm a given parameter set (the paper's published values, or your fit result)
satisfies the qualitative properties.

During a **fitting** job (`de`/`ss`/`pso`/`am`/…), the *same* constraint sets become
**soft weighted penalties added to the objective** (`objective.py:212-213`), so the
optimizer minimizes qualitative violation jointly with quantitative misfit. Same
`.prop` files, two consumers.

**Verification recipe for a constraint-bearing example:** run the fit (penalties active),
then run a second job with `job_type = check` at the fitted (or published) parameters and
confirm `Satisfied M out of M` (or document which the paper itself does not require).

---

## 4. Constraint config keys (`config_schema.py:212-225`)

- **`constraint_scale`** (default `1.0`): global multiplier on every constraint's
  `weight` at load — up/down-weight the whole qualitative block vs. the quantitative
  objective without editing each line. Affects the hinge family only.
- **`qualitative_loss`** (`auto`/`hinge`/`probit`/`logit`, default `auto`): penalty
  **family** override; `auto` keeps each line's authored family. A benchmarking
  convenience — leave `auto` unless comparing losses.
- **`qualitative_scale = fit <param>`** (default none): tie every constraint's scale
  (logit `s` / probit `σ`) to a declared free parameter so the fit **estimates the
  qualitative scale jointly**. Logit/probit only (binding a hinge constraint errors).

There is no separate "hard constraint" key — emulate hardness with a large `weight`/`min`.

---

## 5. What to author, end to end

1. **Model** `model: <name>.bngl` with observables named for the constraints.
2. **`.prop`** one qualitative statement per line (§1), dotted `suffix.obs` for
   cross-experiment comparisons; `weight` per line, higher for stronger claims.
3. **Conf** attach via `experiment:` — fusion (`data: X.exp, X.prop`) or constraint-only
   (`data: X.prop, t_end: …, n_steps: …`); set `constraint_scale` to balance against the
   quantitative objective.
4. **Verify** with `check_conf.py` (tier-1, reports the constraint sets + native-only
   flag), a bounded fit (finite objective), and a `job_type = check` run (satisfaction).
   Do **not** run `petab_roundtrip.py` on a constraint-bearing job — it will (correctly)
   refuse.
5. **Register** as native-only: manifest `system` notes the BPSL constraints; README
   coverage row marks it not PEtab-exportable; **skip** the PEtab test assertion for this
   example (add, instead, a `check`-satisfaction or finite-objective assertion).
