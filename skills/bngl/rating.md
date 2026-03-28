# BNGL Model Rating System

## 0. Purpose

A numeric rating system for individual `.bngl` files in the model library. Ratings communicate trust, annotation depth, and formatting compliance in a single number. Ratings are assigned per-file, not per-folder.

The system is designed to be **fully automatable** through level 5.4. A grading script reads the `.bngl` file, its folder's `metadata.yaml`, and any reference data, then computes the rating.

---

## 1. Rating Scale Overview

| Rating | Name | What it means |
|---|---|---|
| 1 | Registered | File exists with valid metadata and a contactable point of contact |
| 2 | Runnable | Parses, executes, produces output |
| 3 | Reproducible | Output matches committed reference data |
| 4.x | Annotation-compliant | Meets declared annotation standard (see §3) |
| 5.x | Fully compliant | Meets declared annotation + formatting standard (see §4) |
| 5.5+ | Bonus features | Comprehensive + specific enrichments (see §5) |

The decimal at levels 4 and 5 encodes the **documentation level**:

| Decimal | Documentation level |
|---|---|
| .0 | `none` |
| .1 | `minimal` |
| .2 | `standard` |
| .3 | `extra` |
| .4 | `comprehensive` |

---

## 2. Trust Ladder (Ratings 1–3)

These levels are independent of documentation level. Every file climbs the same ladder.

### Rating 1 — Registered

All of the following:
- A `.bngl` file is present.
- A `metadata.yaml` exists in the same folder with all required fields passing validation:
  - `id` matches the folder name.
  - `created` is a valid date, not in the future.
  - `point_of_contact` includes `name` and a valid `email`.
  - `source.tags` is a non-empty list using vocabulary from `skill.md` §5.2.
  - `files` manifest lists the `.bngl` file with a `role`.
  - Exactly one file has `role: primary`.
- The point of contact's email has been verified (confirmation sent and acknowledged).

### Rating 2 — Runnable

Rating 1, plus all of the following:
- The file parses as legal BNGL with `begin model` / `end model`.
- All provided actions (uncommented) execute without error.
- At least one simulation completes successfully.
- At least one output file (`.gdat` or `.scan`) is generated and non-empty.

### Rating 3 — Reproducible

Rating 2, plus all of the following:
- Reference output files are committed in the `reference/` subdirectory of the model folder for all uncommented simulation protocols:
  - `.gdat` and/or `.scan` files from `simulate` / `parameter_scan` commands.
  - `.net` file if any generate-first simulation is present.
  - `.xml` file if any network-free simulation is present.
  - Scan subdirectories (e.g., `*_scan/`) containing per-point `.gdat`/`.cdat`/`.net` files.
- The `reference/` directory uses flat organization — filenames encode the parent `.bngl` file via the naming convention.
- Re-running the model produces output that matches the reference data within tolerance.

**Comparison rules:**
- `.net` and `.xml`: exact structural match (ignoring comment lines and whitespace).
- `.gdat` and `.scan`: numeric comparison passes if, for each value, either:
  - relative error <= 1e-3, OR
  - absolute error <= `atol`, where `atol = 1e-6 * max(|column|)`.

  The absolute-error fallback handles near-zero trajectories where relative error is undefined or inflated.

Reference data for commented-out protocols is allowed and encouraged but not required.

---

## 3. Annotation Compliance (Rating 4.x)

A file reaches rating 4.x when it satisfies rating 3 **and** meets the annotation requirements for its documentation level. The documentation level is declared in `metadata.yaml` via the `documentation_target` field on the file entry (see `skill.md` §5.1).

**Auto-bump rule:** The grader computes the actual annotation level based on file contents. If the detected level exceeds the declared level, the rating is bumped to match reality. The declared level is a floor, not a ceiling.

### 4.0 — Annotation: `none`

No annotation requirements beyond what rating 3 already provides. Intended for test files and feature demos.

### 4.1 — Annotation: `minimal`

4.0, plus:
- `#@title:` tag present.
- `#@description:` tag present.

### 4.2 — Annotation: `standard`

4.1, plus:
- `#@keyword:` tag present.
- `#@reference:` tag present.
- Units comments on all parameters (same-line `# <unit>`).
- Units comments on all seed species (same-line `# <unit>`).
- `begin molecule types` block present with all molecules typed.
- Comment-to-code ratio >= T1 (threshold TBD).
- If `source.tags` includes `literature`: at least one `#@figure:` tag linking a simulation protocol to a paper figure, and a verification notebook (`.ipynb`) in the folder.

### 4.3 — Annotation: `extra`

4.2, plus:
- `#@protocol:` tags on simulation actions.
- `#@id:` tag on each molecule type.
- `#@id:` tag on components and internal states.
- Comment-to-code ratio >= T2 (where T2 > T1; threshold TBD).
- At least one `#@cite:` tag in the body of the file, tying a specific model element (parameter, rule, etc.) to a source in `#@reference:`.

### 4.4 — Annotation: `comprehensive`

4.3, plus:
- Comment-to-code ratio >= T3 (where T3 > T2; threshold TBD).
- Entity grounding: molecule type names use standard nomenclature (HGNC for human, MGI for mouse, etc., per `skill.md` §4.2).

---

## 4. Formatting Compliance (Rating 5.x)

A file reaches rating 5.x when it satisfies rating 4.x **and** meets the formatting requirements for its documentation level. The decimal carries forward from the annotation level.

### 5.0 — Formatting: `none`

No formatting requirements. The file reached 4.0 and there is nothing more to check.

### 5.1 — Formatting: `minimal`

4.1, plus:
- Correct indentation: 2-space indent inside `begin`/`end` blocks, no tabs.

### 5.2 — Formatting: `standard`

4.2, plus:
- Full compliance with formatting conventions in `skill.md`:
  - Correct block order (§2).
  - 2-space indentation, no tabs (§3.1).
  - Maximum 100-character line length (§3.1).
  - Comment placement rules (§3.3).
  - Entity-level annotation attachment (§3.4).
  - Naming conventions (§4).
  - Suffix tags on all simulate commands (§1.1).
  - State management (`saveConcentrations`/`resetConcentrations`) when multiple simulations are present (§1.1).
  - Preferred-method convention (one uncommented, rest commented) (§1.1).

### 5.3 — Formatting: `extra`

4.3, plus:
- Same full formatting compliance as 5.2 (all `skill.md` formatting rules).

### 5.4 — Formatting: `comprehensive`

4.4, plus:
- Same full formatting compliance as 5.2 (all `skill.md` formatting rules).

---

## 5. Bonus Features (Rating 5.5+)

Ratings 5.5 and above are reserved for `comprehensive`-level models that include specific enrichment features beyond full compliance. These will be defined as the library matures. Candidate features include:

- Variant files exploring alternate parameterizations or mechanisms.
- Behavioral demonstration actions (bistability switching, parameter scans, etc.).
- Multiple verified simulation modes (ODE + SSA + NFsim, each with reference data).

Each bonus feature, once defined, will be assigned a specific 5.x number.

---

## 6. Comment-to-Code Ratio

The comment-to-code ratio is computed as:

```
ratio = (number of comment lines) / (number of BNGL code lines)
```

**Counting rules:**
- **Comment lines:** Count each `#` sign that begins a comment (after optional leading whitespace). Multi-line continuations (`#  ...`) each count as one comment line.
- **BNGL code lines:** Count non-blank, non-comment lines inside `begin model` / `end model` and any actions block. Line continuations (`\` at end of line) count the continuation group as one logical line, not multiple.
- Blank lines are excluded from both counts.

**Thresholds (TBD — to be calibrated against existing models):**
- T1 (`standard`): TBD
- T2 (`extra`): TBD, where T2 > T1
- T3 (`comprehensive`): TBD, where T3 > T2

---

## 7. Documentation Level Declaration

The documentation level is declared per-file in `metadata.yaml`:

```yaml
files:
  - name: my_model.bngl
    role: primary
    documentation_target: standard   # none | minimal | standard | extra | comprehensive
```

If `documentation_target` is omitted, the grader detects the level from file contents (auto-detection). The auto-bump rule (§3) ensures that models are never penalized for exceeding their declared level.

---

## 8. Rating Display

Ratings are displayed as a single number:
- `1`, `2`, `3` for the trust ladder.
- `4.0`–`4.4` for annotation compliance.
- `5.0`–`5.4` for full compliance.
- `5.5+` for bonus features.

When a model has passed a level but has scattered features from higher levels without meeting the next full level, the rating is **not** inflated. The rating reflects the highest *fully satisfied* level.

The rating for each `.bngl` file is recorded in `metadata.yaml`:

```yaml
files:
  - name: my_model.bngl
    role: primary
    documentation_target: standard
    rating: 5.2
```

The `rating` field is computed by the grader and should not be manually edited.
