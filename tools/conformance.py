#!/usr/bin/env python3
"""Conformance validator for the curated collection in `models/`.

This is a **regression guard**, not a cleanup tool. Every ERROR rule below passes over the
whole collection today; the job of the validator is to keep it that way, so that a folder
added or edited six months from now cannot quietly drop a required tag, invent a source tag,
or nest its `reference/` output. It is deliberately structural — it reads `metadata.yaml`,
the `.bngl` header comments, the `README.md` Models table, and the shape of `reference/`.
It never runs BioNetGen and never judges whether a model is *correct*; that is what each
folder's `verify_<author><year>` artifact is for.

Rules are transcribed from the house style, and each finding names its section:

  `skills/bngl/skill.md`  §1.3 item 4  simulation intent declared in the header
                          §4.8         primary `.bngl` stem matches the folder name
                          §5.2         source-tag vocabulary
                          §5.3         file-role vocabulary
                          §5.4         one `metadata.yaml` per folder; `reference/` is flat
                          §5.6         `scale:` vocabulary, required on every `.bngl` entry
                          §6.1.1       the four required header keys
                          §7.1         `#@runtime_expectation:` when `scale:` is heavy
                          §9.4a        folder-metadata lint rules
  `skills/curate-model/SKILL.md` step 9 the folder has a row in the README Models table

Severity follows §9: ERROR fails CI, WARN is reported and does not. Known pre-existing ERROR
findings are grandfathered in `conformance_baseline.yaml` — see `apply_baseline` for why a
*stale* baseline entry is itself an error.

Usage:

    uv run python tools/conformance.py            # errors only; exit 1 if any
    uv run python tools/conformance.py -w         # errors and warnings
    uv run python tools/conformance.py --json     # machine-readable
    uv run python tools/conformance.py --write-baseline   # accept current ERRORs (deliberate)
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
BASELINE_PATH = Path(__file__).resolve().parent / "conformance_baseline.yaml"

ERROR = "error"
WARN = "warn"

# skills/bngl/skill.md §5.2. The table is open ("New tags MAY be added as needed"), so adding
# a tag here and to §5.2 in the same commit is the supported way to extend it.
SOURCE_TAGS = frozenset(
    {
        "literature",
        "new_model",
        "classic",
        "originally_rule_based",
        "reformulation",
        "extension",
        "bug_fix",
        "refinement",
        "restructuration",
        "translation",
        "formulated_for_testing_purposes",
        "ai_assisted",
        "adapted_from_ode",
        "adapted_from_sbml",
        "adapted_from_antimony",
        "adapted_from_bngl",
        "adapted_from_python",
        "biomodels_database",
        "benchmark_collection",
        "network_free",
    }
)

# skills/bngl/skill.md §5.3.
FILE_ROLES = frozenset(
    {"primary", "variant", "related", "verification", "reference", "output", "figure", "tooling"}
)

# skills/bngl/skill.md §5.6, in increasing cost. Ordering matters: a folder's README scale is
# the maximum over its files.
SCALES = ("trivial", "minutes", "hours", "cluster")

REQUIRED_METADATA_KEYS = ("id", "created", "point_of_contact", "source", "files")

# skills/bngl/skill.md §6.1.1 minimal mode.
REQUIRED_HEADER_TAGS = ("#@title:", "#@description:", "#@keyword:", "#@reference:")

# skills/bngl/skill.md §1.3 item 4 — the header must say whether the model is read as molecule
# counts or as concentrations converted through NA and a volume. There is no dedicated tag for
# this (§6.5 keeps the required header set to four keys), so it is declared in prose inside
# `#@description:` or `#@note:` and detected by the vocabulary the collection actually uses.
# Matched against the header comment block only, so a `# dimensionless` units comment on some
# parameter three hundred lines down cannot satisfy it.
SIMULATION_INTENT = re.compile(
    r"population[-\s]based|concentration[-\s]based|molecule[-\s]counts?|"
    r"particle[-\s]counts?|copy numbers?|dimensionless",
    re.IGNORECASE,
)

# A committed BioNetGen `parameter_scan()` output directory. BNG names it `<prefix>_scan`, and
# several models suffix the scanned parameter onto that (`..._scan_pp1`,
# `..._phosphatase_scan_ABC`).
SCAN_DIR = re.compile(r"_scan(_|$)")

# `| `folder/`<br>`file.bngl` | trivial | description | refs |` in the README Models table.
README_ROW = re.compile(
    r"^\|\s*`(?P<folder>[a-z0-9_]+)/`.*?\|\s*(?P<scale>trivial|minutes|hours|cluster)\s*\|"
)

# Deliverables a curated folder is expected to declare in its manifest (§9.4a WARN): the models,
# the verification artifact and its figures, and the helper scripts curate-model asks to be
# committed beside them (`independent_`, `digitize_`, `extract_`, `run_`, `generator/build_`).
MANIFEST_SUFFIXES = (".bngl", ".ipynb", ".png", ".pdf", ".svg", ".py")


@dataclass(frozen=True, order=True)
class Finding:
    """One rule violation, anchored to a repo-relative path."""

    rule: str
    severity: str
    where: str
    message: str

    @property
    def key(self) -> str:
        """Stable identity used to match against the baseline."""
        return f"{self.rule} {self.where}"

    def __str__(self) -> str:
        return f"{self.where}: {self.message}"


def header_text(text: str) -> str:
    """The leading comment block of a `.bngl` file.

    Everything from the top of the file down to the first line that is neither blank, nor a
    comment, nor `begin model`. This covers both shapes in the collection: files that open with
    `begin model` and comment underneath, and files that carry a prose banner above it.
    """
    out: list[str] = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("begin model"):
            out.append(line)
            continue
        break
    return "\n".join(out)


def model_bngl_files(folder: Path) -> list[Path]:
    """Every `.bngl` in a model folder that is a model, not committed reference output."""
    return sorted(
        p for p in folder.rglob("*.bngl") if "reference" not in p.relative_to(folder).parts
    )


def reaction_count(net_path: Path) -> int | None:
    """Number of reactions in a committed `.net`, or None if it has no reactions block."""
    try:
        lines = net_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return None
    count, inside = 0, False
    for line in lines:
        stripped = line.strip()
        if stripped == "begin reactions":
            inside = True
        elif stripped == "end reactions":
            return count
        elif inside and stripped and not stripped.startswith("#"):
            count += 1
    return count if inside else None


def is_network_free(text: str) -> bool:
    """True when the file has an *active* network-free protocol.

    §5.6.1 step 3 floors those at `minutes`, never `trivial`."""
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if 'method=>"nf"' in stripped.replace(" ", "") or "method=>'nf'" in stripped.replace(
            " ", ""
        ):
            return True
    return False


def parse_readme(readme: str) -> dict[str, str]:
    """Map each folder named in the README Models table to its declared Scale."""
    rows: dict[str, str] = {}
    for line in readme.splitlines():
        match = README_ROW.match(line)
        if match:
            rows[match.group("folder")] = match.group("scale")
    return rows


def check_folder(
    folder: Path, readme_rows: dict[str, str], repo_root: Path = REPO_ROOT
) -> list[Finding]:
    """Every rule that applies to one `models/<name>/` folder.

    Paths in findings are relative to `repo_root`, which is a parameter rather than the module
    global so the drift fixtures in tests/test_conformance.py can run the rules over a
    synthetic one-folder collection in a tmpdir."""
    findings: list[Finding] = []
    name = folder.name
    rel = folder.relative_to(repo_root).as_posix()

    def add(rule: str, severity: str, where: str, message: str) -> None:
        findings.append(Finding(rule, severity, where, message))

    # --- metadata.yaml exists and parses (§9.4a) ------------------------------------------
    metadata_path = folder / "metadata.yaml"
    if not metadata_path.exists():
        add("metadata-present", ERROR, rel, "no metadata.yaml (skills/bngl/skill.md §9.4a)")
        return findings
    metadata_rel = metadata_path.relative_to(repo_root).as_posix()
    try:
        metadata = yaml.safe_load(metadata_path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        add("metadata-present", ERROR, metadata_rel, f"not valid YAML: {exc}")
        return findings
    if not isinstance(metadata, dict):
        add("metadata-present", ERROR, metadata_rel, "top level is not a mapping (§5.1)")
        return findings

    # --- required keys, id, point of contact, source tags (§9.4a) --------------------------
    for key in REQUIRED_METADATA_KEYS:
        if key not in metadata:
            add("metadata-schema", ERROR, metadata_rel, f"missing required key `{key}` (§9.4a)")

    if metadata.get("id") != name:
        add(
            "metadata-id",
            ERROR,
            metadata_rel,
            f"id is `{metadata.get('id')}`, but the folder is `{name}` (§5.4 rule 1)",
        )

    if "rating" in metadata:
        add(
            "no-rating",
            ERROR,
            metadata_rel,
            "folder-level `rating:` was removed 2026-08-04 (§9.4a)",
        )

    poc = metadata.get("point_of_contact")
    if not isinstance(poc, dict):
        add(
            "point-of-contact",
            ERROR,
            metadata_rel,
            "point_of_contact is not a mapping (§5.4 rule 2)",
        )
    else:
        for field in ("name", "email"):
            if not poc.get(field):
                add(
                    "point-of-contact",
                    ERROR,
                    metadata_rel,
                    f"point_of_contact is missing `{field}` (§5.4 rule 2)",
                )

    source = metadata.get("source")
    tags = source.get("tags") if isinstance(source, dict) else None
    if not tags or not isinstance(tags, list):
        add("source-tags", ERROR, metadata_rel, "source.tags is missing or empty (§9.4a)")
    else:
        for tag in tags:
            if tag not in SOURCE_TAGS:
                add(
                    "source-tags",
                    ERROR,
                    metadata_rel,
                    f"source tag `{tag}` is not in the §5.2 vocabulary — add it to §5.2 and to "
                    "SOURCE_TAGS in this file if it is a real new tag",
                )

    # --- the files manifest (§5.3, §5.6, §9.4a) --------------------------------------------
    entries = metadata.get("files")
    if not isinstance(entries, list):
        add("metadata-schema", ERROR, metadata_rel, "files is not a list (§5.1)")
        entries = []

    primaries: list[str] = []
    for entry in entries:
        if not isinstance(entry, dict):
            add("metadata-schema", ERROR, metadata_rel, f"files entry is not a mapping: {entry!r}")
            continue
        entry_name = str(entry.get("name", ""))
        where = f"{metadata_rel}:{entry_name}" if entry_name else metadata_rel

        role = entry.get("role")
        if role not in FILE_ROLES:
            add("file-roles", ERROR, where, f"role `{role}` is not in the §5.3 vocabulary")
        if role == "primary":
            primaries.append(entry_name)

        if "rating" in entry:
            add(
                "no-rating",
                ERROR,
                where,
                "`rating:` was removed 2026-08-04; use documentation_target (§9.4a)",
            )

        # Every manifest entry must resolve. A trailing "/" marks a directory entry (§5.4 rule 4
        # lets a scan directory be listed as one entry rather than one per point).
        if entry_name and not (folder / entry_name.rstrip("/")).exists():
            add("manifest-paths", ERROR, where, "listed in the manifest but not present on disk")

        if entry_name.endswith(".bngl"):
            scale = entry.get("scale")
            if scale is None:
                add(
                    "bngl-scale",
                    ERROR,
                    where,
                    "`.bngl` entry does not declare `scale:` (§5.6, §9.4a)",
                )
            elif scale not in SCALES:
                add("bngl-scale", ERROR, where, f"scale `{scale}` is not in the §5.6 vocabulary")
            if not entry.get("documentation_target"):
                add(
                    "documentation-target",
                    WARN,
                    where,
                    "`.bngl` entry has no `documentation_target` (§9.4a)",
                )

            bngl_path = folder / entry_name
            if bngl_path.exists():
                text = bngl_path.read_text(encoding="utf-8", errors="replace")
                if scale in ("hours", "cluster") and "#@runtime_expectation:" not in text:
                    add(
                        "runtime-expectation",
                        WARN,
                        where,
                        f"scale `{scale}` with no `#@runtime_expectation:` in the .bngl (§7.1)",
                    )
                if scale == "trivial":
                    stem = Path(entry_name).stem
                    net = folder / "reference" / f"{stem}.net"
                    reactions = reaction_count(net) if net.exists() else None
                    if reactions is not None and reactions >= 1000:
                        add(
                            "scale-evidence",
                            WARN,
                            where,
                            f"scale `trivial` but the committed .net has {reactions} "
                            "reactions (§5.6.1)",
                        )
                    if is_network_free(text):
                        add(
                            "scale-evidence",
                            WARN,
                            where,
                            "scale `trivial` on a network-free file; §5.6.1 step 3 floors "
                            "those at `minutes`",
                        )

    # §5.3: exactly one primary, and its stem matches the folder (§4.8).
    if len(primaries) != 1:
        add(
            "primary-file",
            ERROR,
            metadata_rel,
            f"{len(primaries)} entries have `role: primary`; §5.3 requires exactly one",
        )
    elif Path(primaries[0]).stem != name:
        add(
            "primary-file",
            ERROR,
            metadata_rel,
            f"primary is `{primaries[0]}`, whose stem does not match the folder name (§4.8)",
        )

    # §9.4a WARN: deliverables present on disk but absent from the manifest. A directory entry
    # covers everything under it — §5.4 rule 4 lets a scan directory be one entry, and
    # lambda_switch_arkin1998 lists its whole `generator/` the same way.
    listed = {str(e.get("name", "")).rstrip("/") for e in entries if isinstance(e, dict)}
    listed_dirs = {n for n in listed if n and (folder / n).is_dir()}
    for path in sorted(folder.rglob("*")):
        if not path.is_file() or path.suffix not in MANIFEST_SUFFIXES:
            continue
        relative = path.relative_to(folder).as_posix()
        if relative.startswith("reference/"):
            continue  # reference output is covered by its own directory entries
        if any(relative.startswith(f"{d}/") for d in listed_dirs):
            continue
        if relative not in listed:
            add(
                "manifest-coverage",
                WARN,
                f"{rel}/{relative}",
                "present in the folder but not listed in the metadata.yaml files manifest (§9.4a)",
            )

    # --- .bngl header contents (§6.1.1, §1.3 item 4) ---------------------------------------
    for bngl in model_bngl_files(folder):
        bngl_rel = bngl.relative_to(repo_root).as_posix()
        text = bngl.read_text(encoding="utf-8", errors="replace")
        missing = [tag for tag in REQUIRED_HEADER_TAGS if tag not in text]
        if missing:
            add(
                "header-tags",
                ERROR,
                bngl_rel,
                f"header is missing {', '.join(f'`{t}`' for t in missing)} (§6.1.1, §9.4)",
            )
        if not SIMULATION_INTENT.search(header_text(text)):
            add(
                "simulation-intent",
                ERROR,
                bngl_rel,
                "header does not declare the simulation intent — population-based (molecule "
                "counts) or concentration-based with conversion through NA and a volume "
                "(§1.3 item 4)",
            )

    # --- the verification figure (curate-model completion criteria) ------------------------
    if not list(folder.glob("verify_*.png")):
        add(
            "verification-figure",
            ERROR,
            rel,
            "no `verify_*.png`; curate-model requires the figure in both artifact shapes",
        )

    # --- reference/ layout (§5.4 rule 4) ----------------------------------------------------
    reference = folder / "reference"
    if not reference.is_dir():
        add("reference-layout", ERROR, rel, "no `reference/` directory (§5.4 rule 4)")
    else:
        for sub in sorted(p for p in reference.iterdir() if p.is_dir()):
            sub_rel = sub.relative_to(repo_root).as_posix()
            if not SCAN_DIR.search(sub.name):
                add(
                    "reference-layout",
                    ERROR,
                    sub_rel,
                    "`reference/` is flat apart from BioNetGen scan output directories; filenames "
                    "encode the parent .bngl rather than a directory doing it (§5.4 rule 4)",
                )
            elif any(child.is_dir() for child in sub.iterdir()):
                add(
                    "reference-layout",
                    ERROR,
                    sub_rel,
                    "a scan output directory holds per-point files, not further "
                    "directories (§5.4 rule 4)",
                )

    # --- the README Models table (curate-model workflow step 9) -----------------------------
    # The README row is the maximum over *every* scaled entry, not only the `.bngl` ones: §5.6
    # allows `scale:` on a driver script or notebook, and four folders (arkin1998, chylek2014,
    # creamer2012, dolan2015) are `hours` because of the campaign rather than the model.
    declared = [e.get("scale") for e in entries if isinstance(e, dict) and e.get("scale") in SCALES]
    if name not in readme_rows:
        add("readme-row", ERROR, "README.md", f"`{name}/` has no row in the Models table")
    elif declared:
        expected = max(declared, key=SCALES.index)
        if readme_rows[name] != expected:
            add(
                "readme-scale",
                ERROR,
                "README.md",
                f"`{name}/` is listed as `{readme_rows[name]}` but the maximum `scale:` over its "
                f".bngl entries is `{expected}` (§5.6)",
            )

    return findings


def run(repo_root: Path = REPO_ROOT) -> list[Finding]:
    """Every finding over the whole collection, sorted."""
    models = repo_root / "models"
    readme_rows = parse_readme((repo_root / "README.md").read_text(encoding="utf-8"))
    findings: list[Finding] = []
    for folder in sorted(p for p in models.iterdir() if p.is_dir()):
        findings.extend(check_folder(folder, readme_rows, repo_root))
    return sorted(findings)


def load_baseline(path: Path = BASELINE_PATH) -> set[str]:
    """The grandfathered ERROR keys, as `"<rule> <path>"` strings."""
    if not path.exists():
        return set()
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    keys: set[str] = set()
    for rule, wheres in data.items():
        for where in wheres or []:
            keys.add(f"{rule} {where}")
    return keys


def apply_baseline(
    findings: list[Finding], baseline: set[str]
) -> tuple[list[Finding], list[Finding], list[str]]:
    """Split findings against the baseline.

    Returns `(reported, suppressed, stale)`. A *stale* key — one in the baseline that no
    longer fires — is itself reported as a failure. That is what keeps the baseline shrinking:
    fixing a folder without deleting its baseline line breaks the build, so the file cannot
    quietly outlive the debt it records.
    """
    live = {f.key for f in findings}
    reported = [f for f in findings if f.severity != ERROR or f.key not in baseline]
    suppressed = [f for f in findings if f.severity == ERROR and f.key in baseline]
    stale = sorted(baseline - live)
    return reported, suppressed, stale


def format_report(findings: list[Finding], *, show_warnings: bool) -> str:
    """Findings grouped by rule, most-populous first."""
    wanted = [f for f in findings if show_warnings or f.severity == ERROR]
    if not wanted:
        return ""
    by_rule: dict[tuple[str, str], list[Finding]] = defaultdict(list)
    for finding in wanted:
        by_rule[(finding.severity, finding.rule)].append(finding)

    lines: list[str] = []
    for (severity, rule), group in sorted(
        by_rule.items(), key=lambda kv: (kv[0][0] != ERROR, -len(kv[1]), kv[0][1])
    ):
        lines.append(f"{severity.upper()}  {rule}  ({len(group)})")
        for finding in group:
            lines.append(f"    {finding}")
        lines.append("")
    return "\n".join(lines).rstrip()


def write_baseline(findings: list[Finding], path: Path = BASELINE_PATH) -> int:
    """Rewrite the baseline to accept every current ERROR. Deliberate action, never automatic."""
    grouped: dict[str, list[str]] = defaultdict(list)
    for finding in findings:
        if finding.severity == ERROR:
            grouped[finding.rule].append(finding.where)

    header = (
        "# Grandfathered conformance findings — see tools/conformance.py.\n"
        "#\n"
        "# Each entry is an ERROR that the collection already carried when the validator was\n"
        "# written. They are suppressed so the guard can run green and catch *new* drift; they\n"
        "# are not exemptions. Removing an entry is the point: fix the file, delete its line.\n"
        "# A line that no longer corresponds to a live finding fails the run, so this file\n"
        "# cannot outlive the debt it records.\n"
    )
    body = yaml.safe_dump(
        {rule: sorted(wheres) for rule, wheres in sorted(grouped.items())},
        default_flow_style=False,
        sort_keys=True,
        width=100,
    )
    path.write_text(header + "\n" + body, encoding="utf-8")
    return sum(len(v) for v in grouped.values())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("-w", "--warnings", action="store_true", help="report WARN findings too")
    parser.add_argument("--json", action="store_true", help="emit findings as JSON")
    parser.add_argument(
        "--baseline", type=Path, default=BASELINE_PATH, help="baseline file to apply"
    )
    parser.add_argument(
        "--no-baseline", action="store_true", help="report every finding, ungrandfathered"
    )
    parser.add_argument(
        "--write-baseline",
        action="store_true",
        help="rewrite the baseline to accept every current ERROR (deliberate; review the diff)",
    )
    args = parser.parse_args(argv)

    findings = run()

    if args.write_baseline:
        count = write_baseline(findings, args.baseline)
        print(f"wrote {args.baseline} with {count} grandfathered findings")
        return 0

    baseline = set() if args.no_baseline else load_baseline(args.baseline)
    reported, suppressed, stale = apply_baseline(findings, baseline)

    if args.json:
        print(
            json.dumps(
                {
                    "reported": [f.__dict__ for f in reported],
                    "suppressed": [f.__dict__ for f in suppressed],
                    "stale_baseline": stale,
                },
                indent=2,
            )
        )
    else:
        report = format_report(reported, show_warnings=args.warnings)
        if report:
            print(report)
            print()
        if stale:
            print(f"STALE BASELINE  ({len(stale)})")
            print("    These entries no longer fire. Delete them from the baseline.")
            for key in stale:
                print(f"    {key}")
            print()

        errors = sum(1 for f in reported if f.severity == ERROR)
        warnings = sum(1 for f in reported if f.severity == WARN)
        folders = sum(1 for p in (REPO_ROOT / "models").iterdir() if p.is_dir())
        summary = f"{folders} folders · {errors} errors · {warnings} warnings"
        if suppressed:
            summary += f" · {len(suppressed)} grandfathered"
        if not args.warnings and warnings:
            summary += "  (run with -w to see warnings)"
        print(summary)

    return 1 if any(f.severity == ERROR for f in reported) or stale else 0


if __name__ == "__main__":
    sys.exit(main())
