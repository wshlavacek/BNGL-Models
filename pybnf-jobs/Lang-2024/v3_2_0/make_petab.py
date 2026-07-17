#!/usr/bin/env python
"""Regenerate the committed PEtab v2 bundle for the Lang-2024 v3_2_0 job.

The runnable ``v3_2_0.conf`` is the single source of truth. Its PEtab.v2 compliance
is a *verified round-trip property* (see ``skills/curate-pybnf-job/references/
petab-compliance.md``): ``pybnf.petab.export_job`` transcribes the native edition-2
job into a PEtab v2 problem that lints clean under ``petab.v2`` and re-imports. This
script commits that exported problem as an illustrative ``petab/`` byproduct -- the
same pattern PyBNF's own ``examples/tutorial/12_petab_roundtrip/petab/`` uses -- so
consumers can pick up PEtab v2 directly without running PyBNF.

Do NOT hand-edit anything under ``petab/``: it is generated. Re-run this script after
any change to ``v3_2_0.conf`` / ``v3_2_0.bngl`` / ``v3_2_0.exp`` so the bundle can't
drift from the conf. For Lang-2024/v3_2_0 the export emits 8 observables
(``observableFormula`` = bare function name), 177 estimated parameters and 4800
measurements over the single wildtype condition (no conditions/experiments tables,
since the job does not perturb); the noise is a constant unit sigma, so the objective
is plain sum-of-squares.

Provenance note: the paper's own PEtab files (``paulflang/cell_cycle_petab`` v3.2.0)
are PEtab **v1** (SBML-based); this exported bundle is the first PEtab **v2** rendering
of the problem.

The output is normalized to the repo's file conventions (LF endings, no trailing
whitespace, single final newline) so it is byte-stable under pre-commit and a fresh
run reproduces the committed bytes exactly. Normalization only drops the empty
delimiter of trailing all-empty columns (``noisePlaceholders`` / ``noiseParameters``);
the bundle still lints clean and re-imports (verified).

Requires the PyBNF environment's interpreter (imports ``pybnf``; the lint self-check
also needs the optional ``pybnf[petab]`` extra) and BNGPATH set:
    BNGPATH=... ~/Code/PyBNF/.venv/bin/python make_petab.py
Exit 0 = exported + linted clean (or lint skipped because the petab extra is missing);
1 = a step failed.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _normalize(out_dir: Path) -> None:
    """Rewrite generated text files to the repo's conventions so the committed
    bundle is byte-stable under pre-commit and identical to a fresh regeneration:
    LF endings, no trailing whitespace, exactly one final newline."""
    for f in sorted(out_dir.iterdir()):
        if f.suffix not in (".tsv", ".yaml"):
            continue
        lines = f.read_text().replace("\r\n", "\n").replace("\r", "\n").split("\n")
        lines = [ln.rstrip() for ln in lines]
        while lines and lines[-1] == "":
            lines.pop()
        f.write_text("\n".join(lines) + "\n")


def _lint(problem_yaml: Path) -> tuple[str, str]:
    """Return (status, detail); status in {'clean', 'errors', 'skipped'}."""
    try:
        from petab.v2 import Problem
        from petab.v2.lint import lint_problem
        from pybnf.petab.bngl_model import register_bngl
    except Exception as exc:  # noqa: BLE001 -- petab extra not installed
        return "skipped", f"petab lint unavailable ({exc!r}); install pybnf[petab] to enforce it"

    register_bngl()
    report = lint_problem(Problem.from_yaml(str(problem_yaml)))
    has_errors = report.has_errors() if hasattr(report, "has_errors") else bool(report)
    return ("errors", str(report)) if has_errors else ("clean", "")


def _counts(out_dir: Path) -> str:
    """Human-readable summary of the exported tables (header row excluded)."""

    def _rows(name: str) -> list[str]:
        p = out_dir / name
        if not p.is_file():
            return []
        body = p.read_text().splitlines()
        return body[1:] if body else []

    n_obs = len(_rows("observables.tsv"))
    n_meas = len(_rows("measurements.tsv"))
    n_est = sum(1 for r in _rows("parameters.tsv") if r.split("\t")[1:2] == ["true"])
    parts = [f"{n_obs} observables", f"{n_est} estimated parameters", f"{n_meas} measurements"]
    if (out_dir / "conditions.tsv").is_file():
        parts.append(f"{len(_rows('conditions.tsv'))} conditions")
    return ", ".join(parts)


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--conf",
        default=str(HERE / f"{HERE.name}.conf"),
        help="edition-2 .conf to export (default: <slug>.conf next to this script)",
    )
    ap.add_argument(
        "--out", default=str(HERE / "petab"), help="output bundle directory (default: ./petab)"
    )
    args = ap.parse_args()

    from pybnf.petab import export_job

    conf = Path(args.conf).resolve()
    out_dir = Path(args.out).resolve()
    if not conf.is_file():
        print(f"FAIL  no such conf: {conf}")
        return 1

    # Regenerate from clean so a removed/renamed table can't linger.
    if out_dir.exists():
        shutil.rmtree(out_dir)

    # export_job resolves conf-relative paths (model, .exp) against the CWD.
    home = os.getcwd()
    os.chdir(conf.parent)
    try:
        try:
            export_job(str(conf), str(out_dir))
        except NotImplementedError as exc:
            print(f"FAIL  export: job uses a non-PEtab-exportable feature -- {exc}")
            print("      -> this slug is native-only; do not commit a petab/ bundle for it")
            return 1
    finally:
        os.chdir(home)

    problem_yaml = out_dir / "problem.yaml"
    if not problem_yaml.is_file():
        print(f"FAIL  export: no problem.yaml under {out_dir}")
        return 1

    _normalize(out_dir)

    wrote = sorted(p.name for p in out_dir.iterdir())
    print(
        f"ok    export -> {out_dir.relative_to(HERE) if out_dir.is_relative_to(HERE) else out_dir}"
    )
    print(f"        wrote: {', '.join(wrote)}")
    print(f"        {_counts(out_dir)}")

    status, detail = _lint(problem_yaml)
    if status == "errors":
        print(f"FAIL  lint: petab.v2 reported errors\n{detail}")
        return 1
    if status == "skipped":
        print(f"warn  lint: {detail}")
    else:
        print("ok    lint  -> petab.v2 clean (no errors)")

    print("\nDONE  committed petab/ bundle regenerated from v3_2_0.conf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
