#!/usr/bin/env python
"""Prove an edition-2 PyBNF job setup is PEtab.v2-compliant, by round-trip.

"PEtab.v2-compliant" is not a demand to ship PEtab TSV tables -- it is a *verified
property* of the native ``.conf`` (ADR-0004/0028): the job exports to a PEtab v2
problem, that problem lints clean under ``petab``'s own validator, and re-importing
it reproduces the same runnable job. This script runs all three steps and reports:

  1. EXPORT   ``pybnf.petab.export_job(conf, out/)`` writes problem.yaml + the
              parameters / observables / measurements tables (+ conditions /
              experiments when the job perturbs).
  2. LINT     ``petab.v2.lint.lint_problem(Problem.from_yaml(...))`` after
              ``register_bngl()`` -- must report no errors.
  3. IMPORT   ``pybnf.petab.import_job(problem.yaml, back/, job_type=...)``
              reconstructs a runnable edition-2 conf (the run recipe -- algorithm,
              per-experiment method -- is supplied on import, not recovered).

If EXPORT raises ``NotImplementedError`` the job uses a feature outside the
PEtab-exportable subset (``normalization``, ``cumulative``, ``neg_bin`` /
``lognormal`` noise, ``.prop`` constraints, legacy data linkage, Antimony models --
see references/petab-compliance.md). Rework the conf to stay inside the subset, or
mark the example non-exportable in the manifest and README.

Imports ``pybnf`` (and, for lint, the optional ``petab`` extra), so run it with the
PyBNF environment's interpreter:

    ~/Code/PyBNF/.venv/bin/python \
        skills/curate-pybnf-job/scripts/petab_roundtrip.py \
        ~/Code/PyBNF/examples/real-world/<name>/<name>.conf --job-type de

Exit code 0 = exported + linted clean + re-imported; 1 = a step failed;
3 = exported + imported but lint was skipped (petab extra not installed).
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path


def _run_from(folder: Path):
    """Context helper: PyBNF resolves conf-relative paths against the CWD."""
    class _Chdir:
        def __enter__(self):
            self._home = os.getcwd()
            os.chdir(folder)

        def __exit__(self, *exc):
            os.chdir(self._home)

    return _Chdir()


def _lint(problem_yaml: Path) -> tuple[str, str]:
    """Return (status, detail). status in {'clean', 'errors', 'skipped'}."""
    try:
        from pybnf.petab.bngl_model import register_bngl
        from petab.v2 import Problem
        from petab.v2.lint import lint_problem
    except Exception as exc:  # noqa: BLE001 -- petab extra not installed / import issue
        return "skipped", f"petab lint unavailable ({exc!r}); install pybnf[petab] to enforce it"

    register_bngl()
    report = lint_problem(Problem.from_yaml(str(problem_yaml)))
    # petab's lint report exposes has_errors(); fall back to truthiness defensively.
    has_errors = report.has_errors() if hasattr(report, "has_errors") else bool(report)
    if has_errors:
        return "errors", str(report)
    return "clean", ""


def roundtrip(conf_path: Path, job_type: str, method: str) -> int:
    from pybnf.petab import export_job, import_job

    conf_path = conf_path.resolve()
    folder = conf_path.parent
    workdir = Path(tempfile.mkdtemp(prefix="petab_roundtrip_"))
    petab_out = workdir / "petab"
    back = workdir / "back"

    # 1. EXPORT
    with _run_from(folder):
        try:
            export_job(str(conf_path), str(petab_out))
        except NotImplementedError as exc:
            print(f"FAIL  export: job uses a non-PEtab-exportable feature -- {exc}")
            print("      -> rework the conf to the exportable subset (see references/petab-compliance.md)")
            return 1
        except Exception as exc:  # noqa: BLE001
            print(f"FAIL  export: {exc!r}")
            return 1
    problem_yaml = petab_out / "problem.yaml"
    if not problem_yaml.is_file():
        print(f"FAIL  export: no problem.yaml written under {petab_out}")
        return 1
    wrote = sorted(p.name for p in petab_out.iterdir())
    print(f"ok    export -> {petab_out}")
    print(f"        wrote: {', '.join(wrote)}")

    # 2. LINT
    status, detail = _lint(problem_yaml)
    if status == "errors":
        print(f"FAIL  lint: petab.v2 reported errors\n{detail}")
        return 1
    if status == "skipped":
        print(f"warn  lint: {detail}")
    else:
        print("ok    lint  -> petab.v2 clean (no errors)")

    # 3. IMPORT
    with _run_from(folder):
        try:
            import_job(str(problem_yaml), str(back), job_type=job_type, method=method)
        except Exception as exc:  # noqa: BLE001
            print(f"FAIL  import: re-importing the exported problem failed -- {exc!r}")
            return 1
    imported = sorted(p.name for p in back.iterdir()) if back.is_dir() else []
    print(f"ok    import -> {back}")
    print(f"        wrote: {', '.join(imported)}")

    print()
    if status == "skipped":
        print("EXPORTED + IMPORTED, but LINT SKIPPED (petab extra missing) -- install pybnf[petab] to fully verify")
        return 3
    print("PASS  PEtab.v2 round-trip clean (export -> lint -> import)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("conf", type=Path, help="path to the edition-2 .conf to round-trip")
    ap.add_argument("--job-type", default="de", help="run recipe supplied on import (default: de)")
    ap.add_argument("--method", default="ode", help="simulator supplied on import: ode/ssa/nf (default: ode)")
    args = ap.parse_args()
    if not args.conf.is_file():
        print(f"FAIL  no such conf: {args.conf}")
        return 1
    return roundtrip(args.conf, args.job_type, args.method)


if __name__ == "__main__":
    sys.exit(main())
