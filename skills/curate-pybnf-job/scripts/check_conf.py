#!/usr/bin/env python
"""Backend-free tier-1 check for an edition-2 PyBNF job setup.

Mirrors ``tests/test_real_world_examples.py::test_real_world_conf_is_wellformed``
in the PyBNF repo, so passing this locally means the committed conf will pass the
default (bngsim-free) CI tier. It parses a ``.conf`` exactly the way the test does
-- ``parse.ploop`` then ``config.Configuration`` from *inside* the example folder,
because the conf's ``model:``/``data:`` paths are relative -- and asserts the
edition-2 contract:

  * ``edition == 2``
  * ``job_type`` resolves to a real algorithm (in ``FIT_TYPE_REGISTRY``)
  * at least one experiment bound its data (``exp_data`` non-empty)
  * at least one free parameter, none carrying a ``__FREE`` alias (ADR-0034: free
    params bind to model ids by name under edition 2)

It also reports the resolved simulator (ode / ssa / nf) so you can set the
manifest's ``simulator`` / ``stochastic`` fields correctly -- the test cross-checks
``model.stochastic`` against the manifest, so they must agree.

This imports ``pybnf``, so run it with the PyBNF environment's interpreter, e.g.

    ~/Code/PyBNF/.venv/bin/python \
        skills/curate-pybnf-job/scripts/check_conf.py \
        ~/Code/PyBNF/examples/real-world/<name>/<name>.conf

Exit code 0 = all checks passed; 1 = a check failed; 2 = could not even parse.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def _load_conf(conf_path: Path):
    """Parse a conf into a Configuration, from inside its folder (relative paths)."""
    from pybnf import config, parse  # noqa: WPS433 (import here so --help works w/o pybnf)

    text = conf_path.read_text()
    home = os.getcwd()
    os.chdir(conf_path.parent)
    try:
        raw = parse.ploop(text.splitlines(keepends=True))
        # Pin the noisy/nondeterministic knobs the CI harness also pins, so a
        # check never fails for an unrelated reason (missing verbosity, etc.).
        raw.setdefault("verbosity", 0)
        raw.setdefault("random_seed", 1234)
        return config.Configuration(raw)
    finally:
        os.chdir(home)


def check(conf_path: Path) -> int:
    from pybnf.registry import FIT_TYPE_REGISTRY

    try:
        conf = _load_conf(conf_path)
    except Exception as exc:  # noqa: BLE001 -- surface any parse/load failure verbatim
        print(f"FAIL  parse: {conf_path} did not load: {exc!r}")
        return 2

    ok = True

    def report(passed: bool, label: str, detail: str = "") -> None:
        nonlocal ok
        ok = ok and passed
        tag = "ok  " if passed else "FAIL"
        print(f"{tag}  {label}{(' -- ' + detail) if detail else ''}")

    edition = conf.config.get("edition")
    report(edition == 2, "edition == 2", f"got {edition!r}")

    fit_type = conf.config.get("fit_type")
    report(
        fit_type in FIT_TYPE_REGISTRY,
        "job_type resolves to a real algorithm",
        f"fit_type={fit_type!r}",
    )

    # A job may carry quantitative data (.exp), BPSL constraints (.prop/.con), or both
    # (data fusion). At least one must be present.
    n_constraints = len(getattr(conf, "constraints", ()) or ())
    has_data = bool(conf.exp_data)
    report(
        has_data or n_constraints,
        "an experiment bound data and/or BPSL constraints",
        f"exp_data keys={list(conf.exp_data)}, constraint sets={n_constraints}",
    )

    # `check` jobs evaluate a single point (no search) -> free params are optional;
    # every fitting job_type needs at least one.
    names = [v.name for v in conf.variables]
    if fit_type == "check":
        print(f"info  job_type=check: model checking, no fitting ({len(names)} free params, not required)")
    else:
        report(bool(names), "at least one free parameter", f"{len(names)} declared")
    aliased = [n for n in names if "__FREE" in n]
    report(not aliased, "no free parameter uses a __FREE alias", f"offenders={aliased}" if aliased else "")

    # Resolved simulator -- informational, drives the manifest simulator/stochastic fields.
    stochastic = any(getattr(m, "stochastic", False) for m in conf.models.values())
    print(
        f"info  resolved model.stochastic = {stochastic}  "
        f"-> manifest: stochastic={stochastic} "
        f"(simulator is 'ssa' or 'nf' if stochastic else 'ode'; confirm which from the conf's method:)"
    )
    print(f"info  free parameters: {', '.join(names)}")
    if n_constraints:
        print(
            f"info  BPSL constraints present ({n_constraints} set(s)) -> this job is "
            f"NATIVE-ONLY: not PEtab-exportable. Verify with `job_type = check` "
            f"(satisfaction), not petab_roundtrip.py. See references/bpsl-constraints.md."
        )

    print()
    print("PASS" if ok else "FAILED (fix the FAIL lines above before committing)")
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("conf", type=Path, help="path to the edition-2 .conf to check")
    args = ap.parse_args()
    if not args.conf.is_file():
        print(f"FAIL  no such conf: {args.conf}")
        return 2
    return check(args.conf)


if __name__ == "__main__":
    sys.exit(main())
