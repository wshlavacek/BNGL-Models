#!/usr/bin/env python
"""Score a pybnf-jobs slug against its declared reference objective J*.

The acceptance bar for `curate-pybnf-job` (full rules: references/og-acceptance.md).

    J_paper = -log_likelihood          # Results/information_criteria.txt
    OG      = J_paper - J*             # optimality gap
    solved  = OG < 1.92                # chi^2, alpha = 0.05, 1 dof

PyBNF MINIMIZES a reduced objective (it drops the parameter-independent per-point constants --
1/2 log(2 pi), and, for a log-transformed observable, the change-of-variables Jacobian
sum log(y_obs * ln10)). It then REPORTS the full normalized log-likelihood at the best fit in
`Results/information_criteria.txt`, restoring every dropped constant (ADR-0056), so
`-log_likelihood` is an absolute, transform-agnostic scale comparable against a value another
tool reported. That is why OG is subtracted on that scale and NOT on the reduced one printed in
`sorted_params_final.txt`; this script prints the restored constant so the arithmetic is auditable.

`information_criteria.txt` is written at the END OF A FIT, so a job that has never been run has
none -- whatever its objective. It is genuinely absent only when the RESOLVED objective carries no
per-point log-likelihood (`objective.py:73-96`), which in practice means the legacy edition-1
`objfunc = sos` path (`SumOfSquaresObjective`, value `sum r^2`). Under edition 2, `objective = sos`
resolves to `Gaussian(sigma = 1)` -- a likelihood whose reduced value is `0.5*sum r^2` -- and is
fully scoreable. A constraint-only BPSL job scores no data points and is verified with
`job_type = check` instead of an OG.

Usage:
    score.py [JOBDIR]           # shipped provenance (jstar.txt + information_criteria.txt)
    score.py [JOBDIR] output    # a fresh run's output/ (reads output/Results/)

JOBDIR defaults to the current directory. Exit 0 = solved, 1 = scored but not solved,
2 = could not score (missing or unusable provenance).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

SOLVED_THRESHOLD = 1.92  # 0.5 * chi2_1(0.95)

BADGES = {
    "solved": "SOLVED",
    "objective_validated": "OBJECTIVE VALIDATED (nominal point, not a fit)",
    "regression_anchored": "REGRESSION-ANCHORED (T3 J*, not an optimality claim)",
    "setup_only": "SETUP ONLY",
}


def _die(msg: str, hint: str = "") -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    if hint:
        print(f"       {hint}", file=sys.stderr)
    raise SystemExit(2)


def _read_jstar(job: Path) -> tuple[float, str]:
    """Return (J*, source note). `jstar.txt` is one bare number; `#` lines are provenance."""
    path = job / "jstar.txt"
    if not path.is_file():
        _die(
            f"no {path}",
            "Declare J* before fitting -- see references/og-acceptance.md section 3 "
            "(T1 published, T2 published-parameter, T3 corpus-best).",
        )
    note = ""
    value = None
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith("#"):
            note = note or line.lstrip("# ").strip()
            continue
        if value is None:
            value = float(line.split()[0])
    if value is None:
        _die(f"{path} holds no number")
    return value, note


NON_LIKELIHOOD = ("sos", "sod", "norm_sos", "kl", "wasserstein", "direct_pass")  # edition-1 only


def _conf_objective(job: Path) -> str:
    """The `objective =` token from the slug's conf, or '' if not determinable."""
    for conf in sorted(job.glob("*.conf")):
        if conf.stem.endswith("_check"):
            continue
        for line in conf.read_text().splitlines():
            if line.strip().startswith("objective"):
                _, _, rhs = line.partition("=")
                return rhs.strip()
    return ""


def _read_information_criteria(path: Path, job: Path) -> dict[str, str]:
    if not path.is_file():
        objective = _conf_objective(job)
        hint = (
            "PyBNF writes this at the END OF A FIT, so the usual cause is simply that no fit has "
            "been run yet: run `pybnf -c <slug>.conf`, then re-score with the run's output/ dir. "
            "It is NOT written when the RESOLVED objective carries no per-point log-likelihood "
            "(objective.py:73-96) -- in practice the legacy edition-1 `objfunc = sos` path. Note "
            "edition-2 `objective = sos` resolves to Gaussian(sigma=1) and IS scoreable"
            + (f" (this conf declares `{objective}`)." if objective else ".")
        )
        _die(f"no {path}", hint)
    fields: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 2:
            fields[parts[0]] = parts[1]
    return fields


def _read_reduced_objective(path: Path) -> float | None:
    """Best row of a PyBNF params file: `<name> <objective> ...` on the first data line."""
    if not path.is_file():
        return None
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                return float(parts[1])
            except ValueError:
                return None
    return None


def _nominal_check(job: Path) -> dict:
    """The slug's nominal_check.json, or {} if absent/unreadable."""
    path = job / "nominal_check.json"
    if not path.is_file():
        return {}
    try:
        blob = json.loads(path.read_text())
    except (OSError, ValueError):
        return {}
    return blob if isinstance(blob, dict) else {}


def _report_declared(job: Path, jstar: float, note: str, blob: dict) -> int:
    """No fit artifacts, but the slug declares a nominal OG: surface it, flagged as declared."""
    og = blob.get("OG_nominal")
    status = str(blob.get("status", ""))
    tier = str(blob.get("jstar_tier", ""))
    print(f"slug                          {job.name}")
    print("source                        nominal_check.json (DECLARED, not measured here)")
    if tier:
        print(f"J* provenance tier            {tier}" + (f"  ({note})" if note else ""))
    print(f"reference               J*   = {jstar:.6f}")
    if og is None:
        print("OPTIMALITY GAP  OG           = (not declared)")
        print("=> UNSCORED -- no fit artifacts and no OG_nominal in nominal_check.json")
        return 1
    print(f"OPTIMALITY GAP  OG           = {float(og):.6g}  † at the nominal point, not a fit")
    badge = BADGES.get(status) or (
        BADGES["objective_validated"] if float(og) < SOLVED_THRESHOLD else BADGES["setup_only"]
    )
    print(f"=> {badge} (threshold OG < {SOLVED_THRESHOLD})")
    print(
        "\nThis is a DECLARED value read from nominal_check.json, not one measured from a run. "
        "To measure: `pybnf -c <slug>.conf`, then re-score with the run's output/ dir.",
        file=sys.stderr,
    )
    return 0 if status in ("solved", "objective_validated") else 1


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("jobdir", nargs="?", default=".", help="slug folder (default: cwd)")
    ap.add_argument("rundir", nargs="?", default=None, help="a run's output/; reads its Results/")
    args = ap.parse_args()

    job = Path(args.jobdir).resolve()
    if not job.is_dir():
        _die(f"{job} is not a directory")

    if args.rundir:
        run = Path(args.rundir)
        results = (run if run.is_absolute() else job / run) / "Results"
        ic_path = results / "information_criteria.txt"
        params_path = results / "sorted_params_final.txt"
        provenance = f"run {results}"
    else:
        ic_path = job / "information_criteria.txt"
        params_path = job / "best_fit_params.txt"
        provenance = "shipped provenance"

    jstar, jstar_note = _read_jstar(job)
    blob = _nominal_check(job)

    # A setup-only / objective-validated slug has no fit artifacts but does declare a nominal OG.
    # Surface that rather than dying -- but never let a declared number look like a measured one.
    if not ic_path.is_file() and "OG_nominal" in blob:
        return _report_declared(job, jstar, jstar_note, blob)

    ic = _read_information_criteria(ic_path, job)

    try:
        n = int(ic["n"])
        k = int(ic["k"])
        ln_l = float(ic["log_likelihood"])
    except (KeyError, ValueError):
        _die(f"{ic_path} is missing k / n / log_likelihood")

    j_paper = -ln_l
    og = j_paper - jstar
    solved = og < SOLVED_THRESHOLD

    reduced = _read_reduced_objective(params_path)
    tier = str(blob.get("jstar_tier", ""))
    declared_status = str(blob.get("status", ""))

    print(f"slug                          {job.name}")
    print(f"source                        {provenance}")
    if tier:
        print(f"J* provenance tier            {tier}" + (f"  ({jstar_note})" if jstar_note else ""))
    elif jstar_note:
        print(f"J* source                     {jstar_note}")
    print(f"N scored points           n = {n}")
    print(f"free parameters           k = {k}")
    if reduced is not None:
        restored = j_paper - reduced
        const = n * 0.5 * math.log(2 * math.pi)
        print(f"PyBNF reduced objective      = {reduced:.6f}   (minimized; drops constants)")
        print(
            f"restored constants           = {restored:.6f}   "
            f"(= N/2 log2pi = {const:.4f} for linear; + Jacobian for log-transformed)"
        )
    print(f"J_paper = -log_likelihood    = {j_paper:.6f}   (full normalized NLL; comparable)")
    print(f"reference               J*   = {jstar:.6f}")
    print(f"OPTIMALITY GAP  OG = J-J*    = {og:.6f}")

    if tier == "T3":
        badge = BADGES["regression_anchored"]
    elif declared_status in BADGES:
        badge = BADGES[declared_status]
    else:
        badge = BADGES["solved"] if solved else BADGES["setup_only"]
    print(f"=> {badge} (threshold OG < {SOLVED_THRESHOLD})")

    if not solved and tier != "T3":
        print(
            "\nA large OG is a diagnosis prompt, not a verdict: check the J* scale, then the "
            "forward model at the published parameters, then the budget, then the method. "
            "See references/og-acceptance.md section 6.",
            file=sys.stderr,
        )
    return 0 if solved else 1


if __name__ == "__main__":
    sys.exit(main())
