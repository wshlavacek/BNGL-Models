#!/usr/bin/env python
"""Derive a hand-built slug's reference objective J* at its shipped nominal parameters.

The T2 tier of `references/og-acceptance.md`: for a job reconstructed from a paper, the
model's own parameter values ARE the paper's published best-fit point (that is what
`curate-pybnf-job` step 4 requires of a reconstruction), so PyBNF's objective evaluated
there is the published-parameter objective. That is the anchor an OG is measured against.

It evaluates the objective the same way a real fit does -- through
`likelihood_information_criteria`, which is exactly what writes `information_criteria.txt`
at the end of a run (`algorithms/base.py:1397-1432`) -- so the number is on the comparable
`-log_likelihood` scale and not on PyBNF's reduced one.

WHICH OBJECTIVES CAN BE SCORED. What matters is the RESOLVED objective, not the token in the
conf, and the two editions resolve `sos` differently:

    edition 2  `objective = sos`  -> LikelihoodObjective + Gaussian(sigma = 1)
                                     value = 0.5*sum r^2 (the Gaussian reduced term already
                                     carries the 1/2, `noise/gaussian.py:60`)
                                     supports_pointwise_log_likelihood = True  -> SCOREABLE
    edition 1  `objfunc  = sos`   -> SumOfSquaresObjective (a plain SummationObjective)
                                     value = sum r^2 (no 1/2, `objective.py:1628`)
                                     supports_pointwise_log_likelihood = False -> NOT scoreable

So an edition-2 `sos` job is a Gaussian likelihood and needs no conf change to be scored; only
the legacy edition-1 path is blocked (`objective.py:73-96`). This script therefore gates on the
resolved object, and records the restoration arithmetic auditors should re-run:

    J_paper == reduced + sum log(sigma_i) + (N/2)*log(2*pi)

`sum log(sigma_i)` is the Gaussian normalizer (`noise/gaussian.py:122`) -- exactly zero when
sigma == 1 (edition-2 `sos`), and generally negative for `chi_sq` when the data's _SD are below
1. It is recorded as `sum_log_sigma`. N is also counted independently here and cross-checked
against `ic.n`, so neither is on faith.

KNOWN LIMITATION -- network-free (NFsim) jobs. `Monine-2010/tlbr` and `Kozer-2013/egfr_nf` crash
(`fail_type = 1`) under this single-evaluation harness, although both run to completion under a
full `pybnf` fit and are validated end to end. Treat that as a gap in THIS script's NF path, not
as a defect in those jobs; anchor them from a real run (T3) until it is fixed.

STOCHASTIC MODELS. `reduced` and the information criteria come from two independent simulations,
so for an SSA/NF model they are different realizations and J* is a SINGLE-DRAW estimate. Observed
spread is large -- `Rijal-2025/lacud5_ssa` moved 304.996 -> 360.791 between two draws. Such a slug
is flagged `stochastic` in its record; pin `stochastic_seed` or average replicates before treating
its J* as a stable anchor.

Usage:
    make_jstar.py <slug-dir> [--write]
    make_jstar.py --batch <corpus-dir> [--write] [--skip PATTERN]

Exit 0 = J* derived, 1 = blocked (reason printed), 2 = could not load the job.
"""

from __future__ import annotations

import argparse
import ast
import json
import math
import os
import re
import sys
import tempfile
import traceback
from pathlib import Path

LOG_2PI = math.log(2.0 * math.pi)


def declared_objective(conf_path: Path) -> str:
    """The objective token as WRITTEN in the conf.

    `config.config['objfunc']` reports the schema default ('chi_sq') when the conf uses the
    edition-2 `objective =` spelling, so it cannot be used to report what the job declares.
    """
    for line in conf_path.read_text().splitlines():
        stripped = line.split("#")[0].strip()
        for key in ("objective", "objfunc"):
            if stripped.startswith(key):
                _, sep, rhs = stripped.partition("=")
                if sep:
                    return rhs.strip()
    return ""


_ALLOWED_NODES = (
    ast.Expression, ast.BinOp, ast.UnaryOp, ast.Constant, ast.Name, ast.Load,
    ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.USub, ast.UAdd, ast.Mod,
)
_MATH_NAMES = {k: v for k, v in vars(math).items() if not k.startswith("_")}


def eval_param_expression(expr: str, known: dict[str, float]) -> float:
    """Evaluate a BNGL parameter expression over already-declared parameters.

    BNGL's parameters block is FILE-ORDERED -- a parameter may only reference names declared
    above it -- so a single forward pass resolves every derived value. Restricted to arithmetic
    over known names and math constants/functions; anything else raises.
    """
    tree = ast.parse(expr.replace("^", "**"), mode="eval")
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):  # allow math functions only
            if not isinstance(node.func, ast.Name) or node.func.id not in _MATH_NAMES:
                raise ValueError(f"disallowed call in {expr!r}")
            continue
        if not isinstance(node, _ALLOWED_NODES):
            raise ValueError(f"disallowed syntax in {expr!r}")
    scope = dict(_MATH_NAMES)
    scope.update(known)
    return float(eval(compile(tree, "<bngl>", "eval"), {"__builtins__": {}}, scope))  # noqa: S307


def parse_nominals(bngl_path: Path) -> dict[str, float]:
    """{id: value} for the model's parameters block, resolving derived expressions.

    Handles both BNGL spellings (`id value` and `id=value`) and an optional leading index. A
    parameter whose value is an expression over earlier parameters (e.g. `beta_SRC = 0.0574/
    phen_scale`) is evaluated, because that expression IS its published nominal.
    """
    text = bngl_path.read_text()
    match = re.search(r"begin parameters(.*?)end parameters", text, re.S)
    if not match:
        return {}
    out: dict[str, float] = {}
    for raw in match.group(1).splitlines():
        line = raw.split("#")[0].strip()
        if not line:
            continue
        line = re.sub(r"^\d+\s+", "", line)  # optional leading index
        m = re.match(r"^([A-Za-z_]\w*)\s*=?\s*(.+)$", line)
        if not m:
            continue
        name, value = m.group(1), m.group(2).strip()
        try:
            out[name] = float(value)
            continue
        except ValueError:
            pass
        try:
            out[name] = eval_param_expression(value, out)
        except Exception:  # noqa: BLE001 -- an unresolvable expression is simply not a nominal
            continue
    return out


def noise_nuisances(conf_path: Path) -> set[str]:
    """Free-parameter names that are ESTIMATED NOISE scales, not model parameters.

    These live only in the conf -- `noise_model = <family>, <param> = fit <name>` and
    `objective = chi_sq_dynamic` (implicitly `sigma__FREE`) -- so the model declares no nominal
    for them, by design. A published parameter table fixes the model, not the noise scale, so a
    T2 anchor is undefined when one is free.
    """
    names: set[str] = set()
    for line in conf_path.read_text().splitlines():
        stripped = line.split("#")[0].strip()
        if stripped.startswith("noise_model"):
            for match in re.finditer(r"=\s*fit\s+([A-Za-z_]\w*)", stripped):
                names.add(match.group(1))
                names.add(match.group(1) + "__FREE")
        elif stripped.startswith("objective") and "chi_sq_dynamic" in stripped:
            names.update({"sigma", "sigma__FREE"})
    return names


def count_scored_points(config) -> int:
    """Independent count of scored data points: non-NaN cells in comparable columns.

    Mirrors `SummationObjective.evaluate`'s loop -- skip the independent variable, skip
    `_SD` sidecars, skip NaN. Cross-checked against `ic.n` wherever a likelihood objective
    makes `ic` available, so a drift in either would surface immediately.
    """
    import numpy as np

    def walk(node) -> int:
        # config.exp_data is nested {model: {suffix: Data}}; recurse to the Data leaves.
        if isinstance(node, dict):
            return sum(walk(child) for child in node.values())
        cols = getattr(node, "cols", None)
        data = getattr(node, "data", None)
        if cols is None or data is None:
            return 0
        indvar = min(cols, key=cols.get)

        def is_sidecar(name: str) -> bool:
            # `X_SD` is a sigma sidecar ONLY when `X` is itself a column here. A bare `mRNA_SD`
            # with no `mRNA` column is its own observable and IS scored -- Rijal-2025 ships
            # exactly that, and reading it as a sidecar halves the count.
            return name.endswith("_SD") and name[: -len("_SD")] in cols

        return sum(
            int(np.count_nonzero(~np.isnan(data[:, idx])))
            for name, idx in cols.items()
            if name != indvar and not is_sidecar(name)
        )

    return walk(config.exp_data)


def evaluate_at(config, pset, tag: str):
    """One fresh Result scored at `pset`.

    A measurement-model layer (ADR-0036) materializes observable columns onto the Data in
    place, so a second scoring pass over the SAME Result collides with the columns the first
    one added. A real fit re-simulates for the information criteria for exactly this reason
    (`_compute_information_criteria`); do the same here.
    """
    from pybnf.algorithms import core

    with tempfile.TemporaryDirectory() as sim_dir:
        job = core.Job(
            list(config.models.values()), pset, tag, sim_dir,
            config.config["wall_time_sim"], None, config.config["normalization"],
            config.postprocessing, True,
            stochastic_seed_policy=config.config["stochastic_seed"],
        )
        result = core.run_job(job, debug=True)
    if getattr(result, "failed", False) or result.simdata is None:
        # FailedSimulation carries WHY (`algorithms/core.py:120-136`); a bare "it failed" is not
        # actionable, and for a network-free model the answer is usually the wall clock.
        kinds = {0: "exceeded wall_time_sim", 1: "the simulator crashed", 2: "unknown error"}
        kind = kinds.get(getattr(result, "fail_type", None), "no simulation data returned")
        detail = (getattr(result, "traceback", "") or "").strip().splitlines()
        suffix = f" -- {detail[-1]}" if detail else ""
        raise RuntimeError(
            f"the nominal point failed to simulate: {kind}{suffix}"
            + (
                f". wall_time_sim is {config.config.get('wall_time_sim')}s; a network-free or "
                "dose-scan job may simply need more."
                if getattr(result, "fail_type", None) == 0 else ""
            )
        )
    result.normalize(config.config["normalization"])
    result.postprocess_data(config.postprocessing)
    return result


def derive(slug_dir: Path) -> dict:
    """Return the J* record for one slug, or raise with a reason."""
    from pybnf.objective import likelihood_information_criteria
    from pybnf.parse import load_config
    from pybnf.pset import PSet

    slug = slug_dir.name
    confs = [p for p in sorted(slug_dir.glob("*.conf")) if not p.stem.endswith("_check")]
    if not confs:
        raise RuntimeError("no .conf in the slug folder")
    conf = confs[0]
    if len(confs) > 1:
        named = [p for p in confs if p.stem == slug]
        conf = named[0] if named else conf

    cwd = Path.cwd()
    os.chdir(slug_dir)
    try:
        config = load_config(conf.name)
        free = [v.name for v in config.variables]
        if not free:
            raise RuntimeError("no free parameters")

        models = list(config.models.values())
        nominals: dict[str, float] = {}
        for model in models:
            path = Path(getattr(model, "file_path", "") or "")
            if path.suffix == ".bngl" and path.is_file():
                nominals.update(parse_nominals(path))
        for bngl in sorted(slug_dir.glob("*.bngl")):
            for key, value in parse_nominals(bngl).items():
                nominals.setdefault(key, value)

        # Edition-1 jobs name free params `<id>__FREE` (ADR-0034 retired the alias); the model
        # still declares the bare id, so fall back to the stripped name for those.
        for name in free:
            if name.endswith("__FREE") and name not in nominals:
                base = name[: -len("__FREE")]
                if base in nominals:
                    nominals[name] = nominals[base]

        missing = [n for n in free if n not in nominals]
        if missing:
            nuisance = noise_nuisances(conf)
            if all(n in nuisance for n in missing):
                raise RuntimeError(
                    "the free parameters include estimated NOISE nuisances with no model "
                    "declaration (" + ", ".join(sorted(missing)) + "), so a T2 anchor is "
                    "undefined: the paper's published point fixes the model parameters but not "
                    "an estimated noise scale. Use a T3 corpus-best anchor from a real fit "
                    "instead (references/og-acceptance.md section 3)"
                )
            raise RuntimeError(
                "no bare-numeric nominal in the model for: " + ", ".join(missing)
                + " (an expression-valued or absent declaration -- set the published value "
                  "as a bare nominal, per curate-pybnf-job step 4)"
            )

        objective_name = declared_objective(conf)
        pset = PSet([v.set_value(nominals[v.name]) for v in config.variables])
        n_counted = count_scored_points(config)

        reduced = config.obj.evaluate_multiple(
            evaluate_at(config, pset, "jstar_obj").simdata, config.exp_data, pset,
            config.constraints,
        )
        if reduced is None:
            raise RuntimeError("the objective did not evaluate at the nominal point")

        ic = likelihood_information_criteria(
            config.obj, evaluate_at(config, pset, "jstar_ic").simdata, config.exp_data,
            pset, len(config.variables),
        )

        if ic is None:
            if n_counted == 0:
                raise RuntimeError(
                    "no data points are scored -- this is a constraint-only BPSL job, which has "
                    "no data-fit likelihood and therefore no OG. Verify it with `job_type = check` "
                    "constraint satisfaction instead (references/bpsl-constraints.md)"
                )
            raise RuntimeError(
                f"the resolved objective ({type(config.obj).__name__}, declared "
                f"`{objective_name or '?'}`, edition {config.config.get('edition')}) carries no "
                "per-point log-likelihood (objective.py:73-96), so there is no -lnL scale to put "
                "J* on. Edition-1 `objfunc = sos` is the usual cause; the edition-2 `objective = "
                "sos` spelling resolves to Gaussian(sigma=1) and IS scoreable."
            )

        j_paper = -float(ic.log_likelihood)
        # For a UNIT-sigma Gaussian (edition-2 `sos`) the full NLL is exactly
        #   reduced + (N/2)log(2pi).
        # Any other fixed sigma adds the normalizer sum log(sigma_i) (`noise/gaussian.py:122`),
        # so the residual below IS that sum -- negative when the data's _SD are mostly < 1.
        # It is therefore expected to be non-zero for chi_sq, and is recorded, not flagged.
        sigma1 = float(reduced) + 0.5 * int(ic.n) * LOG_2PI
        residual = j_paper - sigma1
        # Two conditions make that residual NOT the Gaussian normalizer, so it must not be
        # reported as one:
        #  - a constraint-bearing job: `evaluate_multiple` adds the BPSL penalty to `reduced`,
        #    while `likelihood_information_criteria` scores only data points, so the two
        #    quantities measure different things;
        #  - a stochastic model: the two scoring passes are independent realizations, so the
        #    difference carries simulation noise (and J* itself is a single-draw estimate).
        has_constraints = bool(getattr(config, "constraints", None))
        stochastic = bool(any(getattr(m, "stochastic", False) for m in models))
        decomposes = not has_constraints and not stochastic
        return {
            "slug": slug,
            "conf": conf.name,
            "objective": objective_name,
            "objective_class": type(config.obj).__name__,
            "jstar_tier": "T2",
            "jstar_source": "PyBNF objective at the model's shipped nominal (published) parameters",
            "jstar_scale": "-log_likelihood (native, from likelihood_information_criteria)",
            "jstar": j_paper,
            "J_paper": j_paper,
            "reduced_objective": float(reduced),
            "k": int(ic.k),
            "n_scored": int(ic.n),
            "n_counted_here": n_counted,
            "n_count_agrees": bool(int(ic.n) == n_counted),
            "has_constraints": has_constraints,
            "stochastic": stochastic,
            "sigma1_identity": sigma1 if decomposes else None,
            "sigma1_identity_exact": (
                bool(abs(residual) < 1e-6 * max(1.0, abs(j_paper))) if decomposes else None
            ),
            "sum_log_sigma": residual if decomposes else None,
            "sigma1_identity_note": (
                "J_paper - (reduced + (N/2)log(2pi)) = sum log(sigma_i), the Gaussian normalizer. "
                "Zero iff sigma == 1 (edition-2 `sos`); non-zero for chi_sq is expected, not a "
                "defect -- negative when the data's _SD are mostly below 1."
                if decomposes else
                "Not reported: " + (
                    "this job carries BPSL constraints, so `reduced` includes the constraint "
                    "penalty while J_paper scores only data points -- the two are not comparable."
                    if has_constraints else
                    "this model is stochastic, so the two scoring passes are independent "
                    "realizations and J* is a SINGLE-DRAW estimate. Pin `stochastic_seed` in the "
                    "conf, or average replicates, before treating it as a stable anchor."
                )
            ),
            "scoreable_end_to_end": True,
        }
    finally:
        os.chdir(cwd)


def write_out(slug_dir: Path, record: dict) -> None:
    tier, source = record["jstar_tier"], record["jstar_source"]
    lines = [
        f"{record['jstar']!r}",
        f"# tier {tier} -- {source}",
        f"# scale: {record['jstar_scale']}",
        f"# k = {record['k']}, n = {record['n_scored']}, objective = {record['objective']}",
        "# regenerate: skills/curate-pybnf-job/scripts/make_jstar.py <this dir> --write",
    ]
    (slug_dir / "jstar.txt").write_text("\n".join(lines) + "\n")

    path = slug_dir / "nominal_check.json"
    blob = {}
    if path.is_file():
        try:
            blob = json.loads(path.read_text())
        except ValueError:
            blob = {}
    blob.update(record)
    blob["OG_nominal"] = 0.0
    blob.setdefault(
        "interpretation",
        "J* is defined AS the objective at this slug's shipped nominal (published) parameters, "
        "so OG at that point is 0 by construction. It is a T2 anchor: a fit is measured against "
        "it, and a negative OG means PyBNF beat the published parameter set.",
    )
    blob["status"] = "setup_only"
    path.write_text(json.dumps(blob, indent=2) + "\n")


def run_one(slug_dir: Path, write: bool) -> int:
    slug_dir = slug_dir.resolve()
    try:
        record = derive(slug_dir)
    except Exception as exc:  # noqa: BLE001 -- report every failure with its reason
        print(f"BLOCKED  {slug_dir.name}: {exc}")
        if os.environ.get("JSTAR_TRACEBACK"):
            traceback.print_exc()
        return 1
    extra = ""
    if record["has_constraints"]:
        extra = ", +BPSL constraints"
    elif record["stochastic"]:
        extra = ", STOCHASTIC (single draw)"
    elif not record["sigma1_identity_exact"]:
        extra = f", sum_log_sigma={record['sum_log_sigma']:.3f}"
    print(
        f"OK       {slug_dir.name}: J* = {record['jstar']:.6f}  "
        f"(k={record['k']}, n={record['n_scored']}, obj={record['objective']}{extra})"
    )
    if record.get("n_count_agrees") is False:
        print(
            f"  WARNING  n mismatch: ic.n={record['n_scored']} "
            f"counted={record['n_counted_here']}"
        )
    if write:
        write_out(slug_dir, record)
        print(f"  wrote    {slug_dir / 'jstar.txt'}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("target", help="a slug dir, or a corpus dir with --batch")
    ap.add_argument("--batch", action="store_true", help="walk <target>/*/*/ as slugs")
    ap.add_argument("--write", action="store_true", help="write jstar.txt + nominal_check.json")
    ap.add_argument("--skip", default="", help="substring; skip matching slug paths")
    args = ap.parse_args()

    root = Path(args.target).resolve()
    if not args.batch:
        return run_one(root, args.write)

    slugs = [d for d in sorted(root.glob("*/*/")) if any(d.glob("*.conf"))]
    if args.skip:
        slugs = [d for d in slugs if args.skip not in str(d)]
    print(f"{len(slugs)} slugs\n")
    blocked = 0
    for slug in slugs:
        blocked += run_one(slug, args.write)
        sys.stdout.flush()
    print(f"\n{len(slugs) - blocked} derived, {blocked} blocked")
    return 0 if blocked == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
