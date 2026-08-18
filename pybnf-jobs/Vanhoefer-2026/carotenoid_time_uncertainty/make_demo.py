#!/usr/bin/env python
"""Reproduce the Fig 5/6 demonstration fits for a given sigma_t level and copy the results.

The shipped standard.conf / marginal.conf search the full loguniform box [1e-5, 1000] with a
20-start multi-start -- the honest "from scratch" reproduction, which is expensive (the paper
reports ~1500 s per converged start). This script instead runs the CONTROLLED demonstration behind
VALIDATION.md's gate table: it seeds a SINGLE gntr start at theta* (by tightening each parameter's
box to [theta*/100, theta**100], whose log-centre is theta*) so both arms descend from the exact
ground truth, and the standard fit's DRAG away from theta* (or the marginal fit's staying put) is
attributable to the objective, not the start. It then copies each arm's best fit +
information_criteria into demo_results/sigmaT<level>_<arm>/ and scores the pair.

Usage:
    BNGPATH=$HOME/Simulations/BioNetGen-2.9.3 python make_demo.py --level 10
    python make_demo.py --level 0 --iterations 120
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def theta_star():
    lines = [l for l in (HERE / "theta_star_source.txt").read_text().splitlines() if l.strip()]
    names = lines[0].lstrip("#").split()[2:]
    vals = lines[1].split()[2:]
    return {n: float(v) for n, v in zip(names, vals)}


def seeded_conf(base_conf, level, iterations, out_tag):
    """A copy of base_conf seeded at theta* (log-box centre) on the sigmaT<level> data."""
    theta = theta_star()
    src = (HERE / base_conf).read_text()
    for name, value in theta.items():
        src = re.sub(rf"(?m)^loguniform_var = {name} .*$",
                     f"loguniform_var = {name} {value/100:.10g} {value*100:.10g}", src)
    src = (src.replace("population_size = 20", "population_size = 1")
              .replace("max_iterations = 300", f"max_iterations = {iterations}")
              .replace("sigmaT5/", f"sigmaT{level:g}/")
              .replace(f"output/{base_conf[:-5]}", f"output/{out_tag}"))
    path = HERE / f"_demo_{out_tag}.conf"
    path.write_text(src)
    return path


def run(conf_path):
    subprocess.run([sys.executable, "-m", "pybnf", "-c", conf_path.name],
                   cwd=HERE, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=float, required=True, help="sigma_t level (a sigmaT<level>/ dir)")
    ap.add_argument("--iterations", type=int, default=120)
    args = ap.parse_args()
    if "BNGPATH" not in os.environ:
        sys.exit("set BNGPATH to the folder containing BNG2.pl")
    lvl = args.level
    for base, arm in (("standard.conf", "standard"), ("marginal.conf", "marginal")):
        tag = f"demo_sigmaT{lvl:g}_{arm}"
        conf = seeded_conf(base, lvl, args.iterations, tag)
        shutil.rmtree(HERE / "output" / tag, ignore_errors=True)
        run(conf)
        dst = HERE / "demo_results" / f"sigmaT{lvl:g}_{arm}" / "Results"
        dst.mkdir(parents=True, exist_ok=True)
        src = HERE / "output" / tag / "Results"
        for f in ("sorted_params_final.txt", "information_criteria.txt"):
            shutil.copy(src / f, dst / f)
        conf.unlink()
    subprocess.run([sys.executable, str(HERE / "score.py"),
                    "--standard", f"demo_results/sigmaT{lvl:g}_standard",
                    "--marginal", f"demo_results/sigmaT{lvl:g}_marginal",
                    "--injected", str(lvl)], cwd=HERE, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
