#!/usr/bin/env python3
"""Verify linearized Lotka-Volterra near equilibrium."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT / "analytical_odes"))

from _verify_common import parse_gdat, read_numeric_parameters, report_error, require_small_errors, run_bionetgen

MODEL = HERE / "linearized_lotka_volterra.bngl"


def main() -> None:
    p = read_numeric_parameters(MODEL)
    with tempfile.TemporaryDirectory(prefix="linearized_lv_") as tmp:
        columns, data = parse_gdat(run_bionetgen(MODEL, Path(tmp), "linearized_lotka_volterra", REPO_ROOT))

    t = data[:, columns.index("time")]
    omega = np.sqrt(p["alpha"] * p["gamma"])
    exact_p = p["P0"] * np.cos(omega * t) - p["alpha"] * p["Q0"] / omega * np.sin(omega * t)
    exact_q = p["Q0"] * np.cos(omega * t) + p["gamma"] * p["P0"] / omega * np.sin(omega * t)
    prey_star = p["gamma"] / p["delta"]
    predator_star = p["alpha"] / p["beta"]
    exact_prey = prey_star * (1 + exact_p)
    exact_predator = predator_star * (1 + exact_q)

    errors = [
        report_error("Obs_P_vs_python", data[:, columns.index("Obs_P")], exact_p),
        report_error("Obs_Q_vs_python", data[:, columns.index("Obs_Q")], exact_q),
        report_error("Analytical_P_vs_python", data[:, columns.index("Analytical_P")], exact_p),
        report_error("Analytical_Q_vs_python", data[:, columns.index("Analytical_Q")], exact_q),
        report_error("Analytical_Prey_vs_python", data[:, columns.index("Analytical_Prey")], exact_prey),
        report_error("Analytical_Predator_vs_python", data[:, columns.index("Analytical_Predator")], exact_predator),
    ]
    print(f"points: {len(t)}")
    require_small_errors(errors)


if __name__ == "__main__":
    main()
