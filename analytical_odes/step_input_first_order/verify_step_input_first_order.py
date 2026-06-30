#!/usr/bin/env python3
"""Verify first-order system with step input."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT / "analytical_odes"))

from _verify_common import parse_gdat, read_numeric_parameters, report_error, require_small_errors, run_bionetgen

MODEL = HERE / "step_input_first_order.bngl"


def main() -> None:
    p = read_numeric_parameters(MODEL)
    with tempfile.TemporaryDirectory(prefix="step_input_first_order_") as tmp:
        columns, data = parse_gdat(run_bionetgen(MODEL, Path(tmp), "step_input_first_order", REPO_ROOT))

    t = data[:, columns.index("time")]
    x_ss0 = p["J_base"] / p["k"]
    x_ss1 = (p["J_base"] + p["J_step"]) / p["k"]
    x_tau = x_ss0 + (p["X0"] - x_ss0) * np.exp(-p["k"] * p["tau"])
    exact_x = np.where(
        t < p["tau"],
        x_ss0 + (p["X0"] - x_ss0) * np.exp(-p["k"] * t),
        x_ss1 + (x_tau - x_ss1) * np.exp(-p["k"] * (t - p["tau"])),
    )

    errors = [
        report_error("Obs_X_vs_python", data[:, columns.index("Obs_X")], exact_x),
        report_error("Analytical_X_vs_python", data[:, columns.index("Analytical_X")], exact_x),
        report_error("Obs_X_vs_Analytical_X", data[:, columns.index("Obs_X")], data[:, columns.index("Analytical_X")]),
        report_error("PreStep_Steady_X_vs_python", data[:, columns.index("PreStep_Steady_X")], np.full_like(t, x_ss0)),
        report_error("PostStep_Steady_X_vs_python", data[:, columns.index("PostStep_Steady_X")], np.full_like(t, x_ss1)),
    ]
    print(f"points: {len(t)}")
    require_small_errors(errors)


if __name__ == "__main__":
    main()
