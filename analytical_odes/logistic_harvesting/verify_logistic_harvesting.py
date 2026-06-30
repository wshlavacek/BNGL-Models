#!/usr/bin/env python3
"""Verify logistic growth with constant harvesting."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT / "analytical_odes"))

from _verify_common import parse_gdat, read_numeric_parameters, report_error, require_small_errors, run_bionetgen

MODEL = HERE / "logistic_harvesting.bngl"


def main() -> None:
    p = read_numeric_parameters(MODEL)
    with tempfile.TemporaryDirectory(prefix="logistic_harvesting_") as tmp:
        columns, data = parse_gdat(run_bionetgen(MODEL, Path(tmp), "logistic_harvesting", REPO_ROOT))

    t = data[:, columns.index("time")]
    q = p["r"] / p["K"]
    d_root = np.sqrt(p["r"] * p["r"] - 4 * q * p["h"])
    x_plus = (p["r"] + d_root) / (2 * q)
    x_minus = (p["r"] - d_root) / (2 * q)
    ratio = ((p["X0"] - x_plus) / (p["X0"] - x_minus)) * np.exp(-d_root * t)
    exact_x = (x_plus - ratio * x_minus) / (1 - ratio)

    errors = [
        report_error("Obs_X_vs_python", data[:, columns.index("Obs_X")], exact_x),
        report_error("Analytical_X_vs_python", data[:, columns.index("Analytical_X")], exact_x),
        report_error("Obs_X_vs_Analytical_X", data[:, columns.index("Obs_X")], data[:, columns.index("Analytical_X")]),
        report_error("Stable_Root_vs_python", data[:, columns.index("Stable_Root")], np.full_like(t, x_plus)),
        report_error("Unstable_Root_vs_python", data[:, columns.index("Unstable_Root")], np.full_like(t, x_minus)),
    ]
    print(f"points: {len(t)}")
    require_small_errors(errors)


if __name__ == "__main__":
    main()
