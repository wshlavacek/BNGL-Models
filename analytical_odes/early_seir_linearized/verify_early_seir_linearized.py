#!/usr/bin/env python3
"""Verify early linearized SEIR model."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT / "analytical_odes"))

from _verify_common import parse_gdat, read_numeric_parameters, report_error, require_small_errors, run_bionetgen

MODEL = HERE / "early_seir_linearized.bngl"


def main() -> None:
    p = read_numeric_parameters(MODEL)
    with tempfile.TemporaryDirectory(prefix="early_seir_linearized_") as tmp:
        columns, data = parse_gdat(run_bionetgen(MODEL, Path(tmp), "early_seir_linearized", REPO_ROOT))

    t = data[:, columns.index("time")]
    trace = p["sigma"] + p["gamma"]
    gap = np.sqrt((p["sigma"] - p["gamma"]) ** 2 + 4 * p["sigma"] * p["beta_eff"])
    lam1 = (-trace + gap) / 2
    lam2 = (-trace - gap) / 2
    e_dot0 = p["beta_eff"] * p["I0"] - p["sigma"] * p["E0"]
    e1 = (e_dot0 - lam2 * p["E0"]) / (lam1 - lam2)
    e2 = (lam1 * p["E0"] - e_dot0) / (lam1 - lam2)
    i1 = ((lam1 + p["sigma"]) / p["beta_eff"]) * e1
    i2 = ((lam2 + p["sigma"]) / p["beta_eff"]) * e2
    exact_e = e1 * np.exp(lam1 * t) + e2 * np.exp(lam2 * t)
    exact_i = i1 * np.exp(lam1 * t) + i2 * np.exp(lam2 * t)
    exact_r = p["R0_init"] + p["gamma"] * (i1 * (np.exp(lam1 * t) - 1) / lam1 + i2 * (np.exp(lam2 * t) - 1) / lam2)

    errors = [
        report_error("Obs_E_vs_python", data[:, columns.index("Obs_E")], exact_e),
        report_error("Obs_I_vs_python", data[:, columns.index("Obs_I")], exact_i),
        report_error("Obs_R_vs_python", data[:, columns.index("Obs_R")], exact_r),
        report_error("Analytical_E_vs_python", data[:, columns.index("Analytical_E")], exact_e),
        report_error("Analytical_I_vs_python", data[:, columns.index("Analytical_I")], exact_i),
        report_error("Analytical_R_vs_python", data[:, columns.index("Analytical_R")], exact_r),
        report_error("Growth_Eigenvalue_vs_python", data[:, columns.index("Growth_Eigenvalue")], np.full_like(t, lam1)),
    ]
    print(f"points: {len(t)}")
    require_small_errors(errors)


if __name__ == "__main__":
    main()
