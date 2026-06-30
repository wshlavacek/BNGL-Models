#!/usr/bin/env python3
"""Verify damped harmonic oscillator."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT / "analytical_odes"))

from _verify_common import parse_gdat, read_numeric_parameters, report_error, require_small_errors, run_bionetgen

MODEL = HERE / "damped_harmonic_oscillator.bngl"


def main() -> None:
    p = read_numeric_parameters(MODEL)
    with tempfile.TemporaryDirectory(prefix="damped_harmonic_") as tmp:
        columns, data = parse_gdat(run_bionetgen(MODEL, Path(tmp), "damped_harmonic_oscillator", REPO_ROOT))

    t = data[:, columns.index("time")]
    a = p["zeta"] * p["omega0"]
    omega_d = p["omega0"] * np.sqrt(1 - p["zeta"] * p["zeta"])
    x_sin = (p["V0"] + a * p["X0"]) / omega_d
    v_sin = -(p["omega0"] * p["omega0"] * p["X0"] + a * p["V0"]) / omega_d
    decay = np.exp(-a * t)
    exact_x = decay * (p["X0"] * np.cos(omega_d * t) + x_sin * np.sin(omega_d * t))
    exact_v = decay * (p["V0"] * np.cos(omega_d * t) + v_sin * np.sin(omega_d * t))

    errors = [
        report_error("Obs_X_vs_python", data[:, columns.index("Obs_X")], exact_x),
        report_error("Obs_V_vs_python", data[:, columns.index("Obs_V")], exact_v),
        report_error("Analytical_X_vs_python", data[:, columns.index("Analytical_X")], exact_x),
        report_error("Analytical_V_vs_python", data[:, columns.index("Analytical_V")], exact_v),
    ]
    print(f"points: {len(t)}")
    require_small_errors(errors)


if __name__ == "__main__":
    main()
