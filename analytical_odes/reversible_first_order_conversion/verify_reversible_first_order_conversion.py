#!/usr/bin/env python3
"""Verify reversible first-order conversion."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT / "analytical_odes"))

from _verify_common import parse_gdat, read_numeric_parameters, report_error, require_small_errors, run_bionetgen

MODEL = HERE / "reversible_first_order_conversion.bngl"


def main() -> None:
    p = read_numeric_parameters(MODEL)
    with tempfile.TemporaryDirectory(prefix="reversible_first_order_") as tmp:
        columns, data = parse_gdat(run_bionetgen(MODEL, Path(tmp), "reversible_first_order_conversion", REPO_ROOT))

    t = data[:, columns.index("time")]
    total = p["A0"] + p["B0"]
    s_rate = p["kf"] + p["kr"]
    a_eq = p["kr"] * total / s_rate
    b_eq = p["kf"] * total / s_rate
    exact_a = a_eq + (p["A0"] - a_eq) * np.exp(-s_rate * t)
    exact_b = b_eq + (p["B0"] - b_eq) * np.exp(-s_rate * t)
    exact_total = np.full_like(t, total)

    errors = [
        report_error("Obs_A_vs_python", data[:, columns.index("Obs_A")], exact_a),
        report_error("Obs_B_vs_python", data[:, columns.index("Obs_B")], exact_b),
        report_error("Analytical_A_vs_python", data[:, columns.index("Analytical_A")], exact_a),
        report_error("Analytical_B_vs_python", data[:, columns.index("Analytical_B")], exact_b),
        report_error("Conserved_Total_vs_python", data[:, columns.index("Conserved_Total")], exact_total),
    ]
    print(f"points: {len(t)}")
    require_small_errors(errors)


if __name__ == "__main__":
    main()
