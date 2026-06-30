#!/usr/bin/env python3
"""Verify Bateman chain with distinct rates."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT / "analytical_odes"))

from _verify_common import parse_gdat, read_numeric_parameters, report_error, require_small_errors, run_bionetgen

MODEL = HERE / "bateman_chain_distinct_rates.bngl"


def main() -> None:
    p = read_numeric_parameters(MODEL)
    with tempfile.TemporaryDirectory(prefix="bateman_chain_") as tmp:
        columns, data = parse_gdat(run_bionetgen(MODEL, Path(tmp), "bateman_chain_distinct_rates", REPO_ROOT))

    t = data[:, columns.index("time")]
    total = p["A0"] + p["B0"] + p["C0"]
    exact_a = p["A0"] * np.exp(-p["k1"] * t)
    exact_b = p["k1"] * p["A0"] / (p["k2"] - p["k1"]) * (np.exp(-p["k1"] * t) - np.exp(-p["k2"] * t))
    exact_c = total - exact_a - exact_b
    exact_total = np.full_like(t, total)

    errors = [
        report_error("Obs_A_vs_python", data[:, columns.index("Obs_A")], exact_a),
        report_error("Obs_B_vs_python", data[:, columns.index("Obs_B")], exact_b),
        report_error("Obs_C_vs_python", data[:, columns.index("Obs_C")], exact_c),
        report_error("Analytical_A_vs_python", data[:, columns.index("Analytical_A")], exact_a),
        report_error("Analytical_B_vs_python", data[:, columns.index("Analytical_B")], exact_b),
        report_error("Analytical_C_vs_python", data[:, columns.index("Analytical_C")], exact_c),
        report_error("Conserved_Total_vs_python", data[:, columns.index("Conserved_Total")], exact_total),
    ]
    print(f"points: {len(t)}")
    require_small_errors(errors)


if __name__ == "__main__":
    main()
