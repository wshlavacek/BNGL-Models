#!/usr/bin/env python3
"""Verify one-compartment infusion-washout PK."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT / "analytical_odes"))

from _verify_common import parse_gdat, read_numeric_parameters, report_error, require_small_errors, run_bionetgen

MODEL = HERE / "infusion_washout_pk.bngl"


def main() -> None:
    p = read_numeric_parameters(MODEL)
    with tempfile.TemporaryDirectory(prefix="infusion_washout_pk_") as tmp:
        columns, data = parse_gdat(run_bionetgen(MODEL, Path(tmp), "infusion_washout_pk", REPO_ROOT))

    t = data[:, columns.index("time")]
    a_ss = p["R_in"] / p["kel"]
    a_stop = a_ss + (p["A0"] - a_ss) * np.exp(-p["kel"] * p["t_infusion"])
    exact_a = np.where(
        t < p["t_infusion"],
        a_ss + (p["A0"] - a_ss) * np.exp(-p["kel"] * t),
        a_stop * np.exp(-p["kel"] * (t - p["t_infusion"])),
    )
    exact_c = exact_a / p["Vd"]

    errors = [
        report_error("Obs_Drug_vs_python", data[:, columns.index("Obs_Drug")], exact_a),
        report_error("Analytical_Drug_vs_python", data[:, columns.index("Analytical_Drug")], exact_a),
        report_error("Observed_Concentration_vs_python", data[:, columns.index("Observed_Concentration")], exact_c),
        report_error("Analytical_Concentration_vs_python", data[:, columns.index("Analytical_Concentration")], exact_c),
        report_error("Infusion_Stop_Amount_vs_python", data[:, columns.index("Infusion_Stop_Amount")], np.full_like(t, a_stop)),
    ]
    print(f"points: {len(t)}")
    require_small_errors(errors)


if __name__ == "__main__":
    main()
