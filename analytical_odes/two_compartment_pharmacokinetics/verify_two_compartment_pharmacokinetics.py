#!/usr/bin/env python3
"""Verify two-compartment pharmacokinetics against the analytical solution."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "two_compartment_pharmacokinetics.bngl"


def read_numeric_parameters(path: Path) -> dict[str, float]:
    parameters: dict[str, float] = {}
    in_parameters = False

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.split("#", 1)[0].strip()
            if not stripped:
                continue
            if stripped == "begin parameters":
                in_parameters = True
                continue
            if stripped == "end parameters":
                break
            if not in_parameters:
                continue

            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    parameters[parts[0]] = float(parts[1])
                except ValueError:
                    pass

    return parameters


def analytical_solution(
    t: np.ndarray,
    k12: float,
    k21: float,
    ke: float,
    central_0: float,
    peripheral_0: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Closed-form biexponential solution for the two-compartment model."""
    s_rate = k12 + k21 + ke
    d_rate = math.sqrt(s_rate * s_rate - 4 * ke * k21)
    lambda_1 = (-s_rate + d_rate) / 2
    lambda_2 = (-s_rate - d_rate) / 2
    central_dot_0 = -(k12 + ke) * central_0 + k21 * peripheral_0
    coef_1 = (central_dot_0 - lambda_2 * central_0) / (lambda_1 - lambda_2)
    coef_2 = (lambda_1 * central_0 - central_dot_0) / (lambda_1 - lambda_2)
    pcoef_1 = ((lambda_1 + k12 + ke) / k21) * coef_1
    pcoef_2 = ((lambda_2 + k12 + ke) / k21) * coef_2

    central = coef_1 * np.exp(lambda_1 * t) + coef_2 * np.exp(lambda_2 * t)
    peripheral = pcoef_1 * np.exp(lambda_1 * t) + pcoef_2 * np.exp(lambda_2 * t)
    return central, peripheral


def max_scaled_error(left: np.ndarray, right: np.ndarray) -> tuple[float, float]:
    abs_error = np.abs(left - right)
    max_abs_error = float(np.max(abs_error))
    scale = max(float(np.max(np.abs(right))), 1.0)
    return max_abs_error, max_abs_error / scale


def parse_gdat(path: Path) -> tuple[list[str], np.ndarray]:
    columns: list[str] | None = None
    rows: list[list[float]] = []

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                columns = stripped.lstrip("#").split()
            else:
                rows.append([float(value) for value in stripped.split()])

    if columns is None:
        raise ValueError(f"No header found in {path}")
    return columns, np.asarray(rows, dtype=float)


def find_bionetgen() -> str:
    executable = shutil.which("bionetgen")
    if executable is not None:
        return executable

    local_executable = REPO_ROOT / ".venv" / "bin" / "bionetgen"
    if local_executable.exists():
        return str(local_executable)

    raise RuntimeError("BioNetGen executable 'bionetgen' was not found")


def run_bionetgen(output_dir: Path) -> Path:
    result = subprocess.run(
        [find_bionetgen(), "run", "-i", str(MODEL), "-o", str(output_dir)],
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "BioNetGen run failed\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )

    gdat = output_dir / "two_compartment_pharmacokinetics_ode.gdat"
    if not gdat.exists():
        raise FileNotFoundError(f"Expected output not found: {gdat}")
    return gdat


def report_error(label: str, left: np.ndarray, right: np.ndarray) -> float:
    max_abs_error, max_rel_error = max_scaled_error(left, right)
    print(f"{label}_max_abs_error: {max_abs_error:.6g}")
    print(f"{label}_max_rel_error: {max_rel_error:.6g}")
    return max_rel_error


def main() -> None:
    parameters = read_numeric_parameters(MODEL)
    k12 = parameters["k12"]
    k21 = parameters["k21"]
    ke = parameters["ke"]
    central_0 = parameters["Central_0"]
    peripheral_0 = parameters["Peripheral_0"]

    with tempfile.TemporaryDirectory(prefix="two_compartment_pk_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_central = data[:, columns.index("Obs_Central")]
    bng_peripheral = data[:, columns.index("Obs_Peripheral")]
    bng_analytical_central = data[:, columns.index("Analytical_Central")]
    bng_analytical_peripheral = data[:, columns.index("Analytical_Peripheral")]

    exact_central, exact_peripheral = analytical_solution(t, k12, k21, ke, central_0, peripheral_0)

    max_rel_errors = [
        report_error("Obs_Central_vs_python", bng_central, exact_central),
        report_error("Obs_Peripheral_vs_python", bng_peripheral, exact_peripheral),
        report_error("Analytical_Central_vs_python", bng_analytical_central, exact_central),
        report_error("Analytical_Peripheral_vs_python", bng_analytical_peripheral, exact_peripheral),
        report_error("Obs_Central_vs_Analytical_Central", bng_central, bng_analytical_central),
        report_error("Obs_Peripheral_vs_Analytical_Peripheral", bng_peripheral, bng_analytical_peripheral),
    ]

    print(f"points: {len(t)}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
