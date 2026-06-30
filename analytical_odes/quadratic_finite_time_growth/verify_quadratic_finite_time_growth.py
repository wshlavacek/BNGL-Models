#!/usr/bin/env python3
"""Verify quadratic finite-time growth."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "quadratic_finite_time_growth.bngl"


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


def analytical_solution(t: np.ndarray, k: float, x0: float) -> np.ndarray:
    """Closed-form solution for dX/dt = k*X^2."""
    if x0 <= 0:
        raise ValueError("The verified finite-time growth branch expects X0 > 0")
    denominator = 1 - k * x0 * t
    if np.any(denominator <= 0):
        raise ValueError("Analytical solution is singular at or beyond t_blowup")
    return x0 / denominator


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

    gdat = output_dir / "quadratic_finite_time_growth_ode.gdat"
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
    k = parameters["k"]
    x0 = parameters["X0"]

    with tempfile.TemporaryDirectory(prefix="quadratic_growth_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_x = data[:, columns.index("Obs_X")]
    bng_analytical_x = data[:, columns.index("Analytical_X")]
    bng_denominator = data[:, columns.index("Analytical_Denominator")]
    bng_blowup_time = data[:, columns.index("Blowup_Time")]

    exact_x = analytical_solution(t, k, x0)
    exact_denominator = 1 - k * x0 * t
    exact_blowup_time = np.full_like(t, 1 / (k * x0))

    max_rel_errors = [
        report_error("Obs_X_vs_python", bng_x, exact_x),
        report_error("Analytical_X_vs_python", bng_analytical_x, exact_x),
        report_error("Obs_X_vs_Analytical_X", bng_x, bng_analytical_x),
        report_error("Analytical_Denominator_vs_python", bng_denominator, exact_denominator),
        report_error("Blowup_Time_vs_python", bng_blowup_time, exact_blowup_time),
    ]

    print(f"points: {len(t)}")
    print(f"final_Obs_X: {float(bng_x[-1]):.6g}")
    print(f"t_blowup: {float(exact_blowup_time[-1]):.6g}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
