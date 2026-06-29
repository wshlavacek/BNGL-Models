#!/usr/bin/env python3
"""Verify non-Lipschitz square-root growth."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "non_lipschitz_square_root_growth.bngl"


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


def analytical_solution(t: np.ndarray, a: float, x0: float) -> np.ndarray:
    """Unique positive-branch solution for dX/dt = a*sqrt(X), X0 > 0."""
    if x0 <= 0:
        raise ValueError("The verified branch expects X0 > 0")
    return (math.sqrt(x0) + a * t / 2) ** 2


def zero_initial_delayed_solution(t: np.ndarray, a: float, delay_tau: float) -> np.ndarray:
    """One delayed-start solution for the zero-initial-condition problem."""
    return np.where(t < delay_tau, 0.0, (a * (t - delay_tau) / 2) ** 2)


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

    gdat = output_dir / "non_lipschitz_square_root_growth_ode.gdat"
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
    a = parameters["a"]
    x0 = parameters["X0"]
    delay_tau = parameters["delay_tau"]

    with tempfile.TemporaryDirectory(prefix="non_lipschitz_growth_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_x = data[:, columns.index("Obs_X")]
    bng_analytical_x = data[:, columns.index("Analytical_X")]
    bng_stationary = data[:, columns.index("ZeroInitial_Stationary_X")]
    bng_delayed = data[:, columns.index("ZeroInitial_Delayed_X")]

    exact_x = analytical_solution(t, a, x0)
    exact_stationary = np.zeros_like(t)
    exact_delayed = zero_initial_delayed_solution(t, a, delay_tau)

    max_rel_errors = [
        report_error("Obs_X_vs_python", bng_x, exact_x),
        report_error("Analytical_X_vs_python", bng_analytical_x, exact_x),
        report_error("Obs_X_vs_Analytical_X", bng_x, bng_analytical_x),
        report_error("ZeroInitial_Stationary_X_vs_python", bng_stationary, exact_stationary),
        report_error("ZeroInitial_Delayed_X_vs_python", bng_delayed, exact_delayed),
    ]

    print(f"points: {len(t)}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
