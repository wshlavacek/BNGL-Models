#!/usr/bin/env python3
"""Verify the logistic growth model."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "logistic_growth.bngl"


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
    r: float,
    carrying_capacity: float,
    n0: float,
) -> np.ndarray:
    """Closed-form solution for dN/dt = r*N*(1 - N/K)."""
    if carrying_capacity <= 0 or n0 <= 0:
        raise ValueError("The logistic closed form requires K > 0 and N0 > 0")
    ratio = (carrying_capacity - n0) / n0
    return carrying_capacity / (1 + ratio * np.exp(-r * t))


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

    gdat = output_dir / "logistic_growth_ode.gdat"
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
    r = parameters["r"]
    carrying_capacity = parameters["K"]
    n0 = parameters["N0"]

    with tempfile.TemporaryDirectory(prefix="logistic_growth_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_n = data[:, columns.index("Obs_N")]
    bng_analytical_n = data[:, columns.index("Analytical_N")]
    bng_fraction = data[:, columns.index("Fraction_Of_K")]
    bng_inflection_n = data[:, columns.index("Inflection_N")]

    exact_n = analytical_solution(t, r, carrying_capacity, n0)
    exact_fraction = exact_n / carrying_capacity
    exact_inflection = np.full_like(t, carrying_capacity / 2)

    max_rel_errors = [
        report_error("Obs_N_vs_python", bng_n, exact_n),
        report_error("Analytical_N_vs_python", bng_analytical_n, exact_n),
        report_error("Obs_N_vs_Analytical_N", bng_n, bng_analytical_n),
        report_error("Fraction_Of_K_vs_python", bng_fraction, exact_fraction),
        report_error("Inflection_N_vs_python", bng_inflection_n, exact_inflection),
    ]

    print(f"points: {len(t)}")
    print(f"final_Obs_N: {float(bng_n[-1]):.6g}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
