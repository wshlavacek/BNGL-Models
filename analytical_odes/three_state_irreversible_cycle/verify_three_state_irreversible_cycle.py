#!/usr/bin/env python3
"""Verify the three-state irreversible cycle model."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "three_state_irreversible_cycle.bngl"


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
    k: float,
    x1_0: float,
    x2_0: float,
    x3_0: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Closed-form solution for X1 -> X2 -> X3 -> X1 with equal rates."""
    total = x1_0 + x2_0 + x3_0
    x_ss = total / 3
    gamma = 3 * k / 2
    omega = math.sqrt(3) * k / 2

    u1 = x1_0 - x_ss
    u2 = x2_0 - x_ss
    u3 = x3_0 - x_ss
    q1 = (u3 - u2) / math.sqrt(3)
    q2 = (u1 - u3) / math.sqrt(3)
    q3 = (u2 - u1) / math.sqrt(3)

    decay = np.exp(-gamma * t)
    cos_term = np.cos(omega * t)
    sin_term = np.sin(omega * t)

    x1 = x_ss + decay * (u1 * cos_term + q1 * sin_term)
    x2 = x_ss + decay * (u2 * cos_term + q2 * sin_term)
    x3 = x_ss + decay * (u3 * cos_term + q3 * sin_term)
    return x1, x2, x3


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

    gdat = output_dir / "three_state_irreversible_cycle_ode.gdat"
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
    x1_0 = parameters["X1_0"]
    x2_0 = parameters["X2_0"]
    x3_0 = parameters["X3_0"]

    with tempfile.TemporaryDirectory(prefix="three_state_cycle_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_x1 = data[:, columns.index("Obs_X1")]
    bng_x2 = data[:, columns.index("Obs_X2")]
    bng_x3 = data[:, columns.index("Obs_X3")]
    bng_analytical_x1 = data[:, columns.index("Analytical_X1")]
    bng_analytical_x2 = data[:, columns.index("Analytical_X2")]
    bng_analytical_x3 = data[:, columns.index("Analytical_X3")]

    exact_x1, exact_x2, exact_x3 = analytical_solution(t, k, x1_0, x2_0, x3_0)

    max_rel_errors = [
        report_error("Obs_X1_vs_python", bng_x1, exact_x1),
        report_error("Obs_X2_vs_python", bng_x2, exact_x2),
        report_error("Obs_X3_vs_python", bng_x3, exact_x3),
        report_error("Analytical_X1_vs_python", bng_analytical_x1, exact_x1),
        report_error("Analytical_X2_vs_python", bng_analytical_x2, exact_x2),
        report_error("Analytical_X3_vs_python", bng_analytical_x3, exact_x3),
        report_error("Obs_X1_vs_Analytical_X1", bng_x1, bng_analytical_x1),
        report_error("Obs_X2_vs_Analytical_X2", bng_x2, bng_analytical_x2),
        report_error("Obs_X3_vs_Analytical_X3", bng_x3, bng_analytical_x3),
    ]

    print(f"points: {len(t)}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
