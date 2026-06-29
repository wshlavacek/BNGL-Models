#!/usr/bin/env python3
"""Verify deterministic SIS epidemic threshold dynamics."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "sis_epidemic_threshold.bngl"


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
    beta: float,
    gamma: float,
    s0: float,
    i0: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Closed-form SIS solution for beta > gamma and I0 > 0."""
    if beta <= gamma or i0 <= 0:
        raise ValueError("This verifier expects beta > gamma and I0 > 0")

    total = s0 + i0
    r_net = beta - gamma
    i_star = total * r_net / beta
    ratio = (i_star - i0) / i0
    infected = i_star / (1 + ratio * np.exp(-r_net * t))
    susceptible = total - infected
    return susceptible, infected


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

    gdat = output_dir / "sis_epidemic_threshold_ode.gdat"
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
    beta = parameters["beta"]
    gamma = parameters["gamma"]
    s0 = parameters["S0"]
    i0 = parameters["I0"]

    with tempfile.TemporaryDirectory(prefix="sis_epidemic_threshold_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_s = data[:, columns.index("Obs_S")]
    bng_i = data[:, columns.index("Obs_I")]
    bng_analytical_s = data[:, columns.index("Analytical_S")]
    bng_analytical_i = data[:, columns.index("Analytical_I")]

    exact_s, exact_i = analytical_solution(t, beta, gamma, s0, i0)

    max_rel_errors = [
        report_error("Obs_S_vs_python", bng_s, exact_s),
        report_error("Obs_I_vs_python", bng_i, exact_i),
        report_error("Analytical_S_vs_python", bng_analytical_s, exact_s),
        report_error("Analytical_I_vs_python", bng_analytical_i, exact_i),
        report_error("Obs_S_vs_Analytical_S", bng_s, bng_analytical_s),
        report_error("Obs_I_vs_Analytical_I", bng_i, bng_analytical_i),
    ]

    print(f"points: {len(t)}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
