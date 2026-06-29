#!/usr/bin/env python3
"""Verify birth-death-immigration moment dynamics."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "birth_death_immigration_moments.bngl"


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
    immigration: float,
    birth_rate: float,
    death_rate: float,
    mean0: float,
    variance0: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Closed-form mean and variance for birth-death-immigration moments."""
    a_net = birth_rate - death_rate
    if a_net == 0:
        raise ValueError("This verifier expects birth_rate != death_rate")

    c_event = birth_rate + death_rate
    mean_ss = -immigration / a_net
    mean_delta = mean0 - mean_ss
    variance_const = -(immigration + c_event * mean_ss) / (2 * a_net)
    variance_exp1 = -(c_event * mean_delta) / a_net
    variance_exp2 = variance0 - variance_const - variance_exp1

    mean = mean_ss + mean_delta * np.exp(a_net * t)
    variance = (
        variance_const
        + variance_exp1 * np.exp(a_net * t)
        + variance_exp2 * np.exp(2 * a_net * t)
    )
    return mean, variance


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

    gdat = output_dir / "birth_death_immigration_moments_ode.gdat"
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
    immigration = parameters["s"]
    birth_rate = parameters["birth_rate"]
    death_rate = parameters["death_rate"]
    mean0 = parameters["Mean0"]
    variance0 = parameters["Variance0"]

    with tempfile.TemporaryDirectory(prefix="birth_death_immigration_moments_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_mean = data[:, columns.index("Obs_Mean")]
    bng_variance = data[:, columns.index("Obs_Variance")]
    bng_analytical_mean = data[:, columns.index("Analytical_Mean")]
    bng_analytical_variance = data[:, columns.index("Analytical_Variance")]

    exact_mean, exact_variance = analytical_solution(
        t,
        immigration,
        birth_rate,
        death_rate,
        mean0,
        variance0,
    )

    max_rel_errors = [
        report_error("Obs_Mean_vs_python", bng_mean, exact_mean),
        report_error("Obs_Variance_vs_python", bng_variance, exact_variance),
        report_error("Analytical_Mean_vs_python", bng_analytical_mean, exact_mean),
        report_error(
            "Analytical_Variance_vs_python",
            bng_analytical_variance,
            exact_variance,
        ),
        report_error("Obs_Mean_vs_Analytical_Mean", bng_mean, bng_analytical_mean),
        report_error(
            "Obs_Variance_vs_Analytical_Variance",
            bng_variance,
            bng_analytical_variance,
        ),
    ]

    print(f"points: {len(t)}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
