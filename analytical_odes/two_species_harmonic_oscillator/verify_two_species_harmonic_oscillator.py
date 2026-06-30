#!/usr/bin/env python3
"""Verify the two-species harmonic oscillator model."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "two_species_harmonic_oscillator.bngl"


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
    k2: float,
    kd: float,
    k4: float,
    k6: float,
    s1_0: float,
    s2_0: float,
) -> tuple[np.ndarray, np.ndarray, float, float, float, float, float]:
    """Closed-form solution for the constrained 2SHO system."""
    theta = math.sqrt(k2 * kd)
    s1_offset = (k6 - k4) / kd
    s2_offset = k4 / k2 - s1_offset
    s1_delta0 = s1_0 - s1_offset
    s2_delta0 = s2_0 - s2_offset
    s1_sin_coeff = (k2 * s1_delta0 + k2 * s2_delta0) / theta
    s2_sin_coeff = (-(k2 + kd) * s1_delta0 - k2 * s2_delta0) / theta

    cos_term = np.cos(theta * t)
    sin_term = np.sin(theta * t)
    s1 = s1_offset + s1_delta0 * cos_term + s1_sin_coeff * sin_term
    s2 = s2_offset + s2_delta0 * cos_term + s2_sin_coeff * sin_term
    s1_amplitude = math.hypot(s1_delta0, s1_sin_coeff)
    s2_amplitude = math.hypot(s2_delta0, s2_sin_coeff)
    return s1, s2, theta, s1_amplitude, s2_amplitude, s1_offset, s2_offset


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

    gdat = output_dir / "two_species_harmonic_oscillator_ode.gdat"
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
    k2 = parameters["k2"]
    kd = parameters["kd"]
    k4 = parameters["k4"]
    k6 = parameters["k6"]
    s1_0 = parameters["S1_0"]
    s2_0 = parameters["S2_0"]

    with tempfile.TemporaryDirectory(prefix="two_species_harmonic_oscillator_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_s1 = data[:, columns.index("Obs_S1")]
    bng_s2 = data[:, columns.index("Obs_S2")]
    bng_analytical_s1 = data[:, columns.index("Analytical_S1")]
    bng_analytical_s2 = data[:, columns.index("Analytical_S2")]
    bng_frequency = data[:, columns.index("Frequency")]
    bng_amplitude_s1 = data[:, columns.index("Amplitude_S1")]
    bng_amplitude_s2 = data[:, columns.index("Amplitude_S2")]
    bng_offset_s1 = data[:, columns.index("Offset_S1")]
    bng_offset_s2 = data[:, columns.index("Offset_S2")]
    bng_minimum_s1 = data[:, columns.index("Minimum_S1")]
    bng_minimum_s2 = data[:, columns.index("Minimum_S2")]
    bng_det_a = data[:, columns.index("Det_A")]

    (
        exact_s1,
        exact_s2,
        theta,
        s1_amplitude,
        s2_amplitude,
        s1_offset,
        s2_offset,
    ) = analytical_solution(t, k2, kd, k4, k6, s1_0, s2_0)

    expected_frequency = np.full_like(t, theta)
    expected_amplitude_s1 = np.full_like(t, s1_amplitude)
    expected_amplitude_s2 = np.full_like(t, s2_amplitude)
    expected_offset_s1 = np.full_like(t, s1_offset)
    expected_offset_s2 = np.full_like(t, s2_offset)
    expected_minimum_s1 = np.full_like(t, s1_offset - s1_amplitude)
    expected_minimum_s2 = np.full_like(t, s2_offset - s2_amplitude)
    expected_det_a = np.full_like(t, theta * theta)

    max_rel_errors = [
        report_error("Obs_S1_vs_python", bng_s1, exact_s1),
        report_error("Obs_S2_vs_python", bng_s2, exact_s2),
        report_error("Analytical_S1_vs_python", bng_analytical_s1, exact_s1),
        report_error("Analytical_S2_vs_python", bng_analytical_s2, exact_s2),
        report_error("Obs_S1_vs_Analytical_S1", bng_s1, bng_analytical_s1),
        report_error("Obs_S2_vs_Analytical_S2", bng_s2, bng_analytical_s2),
        report_error("Frequency_vs_python", bng_frequency, expected_frequency),
        report_error("Amplitude_S1_vs_python", bng_amplitude_s1, expected_amplitude_s1),
        report_error("Amplitude_S2_vs_python", bng_amplitude_s2, expected_amplitude_s2),
        report_error("Offset_S1_vs_python", bng_offset_s1, expected_offset_s1),
        report_error("Offset_S2_vs_python", bng_offset_s2, expected_offset_s2),
        report_error("Minimum_S1_vs_python", bng_minimum_s1, expected_minimum_s1),
        report_error("Minimum_S2_vs_python", bng_minimum_s2, expected_minimum_s2),
        report_error("Det_A_vs_python", bng_det_a, expected_det_a),
    ]

    print(f"points: {len(t)}")
    print(f"min_Obs_S1: {float(np.min(bng_s1)):.6g}")
    print(f"min_Obs_S2: {float(np.min(bng_s2)):.6g}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
