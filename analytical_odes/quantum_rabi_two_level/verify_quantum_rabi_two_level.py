#!/usr/bin/env python3
"""Verify the quantum two-level Rabi oscillation model."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "quantum_rabi_two_level.bngl"


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
    omega: float,
    u0: float,
    v0: float,
    w0: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Closed-form resonant Rabi solution in Bloch coordinates."""
    u = np.full_like(t, u0, dtype=float)
    v = v0 * np.cos(omega * t) - w0 * np.sin(omega * t)
    w = v0 * np.sin(omega * t) + w0 * np.cos(omega * t)
    p_excited = (1 + w) / 2
    p_ground = (1 - w) / 2
    return u, v, w, p_excited, p_ground


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

    gdat = output_dir / "quantum_rabi_two_level_ode.gdat"
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
    omega = parameters["Omega"]
    u0 = parameters["U0"]
    v0 = parameters["V0"]
    w0 = parameters["W0"]

    with tempfile.TemporaryDirectory(prefix="quantum_rabi_two_level_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_u = data[:, columns.index("Obs_U")]
    bng_v = data[:, columns.index("Obs_V")]
    bng_w = data[:, columns.index("Obs_W")]
    bng_excited = data[:, columns.index("P_Excited")]
    bng_ground = data[:, columns.index("P_Ground")]
    bng_norm2 = data[:, columns.index("Bloch_Norm2")]
    bng_analytical_u = data[:, columns.index("Analytical_U")]
    bng_analytical_v = data[:, columns.index("Analytical_V")]
    bng_analytical_w = data[:, columns.index("Analytical_W")]
    bng_analytical_excited = data[:, columns.index("Analytical_P_Excited")]
    bng_analytical_ground = data[:, columns.index("Analytical_P_Ground")]

    exact_u, exact_v, exact_w, exact_excited, exact_ground = analytical_solution(
        t,
        omega,
        u0,
        v0,
        w0,
    )
    exact_norm2 = exact_u**2 + exact_v**2 + exact_w**2

    max_rel_errors = [
        report_error("Obs_U_vs_python", bng_u, exact_u),
        report_error("Obs_V_vs_python", bng_v, exact_v),
        report_error("Obs_W_vs_python", bng_w, exact_w),
        report_error("P_Excited_vs_python", bng_excited, exact_excited),
        report_error("P_Ground_vs_python", bng_ground, exact_ground),
        report_error("Bloch_Norm2_vs_python", bng_norm2, exact_norm2),
        report_error("Analytical_U_vs_python", bng_analytical_u, exact_u),
        report_error("Analytical_V_vs_python", bng_analytical_v, exact_v),
        report_error("Analytical_W_vs_python", bng_analytical_w, exact_w),
        report_error("Analytical_P_Excited_vs_python", bng_analytical_excited, exact_excited),
        report_error("Analytical_P_Ground_vs_python", bng_analytical_ground, exact_ground),
        report_error("Obs_V_vs_Analytical_V", bng_v, bng_analytical_v),
        report_error("Obs_W_vs_Analytical_W", bng_w, bng_analytical_w),
        report_error("P_Excited_vs_Analytical_P_Excited", bng_excited, bng_analytical_excited),
        report_error("P_Ground_vs_Analytical_P_Ground", bng_ground, bng_analytical_ground),
    ]

    print(f"points: {len(t)}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
