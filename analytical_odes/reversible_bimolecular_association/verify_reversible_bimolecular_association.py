#!/usr/bin/env python3
"""Verify reversible bimolecular association against the analytical solution."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "reversible_bimolecular_association.bngl"


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
    kf: float,
    kr: float,
    a0: float,
    b0: float,
    c0: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Closed-form solution for A + B <-> C."""
    if kf <= 0:
        raise ValueError("The Riccati closed form used here requires kf > 0")

    atot = a0 + c0
    btot = b0 + c0
    s_assoc = kf * (atot + btot) + kr
    d_assoc = math.sqrt(s_assoc * s_assoc - 4 * kf * kf * atot * btot)
    c_low = (s_assoc - d_assoc) / (2 * kf)
    c_high = (s_assoc + d_assoc) / (2 * kf)
    r_init = (c0 - c_low) / (c0 - c_high)
    r_t = r_init * np.exp(-d_assoc * t)

    c_t = (c_low - r_t * c_high) / (1 - r_t)
    a_t = atot - c_t
    b_t = btot - c_t
    return a_t, b_t, c_t


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

    gdat = output_dir / "reversible_bimolecular_association_ode.gdat"
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
    kf = parameters["kf"]
    kr = parameters["kr"]
    a0 = parameters["A0"]
    b0 = parameters["B0"]
    c0 = parameters["C0"]

    with tempfile.TemporaryDirectory(prefix="reversible_assoc_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_a = data[:, columns.index("Obs_A")]
    bng_b = data[:, columns.index("Obs_B")]
    bng_c = data[:, columns.index("Obs_C")]
    bng_analytical_a = data[:, columns.index("Analytical_A")]
    bng_analytical_b = data[:, columns.index("Analytical_B")]
    bng_analytical_c = data[:, columns.index("Analytical_C")]

    exact_a, exact_b, exact_c = analytical_solution(t, kf, kr, a0, b0, c0)

    max_rel_errors = [
        report_error("Obs_A_vs_python", bng_a, exact_a),
        report_error("Obs_B_vs_python", bng_b, exact_b),
        report_error("Obs_C_vs_python", bng_c, exact_c),
        report_error("Analytical_A_vs_python", bng_analytical_a, exact_a),
        report_error("Analytical_B_vs_python", bng_analytical_b, exact_b),
        report_error("Analytical_C_vs_python", bng_analytical_c, exact_c),
        report_error("Obs_A_vs_Analytical_A", bng_a, bng_analytical_a),
        report_error("Obs_B_vs_Analytical_B", bng_b, bng_analytical_b),
        report_error("Obs_C_vs_Analytical_C", bng_c, bng_analytical_c),
    ]

    print(f"points: {len(t)}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
