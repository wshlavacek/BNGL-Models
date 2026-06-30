#!/usr/bin/env python3
"""Verify the N=2 transit absorption pharmacokinetic model."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "transit_absorption_pk_n2.bngl"


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
    dose: float,
    f_bio: float,
    k_transit: float,
    k_abs: float,
    k_elim: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Closed-form solution for n=2 transit absorption after one bolus."""
    effective_dose = f_bio * dose
    transit1 = effective_dose * np.exp(-k_transit * t)
    transit2 = effective_dose * k_transit * t * np.exp(-k_transit * t)

    d_transit_abs = k_transit - k_abs
    absorption = effective_dose * k_transit * k_transit / (d_transit_abs**2) * (
        np.exp(-k_abs * t) - np.exp(-k_transit * t) * (1 + d_transit_abs * t)
    )

    central_scale = effective_dose * k_transit * k_transit * k_abs
    coef_transit_exp = -(
        k_abs + k_elim - 2 * k_transit
    ) / ((k_abs - k_transit) ** 2 * (k_elim - k_transit) ** 2)
    coef_transit_t = 1 / ((k_abs - k_transit) * (k_elim - k_transit))
    coef_abs = 1 / ((k_transit - k_abs) ** 2 * (k_elim - k_abs))
    coef_elim = 1 / ((k_transit - k_elim) ** 2 * (k_abs - k_elim))
    central = central_scale * (
        coef_transit_exp * np.exp(-k_transit * t)
        + coef_transit_t * t * np.exp(-k_transit * t)
        + coef_abs * np.exp(-k_abs * t)
        + coef_elim * np.exp(-k_elim * t)
    )
    return transit1, transit2, absorption, central


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

    gdat = output_dir / "transit_absorption_pk_n2_ode.gdat"
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
    dose = parameters["Dose"]
    f_bio = parameters["F_bio"]
    k_transit = parameters["k_transit"]
    k_abs = parameters["k_abs"]
    k_elim = parameters["k_elim"]

    with tempfile.TemporaryDirectory(prefix="transit_absorption_pk_n2_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_transit1 = data[:, columns.index("Obs_Transit1")]
    bng_transit2 = data[:, columns.index("Obs_Transit2")]
    bng_absorption = data[:, columns.index("Obs_Absorption")]
    bng_central = data[:, columns.index("Obs_Central")]
    bng_analytical_transit1 = data[:, columns.index("Analytical_Transit1")]
    bng_analytical_transit2 = data[:, columns.index("Analytical_Transit2")]
    bng_analytical_absorption = data[:, columns.index("Analytical_Absorption")]
    bng_analytical_central = data[:, columns.index("Analytical_Central")]

    exact_transit1, exact_transit2, exact_absorption, exact_central = analytical_solution(
        t,
        dose,
        f_bio,
        k_transit,
        k_abs,
        k_elim,
    )

    max_rel_errors = [
        report_error("Obs_Transit1_vs_python", bng_transit1, exact_transit1),
        report_error("Obs_Transit2_vs_python", bng_transit2, exact_transit2),
        report_error("Obs_Absorption_vs_python", bng_absorption, exact_absorption),
        report_error("Obs_Central_vs_python", bng_central, exact_central),
        report_error("Analytical_Transit1_vs_python", bng_analytical_transit1, exact_transit1),
        report_error("Analytical_Transit2_vs_python", bng_analytical_transit2, exact_transit2),
        report_error("Analytical_Absorption_vs_python", bng_analytical_absorption, exact_absorption),
        report_error("Analytical_Central_vs_python", bng_analytical_central, exact_central),
        report_error("Obs_Transit1_vs_Analytical_Transit1", bng_transit1, bng_analytical_transit1),
        report_error("Obs_Transit2_vs_Analytical_Transit2", bng_transit2, bng_analytical_transit2),
        report_error("Obs_Absorption_vs_Analytical_Absorption", bng_absorption, bng_analytical_absorption),
        report_error("Obs_Central_vs_Analytical_Central", bng_central, bng_analytical_central),
    ]

    print(f"points: {len(t)}")
    print(f"peak_Obs_Central: {float(np.max(bng_central)):.6g}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
