#!/usr/bin/env python3
"""Verify the piecewise-linear gene switch trajectory."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "piecewise_linear_gene_switch.bngl"


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
    theta_b_low: float,
    kappa_a: float,
    gamma_a: float,
    gamma_b: float,
    a0: float,
    b0: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Piecewise solution for the default one-switch trajectory."""
    a_on = kappa_a / gamma_a
    t_switch = math.log(b0 / theta_b_low) / gamma_b
    a_switch = a_on + (a0 - a_on) * math.exp(-gamma_a * t_switch)

    before = t < t_switch
    after_time = t - t_switch
    a = np.where(
        before,
        a_on + (a0 - a_on) * np.exp(-gamma_a * t),
        a_switch * np.exp(-gamma_a * after_time),
    )
    b = np.where(
        before,
        b0 * np.exp(-gamma_b * t),
        theta_b_low * np.exp(-gamma_b * after_time),
    )
    return a, b, t_switch, a_switch


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

    gdat = output_dir / "piecewise_linear_gene_switch_ode.gdat"
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
    theta_b_low = parameters["theta_B_low"]
    kappa_a = parameters["kappa_A"]
    gamma_a = parameters["gamma_A"]
    gamma_b = parameters["gamma_B"]
    a0 = parameters["A0"]
    b0 = parameters["B0"]

    with tempfile.TemporaryDirectory(prefix="piecewise_gene_switch_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_a = data[:, columns.index("Obs_A")]
    bng_b = data[:, columns.index("Obs_B")]
    bng_analytical_a = data[:, columns.index("Analytical_A")]
    bng_analytical_b = data[:, columns.index("Analytical_B")]
    bng_switch_time = data[:, columns.index("Switch_Time")]
    bng_switch_a = data[:, columns.index("Switch_A")]

    exact_a, exact_b, t_switch, a_switch = analytical_solution(
        t,
        theta_b_low,
        kappa_a,
        gamma_a,
        gamma_b,
        a0,
        b0,
    )
    exact_switch_time = np.full_like(t, t_switch)
    exact_switch_a = np.full_like(t, a_switch)

    max_rel_errors = [
        report_error("Obs_A_vs_python", bng_a, exact_a),
        report_error("Obs_B_vs_python", bng_b, exact_b),
        report_error("Analytical_A_vs_python", bng_analytical_a, exact_a),
        report_error("Analytical_B_vs_python", bng_analytical_b, exact_b),
        report_error("Obs_A_vs_Analytical_A", bng_a, bng_analytical_a),
        report_error("Obs_B_vs_Analytical_B", bng_b, bng_analytical_b),
        report_error("Switch_Time_vs_python", bng_switch_time, exact_switch_time),
        report_error("Switch_A_vs_python", bng_switch_a, exact_switch_a),
    ]

    print(f"points: {len(t)}")
    print(f"switch_time: {t_switch:.6g}")
    print(f"switch_A: {a_switch:.6g}")
    print(f"final_Obs_A: {float(bng_a[-1]):.6g}")
    print(f"final_Obs_B: {float(bng_b[-1]):.6g}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
