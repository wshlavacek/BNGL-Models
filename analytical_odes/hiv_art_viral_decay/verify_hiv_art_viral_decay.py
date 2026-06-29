#!/usr/bin/env python3
"""Verify HIV ART viral decay BNGL output against the analytical solution."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
MODEL = HERE / "hiv_art_viral_decay.bngl"


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


def analytical_solution(t: np.ndarray, p_pre: float, p_art: float, c: float) -> np.ndarray:
    """Closed-form solution for dV/dt = P_art - c*V and V(0) = P_pre/c."""
    post_treatment_steady_state = p_art / c
    initial_value = p_pre / c
    return post_treatment_steady_state + (
        initial_value - post_treatment_steady_state
    ) * np.exp(-c * t)


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


def run_bionetgen(output_dir: Path) -> Path:
    executable = shutil.which("bionetgen")
    if executable is None:
        raise RuntimeError("BioNetGen executable 'bionetgen' was not found on PATH")

    result = subprocess.run(
        [executable, "run", "-i", str(MODEL), "-o", str(output_dir)],
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

    gdat = output_dir / "hiv_art_viral_decay_ode.gdat"
    if not gdat.exists():
        raise FileNotFoundError(f"Expected output not found: {gdat}")
    return gdat


def main() -> None:
    parameters = read_numeric_parameters(MODEL)
    p_pre = parameters["P_pre"]
    p_art = parameters["P_art"]
    c = parameters["c"]

    with tempfile.TemporaryDirectory(prefix="hiv_art_decay_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    time_idx = columns.index("time")
    viral_idx = columns.index("Obs_V")
    bng_analytical_idx = columns.index("Analytical_V")

    t = data[:, time_idx]
    bng_v = data[:, viral_idx]
    bng_analytical_v = data[:, bng_analytical_idx]
    exact_v = analytical_solution(t, p_pre, p_art, c)

    max_abs_error, max_rel_error = max_scaled_error(bng_v, exact_v)
    max_func_abs_error, max_func_rel_error = max_scaled_error(bng_analytical_v, exact_v)
    max_bng_abs_error, max_bng_rel_error = max_scaled_error(bng_v, bng_analytical_v)

    print(f"points: {len(t)}")
    print(f"Obs_V_vs_python_max_abs_error: {max_abs_error:.6g}")
    print(f"Obs_V_vs_python_max_rel_error: {max_rel_error:.6g}")
    print(f"Analytical_V_vs_python_max_abs_error: {max_func_abs_error:.6g}")
    print(f"Analytical_V_vs_python_max_rel_error: {max_func_rel_error:.6g}")
    print(f"Obs_V_vs_Analytical_V_max_abs_error: {max_bng_abs_error:.6g}")
    print(f"Obs_V_vs_Analytical_V_max_rel_error: {max_bng_rel_error:.6g}")

    max_reported_rel_error = max(max_rel_error, max_func_rel_error, max_bng_rel_error)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
