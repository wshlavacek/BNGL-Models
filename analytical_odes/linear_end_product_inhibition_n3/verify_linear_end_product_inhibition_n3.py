#!/usr/bin/env python3
"""Verify the N=3 linear end-product inhibition model."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "linear_end_product_inhibition_n3.bngl"


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
    x0: float,
    k0: float,
    k: float,
    h: float,
    x1_0: float,
    x2_0: float,
    x3_0: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Closed-form solution for the equal-rate N=3 feedback pathway."""
    if k <= 0 or h <= 0:
        raise ValueError("The closed form used here requires k > 0 and h > 0")

    j = k0 * x0
    x_ss = j / (k + h)
    alpha = (h * k * k) ** (1 / 3)
    omega = math.sqrt(3) * alpha / 2

    u1 = x1_0 - x_ss
    u2 = x2_0 - x_ss
    u3 = x3_0 - x_ss

    a_mode = (k * k * u1 - alpha * k * u2 + alpha * alpha * u3) / (
        3 * alpha * alpha
    )
    b_mode = u3 - a_mode
    c_mode = (2 * k * u2 / alpha - u3 + 3 * a_mode) / math.sqrt(3)

    exp_alpha = np.exp(alpha * t / 2)
    cos_term = np.cos(omega * t)
    sin_term = np.sin(omega * t)
    trig_sum = b_mode * cos_term + c_mode * sin_term
    trig_deriv = -b_mode * sin_term + c_mode * cos_term

    w = a_mode * np.exp(-alpha * t) + exp_alpha * trig_sum
    wp = (
        -alpha * a_mode * np.exp(-alpha * t)
        + exp_alpha * ((alpha / 2) * trig_sum + omega * trig_deriv)
    )
    wpp = (
        alpha * alpha * a_mode * np.exp(-alpha * t)
        + exp_alpha * ((-alpha * alpha / 2) * trig_sum + alpha * omega * trig_deriv)
    )

    exp_decay = np.exp(-k * t)
    x1 = x_ss + exp_decay * wpp / (k * k)
    x2 = x_ss + exp_decay * wp / k
    x3 = x_ss + exp_decay * w
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

    gdat = output_dir / "linear_end_product_inhibition_n3_ode.gdat"
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
    x0 = parameters["x0"]
    k0 = parameters["k0"]
    k = parameters["k"]
    h = parameters["h"]
    x1_0 = parameters["X1_0"]
    x2_0 = parameters["X2_0"]
    x3_0 = parameters["X3_0"]

    with tempfile.TemporaryDirectory(prefix="linear_feedback_n3_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_x1 = data[:, columns.index("Obs_X1")]
    bng_x2 = data[:, columns.index("Obs_X2")]
    bng_x3 = data[:, columns.index("Obs_X3")]
    bng_analytical_x1 = data[:, columns.index("Analytical_X1")]
    bng_analytical_x2 = data[:, columns.index("Analytical_X2")]
    bng_analytical_x3 = data[:, columns.index("Analytical_X3")]

    exact_x1, exact_x2, exact_x3 = analytical_solution(t, x0, k0, k, h, x1_0, x2_0, x3_0)

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
