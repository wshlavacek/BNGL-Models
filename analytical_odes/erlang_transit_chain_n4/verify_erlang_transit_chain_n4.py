#!/usr/bin/env python3
"""Verify the N=4 Erlang transit-chain delay model."""

from __future__ import annotations

import math
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
MODEL = HERE / "erlang_transit_chain_n4.bngl"


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
    dose: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Closed-form solution for a four-stage equal-rate transit chain."""
    tau = k * t
    decay = np.exp(-tau)
    t0 = dose * decay
    t1 = dose * decay * tau
    t2 = dose * decay * tau**2 / 2
    t3 = dose * decay * tau**3 / 6
    product = dose * (1 - decay * (1 + tau + tau**2 / 2 + tau**3 / 6))
    return t0, t1, t2, t3, product


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

    gdat = output_dir / "erlang_transit_chain_n4_ode.gdat"
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
    dose = parameters["Dose"]

    with tempfile.TemporaryDirectory(prefix="erlang_transit_chain_n4_") as tmp:
        gdat = run_bionetgen(Path(tmp))
        columns, data = parse_gdat(gdat)

    t = data[:, columns.index("time")]
    bng_t0 = data[:, columns.index("Obs_T0")]
    bng_t1 = data[:, columns.index("Obs_T1")]
    bng_t2 = data[:, columns.index("Obs_T2")]
    bng_t3 = data[:, columns.index("Obs_T3")]
    bng_product = data[:, columns.index("Obs_Product")]
    bng_analytical_t0 = data[:, columns.index("Analytical_T0")]
    bng_analytical_t1 = data[:, columns.index("Analytical_T1")]
    bng_analytical_t2 = data[:, columns.index("Analytical_T2")]
    bng_analytical_t3 = data[:, columns.index("Analytical_T3")]
    bng_analytical_product = data[:, columns.index("Analytical_Product")]

    exact_t0, exact_t1, exact_t2, exact_t3, exact_product = analytical_solution(t, k, dose)

    max_rel_errors = [
        report_error("Obs_T0_vs_python", bng_t0, exact_t0),
        report_error("Obs_T1_vs_python", bng_t1, exact_t1),
        report_error("Obs_T2_vs_python", bng_t2, exact_t2),
        report_error("Obs_T3_vs_python", bng_t3, exact_t3),
        report_error("Obs_Product_vs_python", bng_product, exact_product),
        report_error("Analytical_T0_vs_python", bng_analytical_t0, exact_t0),
        report_error("Analytical_T1_vs_python", bng_analytical_t1, exact_t1),
        report_error("Analytical_T2_vs_python", bng_analytical_t2, exact_t2),
        report_error("Analytical_T3_vs_python", bng_analytical_t3, exact_t3),
        report_error("Analytical_Product_vs_python", bng_analytical_product, exact_product),
        report_error("Obs_T0_vs_Analytical_T0", bng_t0, bng_analytical_t0),
        report_error("Obs_T1_vs_Analytical_T1", bng_t1, bng_analytical_t1),
        report_error("Obs_T2_vs_Analytical_T2", bng_t2, bng_analytical_t2),
        report_error("Obs_T3_vs_Analytical_T3", bng_t3, bng_analytical_t3),
        report_error("Obs_Product_vs_Analytical_Product", bng_product, bng_analytical_product),
    ]

    print(f"points: {len(t)}")
    max_reported_rel_error = max(max_rel_errors)
    if not math.isfinite(max_reported_rel_error) or max_reported_rel_error > 1.0e-6:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
