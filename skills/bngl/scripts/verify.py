#!/usr/bin/env python3
"""verify.py — Post-modification verification harness for BNGL models.

Compares simulation output from an original and modified model to verify
that dynamics are preserved after reformulation (unit conversion, format
conversion, compartment restructuring, etc.).

See skills/bngl/skill.md §8.7 for requirements.

Usage:
    python verify.py --original path/to/original.bngl \
                     --modified path/to/modified.bngl \
                     --mode ode --tolerance 0.01

    # With pre-computed output:
    python verify.py --original path/to/original_ode.gdat \
                     --modified path/to/modified_ode.gdat \
                     --tolerance 0.01

    # With explicit column mapping:
    python verify.py --original orig.bngl --modified mod.bngl \
                     --column-map "R1=R1_conc,R2=R2_conc"

    # SSA mode with ensemble comparison:
    python verify.py --original orig.bngl --modified mod.bngl \
                     --mode ssa --replicates 50
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

# ---------------------------------------------------------------------------
# Data I/O
# ---------------------------------------------------------------------------

def parse_data_file(path):
    """Parse a .gdat or .scan file (whitespace-delimited, # header row).

    Returns:
        columns: list of column names (including the independent variable)
        data: numpy array of shape (n_rows, n_cols)
    """
    columns = None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                # Header row: strip leading # and parse column names
                header = line.lstrip("#").strip()
                columns = header.split()
                continue
            values = line.split()
            rows.append([float(v) for v in values])
    if columns is None:
        raise ValueError(f"No header row found in {path}")
    if not rows:
        raise ValueError(f"No data rows found in {path}")
    data = np.array(rows)
    if data.shape[1] != len(columns):
        raise ValueError(
            f"Column count mismatch in {path}: header has {len(columns)} "
            f"columns but data has {data.shape[1]}"
        )
    return columns, data


# ---------------------------------------------------------------------------
# Simulation runner
# ---------------------------------------------------------------------------

def find_output_file(output_dir, model_basename, suffix, output_type):
    """Find the simulation output file in the output directory.

    BioNetGen produces files like: modelname_suffix.gdat or modelname_suffix.scan
    """
    expected = f"{model_basename}_{suffix}.{output_type}"
    candidate = os.path.join(output_dir, expected)
    if os.path.isfile(candidate):
        return candidate
    # Fallback: look for any matching file
    for fname in os.listdir(output_dir):
        if fname.endswith(f".{output_type}"):
            return os.path.join(output_dir, fname)
    return None


def run_bngl(model_path, output_dir, suffix="ode"):
    """Run a BNGL model via PyBioNetGen.

    Returns the path to the output directory.
    """
    cmd = ["bionetgen", "run", "-i", model_path, "-o", output_dir]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"bionetgen run failed for {model_path}:\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )
    return output_dir


def run_sbml(model_path, output_dir, suffix="ode"):
    """Stub: run an SBML model. Not yet implemented."""
    raise NotImplementedError(
        "SBML simulation is not yet implemented. "
        "Please provide a pre-computed .gdat or .scan file instead."
    )


def run_antimony(model_path, output_dir, suffix="ode"):
    """Stub: run an Antimony model. Not yet implemented."""
    raise NotImplementedError(
        "Antimony simulation is not yet implemented. "
        "Please provide a pre-computed .gdat or .scan file instead."
    )


RUNNERS = {
    "bngl": run_bngl,
    "sbml": run_sbml,
    "antimony": run_antimony,
}

FORMAT_EXTENSIONS = {
    ".bngl": "bngl",
    ".xml": "sbml",
    ".ant": "antimony",
}


def is_data_file(path):
    """Check if path is a pre-computed output file (.gdat or .scan)."""
    return path.endswith(".gdat") or path.endswith(".scan")


def detect_format(path):
    """Detect model format from file extension."""
    _, ext = os.path.splitext(path)
    fmt = FORMAT_EXTENSIONS.get(ext)
    if fmt is None:
        raise ValueError(
            f"Cannot detect format for {path}. "
            f"Supported extensions: {list(FORMAT_EXTENSIONS.keys())}"
        )
    return fmt


def simulate_model(model_path, fmt, output_type, sim_suffix="ode"):
    """Run simulation and return (columns, data) from the output file.

    If model_path is a .gdat or .scan file, parse it directly.
    Otherwise, run the appropriate simulator.
    """
    if is_data_file(model_path):
        return parse_data_file(model_path)

    runner = RUNNERS[fmt]
    model_basename = os.path.splitext(os.path.basename(model_path))[0]
    tmpdir = tempfile.mkdtemp(prefix="verify_")
    try:
        runner(model_path, tmpdir, suffix=sim_suffix)
        out_file = find_output_file(tmpdir, model_basename, sim_suffix, output_type)
        if out_file is None:
            # List what we got for debugging
            contents = os.listdir(tmpdir)
            raise FileNotFoundError(
                f"No .{output_type} file found in {tmpdir} after running "
                f"{model_path}. Directory contains: {contents}"
            )
        columns, data = parse_data_file(out_file)
        return columns, data
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Unit reconciliation: dimensional conversion functions
# ---------------------------------------------------------------------------

DIMENSIONAL_HEADER_RE = re.compile(
    r"^\s*#\s*Dimensional\s+conversion", re.IGNORECASE
)
COMMENTED_FUNC_RE = re.compile(
    r"^#\s+(\w+\(\))\s*=\s*(.+)$"
)


def detect_dimensional_functions(model_path):
    """Detect commented-out dimensional conversion functions in a BNGL model.

    Returns True if the model has commented-out dimensional conversion
    functions that should be uncommented for comparison.
    """
    in_functions = False
    in_dimensional = False
    with open(model_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("begin functions"):
                in_functions = True
                continue
            if stripped.startswith("end functions"):
                break
            if in_functions and DIMENSIONAL_HEADER_RE.match(stripped):
                in_dimensional = True
                continue
            if in_dimensional and COMMENTED_FUNC_RE.match(stripped):
                return True
    return False


def uncomment_dimensional_functions(model_path, tmp_dir):
    """Create a copy of the model with dimensional conversion functions uncommented.

    Returns the path to the temporary model file.
    """
    tmp_model = os.path.join(tmp_dir, os.path.basename(model_path))
    in_functions = False
    in_dimensional = False
    with open(model_path) as fin, open(tmp_model, "w") as fout:
        for line in fin:
            stripped = line.strip()
            if stripped.startswith("begin functions"):
                in_functions = True
                fout.write(line)
                continue
            if stripped.startswith("end functions"):
                in_functions = False
                in_dimensional = False
                fout.write(line)
                continue
            if in_functions and DIMENSIONAL_HEADER_RE.match(stripped):
                in_dimensional = True
                # Uncomment the header line itself
                fout.write(line.replace("# ", "  # ", 1))
                continue
            if in_dimensional and COMMENTED_FUNC_RE.match(stripped):
                # Uncomment: remove leading "# " to produce "  func() = ..."
                uncommented = re.sub(r"^#\s+", "  ", line)
                fout.write(uncommented)
                continue
            fout.write(line)
    return tmp_model


def get_dimensional_function_names(model_path):
    """Extract the names of dimensional conversion functions from a model.

    Returns a list of function names (without parentheses), e.g. ["R1_conc", "R2_conc"].
    """
    names = []
    in_functions = False
    in_dimensional = False
    with open(model_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("begin functions"):
                in_functions = True
                continue
            if stripped.startswith("end functions"):
                break
            if in_functions and DIMENSIONAL_HEADER_RE.match(stripped):
                in_dimensional = True
                continue
            if in_dimensional:
                m = COMMENTED_FUNC_RE.match(stripped)
                if m:
                    # e.g. "R1_conc()" -> "R1_conc"
                    names.append(m.group(1).rstrip("()"))
    return names


# ---------------------------------------------------------------------------
# Column mapping
# ---------------------------------------------------------------------------

def parse_column_map(map_str):
    """Parse a column map string like 'R1=R1_conc,R2=R2_conc'.

    Returns a dict mapping original column names to modified column names.
    """
    mapping = {}
    for pair in map_str.split(","):
        pair = pair.strip()
        if "=" not in pair:
            raise ValueError(f"Invalid column map entry: {pair!r} (expected 'orig=modified')")
        orig, mod = pair.split("=", 1)
        mapping[orig.strip()] = mod.strip()
    return mapping


def auto_match_columns(orig_cols, mod_cols):
    """Auto-match columns by exact name, suffix, or substring.

    Returns a dict mapping original column names to modified column names.
    Skips the first column (independent variable, usually 'time').
    """
    mapping = {}
    orig_data_cols = orig_cols[1:]  # skip independent variable
    mod_data_cols = mod_cols[1:]

    for orig_name in orig_data_cols:
        # 1. Exact match
        if orig_name in mod_data_cols:
            mapping[orig_name] = orig_name
            continue

        # 2. Suffix/prefix match: check if one name is a suffix/prefix of the other
        #    separated by underscores.
        candidates = []
        orig_parts = orig_name.split("_")
        for mod_name in mod_data_cols:
            mod_parts = mod_name.split("_")
            # e.g., "Obs_Tot_R1" ends with "R1", or "R1_conc" starts with "R1"
            if mod_name.endswith(f"_{orig_name}") or mod_name.startswith(f"{orig_name}_"):
                candidates.append(mod_name)
            elif orig_name.endswith(f"_{mod_name}") or orig_name.startswith(f"{mod_name}_"):
                candidates.append(mod_name)
            # 3. Shared trailing part: e.g., Obs_R1 and Obs_Tot_R1 both end with R1
            elif orig_parts[-1] == mod_parts[-1] and len(orig_parts[-1]) > 0:
                candidates.append(mod_name)

        if len(candidates) == 1:
            mapping[orig_name] = candidates[0]
        elif len(candidates) > 1:
            # Prefer exact suffix match (Obs_Tot_R1 for R1)
            for c in candidates:
                if c.endswith(f"_{orig_name}"):
                    mapping[orig_name] = c
                    break
            else:
                # Prefer the candidate sharing the most trailing parts
                def shared_tail(a, b):
                    pa, pb = a.split("_")[::-1], b.split("_")[::-1]
                    count = 0
                    for x, y in zip(pa, pb):
                        if x == y:
                            count += 1
                        else:
                            break
                    return count
                candidates.sort(key=lambda c: shared_tail(orig_name, c), reverse=True)
                mapping[orig_name] = candidates[0]
        # If no match found, leave unmapped (will be reported)

    return mapping


def resolve_column_mapping(orig_cols, mod_cols, explicit_map_str=None):
    """Build the final column mapping, with explicit overriding auto.

    Returns a dict mapping original column names to modified column names.
    """
    auto_map = auto_match_columns(orig_cols, mod_cols)

    if explicit_map_str:
        explicit_map = parse_column_map(explicit_map_str)
        auto_map.update(explicit_map)

    return auto_map


# ---------------------------------------------------------------------------
# Comparison metrics
# ---------------------------------------------------------------------------

def pointwise_errors(a, b, atol):
    """Compute pointwise absolute error and scale for mixed tolerance.

    The mixed tolerance criterion at each point is:
        |a_i - b_i| <= atol + rtol * |a_i|

    This is equivalent to numpy.allclose's criterion and naturally handles
    near-zero regions: when |a_i| is small, the absolute tolerance dominates,
    so we don't get spurious failures from dividing by ~0.

    Returns:
        abs_err: array of |a_i - b_i|
        scale: array of atol + |a_i| (the denominator for normalized error)
    """
    abs_err = np.abs(a - b)
    scale = atol + np.abs(a)
    return abs_err, scale


def normalized_max_error(a, b, atol):
    """Max of |a_i - b_i| / (atol + |a_i|).

    A value <= rtol means every point satisfies |a-b| <= atol + rtol*|a|.
    """
    abs_err, scale = pointwise_errors(a, b, atol)
    return np.max(abs_err / scale)


def normalized_l2_error(a, b, atol):
    """RMS of |a_i - b_i| / (atol + |a_i|).

    An aggregate measure of trajectory agreement that is well-behaved
    near zero thanks to the absolute tolerance floor.
    """
    abs_err, scale = pointwise_errors(a, b, atol)
    return np.sqrt(np.mean((abs_err / scale) ** 2))


def interpolate_to_common_grid(t_a, y_a, t_b, y_b):
    """Interpolate two trajectories onto a common time grid.

    Uses the denser grid as the common grid. Only compares over the
    overlapping time range.
    """
    t_min = max(t_a[0], t_b[0])
    t_max = min(t_a[-1], t_b[-1])

    if t_min >= t_max:
        raise ValueError(
            f"No overlapping time range: original [{t_a[0]}, {t_a[-1]}], "
            f"modified [{t_b[0]}, {t_b[-1]}]"
        )

    # Use the denser grid within the overlap
    t_common_a = t_a[(t_a >= t_min) & (t_a <= t_max)]
    t_common_b = t_b[(t_b >= t_min) & (t_b <= t_max)]
    t_common = t_common_a if len(t_common_a) >= len(t_common_b) else t_common_b

    y_a_interp = np.interp(t_common, t_a, y_a)
    y_b_interp = np.interp(t_common, t_b, y_b)

    return t_common, y_a_interp, y_b_interp


# ---------------------------------------------------------------------------
# ODE comparison
# ---------------------------------------------------------------------------

def compare_ode(orig_cols, orig_data, mod_cols, mod_data, column_map, rtol, atol):
    """Compare ODE trajectories column by column.

    Uses mixed tolerance: a point passes if |a-b| <= atol + rtol*|a|.
    Reports normalized_max (must be <= rtol for pass) and normalized_l2
    (aggregate measure, informational).

    Returns a list of result dicts, one per mapped column.
    """
    results = []
    orig_time = orig_data[:, 0]
    mod_time = mod_data[:, 0]

    for orig_name, mod_name in column_map.items():
        if orig_name not in orig_cols:
            results.append({
                "original_col": orig_name,
                "modified_col": mod_name,
                "status": "ERROR",
                "message": f"Column {orig_name!r} not found in original data",
            })
            continue
        if mod_name not in mod_cols:
            results.append({
                "original_col": orig_name,
                "modified_col": mod_name,
                "status": "ERROR",
                "message": f"Column {mod_name!r} not found in modified data",
            })
            continue

        orig_idx = orig_cols.index(orig_name)
        mod_idx = mod_cols.index(mod_name)

        orig_y = orig_data[:, orig_idx]
        mod_y = mod_data[:, mod_idx]

        # Interpolate to common grid if time vectors differ
        if len(orig_time) != len(mod_time) or not np.allclose(orig_time, mod_time):
            _, orig_y, mod_y = interpolate_to_common_grid(
                orig_time, orig_y, mod_time, mod_y
            )

        norm_max = normalized_max_error(orig_y, mod_y, atol)
        norm_l2 = normalized_l2_error(orig_y, mod_y, atol)
        passed = norm_max <= rtol

        results.append({
            "original_col": orig_name,
            "modified_col": mod_name,
            "status": "PASS" if passed else "FAIL",
            "norm_max": norm_max,
            "norm_l2": norm_l2,
            "rtol": rtol,
            "atol": atol,
        })

    return results


# ---------------------------------------------------------------------------
# SSA comparison (basic v1: ensemble means)
# ---------------------------------------------------------------------------

def run_ssa_ensemble(model_path, fmt, output_type, replicates, sim_suffix="ssa"):
    """Run SSA replicates and collect per-column means and variances.

    Returns (columns, mean_data, var_data) where mean_data and var_data
    are arrays of shape (n_timesteps, n_cols).
    """
    if is_data_file(model_path):
        raise ValueError("SSA ensemble mode requires a model file, not a pre-computed data file")

    runner = RUNNERS[fmt]
    model_basename = os.path.splitext(os.path.basename(model_path))[0]
    all_data = []
    columns = None

    for i in range(replicates):
        tmpdir = tempfile.mkdtemp(prefix=f"verify_ssa_{i}_")
        try:
            runner(model_path, tmpdir, suffix=sim_suffix)
            out_file = find_output_file(tmpdir, model_basename, sim_suffix, output_type)
            if out_file is None:
                raise FileNotFoundError(
                    f"No .{output_type} file found after SSA replicate {i}"
                )
            cols, data = parse_data_file(out_file)
            if columns is None:
                columns = cols
            all_data.append(data)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    # Stack and compute statistics
    # All replicates should have the same shape; use the minimum row count
    min_rows = min(d.shape[0] for d in all_data)
    stacked = np.stack([d[:min_rows, :] for d in all_data], axis=0)  # (replicates, rows, cols)
    mean_data = np.mean(stacked, axis=0)
    var_data = np.var(stacked, axis=0)

    return columns, mean_data, var_data


def compare_ssa(orig_cols, orig_mean, mod_cols, mod_mean, column_map, rtol, atol):
    """Compare SSA ensemble means (same logic as ODE comparison on means)."""
    return compare_ode(orig_cols, orig_mean, mod_cols, mod_mean, column_map, rtol, atol)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def generate_plots(orig_cols, orig_data, mod_cols, mod_data, column_map, plot_dir):
    """Generate comparison plots for each mapped column pair."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(plot_dir, exist_ok=True)
    orig_time = orig_data[:, 0]
    mod_time = mod_data[:, 0]

    for orig_name, mod_name in column_map.items():
        if orig_name not in orig_cols or mod_name not in mod_cols:
            continue

        orig_idx = orig_cols.index(orig_name)
        mod_idx = mod_cols.index(mod_name)

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(orig_time, orig_data[:, orig_idx], label=f"Original: {orig_name}", linewidth=2)
        ax.plot(mod_time, mod_data[:, mod_idx], label=f"Modified: {mod_name}",
                linewidth=2, linestyle="--")
        ax.set_xlabel(orig_cols[0])
        ax.set_ylabel("Value")
        ax.set_title(f"Verification: {orig_name} vs {mod_name}")
        ax.legend()
        ax.grid(True, alpha=0.3)

        safe_name = re.sub(r"[^\w]", "_", f"{orig_name}_vs_{mod_name}")
        fig.savefig(os.path.join(plot_dir, f"{safe_name}.png"), dpi=150, bbox_inches="tight")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def print_report(results, mode):
    """Print a pass/fail report to stdout."""
    print(f"\n{'='*60}")
    print(f"Verification Report (mode: {mode})")
    print(f"{'='*60}\n")

    all_passed = True
    for r in results:
        status = r["status"]
        orig = r["original_col"]
        mod = r["modified_col"]

        if status == "ERROR":
            print(f"  ERROR  {orig} -> {mod}: {r['message']}")
            all_passed = False
        else:
            norm_max = r["norm_max"]
            norm_l2 = r["norm_l2"]
            rtol = r["rtol"]
            atol = r["atol"]
            marker = "PASS" if status == "PASS" else "FAIL"
            print(f"  {marker}  {orig} -> {mod}")
            print(f"         norm_max={norm_max:.6e}  norm_l2={norm_l2:.6e}  "
                  f"rtol={rtol:.2e}  atol={atol:.2e}")
            if status == "FAIL":
                all_passed = False

    print(f"\n{'─'*60}")
    overall = "PASS" if all_passed else "FAIL"
    print(f"Overall: {overall}")
    print(f"{'─'*60}\n")

    return all_passed


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(
        description="Verify that a modified BNGL model reproduces the dynamics of the original.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--original", required=True,
                   help="Path to original model (.bngl) or pre-computed output (.gdat/.scan)")
    p.add_argument("--modified", required=True,
                   help="Path to modified model (.bngl) or pre-computed output (.gdat/.scan)")
    p.add_argument("--mode", choices=["ode", "ssa"], default="ode",
                   help="Comparison mode (default: ode)")
    p.add_argument("--tolerance", type=float, default=0.01,
                   help="Relative tolerance (rtol) for trajectory agreement (default: 0.01)")
    p.add_argument("--atol", type=float, default=1e-8,
                   help="Absolute tolerance floor for near-zero values (default: 1e-8). "
                        "Error criterion: |a-b| <= atol + rtol*|a|")
    p.add_argument("--format", choices=["bngl", "sbml", "antimony"], default=None,
                   help="Model format (auto-detected from extension if omitted)")
    p.add_argument("--output-type", choices=["gdat", "scan"], default="gdat",
                   help="Output file type to compare (default: gdat)")
    p.add_argument("--column-map", default=None,
                   help='Explicit column mapping, e.g. "R1=R1_conc,R2=R2_conc"')
    p.add_argument("--replicates", type=int, default=100,
                   help="Number of SSA replicates (default: 100, only used in ssa mode)")
    p.add_argument("--plot", action="store_true",
                   help="Generate comparison plots")
    p.add_argument("--plot-dir", default="./verify_plots",
                   help="Directory for comparison plots (default: ./verify_plots)")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    orig_path = args.original
    mod_path = args.modified

    # Detect formats
    orig_is_data = is_data_file(orig_path)
    mod_is_data = is_data_file(mod_path)

    if args.format:
        orig_fmt = mod_fmt = args.format
    else:
        orig_fmt = None if orig_is_data else detect_format(orig_path)
        mod_fmt = None if mod_is_data else detect_format(mod_path)

    # --- Unit reconciliation for the modified model ---
    mod_run_path = mod_path
    tmp_reconcile_dir = None
    if not mod_is_data and mod_fmt == "bngl" and detect_dimensional_functions(mod_path):
        print("Detected dimensional conversion functions in modified model — "
              "uncommenting for comparison.")
        tmp_reconcile_dir = tempfile.mkdtemp(prefix="verify_dim_")
        mod_run_path = uncomment_dimensional_functions(mod_path, tmp_reconcile_dir)

    try:
        if args.mode == "ode":
            # --- ODE mode ---
            sim_suffix = "ode"
            orig_cols, orig_data = simulate_model(
                orig_path, orig_fmt, args.output_type, sim_suffix
            )
            mod_cols, mod_data = simulate_model(
                mod_run_path, mod_fmt, args.output_type, sim_suffix
            )

            column_map = resolve_column_mapping(
                orig_cols, mod_cols, args.column_map
            )
            if not column_map:
                print("ERROR: No columns could be matched between original and modified output.")
                print(f"  Original columns: {orig_cols}")
                print(f"  Modified columns: {mod_cols}")
                sys.exit(1)

            results = compare_ode(
                orig_cols, orig_data, mod_cols, mod_data,
                column_map, args.tolerance, args.atol,
            )

            if args.plot:
                generate_plots(
                    orig_cols, orig_data, mod_cols, mod_data, column_map, args.plot_dir
                )
                print(f"Plots saved to {args.plot_dir}/")

        elif args.mode == "ssa":
            # --- SSA mode ---
            sim_suffix = "ssa"
            if orig_is_data or mod_is_data:
                # If pre-computed, treat as single-replicate comparison
                orig_cols, orig_data = simulate_model(
                    orig_path, orig_fmt, args.output_type, sim_suffix
                )
                mod_cols, mod_data = simulate_model(
                    mod_run_path, mod_fmt, args.output_type, sim_suffix
                )
                orig_mean = orig_data
                mod_mean = mod_data
            else:
                print(f"Running {args.replicates} SSA replicates for original model...")
                orig_cols, orig_mean, _orig_var = run_ssa_ensemble(
                    orig_path, orig_fmt, args.output_type, args.replicates, sim_suffix
                )
                print(f"Running {args.replicates} SSA replicates for modified model...")
                mod_cols, mod_mean, _mod_var = run_ssa_ensemble(
                    mod_run_path, mod_fmt, args.output_type, args.replicates, sim_suffix
                )

            column_map = resolve_column_mapping(
                orig_cols, mod_cols, args.column_map
            )
            if not column_map:
                print("ERROR: No columns could be matched between original and modified output.")
                sys.exit(1)

            results = compare_ssa(
                orig_cols, orig_mean, mod_cols, mod_mean,
                column_map, args.tolerance, args.atol,
            )

            if args.plot:
                generate_plots(
                    orig_cols, orig_mean, mod_cols, mod_mean, column_map, args.plot_dir
                )
                print(f"Plots saved to {args.plot_dir}/")

    finally:
        if tmp_reconcile_dir:
            shutil.rmtree(tmp_reconcile_dir, ignore_errors=True)

    passed = print_report(results, args.mode)
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
