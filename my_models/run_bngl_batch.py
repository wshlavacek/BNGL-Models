#!/usr/bin/env python3
"""Batch-run BNGL models with BNG2.pl and collect their artifacts."""

import argparse
import concurrent.futures
import os
import shutil
import subprocess
import sys
from pathlib import Path

DEFAULT_BNGL_EXE = Path("~/Simulations/BioNetGen-2.9.3/BNG2.pl").expanduser()
DEFAULT_TIMEOUT = 120
DEFAULT_WORKERS = 12
DEFAULT_OUTPUT_FILE = Path("bngl_results.txt")
DEFAULT_REFERENCE_DIR = Path("reference")
ARTIFACT_EXTENSIONS = {
    ".cdat",
    ".c",
    ".gdat",
    ".graphml",
    ".m",
    ".net",
    ".scan",
    ".species",
    ".tex",
    ".xml",
}
PROTECTED_DIRS = {"nf", "ode", "ssa", "reference"}


def resolve_bngl_exe(cli_path):
    """Resolve BNG2.pl. An explicit --bngl-exe or $BNG2_PL is required to exist;
    only the built-in default falls back silently when missing."""
    if cli_path is not None:
        path = Path(cli_path).expanduser()
        if not path.is_file():
            sys.exit(f"--bngl-exe path does not exist: {path}")
        return path
    env = os.environ.get("BNG2_PL")
    if env:
        path = Path(env).expanduser()
        if not path.is_file():
            sys.exit(f"$BNG2_PL points to a non-existent path: {path}")
        return path
    if DEFAULT_BNGL_EXE.is_file():
        return DEFAULT_BNGL_EXE
    sys.exit(
        f"BNG2.pl not found at default location ({DEFAULT_BNGL_EXE}). "
        f"Pass --bngl-exe or set $BNG2_PL."
    )


def run_bngl(file_path, bngl_exe, timeout, log_dir):
    """Run a single BNGL file. Returns (file_path, status). Writes a log on failure."""
    try:
        result = subprocess.run(
            ["perl", str(bngl_exe), str(file_path)],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if result.returncode == 0:
            return file_path, "complete"
        _write_log(
            log_dir,
            file_path,
            "crash",
            returncode=result.returncode,
            stdout=result.stdout,
            stderr=result.stderr,
        )
        return file_path, "crash"
    except subprocess.TimeoutExpired as e:
        _write_log(log_dir, file_path, "timeout", stdout=e.stdout or "", stderr=e.stderr or "")
        return file_path, "timeout"
    except Exception as e:
        return file_path, f"error: {e}"


def _write_log(log_dir, file_path, status, returncode=None, stdout="", stderr=""):
    log_path = log_dir / f"{file_path.stem}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    header = [f"# {file_path}", f"# status={status}"]
    if returncode is not None:
        header.append(f"# returncode={returncode}")
    log_path.write_text(
        "\n".join(header) + f"\n\n--- STDOUT ---\n{stdout}\n\n--- STDERR ---\n{stderr}\n"
    )


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bngl-exe", help="Path to BNG2.pl (overrides $BNG2_PL and default)")
    p.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        help=f"Per-model timeout in seconds (default: {DEFAULT_TIMEOUT})",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"Parallel workers (default: {DEFAULT_WORKERS})",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_FILE,
        help=f"Results summary file (default: {DEFAULT_OUTPUT_FILE})",
    )
    p.add_argument(
        "--reference-dir",
        type=Path,
        default=DEFAULT_REFERENCE_DIR,
        help=f"Where to collect artifacts (default: {DEFAULT_REFERENCE_DIR})",
    )
    p.add_argument(
        "--clean",
        action="store_true",
        help="Remove an existing reference dir before running (default: refuse)",
    )
    p.add_argument(
        "--target-dir",
        type=Path,
        help="Directory to scan for .bngl files (default: this script's directory)",
    )
    return p.parse_args()


def write_results(output_file, results):
    with open(output_file, "w") as f:
        for category in ("complete", "crash", "timeout", "error"):
            f.write(f"=== {category.upper()} ===\n")
            f.write("\n".join(results[category]) + "\n\n")


def move_artifacts(reference_dir):
    """Sweep stray BNG2.pl artifacts (top-level dirs and tagged files) into reference/."""
    print(f"\nMoving simulation artifacts to {reference_dir}/...")
    reference_dir.mkdir(exist_ok=True)

    for item in list(Path(".").iterdir()):
        if item.is_dir() and item.name not in PROTECTED_DIRS and not item.name.startswith("."):
            dest = reference_dir / item.name
            print(f"  dir  {item.name} -> {dest}")
            shutil.move(str(item), str(dest))

    for path in list(Path(".").rglob("*")):
        if not path.is_file() or path.suffix not in ARTIFACT_EXTENSIONS:
            continue
        if reference_dir in path.parents:
            continue
        rel_path = path.relative_to(".")
        dest = reference_dir / rel_path
        dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            print(f"  file {rel_path} -> {dest}")
            shutil.move(str(path), str(dest))
        except FileNotFoundError:
            continue


def main():
    args = parse_args()

    # Anchor to a known directory so relative paths and artifact sweeping behave.
    script_dir = Path(__file__).resolve().parent
    target_dir = (args.target_dir or script_dir).resolve()
    if not target_dir.is_dir():
        sys.exit(f"Target dir does not exist: {target_dir}")
    os.chdir(target_dir)

    # Validate dependencies up front, not 12x inside workers.
    if shutil.which("perl") is None:
        sys.exit("'perl' not found on PATH.")
    bngl_exe = resolve_bngl_exe(args.bngl_exe)

    print(f"BNG2.pl:    {bngl_exe}")
    print(f"Working in: {target_dir}")

    reference_dir = args.reference_dir
    if reference_dir.exists():
        if args.clean:
            print(f"Removing existing {reference_dir}/...")
            shutil.rmtree(reference_dir)
        else:
            sys.exit(
                f"{reference_dir}/ already exists. Pass --clean to remove it, "
                f"or --reference-dir to use a different path."
            )
    log_dir = reference_dir / "_logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    bngl_files = sorted(Path(".").rglob("*.bngl"))
    if not bngl_files:
        sys.exit(f"No .bngl files found under {target_dir}")
    print(f"Found {len(bngl_files)} BNGL files.\n")

    results = {"complete": [], "crash": [], "timeout": [], "error": []}
    total = len(bngl_files)
    interrupted = False

    executor = concurrent.futures.ProcessPoolExecutor(max_workers=args.workers)
    try:
        future_to_file = {
            executor.submit(run_bngl, f, bngl_exe, args.timeout, log_dir): f for f in bngl_files
        }
        done = 0
        for future in concurrent.futures.as_completed(future_to_file):
            file_path, status = future.result()
            done += 1
            if status in results:
                results[status].append(str(file_path))
            else:
                results["error"].append(f"{file_path} ({status})")
            print(f"[{done}/{total}] {file_path}: {status}")
    except KeyboardInterrupt:
        interrupted = True
        print("\nInterrupted — cancelling pending jobs and writing partial results...")
    finally:
        executor.shutdown(wait=not interrupted, cancel_futures=True)

    write_results(args.output, results)
    print(f"\nResults written to {args.output}")
    for k in ("complete", "crash", "timeout", "error"):
        print(f"  {k.capitalize():9s} {len(results[k])}")

    move_artifacts(reference_dir)
    print("Done." if not interrupted else "Done (partial).")
    sys.exit(1 if interrupted else 0)


if __name__ == "__main__":
    main()
