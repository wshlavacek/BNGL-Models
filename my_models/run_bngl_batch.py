import os
import subprocess
import concurrent.futures
import shutil
from pathlib import Path

# Configuration
BNGL_EXE = Path("~/Simulations/BioNetGen-2.9.3/BNG2.pl").expanduser()
TIMEOUT = 120  # 2 minutes
MAX_WORKERS = 12
OUTPUT_FILE = "bngl_results.txt"
REFERENCE_DIR = Path("reference")
ARTIFACT_EXTENSIONS = {".gdat", ".cdat", ".xml", ".scan", ".net", ".species", ".graphml", ".m", ".tex", ".c"}
PROTECTED_DIRS = {"nf", "ode", "ssa", "reference"}

def run_bngl(file_path):
    """Runs a single BNGL file and returns the status."""
    try:
        # BNG2.pl is a Perl script, so we call it via 'perl' or directly if executable
        result = subprocess.run(
            ["perl", str(BNGL_EXE), str(file_path)],
            capture_output=True,
            text=True,
            timeout=TIMEOUT
        )
        
        if result.returncode == 0:
            return (file_path, "complete")
        else:
            # If return code is non-zero, it's a crash/error
            return (file_path, "crash")
            
    except subprocess.TimeoutExpired:
        return (file_path, "timeout")
    except Exception as e:
        return (file_path, f"error: {str(e)}")

def main():
    # 1. Handle reference directory at start
    if REFERENCE_DIR.exists() and REFERENCE_DIR.is_dir():
        print(f"Removing existing {REFERENCE_DIR} folder...")
        shutil.rmtree(REFERENCE_DIR)

    # Find all .bngl files in the current directory and subdirectories
    bngl_files = list(Path(".").rglob("*.bngl"))
    print(f"Found {len(bngl_files)} BNGL files.")

    results = {"complete": [], "crash": [], "timeout": [], "error": []}

    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_file = {executor.submit(run_bngl, f): f for f in bngl_files}
        
        for future in concurrent.futures.as_completed(future_to_file):
            file_path, status = future.result()
            if status == "complete":
                results["complete"].append(str(file_path))
            elif status == "crash":
                results["crash"].append(str(file_path))
            elif status == "timeout":
                results["timeout"].append(str(file_path))
            else:
                results["error"].append(f"{file_path} ({status})")
            
            print(f"Processed {file_path}: {status}")

    # Write results to file
    with open(OUTPUT_FILE, "w") as f:
        f.write("=== COMPLETE ===\n")
        f.write("\n".join(results["complete"]) + "\n\n")
        f.write("=== CRASHED ===\n")
        f.write("\n".join(results["crash"]) + "\n\n")
        f.write("=== TIMEOUT ===\n")
        f.write("\n".join(results["timeout"]) + "\n\n")
        f.write("=== ERRORS ===\n")
        f.write("\n".join(results["error"]) + "\n")

    print(f"\nResults written to {OUTPUT_FILE}")
    print(f"Complete: {len(results['complete'])}")
    print(f"Crashed: {len(results['crash'])}")
    print(f"Timeout: {len(results['timeout'])}")
    print(f"Errors: {len(results['error'])}")

    # 2. Move artifacts to reference folder
    print("\nMoving simulation artifacts to reference/...")
    REFERENCE_DIR.mkdir(exist_ok=True)

    # First, move directories (excluding protected ones and the reference dir itself)
    # We do this first so that files inside these dirs are moved as part of the directory
    for item in list(Path(".").iterdir()):
        if item.is_dir() and item.name not in PROTECTED_DIRS:
            dest_path = REFERENCE_DIR / item.name
            print(f"Moving directory {item.name} -> {dest_path}")
            shutil.move(str(item), str(dest_path))

    # Then, move any remaining files with specific extensions that weren't in the moved directories
    # Use a list to avoid issues with modifying the directory structure while iterating
    all_files = list(Path(".").rglob("*"))
    for path in all_files:
        if path.is_file() and path.suffix in ARTIFACT_EXTENSIONS:
            # Check if it's already inside reference/ to avoid moving it into itself
            if REFERENCE_DIR in path.parents:
                continue
                
            try:
                rel_path = path.relative_to(".")
                dest_path = REFERENCE_DIR / rel_path
                dest_path.parent.mkdir(parents=True, exist_ok=True)
                print(f"Moving file {rel_path} -> {dest_path}")
                shutil.move(str(path), str(dest_path))
            except FileNotFoundError:
                # This happens if the file was already moved as part of a directory move
                continue

    print("Artifacts moved successfully.")

if __name__ == "__main__":
    main()
