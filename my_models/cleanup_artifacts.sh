#!/bin/bash

# Cleanup script to remove BioNetGen simulation artifacts
echo "Cleaning up simulation artifacts..."

# Define extensions to be removed
EXTENSIONS=("gdat" "cdat" "xml" "scan" "net" "species" "graphml" "m" "tex" "c")

for ext in "${EXTENSIONS[@]}"; do
    echo "Removing *.$ext files..."
    find . -type f -name "*.$ext" -delete
done

echo "Removing artifact directories..."
# Find all directories, excluding the current directory and protected ones
find . -maxdepth 1 -type d ! -name "." ! -name "nf" ! -name "ode" ! -name "ssa" -exec rm -rf {} +

echo "Cleanup complete."
