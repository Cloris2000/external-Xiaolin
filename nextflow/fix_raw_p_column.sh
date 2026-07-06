#!/bin/bash
# Fix the P column in all .raw_p files.
# Bug: P was computed as 10^(-CHISQ) [col 12] instead of 10^(-LOG10P) [col 13].
# Fix: recompute P = 10^(-LOG10P) using the correct column 13.
#
#SBATCH --job-name=fix_raw_p
#SBATCH --partition=mediumtmp
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=logs/fix_raw_p_%j.out
#SBATCH --error=logs/fix_raw_p_%j.err

set -euo pipefail

RESULTS=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results

# Find all raw_p files
mapfile -t files < <(find "$RESULTS" -name "*.raw_p")
total=${#files[@]}
echo "Found $total .raw_p files to fix"

fix_file() {
    local fpath="$1"
    # Verify column 13 is LOG10P (sanity check on first file)
    local col13_name
    col13_name=$(awk 'NR==1{print $13; exit}' "$fpath")
    if [ "$col13_name" != "LOG10P" ]; then
        echo "SKIP (unexpected col13='$col13_name'): $fpath"
        return
    fi
    # Replace last column (P) with 10^(-LOG10P) = 10^(-col13)
    local tmpfile="${fpath}.tmp"
    awk 'NR==1{print} NR>1{$NF=10^(-$13); print}' "$fpath" > "$tmpfile" && mv "$tmpfile" "$fpath"
    echo "FIXED: $fpath"
}

export -f fix_file

# Run in parallel using 4 cores
printf '%s\n' "${files[@]}" | xargs -P 4 -I{} bash -c 'fix_file "$@"' _ {}

echo "Done. Fixed $total files."
