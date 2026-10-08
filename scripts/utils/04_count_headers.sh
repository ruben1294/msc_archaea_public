#!/usr/bin/env bash

# Script: check_protein_headers.sh
# Description: Extract and count protein headers in FASTA file

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

load_project_config
require_vars RESULTS_DIR GENE_PREDICTION_SUBDIR

# Configuration
INPUT_DIR="${1:-${RESULTS_DIR}/${GENE_PREDICTION_SUBDIR}}"
INPUT_FILE="$INPUT_DIR/CMRA_all_proteins_final.faa"  # change to your file name
OUTPUT_FILE="$INPUT_DIR/header_counts.txt"
DUPLICATES_FILE="$INPUT_DIR/duplicate_headers.txt"

echo "=========================================="
echo "Protein Header Analysis"
echo "=========================================="
echo "Input file: $INPUT_FILE"
echo ""

# Check if file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "ERROR: File not found: $INPUT_FILE"
    exit 1
fi

# Extract headers (between > and first space)
echo "Extracting headers..."
grep "^>" "$INPUT_FILE" | awk '{print substr($1, 2)}' | sort > all_headers.tmp

# Count total headers
total_headers=$(wc -l < all_headers.tmp)
echo "Total headers found: $total_headers"

# Count unique headers
unique_headers=$(sort -u all_headers.tmp | wc -l)
echo "Unique headers: $unique_headers"

# Calculate duplicates
duplicates=$((total_headers - unique_headers))
echo "Duplicate headers: $duplicates"
echo ""

# Generate full count report
echo "Generating header count report..."
sort all_headers.tmp | uniq -c | sort -rn > "$OUTPUT_FILE"

echo "[OK] Full report saved to: $OUTPUT_FILE"
echo ""

# Extract only duplicated headers (count > 1)
if [[ $duplicates -gt 0 ]]; then
    echo "=========================================="
    echo "DUPLICATES DETECTED!"
    echo "=========================================="
    
    awk '$1 > 1' "$OUTPUT_FILE" > "$DUPLICATES_FILE"
    
    echo "Duplicated headers saved to: $DUPLICATES_FILE"
    echo ""
    echo "Top 10 most duplicated headers:"
    head -10 "$DUPLICATES_FILE"
    
    # Show statistics
    max_duplicates=$(head -1 "$DUPLICATES_FILE" | awk '{print $1}')
    n_duplicated=$(wc -l < "$DUPLICATES_FILE")
    
    echo ""
    echo "Statistics:"
    echo "  - Headers with duplicates: $n_duplicated"
    echo "  - Maximum occurrences: $max_duplicates"
else
    echo "[OK] No duplicates found! All headers are unique."
fi

# Cleanup
rm all_headers.tmp

echo ""
echo "=========================================="
echo "Analysis complete!"
echo "=========================================="
