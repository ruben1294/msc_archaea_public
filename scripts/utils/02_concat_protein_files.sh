#!/usr/bin/env bash

# Script: 02_concat_protein_files.sh
# Description: Concatenate all .faa files into a single file

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

load_project_config
require_vars RESULTS_DIR GENE_PREDICTION_SUBDIR

# Configuration
INPUT_DIR="${1:-${RESULTS_DIR}/${GENE_PREDICTION_SUBDIR}}"
OUTPUT_PREFIX="$INPUT_DIR/CMRA_all_proteins_prodigal_prefixed"
SEQUENCES_PER_FILE=10000
FINAL_OUTPUT="$INPUT_DIR/CMRA_all_proteins_final.faa"

# Count only original .faa files (in subdirectories, excluding generated files)
n_files=$(find "$INPUT_DIR" -mindepth 2 -maxdepth 2 -name "*.faa" -type f \
    ! -name "*prefixed.faa" \
    ! -name "CMRA_all_proteins*.faa" | wc -l)

echo "Found $n_files original .faa files"
echo "Output prefix: $OUTPUT_PREFIX"
echo "Sequences per file: $SEQUENCES_PER_FILE"

# Remove old output files if they exist
rm -f "$OUTPUT_PREFIX"_part*.faa
rm -f "$FINAL_OUTPUT"
find "$INPUT_DIR" -mindepth 2 -name "*.prefixed.faa" -type f -delete

# Concatenate all .faa files (including those in subdirectories)
# First, create prefixed versions with MAG names in headers
find "$INPUT_DIR" -mindepth 2 -maxdepth 2 -name "*.faa" -type f \
    ! -name "*prefixed.faa" \
    ! -name "CMRA_all_proteins*.faa" | while read -r faa; do
    mag=$(basename "$(dirname "$faa")")
    awk -v mag="$mag" '
        /^>/ {print ">" mag "|" substr($0,2); next}
        {print}
    ' "$faa" > "${faa%.faa}.prefixed.faa"
done

# Concatenate and split into multiple files based on sequence count
# Create a temporary file to hold all sequences
temp_file=$(mktemp)

find "$INPUT_DIR" -name "*.prefixed.faa" -type f | sort | while read -r prefixed_faa; do
    cat "$prefixed_faa" >> "$temp_file"
done

# Split the sequences into multiple output files
seq_count=0
file_count=1
current_output="${OUTPUT_PREFIX}_part${file_count}.faa"
: > "$current_output"  # Initialize empty file

while IFS= read -r line; do
    if [[ $line =~ ^">" ]]; then
        # If we've reached the sequence limit, start a new file BEFORE writing the header
        if (( seq_count > 0 && seq_count % SEQUENCES_PER_FILE == 0 )); then
            ((++file_count))
            current_output="${OUTPUT_PREFIX}_part${file_count}.faa"
            : > "$current_output"  # Initialize new empty file
        fi
        ((++seq_count))
    fi
    
    echo "$line" >> "$current_output"
done < "$temp_file"

rm "$temp_file"

echo "[OK] Done!"
echo "Total output files: $file_count"

# Count sequences in each output file and display summary
total_seqs=0
for i in $(seq 1 $file_count); do
    output_file="${OUTPUT_PREFIX}_part${i}.faa"
    n_seqs=$(grep -c "^>" "$output_file" || echo 0)
    file_size=$(du -h "$output_file" | cut -f1)
    echo "  ${output_file##*/}: $n_seqs sequences ($file_size)"
    ((total_seqs += n_seqs))
done

echo "Total sequences across all files: $total_seqs"

# Verification
find "$INPUT_DIR" -mindepth 2 -type f -name "*.prefixed.faa" \
! -name "${OUTPUT_PREFIX##*/}_part*.faa" \
-exec grep '^>' {} \; \
| sort > individual_headers.txt

for i in $(seq 1 $file_count); do
    grep '^>' "${OUTPUT_PREFIX}_part${i}.faa" || true
done | sort > concatenated_headers.txt

diff -s individual_headers.txt concatenated_headers.txt && echo "[OK] Headers match!" || echo "[FAIL] Headers do not match!"

# === Concatenate all parts into a single file with modified headers ===
echo ""
echo "=========================================="
echo "Concatenating all parts into final file with modified headers"
echo "=========================================="

# Remove final output if it exists
[[ -f "$FINAL_OUTPUT" ]] && rm "$FINAL_OUTPUT"

# Concatenate all part files and transform headers
for i in $(seq 1 $file_count); do
    part_file="${OUTPUT_PREFIX}_part${i}.faa"
    echo "Processing part ${i}/${file_count}..."
    
    awk '
        /^>/ {
            # Extract protein name (between > and first space)
            if (match($0, /^>([^ ]+)(.*)/, arr)) {
                protein_name = arr[1]
                description = arr[2]
                
                # Step 1: Replace all underscores with dots
                gsub(/_/, ".", protein_name)
                
                # Step 2: Replace all pipes with underscores
                gsub(/\|/, "_", protein_name)
                
                # Reconstruct header
                print ">" protein_name description
            } else {
                # If no space found, process entire header
                header = substr($0, 2)  # Remove >
                gsub(/_/, ".", header)
                gsub(/\|/, "_", header)
                print ">" header
            }
            next
        }
        {print}
    ' "$part_file" >> "$FINAL_OUTPUT"
done

# Count sequences in final file
final_seqs=$(grep -c "^>" "$FINAL_OUTPUT")
final_size=$(du -h "$FINAL_OUTPUT" | cut -f1)

echo ""
echo "[OK] Final concatenated file created!"
echo "File: $FINAL_OUTPUT"
echo "Total sequences: $final_seqs"
echo "File size: $final_size"

# Verify sequence count matches
if [[ $final_seqs -eq $total_seqs ]]; then
    echo "[OK] Sequence count matches!"
else
    echo "[FAIL] WARNING: Sequence count mismatch! (Expected: $total_seqs, Got: $final_seqs)"
fi

# Show example of header transformation
echo ""
echo "Example header transformations:"
echo "Original headers (from part files):"
head -5 "${OUTPUT_PREFIX}_part1.faa" | grep "^>" | head -3

echo ""
echo "Transformed headers (in final file):"
head -15 "$FINAL_OUTPUT" | grep "^>" | head -3

echo ""
echo "=========================================="
echo "Pipeline completed!"
echo "=========================================="
