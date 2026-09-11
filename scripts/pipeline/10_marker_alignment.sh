#!/usr/bin/env bash
set -euo pipefail

# Step 10. Aligns each archaeal marker from the ar53 set with MAFFT, running
# three strategies over the same input so the downstream tree can be compared
# across them: genafpair (E-INS-i, tolerant of unalignable regions between
# conserved blocks), localpair (L-INS-i) and globalpair (G-INS-i).
#
# Each strategy writes to its own directory, and steps 11 to 13 carry the three
# branches forward independently.
#
# Markers with fewer than two sequences cannot be aligned; they are copied
# through unchanged so the file set stays complete for the concatenation in
# step 12.

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

trap 'log_error "Failed at line $LINENO with exit code $?"; exit 1' ERR

log_job_header "MAFFT marker alignment, three strategies"

load_project_config
require_vars PROJECT_DIR

# MAFFT_PATH in server.env points at a specific build; without it the script
# falls back to whatever mafft is on PATH.
if [[ -n "${MAFFT_PATH:-}" && -x "${MAFFT_PATH}" ]]; then
    MAFFT_EXEC="$MAFFT_PATH"
elif command -v mafft > /dev/null 2>&1; then
    MAFFT_EXEC="$(command -v mafft)"
else
    log_error "MAFFT not found. Set MAFFT_PATH in $CONFIG_FILE or put mafft on PATH."
    log_error "PATH is: $PATH"
    exit 1
fi

mafft_version=$("$MAFFT_EXEC" --version 2>&1 | head -n 1 || true)
log_info "MAFFT: $MAFFT_EXEC ($mafft_version)"

WORK_DIR="${PROJECT_DIR}/${MARKER_ALIGNMENT_SUBDIR:-10_marker_alignment}"
INPUT_DIR="${WORK_DIR}/ar53_filtered"

OUTPUT_DIR_GENAFPAIR="${WORK_DIR}/genafpair"
OUTPUT_DIR_LOCALPAIR="${WORK_DIR}/localpair"
OUTPUT_DIR_GLOBALPAIR="${WORK_DIR}/globalpair"

mkdir -p "$OUTPUT_DIR_GENAFPAIR" "$OUTPUT_DIR_LOCALPAIR" "$OUTPUT_DIR_GLOBALPAIR"

if [[ ! -d "$INPUT_DIR" ]]; then
    log_error "Input directory not found: $INPUT_DIR"
    log_error "Run utils/07_filter_ar53_markers.sh first."
    exit 1
fi

shopt -s nullglob
fasta_files=("$INPUT_DIR"/*.fa)
total_files=${#fasta_files[@]}
shopt -u nullglob

if [[ $total_files -eq 0 ]]; then
    log_error "No .fa files in $INPUT_DIR"
    exit 1
fi

log_info "Found $total_files marker files"
log_info "Input:      $INPUT_DIR"
log_info "GenAFpair:  $OUTPUT_DIR_GENAFPAIR"
log_info "Localpair:  $OUTPUT_DIR_LOCALPAIR"
log_info "Globalpair: $OUTPUT_DIR_GLOBALPAIR"

total_processed=0
total_success=0
total_failed=0

genafpair_success=0;  genafpair_failed=0
localpair_success=0;  localpair_failed=0
globalpair_success=0; globalpair_failed=0

threads="${SLURM_CPUS_PER_TASK:-24}"

# Runs one MAFFT strategy over one marker and reports how long it took.
# Arguments: flag, input file, output file, label.
align_marker() {
    local flag="$1" infile="$2" outfile="$3" label="$4"
    local start_time end_time elapsed

    start_time=$(date +%s)
    if "$MAFFT_EXEC" "$flag" --maxiterate 1000 --thread "$threads" "$infile" > "$outfile"; then
        end_time=$(date +%s)
        elapsed=$((end_time - start_time))
        log_info "    [OK] $label (${elapsed}s)"
        return 0
    fi
    log_error "    [FAIL] $label"
    return 1
}

for fasta_file in "${fasta_files[@]}"; do
    basename_file=$(basename "$fasta_file")
    stem="${basename_file%.fa}"

    seq_count=$(grep -c "^>" "$fasta_file" || echo "0")

    log_info "Processing [$((total_processed + 1))/$total_files]: $basename_file ($seq_count sequences)"

    if [[ $seq_count -lt 2 ]]; then
        log_warn "  Fewer than 2 sequences, copying through unaligned"
        cp "$fasta_file" "$OUTPUT_DIR_GENAFPAIR/${stem}_aligned.fa"
        cp "$fasta_file" "$OUTPUT_DIR_LOCALPAIR/${stem}_aligned.fa"
        cp "$fasta_file" "$OUTPUT_DIR_GLOBALPAIR/${stem}_aligned.fa"

        ((++genafpair_success))
        ((++localpair_success))
        ((++globalpair_success))
        ((total_success += 3))
        ((++total_processed))
        continue
    fi

    log_info "  [1/3] GenAFpair"
    if align_marker --genafpair "$fasta_file" "$OUTPUT_DIR_GENAFPAIR/${stem}_aligned.fa" GenAFpair; then
        ((++genafpair_success)); ((++total_success))
    else
        ((++genafpair_failed)); ((++total_failed))
    fi

    log_info "  [2/3] Localpair"
    if align_marker --localpair "$fasta_file" "$OUTPUT_DIR_LOCALPAIR/${stem}_aligned.fa" Localpair; then
        ((++localpair_success)); ((++total_success))
    else
        ((++localpair_failed)); ((++total_failed))
    fi

    log_info "  [3/3] Globalpair"
    if align_marker --globalpair "$fasta_file" "$OUTPUT_DIR_GLOBALPAIR/${stem}_aligned.fa" Globalpair; then
        ((++globalpair_success)); ((++total_success))
    else
        ((++globalpair_failed)); ((++total_failed))
    fi

    ((++total_processed))
done

log_info "Alignment finished"
log_info "  Markers processed:     $total_processed"
log_info "  Alignments attempted:  $((total_processed * 3))"
log_info "  GenAFpair:             $genafpair_success ok, $genafpair_failed failed"
log_info "  Localpair:             $localpair_success ok, $localpair_failed failed"
log_info "  Globalpair:            $globalpair_success ok, $globalpair_failed failed"
log_info "  Total:                 $total_success ok, $total_failed failed"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"

if [[ $total_failed -gt 0 ]]; then
    log_error "$total_failed alignments failed."
    exit 1
fi
