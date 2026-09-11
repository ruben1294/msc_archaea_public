#!/usr/bin/env bash
set -euo pipefail

# Step 11. Trims poorly aligned columns out of every marker alignment with
# trimAl, keeping the three MAFFT strategies from step 10 in separate branches.
#
# The heuristic is -automated1, which picks between gappyout and strict based on
# the alignment's own gap and similarity distribution. trimAl's authors
# recommend it for alignments headed into maximum likelihood tree inference,
# which is what step 13 does.
#
# Each branch gets a trimmed/ directory, an HTML report per marker, and a TSV of
# per-marker before and after counts. The TSV is what tells you whether a branch
# lost so many columns that its tree is not worth inferring.

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

trap 'log_error "Failed at line $LINENO with exit code $?"; exit 1' ERR

log_job_header "trimAl column trimming, three branches"

load_project_config
require_vars PROJECT_DIR

# TRIMAL_PATH in server.env points at a specific build; without it the script
# falls back to whatever trimal is on PATH.
if [[ -n "${TRIMAL_PATH:-}" && -x "${TRIMAL_PATH}" ]]; then
    TRIMAL_EXEC="$TRIMAL_PATH"
elif command -v trimal > /dev/null 2>&1; then
    TRIMAL_EXEC="$(command -v trimal)"
else
    log_error "trimAl not found. Set TRIMAL_PATH in $CONFIG_FILE or put trimal on PATH."
    log_error "PATH is: $PATH"
    exit 1
fi

trimal_version=$("$TRIMAL_EXEC" --version 2>&1 | head -n 1 || true)
log_info "trimAl: $TRIMAL_EXEC ($trimal_version)"

WORK_DIR="${PROJECT_DIR}/${MARKER_ALIGNMENT_SUBDIR:-10_marker_alignment}"
OUTPUT_BASE_DIR="${PROJECT_DIR}/${TRIMAL_SUBDIR:-11_remove_uninformative_seqs}"

METHODS=(genafpair localpair globalpair)

for method in "${METHODS[@]}"; do
    if [[ ! -d "${WORK_DIR}/${method}" ]]; then
        log_error "Input directory not found: ${WORK_DIR}/${method}"
        log_error "Run 10_marker_alignment.sh first."
        exit 1
    fi
    mkdir -p "${OUTPUT_BASE_DIR}/${method}/trimmed" "${OUTPUT_BASE_DIR}/${method}/reports"
done

log_info "Input:  $WORK_DIR"
log_info "Output: $OUTPUT_BASE_DIR"

total_processed=0
total_success=0
total_failed=0

# process_method reaches these through namerefs and the summary reads them back
# through ${!var}, neither of which shellcheck can follow.
# shellcheck disable=SC2034
{
    genafpair_success=0;  genafpair_failed=0
    localpair_success=0;  localpair_failed=0
    globalpair_success=0; globalpair_failed=0
}

DATE_TAG=$(date '+%Y%m%d')
JOB_TAG="${SLURM_JOB_ID:-local}"

STATS_HEADER="marker\tseqs_before\tseqs_after\tresidues_before\tresidues_after\talignment_length_before\talignment_length_after\tpercent_residues_retained"

# Trims one branch. The two counter arguments are variable names; namerefs let
# the function update the caller's counters without eval.
process_method() {
    local method_name="$1"
    local input_dir="$2"
    local output_dir="$3"
    local stats_file="$4"
    local -n success_ref="$5"
    local -n failed_ref="$6"

    log_info "Trimming $method_name alignments"

    shopt -s nullglob
    local input_files=("$input_dir"/*.fa)
    shopt -u nullglob

    if [[ ${#input_files[@]} -eq 0 ]]; then
        log_error "No .fa alignments in $input_dir"
        return 1
    fi

    local total_files=${#input_files[@]}
    log_info "  $total_files alignment files"

    local method_processed=0

    for input_file in "${input_files[@]}"; do
        [[ -f "$input_file" ]] || continue

        local file_name marker_name output_file report_file
        file_name=$(basename "$input_file")
        marker_name="${file_name%.*}"
        output_file="${output_dir}/trimmed/${marker_name}_trimmed.fa"
        report_file="${output_dir}/reports/${marker_name}_trimal_report.html"

        log_info "  [$((method_processed + 1))/$total_files] $marker_name"

        if "$TRIMAL_EXEC" \
            -in "$input_file" \
            -out "$output_file" \
            -fasta \
            -automated1 \
            -htmlout "$report_file"; then

            if [[ ! -f "$output_file" ]]; then
                log_error "    [FAIL] trimAl returned 0 but wrote nothing to $output_file"
                ((++failed_ref))
                ((++total_failed))
                ((++method_processed))
                ((++total_processed))
                continue
            fi

            local seqs_before seqs_after residues_before residues_after
            local align_len_before align_len_after percent_retained
            seqs_before=$(grep -c "^>" "$input_file" || echo "0")
            seqs_after=$(grep -c "^>" "$output_file" || echo "0")

            # Residue counts drop gaps and newlines, so they measure real
            # sequence rather than alignment width.
            residues_before=$(grep -v "^>" "$input_file" | tr -d '\n-' | wc -c)
            residues_after=$(grep -v "^>" "$output_file" | tr -d '\n-' | wc -c)

            # Alignment width comes from the first sequence line, which is the
            # same for every record in an aligned FASTA.
            align_len_before=$(awk '/^>/ {next} {print length($0); exit}' "$input_file")
            align_len_after=$(awk '/^>/ {next} {print length($0); exit}' "$output_file")

            if [[ $residues_before -gt 0 ]]; then
                percent_retained=$(awk "BEGIN {printf \"%.2f\", ($residues_after/$residues_before)*100}")
            else
                percent_retained="0.00"
            fi

            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$marker_name" "$seqs_before" "$seqs_after" \
                "$residues_before" "$residues_after" \
                "$align_len_before" "$align_len_after" "$percent_retained" \
                >> "$stats_file"

            log_info "    [OK] ${seqs_before} seqs / ${align_len_before} cols to ${seqs_after} seqs / ${align_len_after} cols (${percent_retained}% residues kept)"
            ((++success_ref))
            ((++total_success))
        else
            log_error "    [FAIL] $marker_name"
            ((++failed_ref))
            ((++total_failed))
        fi

        ((++method_processed))
        ((++total_processed))
    done

    log_info "  $method_name: $method_processed files processed"
}

for method in "${METHODS[@]}"; do
    stats_file="${OUTPUT_BASE_DIR}/${method}/${DATE_TAG}_${JOB_TAG}_trim_statistics.tsv"
    printf '%b\n' "$STATS_HEADER" > "$stats_file"
    process_method "$method" \
        "${WORK_DIR}/${method}" \
        "${OUTPUT_BASE_DIR}/${method}" \
        "$stats_file" \
        "${method}_success" \
        "${method}_failed"
done

log_info "Summary"
log_info "  Alignments processed: $total_processed"
log_info "  Successful:           $total_success"
log_info "  Failed:               $total_failed"

for method in "${METHODS[@]}"; do
    trimmed_count=$(find "${OUTPUT_BASE_DIR}/${method}/trimmed" -name "*_trimmed.fa" -type f 2>/dev/null | wc -l)
    reports_count=$(find "${OUTPUT_BASE_DIR}/${method}/reports" -name "*_report.html" -type f 2>/dev/null | wc -l)
    success_var="${method}_success"
    failed_var="${method}_failed"
    log_info "  ${method}: ${!success_var} ok, ${!failed_var} failed (${trimmed_count} trimmed, ${reports_count} reports)"
done

log_info "Output layout, one block per branch:"
log_info "  ${OUTPUT_BASE_DIR}/<branch>/trimmed/          trimmed alignments"
log_info "  ${OUTPUT_BASE_DIR}/<branch>/reports/          trimAl HTML reports"
log_info "  ${OUTPUT_BASE_DIR}/<branch>/*_statistics.tsv  per-marker counts"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"

if [[ $total_failed -gt 0 ]]; then
    log_error "$total_failed alignments failed."
    exit 1
fi
