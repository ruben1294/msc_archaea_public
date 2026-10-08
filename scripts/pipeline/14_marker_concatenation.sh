#!/usr/bin/env bash
set -euo pipefail

# Step 14. Drops markers with too few taxa, then concatenates what survives into
# a supermatrix per MAFFT branch with FASconCAT-G, writing the partition file
# that step 15 hands to IQ-TREE.
#
# MIN_SEQUENCES is 5. A marker present in fewer genomes contributes almost no
# resolution and a lot of missing data, and the missing-data fraction reported
# at the end is the number to watch when comparing the three branches.
#
# FASconCAT-G is interactive: it reads single-letter commands from stdin, so the
# heredoc below sends "l" to turn on the partition file and "s" to start. It
# also only accepts .fas and .fasta, and writes its output into the current
# working directory, which is why each run happens inside a subshell that cds
# into that branch's output directory first.
#
# utils/install_FASconCAT-G.sh installs it and writes the fasconcat_config.env
# that this script sources.

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

trap 'log_error "Failed at line $LINENO with exit code $?"; exit 1' ERR

log_job_header "FASconCAT-G marker concatenation, three branches"

load_project_config
require_commands perl
require_vars PROJECT_DIR RESULTS_DIR

WORK_DIR="${RESULTS_DIR}/${MARKER_CONCAT_SUBDIR:-14_marker_concatenation}"
TRIMMED_DIR="${RESULTS_DIR}/${TRIMAL_SUBDIR:-13_remove_uninformative_seqs}"
FASCONCAT_CONFIG="${WORK_DIR}/fasconcat_config.env"

if [[ ! -f "$FASCONCAT_CONFIG" ]]; then
    log_error "FASconCAT-G configuration not found: $FASCONCAT_CONFIG"
    log_error "Run utils/install_FASconCAT-G.sh first."
    exit 1
fi
# shellcheck source=/dev/null
source "$FASCONCAT_CONFIG"

if [[ ! -f "${FASCONCAT_EXEC:-}" ]]; then
    log_error "FASconCAT-G executable not found: ${FASCONCAT_EXEC:-unset}"
    log_error "Run utils/install_FASconCAT-G.sh first."
    exit 1
fi

log_info "FASconCAT-G: $FASCONCAT_EXEC"

MIN_SEQUENCES=5
METHODS=(genafpair localpair globalpair)

log_info "Minimum sequences per marker: $MIN_SEQUENCES"

for method in "${METHODS[@]}"; do
    if [[ ! -d "${TRIMMED_DIR}/${method}/trimmed" ]]; then
        log_error "Input directory not found: ${TRIMMED_DIR}/${method}/trimmed"
        log_error "Run 13_remove_uninformative_seqs.sh first."
        exit 1
    fi
done

# ---- step 0: clear previous runs ----

for method in "${METHODS[@]}"; do
    if [[ -d "${WORK_DIR}/${method}" ]]; then
        log_info "Removing previous ${method} output"
        rm -rf "${WORK_DIR:?}/${method}"
    fi
    mkdir -p "${WORK_DIR}/${method}"/{input,output,logs}
done

# ---- step 1: filter markers by taxon count ----

# Sets COPIED_COUNT rather than echoing it. The log helpers write to stdout, so
# a command substitution here would capture the log lines along with the number.
COPIED_COUNT=0
filter_and_copy_markers() {
    local method_name="$1" input_dir="$2" output_dir="$3"
    local total_files=0 copied_files=0 skipped_files=0
    local marker_list="${output_dir}/logs/marker_list.txt"

    log_info "Filtering $method_name markers"
    printf 'marker_name\tsequence_count\tstatus\n' > "$marker_list"

    shopt -s nullglob
    for fasta_file in "$input_dir"/*.fa; do
        [[ -f "$fasta_file" ]] || continue
        ((++total_files))

        local file_name marker_name seq_count
        file_name=$(basename "$fasta_file")
        marker_name="${file_name%_trimmed.fa}"
        seq_count=$(grep -c "^>" "$fasta_file" || echo "0")

        if [[ $seq_count -ge $MIN_SEQUENCES ]]; then
            cp "$fasta_file" "${output_dir}/input/${file_name}"
            printf '%s\t%s\tCOPIED\n' "$marker_name" "$seq_count" >> "$marker_list"
            ((++copied_files))
        else
            printf '%s\t%s\tSKIPPED\n' "$marker_name" "$seq_count" >> "$marker_list"
            ((++skipped_files))
        fi
    done
    shopt -u nullglob

    log_info "  Total markers:                $total_files"
    log_info "  Kept (>= $MIN_SEQUENCES seqs):            $copied_files"
    log_info "  Dropped (< $MIN_SEQUENCES seqs):          $skipped_files"
    log_info "  Marker list:                  $marker_list"

    COPIED_COUNT=$copied_files
}

declare -A copied_per_method
for method in "${METHODS[@]}"; do
    filter_and_copy_markers "$method" \
        "${TRIMMED_DIR}/${method}/trimmed" \
        "${WORK_DIR}/${method}"
    copied_per_method[$method]=$COPIED_COUNT
done

# ---- step 2: check which taxa are present in which markers ----

# A taxon missing from many markers becomes a mostly-gap row in the supermatrix,
# so this reports the presence counts before concatenation rather than leaving
# them to be discovered in the tree.
verify_headers() {
    local method_name="$1" input_dir="$2" output_dir="$3"
    local headers_file="${output_dir}/logs/all_headers.txt"
    local unique_headers="${output_dir}/logs/unique_headers.txt"
    local taxa_report="${output_dir}/logs/taxa_report.txt"

    log_info "Checking $method_name headers"

    : > "$headers_file"
    shopt -s nullglob
    local marker_files=("$input_dir"/*.fa)
    shopt -u nullglob

    for fasta_file in "${marker_files[@]}"; do
        grep "^>" "$fasta_file" >> "$headers_file"
    done

    sort "$headers_file" | uniq > "$unique_headers"

    local total_headers unique_count file_count
    total_headers=$(wc -l < "$headers_file")
    unique_count=$(wc -l < "$unique_headers")
    file_count=${#marker_files[@]}

    log_info "  Headers: $total_headers across $file_count markers"
    log_info "  Unique taxa: $unique_count"

    {
        echo "Taxon presence across markers, $method_name"
        echo ""
    } > "$taxa_report"

    while IFS= read -r header; do
        local header_clean count
        header_clean="${header#>}"
        count=$(grep -c "^${header}$" "$headers_file" || echo "0")
        printf "%-50s %3d/%3d markers\n" "$header_clean" "$count" "$file_count" >> "$taxa_report"

        if [[ $count -ne $file_count ]]; then
            log_warn "  $header_clean present in $count/$file_count markers"
        fi
    done < "$unique_headers"

    log_info "  Taxa report: $taxa_report"
}

for method in "${METHODS[@]}"; do
    verify_headers "$method" "${WORK_DIR}/${method}/input" "${WORK_DIR}/${method}"
done

# ---- step 3: concatenate ----

run_fasconcat() {
    local method_name="$1" input_dir="$2" output_dir="$3"
    local log_file="${output_dir}/logs/fasconcat.log"
    local out_dir="${output_dir}/output"
    local supermatrix_name="supermatrix_${method_name}"

    log_info "Concatenating $method_name"

    shopt -s nullglob
    local input_files=("$input_dir"/*.fa)
    shopt -u nullglob

    if [[ ${#input_files[@]} -eq 0 ]]; then
        log_error "  No markers survived filtering for $method_name"
        return 1
    fi

    log_info "  Input markers: ${#input_files[@]}"

    # FASconCAT-G only reads .fas and .fasta.
    for fa_file in "${input_files[@]}"; do
        local stem
        stem=$(basename "$fa_file" .fa)
        cp "$fa_file" "${out_dir}/${stem}.fas"
    done

    if ! (cd "$out_dir" && LC_ALL=C perl "$FASCONCAT_EXEC" ./*.fas <<'FCC'
l
s
FCC
        ) > "$log_file" 2>&1; then
        log_error "  [FAIL] FASconCAT-G, see $log_file"
        tail -30 "$log_file" | while read -r line; do log_error "    $line"; done
        return 1
    fi

    log_info "  [OK] FASconCAT-G finished"

    shopt -s nullglob
    local produced=("$out_dir"/FcC_supermatrix*.fas)
    local partitions=("$out_dir"/FcC_supermatrix*_partition.txt)
    shopt -u nullglob

    if [[ ${#produced[@]} -eq 0 ]]; then
        log_error "  FASconCAT-G wrote no supermatrix into $out_dir"
        return 1
    fi

    local fasta_out="${produced[0]}"
    local taxa_count aln_length
    taxa_count=$(grep -c "^>" "$fasta_out")
    aln_length=$(head -2 "$fasta_out" | tail -1 | tr -d '\n' | wc -c)
    log_info "  Taxa:   $taxa_count"
    log_info "  Length: $aln_length positions"

    # Rename off the FASconCAT-G default so the three branches stay apart.
    mv "$fasta_out" "${out_dir}/FcC_${supermatrix_name}.fas"
    log_info "  Supermatrix: ${out_dir}/FcC_${supermatrix_name}.fas"

    if [[ ${#partitions[@]} -gt 0 ]]; then
        mv "${partitions[0]}" "${out_dir}/FcC_${supermatrix_name}_partition.txt"
        log_info "  Partitions:  ${out_dir}/FcC_${supermatrix_name}_partition.txt"
    else
        log_warn "  No partition file produced"
    fi

    # Only the per-marker copies go; the supermatrix has already been renamed.
    rm -f "${out_dir}"/*_trimmed.fas
}

declare -A method_status
for method in "${METHODS[@]}"; do
    if run_fasconcat "$method" "${WORK_DIR}/${method}/input" "${WORK_DIR}/${method}"; then
        method_status[$method]=0
    else
        method_status[$method]=1
    fi
done

# ---- step 4: quality of each supermatrix ----

assess_quality() {
    local method_name="$1" output_dir="$2"
    local out_dir="${output_dir}/output"
    local fasta_out="${out_dir}/FcC_supermatrix_${method_name}.fas"
    local nexus_part="${out_dir}/FcC_supermatrix_${method_name}_partition.txt"
    local quality_report="${output_dir}/logs/quality_report.txt"

    if [[ ! -f "$fasta_out" ]]; then
        log_warn "  No supermatrix for $method_name, skipping assessment"
        return 1
    fi

    local taxa_count aln_length total_chars gap_chars missing_chars
    local total_missing missing_percent
    taxa_count=$(grep -c "^>" "$fasta_out")
    aln_length=$(head -2 "$fasta_out" | tail -1 | tr -d '\n' | wc -c)
    total_chars=$(grep -v "^>" "$fasta_out" | tr -d '\n' | wc -c)
    gap_chars=$(grep -v "^>" "$fasta_out" | tr -dc '-' | wc -c)
    missing_chars=$(grep -v "^>" "$fasta_out" | tr -dc '?' | wc -c)
    total_missing=$((gap_chars + missing_chars))

    if [[ $total_chars -gt 0 ]]; then
        missing_percent=$(awk "BEGIN {printf \"%.2f\", ($total_missing/$total_chars)*100}")
    else
        missing_percent="0.00"
    fi

    {
        echo "Quality assessment, $method_name"
        echo "Generated: $(date '+%Y-%m-%d %H:%M:%S')"
        echo ""
        echo "Taxa:                    $taxa_count"
        echo "Alignment length:        $aln_length"
        echo "Total characters:        $total_chars"
        echo "Gap characters (-):      $gap_chars"
        echo "Missing characters (?):  $missing_chars"
        echo "Missing data:            ${missing_percent}%"
    } > "$quality_report"

    if [[ -f "$nexus_part" ]]; then
        local partition_count
        partition_count=$(wc -l < "$nexus_part")
        {
            echo "Partitions:              $partition_count"
            echo ""
            echo "First partitions:"
            head -10 "$nexus_part"
            [[ $partition_count -gt 10 ]] && echo "(showing 10 of $partition_count)"
        } >> "$quality_report"
    fi

    log_info "  $method_name: $taxa_count taxa, $aln_length positions, ${missing_percent}% missing"
    log_info "    Report: $quality_report"
}

log_info "Quality assessment"
for method in "${METHODS[@]}"; do
    assess_quality "$method" "${WORK_DIR}/${method}" || true
done

# ---- summary ----

log_info "Summary"
failed_any=0
for method in "${METHODS[@]}"; do
    if [[ ${method_status[$method]} -eq 0 ]]; then
        log_info "  ${method}: ${copied_per_method[$method]} markers concatenated"
    else
        log_error "  ${method}: failed"
        failed_any=1
    fi
done

log_info "Output layout, one block per branch:"
log_info "  ${WORK_DIR}/<branch>/input/   filtered markers"
log_info "  ${WORK_DIR}/<branch>/output/  supermatrix and partitions"
log_info "  ${WORK_DIR}/<branch>/logs/    reports and logs"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"

if [[ $failed_any -ne 0 ]]; then
    log_error "At least one branch failed."
    exit 1
fi
