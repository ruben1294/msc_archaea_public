#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 3. Calls genes on every filtered MAG with Prodigal and writes one
# subdirectory per MAG under 04_gene_prediction/, each holding .faa proteins,
# .fna nucleotide CDS and a .gff. The .faa files feed both eggNOG-mapper
# (step 4) and DRAM (step 5).
#
# Prodigal runs in single-genome mode rather than meta mode: the input is a
# binned genome, and meta mode trades per-genome training for generality, which
# costs accuracy here. Translation table 11 is the bacterial and archaeal code.
#
# Hyatt et al. (2010) BMC Bioinformatics 11:119.

load_project_config

INPUT_DIR="${PROJECT_DIR}/${FILTERED_SUBDIR}"
OUTPUT_DIR="${PROJECT_DIR}/${GENE_PREDICTION_SUBDIR}"
LOG_DIR="${PROJECT_DIR}/${LOGS_SUBDIR}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
MAIN_LOG="${LOG_DIR}/03_gene_pred_${TIMESTAMP}.log"
SUMMARY_FILE="${OUTPUT_DIR}/03_gene_pred_summary_${TIMESTAMP}.tsv"

PRODIGAL_MODE="single"
TRANSLATION_TABLE=11

readonly GREEN='\033[0;32m'
readonly BLUE='\033[0;34m'
readonly RED='\033[0;31m'
readonly NC='\033[0m'

log_message() {
    echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$MAIN_LOG"
}

log_colored() {
    local color="$1"
    local message="$2"
    echo -e "${color}[$(date '+%Y-%m-%d %H:%M:%S')] ${message}${NC}" | tee -a "$MAIN_LOG"
}

check_dependencies() {
    log_message "Checking dependencies..."

    if ! command -v prodigal &> /dev/null; then
        log_colored "$RED" "ERROR: Prodigal not found in PATH"
        log_message "Install Prodigal or load the module that provides it"
        return 1
    fi

    # grep exits 1 when the version string does not match, and pipefail would
    # take the whole script down with it, so absorb the failure and fall back.
    local prodigal_version
    prodigal_version=$(prodigal -v 2>&1 | grep -oP 'Prodigal V\d+\.\d+\.\d+' || true)
    [[ -z "$prodigal_version" ]] && prodigal_version="unknown"

    log_message "Found: $prodigal_version"
    return 0
}

validate_prodigal_output() {
    local mag_name="$1"
    local mag_output_dir="$2"
    local required_files=("${mag_output_dir}/${mag_name}.faa"
                          "${mag_output_dir}/${mag_name}.fna"
                          "${mag_output_dir}/${mag_name}.gff")

    for file in "${required_files[@]}"; do
        if [[ ! -f "$file" ]]; then
            log_colored "$RED" "ERROR: missing output file: $file"
            return 1
        fi
    done
    return 0
}

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

if [[ ! -w "$LOG_DIR" ]]; then
    echo "Log directory is not writable: $LOG_DIR" >&2
    exit 1
fi

log_colored "$BLUE" "Prodigal gene prediction"
log_message ""
log_message "Configuration:"
log_message "  Input directory:    $INPUT_DIR"
log_message "  Output directory:   $OUTPUT_DIR"
log_message "  Log directory:      $LOG_DIR"
log_message "  Prodigal mode:      $PRODIGAL_MODE"
log_message "  Translation table:  $TRANSLATION_TABLE (bacterial and archaeal)"
log_message "  Output formats:     .faa proteins, .fna nucleotides, .gff coordinates"
log_message ""

if ! check_dependencies; then
    exit 1
fi
log_message ""

if [[ ! -d "$INPUT_DIR" ]]; then
    log_colored "$RED" "ERROR: input directory not found: $INPUT_DIR"
    exit 1
fi

mag_count=$(find "$INPUT_DIR" -maxdepth 1 -name "*.fa" -type f | wc -l)
if [[ $mag_count -eq 0 ]]; then
    log_colored "$RED" "ERROR: no *.fa files in $INPUT_DIR"
    exit 1
fi

log_message "Found $mag_count MAGs to process"
log_message ""

processed=0
failed=0
failed_mags=()
total_genes=0

for mag_file in "$INPUT_DIR"/*.fa; do
    [[ -f "$mag_file" ]] || continue

    mag_name=$(basename "$mag_file" .fa)
    mag_output_dir="${OUTPUT_DIR}/${mag_name}"
    mkdir -p "$mag_output_dir"

    log_colored "$BLUE" "Processing: $mag_name"

    mag_log="${mag_output_dir}/${mag_name}_prod.log"

    if prodigal \
        -i "$mag_file" \
        -a "${mag_output_dir}/${mag_name}.faa" \
        -d "${mag_output_dir}/${mag_name}.fna" \
        -f gff \
        -o "${mag_output_dir}/${mag_name}.gff" \
        -p "$PRODIGAL_MODE" \
        -g "$TRANSLATION_TABLE" \
        > "$mag_log" 2>&1; then

        if ! validate_prodigal_output "$mag_name" "$mag_output_dir"; then
            # Pre-increment, because the post-increment form evaluates to the old
            # value and a zero there makes set -e abort the loop on the first MAG.
            ((++failed))
            failed_mags+=("$mag_name")
            log_colored "$RED" "  [FAIL] output validation"
            log_message "Progress: $((processed + failed))/$mag_count"
            continue
        fi

        # grep -c would exit 1 on an empty file and trip pipefail; counting the
        # lines separately keeps the failure local to this MAG.
        # shellcheck disable=SC2126
        n_genes=$(grep "^>" "${mag_output_dir}/${mag_name}.faa" 2>/dev/null | wc -l)

        if [[ $n_genes -eq 0 ]]; then
            log_colored "$RED" "  [FAIL] no genes predicted"
            ((++failed))
            failed_mags+=("$mag_name")
        else
            log_colored "$GREEN" "  [OK] $n_genes genes"
            ((++processed))
            total_genes=$((total_genes + n_genes))
            cat "$mag_log" >> "$MAIN_LOG"
        fi

    else
        exit_code=$?
        log_colored "$RED" "  [FAIL] Prodigal exited $exit_code"
        ((++failed))
        failed_mags+=("$mag_name")
        [[ -f "$mag_log" ]] && cat "$mag_log" >> "$MAIN_LOG"
    fi

    log_message "Progress: $((processed + failed))/$mag_count (ok: $processed, failed: $failed)"
    log_message ""
done

log_message ""
log_colored "$BLUE" "Gene prediction finished"
log_message ""
log_message "Statistics:"
log_message "  Total MAGs:             $mag_count"
log_message "  Successfully processed: $processed"
log_message "  Failed:                 $failed"
log_message "  Total genes predicted:  $total_genes"

if [[ $processed -gt 0 ]]; then
    log_message "  Average genes per MAG:  $((total_genes / processed))"
fi

if [[ $failed -gt 0 ]]; then
    log_message ""
    log_message "Failed MAGs:"
    for mag in "${failed_mags[@]}"; do
        log_message "  $mag"
    done
fi

log_message ""
log_message "Output locations:"
log_message "  Gene predictions: $OUTPUT_DIR"
log_message "  Main log:         $MAIN_LOG"
log_message ""

log_message "Generating summary table..."
{
    echo -e "MAG_Name\tTotal_Genes\tFasta_File\tProtein_File\tNucleotide_File\tGFF_File"
    for mag_dir in "$OUTPUT_DIR"/*/; do
        [[ -d "$mag_dir" ]] || continue

        mag_name=$(basename "$mag_dir")
        faa_file="${mag_dir}/${mag_name}.faa"
        fna_file="${mag_dir}/${mag_name}.fna"
        gff_file="${mag_dir}/${mag_name}.gff"
        original_fa="${INPUT_DIR}/${mag_name}.fa"

        # Only MAGs with the full set of outputs make it into the table.
        if [[ -f "$faa_file" && -f "$fna_file" && -f "$gff_file" ]]; then
            n_genes=$(grep -c "^>" "$faa_file" 2>/dev/null || echo "0")
            echo -e "${mag_name}\t${n_genes}\t${original_fa}\t${faa_file}\t${fna_file}\t${gff_file}"
        fi
    done
} > "$SUMMARY_FILE"

log_colored "$GREEN" "Summary table: $SUMMARY_FILE"
log_message ""

if [[ $failed -eq 0 ]]; then
    log_colored "$GREEN" "All $processed MAGs completed."
    exit 0
else
    log_colored "$RED" "$failed of $mag_count MAGs failed. See $MAIN_LOG"
    exit 1
fi
