#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 6. Annotates the Prodigal proteins with eggNOG-mapper v2, one output
# directory per MAG under 06_functional_annotation/, and writes a summary table
# with the per-MAG annotation rate.
#
# tax_scope is archaea, so orthology is resolved against the archaeal subset of
# eggNOG rather than the full tree.
#
# The seed ortholog e-value is 1e-3, looser than the 1e-5 default. For scale:
# 1e-3 permissive, 1e-5 the eggNOG default, 1e-10 strict. The looser cutoff is
# paired with a minimum seed ortholog bit score of 60, which carries most of the
# filtering.

load_project_config
require_vars RESULTS_DIR GENE_PREDICTION_SUBDIR FUNCTIONAL_ANNOTATION_SUBDIR

INPUT_DIR="${RESULTS_DIR}/${GENE_PREDICTION_SUBDIR}"
OUTPUT_DIR="${RESULTS_DIR}/${FUNCTIONAL_ANNOTATION_SUBDIR}"
LOG_DIR="${PROJECT_DIR}/${LOGS_SUBDIR}"
TEMP_DIR="${PROJECT_DIR}/${TEMP_SUBDIR}"
TIMESTAMP=$(date +%Y%m%d)
MAIN_LOG="${LOG_DIR}/05_func_annot_egg_${TIMESTAMP}.log"
SUMMARY_FILE="${OUTPUT_DIR}/05_func_annot_egg_summ_${TIMESTAMP}.tsv"

CPU_CORES="${SLURM_CPUS_PER_TASK:-4}"
EVALUE_THRESHOLD="1e-3"
MIN_ALIGNMENT_SCORE=60
MODE="diamond"          # diamond, mmseqs2 or hmmer

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$MAIN_LOG"
}

log "eggNOG-mapper v2 functional annotation"
log ""
log "Input:        $INPUT_DIR"
log "Output:       $OUTPUT_DIR"
log "Log file:     $MAIN_LOG"
log "Host:         $(hostname)"
log "SLURM job ID: ${SLURM_JOB_ID:-NA}"
log "CPUs:         $CPU_CORES"
log "E-value:      $EVALUE_THRESHOLD"
log "Min score:    $MIN_ALIGNMENT_SCORE"
log ""

if ! command -v emapper.py &> /dev/null; then
    log "ERROR: eggNOG-mapper not found in PATH"
    log "Load the module or activate the conda environment that provides it"
    exit 1
fi

emapper_version=$(emapper.py --version 2>&1 | head -1 || echo "unknown")
log "Found: $emapper_version"
log ""

# An empty EGGNOG_DATA_DIR makes emapper.py fall back to its own default
# location, so this warns instead of exiting.
if [[ ! -d "${EGGNOG_DATA_DIR:-}" ]]; then
    log "WARNING: eggNOG database not found at '${EGGNOG_DATA_DIR:-unset}'"
    log "Falling back to the emapper.py default location"
    EGGNOG_DATA_DIR=""
fi

# Step 05 writes one directory per MAG, so the .faa files sit exactly two levels
# below INPUT_DIR.
log "Searching for protein files..."
faa_files=()
while IFS= read -r faa_file; do
    faa_files+=("$faa_file")
done < <(find "$INPUT_DIR" -mindepth 2 -maxdepth 2 -name "*.faa" -type f | sort)

mag_count=${#faa_files[@]}
if [[ $mag_count -eq 0 ]]; then
    log "ERROR: no .faa files under $INPUT_DIR"
    exit 1
fi

log "Found $mag_count .faa files to annotate"
log ""

processed=0
failed=0
declare -a failed_mags
total_annotations=0
counter=0

for faa_file in "${faa_files[@]}"; do
    ((++counter))

    # Path shape is <INPUT_DIR>/<MAG_NAME>/<MAG_NAME>.faa
    mag_name=$(basename "$(dirname "$faa_file")")
    mag_output_dir="${OUTPUT_DIR}/${mag_name}"
    mkdir -p "$mag_output_dir"

    log "[$counter/$mag_count] Annotating $mag_name..."

    n_proteins=$(grep -c "^>" "$faa_file" 2>/dev/null || echo "0")
    log "  Proteins: $n_proteins"

    if [[ $n_proteins -eq 0 ]]; then
        log "  [SKIP] no proteins"
        ((++failed))
        failed_mags+=("$mag_name")
        continue
    fi

    output_prefix="${mag_output_dir}/${mag_name}"
    emapper_log="${mag_output_dir}/${mag_name}_egg.log"

    emapper_cmd=(
        "emapper.py"
        "-i" "$faa_file"
        "-o" "${mag_name}"
        "--output_dir" "$mag_output_dir"
        "-m" "$MODE"
        "--itype" "proteins"
        "--tax_scope" "archaea"
        "--cpu" "$CPU_CORES"
        "--seed_ortholog_score" "$MIN_ALIGNMENT_SCORE"
        "--target_orthologs" "all"
        "--seed_ortholog_evalue" "$EVALUE_THRESHOLD"
        "--pfam_realign" "realign"
        "--override"
        "--temp_dir" "$TEMP_DIR"
    )

    if [[ -n "$EGGNOG_DATA_DIR" ]]; then
        emapper_cmd+=("--data_dir" "$EGGNOG_DATA_DIR")
    fi

    if "${emapper_cmd[@]}" >> "$emapper_log" 2>&1; then

        annot_file="${output_prefix}.emapper.annotations"

        if [[ -f "$annot_file" ]]; then
            # The emapper output carries # comment lines above and below the
            # table. grep -c would exit 1 on a file that is all comments and
            # pipefail would take the script with it.
            # shellcheck disable=SC2126
            n_annotated=$(grep -v "^#" "$annot_file" | wc -l)

            if [[ $n_annotated -gt 0 ]]; then
                pct_annotated=$(awk "BEGIN {printf \"%.1f\", ($n_annotated/$n_proteins)*100}")
                log "  [OK] $n_annotated annotated (${pct_annotated}%)"
                ((++processed))
                ((total_annotations+=n_annotated))
            else
                log "  [FAIL] no annotations generated"
                ((++failed))
                failed_mags+=("$mag_name")
            fi
        else
            log "  [FAIL] output file missing"
            ((++failed))
            failed_mags+=("$mag_name")
        fi
    else
        exit_code=$?
        log "  [FAIL] eggNOG-mapper exited $exit_code"
        ((++failed))
        failed_mags+=("$mag_name")
    fi

    log ""
done

log "Summary"
log "  Total MAGs:          $mag_count"
log "  Successful:          $processed"
log "  Failed:              $failed"
log "  Total annotations:   $total_annotations"

if [[ $processed -gt 0 ]]; then
    log "  Avg annotations/MAG: $((total_annotations / processed))"
fi

if [[ $failed -gt 0 ]]; then
    log ""
    log "Failed MAGs:"
    for mag in "${failed_mags[@]}"; do
        log "  $mag"
    done
fi

log ""
log "Output: $OUTPUT_DIR"
log "Log:    $MAIN_LOG"
log ""
log "Generating summary table..."

{
    echo -e "MAG_Name\tTotal_Proteins\tAnnotated_Proteins\tPercent_Annotated\tAnnotation_File\tSeed_Orthologs_File"

    for mag_dir in "$OUTPUT_DIR"/*/; do
        [[ -d "$mag_dir" ]] || continue

        mag_name=$(basename "$mag_dir")
        annot_file="${mag_dir}/${mag_name}.emapper.annotations"
        seed_file="${mag_dir}/${mag_name}.emapper.seed_orthologs"

        if [[ -f "$annot_file" ]]; then
            # shellcheck disable=SC2126
            n_annotated=$(grep -v "^#" "$annot_file" | wc -l)

            faa_file="${INPUT_DIR}/${mag_name}/${mag_name}.faa"
            n_proteins=$(grep -c "^>" "$faa_file" 2>/dev/null || echo "0")

            if [[ $n_proteins -gt 0 ]]; then
                pct=$(awk "BEGIN {printf \"%.2f\", ($n_annotated/$n_proteins)*100}")
            else
                pct="0.00"
            fi

            # Paths go in relative to the project root so the table survives a
            # move to another machine.
            rel_annot_file="${annot_file/#${PROJECT_DIR}\//}"
            rel_seed_file="${seed_file/#${PROJECT_DIR}\//}"
            echo -e "${mag_name}\t${n_proteins}\t${n_annotated}\t${pct}\t${rel_annot_file}\t${rel_seed_file}"
        fi
    done
} > "$SUMMARY_FILE"

log "Summary table: $SUMMARY_FILE"
log ""

if [[ $failed -eq 0 ]]; then
    log "All $processed MAGs annotated."
    exit 0
else
    log "$failed of $mag_count MAGs failed. See $MAIN_LOG"
    exit 1
fi
