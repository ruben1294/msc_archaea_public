#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 5. Annotates each MAG with DRAM v1 and, at the end, runs DRAM distill
# over the whole set to produce metabolism_summary.xlsx and product.html.
#
# Step 8 runs DRAM2 through Nextflow over the same MAGs and supersedes this for
# the thesis figures. Step 5 stays because it produced the first pass and
# because DRAM v1 and DRAM2 disagree on enough trait calls that the comparison
# is worth keeping reproducible.
#
# Two modes, selected through the DRAM_MODE environment variable:
#   annotate        starts from .fa contigs and calls genes itself
#   annotate_genes  starts from the .faa proteins that step 3 already predicted
#
# In annotate mode DRAM runs Prodigal in meta mode, unlike step 3 which uses
# single mode. Gene sets from the two paths are therefore not identical.

DRAM_MODE="${DRAM_MODE:-annotate}"

load_project_config

if ! command -v DRAM.py &> /dev/null; then
    echo "DRAM.py not found in PATH" >&2
    exit 1
fi

echo "DRAM version: $(DRAM-setup.py version 2>&1)"
echo ""

case "$DRAM_MODE" in
    annotate)
        INPUT_DIR="${PROJECT_DIR}/${FILTERED_SUBDIR}"
        INPUT_EXTENSION="fa"
        INPUT_TYPE="genomic DNA"
        ;;
    annotate_genes)
        INPUT_DIR="${PROJECT_DIR}/${GENE_PREDICTION_SUBDIR}"
        INPUT_EXTENSION="faa"
        INPUT_TYPE="predicted proteins"
        ;;
    *)
        echo "DRAM_MODE must be 'annotate' or 'annotate_genes', got '$DRAM_MODE'" >&2
        exit 1
        ;;
esac

OUTPUT_DIR="${PROJECT_DIR}/${FUNCTIONAL_ANNOTATION_SUBDIR}/dram"
LOG_DIR="${PROJECT_DIR}/${LOGS_SUBDIR}"
TEMP_DIR="${PROJECT_DIR}/${TEMP_SUBDIR}/dram"
TIMESTAMP=$(date +%Y%m%d)
MAIN_LOG="${LOG_DIR}/05_func_annot_dram_${DRAM_MODE}_${TIMESTAMP}.log"
SUMMARY_FILE="${OUTPUT_DIR}/05_func_annot_dram_summ_${TIMESTAMP}.tsv"

DRAM_DB_DIR="${DRAM_DB_DIR:-}"
CPU_CORES="${SLURM_CPUS_PER_TASK:-4}"
MIN_CONTIG_SIZE=2500
PRODIGAL_MODE="meta"
USE_UNIREF=true
RUN_DISTILL=true

mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$TEMP_DIR"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$MAIN_LOG"
}

log "DRAM v1 functional annotation"
log ""
log "Mode:            $DRAM_MODE"
log "Input type:      $INPUT_TYPE"
log "Input:           $INPUT_DIR"
log "Output:          $OUTPUT_DIR"
log "Log file:        $MAIN_LOG"
log "Host:            $(hostname)"
log "SLURM job ID:    ${SLURM_JOB_ID:-NA}"
log "CPUs:            $CPU_CORES"
log "Min contig size: ${MIN_CONTIG_SIZE} bp"
if [[ "$DRAM_MODE" == "annotate" ]]; then
    log "Prodigal mode:   $PRODIGAL_MODE"
fi
log "Use UniRef90:    $USE_UNIREF"
log ""

# An empty DRAM_DB_DIR makes DRAM fall back to the location in its own CONFIG.
if [[ -n "$DRAM_DB_DIR" && ! -d "$DRAM_DB_DIR" ]]; then
    log "WARNING: DRAM database not found at $DRAM_DB_DIR"
    log "Falling back to the location in DRAM's own CONFIG"
    DRAM_DB_DIR=""
fi

# In annotate mode the .fa files sit one level down; in annotate_genes mode the
# .faa files sit two levels down, one directory per MAG from step 3.
log "Searching for .${INPUT_EXTENSION} files..."
input_files=()
if [[ "$DRAM_MODE" == "annotate" ]]; then
    while IFS= read -r input_file; do
        input_files+=("$input_file")
    done < <(find "$INPUT_DIR" -mindepth 1 -maxdepth 2 -name "*.${INPUT_EXTENSION}" -type f | sort)
else
    while IFS= read -r input_file; do
        input_files+=("$input_file")
    done < <(find "$INPUT_DIR" -mindepth 2 -maxdepth 2 -name "*.${INPUT_EXTENSION}" -type f | sort)
fi

mag_count=${#input_files[@]}
if [[ $mag_count -eq 0 ]]; then
    log "ERROR: no .${INPUT_EXTENSION} files under $INPUT_DIR"
    exit 1
fi

log "Found $mag_count .${INPUT_EXTENSION} files to annotate"
log ""

processed=0
failed=0
declare -a failed_mags
total_genes=0
counter=0

for input_file in "${input_files[@]}"; do
    ((++counter))

    if [[ "$DRAM_MODE" == "annotate" ]]; then
        mag_name=$(basename "$input_file" .fa)
    else
        mag_name=$(basename "$(dirname "$input_file")")
    fi

    mag_output_dir="${OUTPUT_DIR}/${mag_name}"
    mkdir -p "$mag_output_dir"

    log "[$counter/$mag_count] Processing $mag_name..."

    n_seqs=$(grep -c "^>" "$input_file" 2>/dev/null || echo "0")
    if [[ "$DRAM_MODE" == "annotate" ]]; then
        log "  Contigs: $n_seqs"
    else
        log "  Proteins: $n_seqs"
    fi

    if [[ $n_seqs -eq 0 ]]; then
        log "  [SKIP] no sequences"
        ((++failed))
        failed_mags+=("$mag_name")
        continue
    fi

    dram_log="${mag_output_dir}/${mag_name}_dram.log"

    dram_cmd=(
        "DRAM.py" "${DRAM_MODE}"
        "-i" "$input_file"
        "-o" "$mag_output_dir"
        "--threads" "$CPU_CORES"
        "--verbose"
    )

    if [[ "$DRAM_MODE" == "annotate" ]]; then
        dram_cmd+=(
            "--min_contig_size" "$MIN_CONTIG_SIZE"
            "--prodigal_mode" "$PRODIGAL_MODE"
        )
    fi

    if [[ "$USE_UNIREF" == true ]]; then
        dram_cmd+=("--use_uniref")
    fi

    if [[ -n "$DRAM_DB_DIR" ]]; then
        dram_cmd+=("--config_loc" "${DRAM_DB_DIR}/CONFIG")
    fi

    log "  Running DRAM ${DRAM_MODE}..."
    if "${dram_cmd[@]}" >> "$dram_log" 2>&1; then
        annot_file="${mag_output_dir}/annotations.tsv"

        if [[ -f "$annot_file" ]]; then
            n_annotated=$(tail -n +2 "$annot_file" | wc -l)

            if [[ $n_annotated -gt 0 ]]; then
                log "  [OK] $n_annotated genes annotated"
                ((++processed))
                ((total_genes+=n_annotated))
            else
                log "  [FAIL] no genes annotated"
                ((++failed))
                failed_mags+=("$mag_name")
            fi
        else
            log "  [FAIL] annotations.tsv missing"
            ((++failed))
            failed_mags+=("$mag_name")
        fi
    else
        exit_code=$?
        log "  [FAIL] DRAM exited $exit_code"
        ((++failed))
        failed_mags+=("$mag_name")
    fi

    log ""
done

if [[ "$RUN_DISTILL" == true ]] && [[ $processed -gt 0 ]]; then
    log "Running DRAM distill over $processed MAGs"

    distill_output="${OUTPUT_DIR}/distillate"
    mkdir -p "$distill_output"

    annot_files_str=""
    for mag_dir in "$OUTPUT_DIR"/*/; do
        [[ -d "$mag_dir" ]] || continue
        annot_file="${mag_dir}/annotations.tsv"
        if [[ -f "$annot_file" ]]; then
            annot_files_str+="${annot_file} "
        fi
    done

    if [[ -n "$annot_files_str" ]]; then
        distill_log="${distill_output}/distill.log"

        # -i takes several paths, so this one expansion stays unquoted on purpose.
        # shellcheck disable=SC2086
        if DRAM.py distill \
            -i $annot_files_str \
            -o "$distill_output" \
            --rrna_path "${OUTPUT_DIR}/*/rrnas.tsv" \
            --trna_path "${OUTPUT_DIR}/*/trnas.tsv" \
            >> "$distill_log" 2>&1; then
            log "  [OK] distillation finished"
            log "  Product:  ${distill_output}/product.html"
            log "  Summary:  ${distill_output}/metabolism_summary.xlsx"
        else
            # Annotation is the deliverable; a failed distill does not invalidate it.
            log "  [FAIL] distillation failed, continuing"
        fi
    fi
    log ""
fi

log "Summary"
log "  Total MAGs:            $mag_count"
log "  Successful:            $processed"
log "  Failed:                $failed"
log "  Total annotated genes: $total_genes"

if [[ $processed -gt 0 ]]; then
    log "  Avg genes/MAG:         $((total_genes / processed))"
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
    echo -e "MAG_Name\tTotal_Genes\tAnnotated_Genes\tAnnotation_File\tGFF_File\tGenBank_File"

    for mag_dir in "$OUTPUT_DIR"/*/; do
        [[ -d "$mag_dir" ]] || continue
        mag_name=$(basename "$mag_dir")

        # distillate holds the cross-MAG summary, not a genome.
        [[ "$mag_name" == "distillate" ]] && continue

        annot_file="${mag_dir}/annotations.tsv"
        gff_file="${mag_dir}/genes.gff"
        gbk_file="${mag_dir}/genes.gbk"

        if [[ -f "$annot_file" ]]; then
            n_genes=$(tail -n +2 "$annot_file" | wc -l)

            rel_annot="${annot_file/#${PROJECT_DIR}\//}"
            rel_gff="${gff_file/#${PROJECT_DIR}\//}"
            rel_gbk="${gbk_file/#${PROJECT_DIR}\//}"

            echo -e "${mag_name}\t${n_genes}\t${n_genes}\t${rel_annot}\t${rel_gff}\t${rel_gbk}"
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
