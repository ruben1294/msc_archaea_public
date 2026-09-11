#!/usr/bin/env bash
set -euo pipefail

# Step 7. Turns the KO table from step 6 into pathway completeness matrices with
# KEGG-Decoder, then optionally extends them with KEGG-Expander and merges both
# into a single heatmap.
#
# Decoder runs twice over the same input. Static writes function_heatmap.svg,
# one row per distinct identifier before the underscore in the FASTA header.
# Interactive writes the same matrix as function_heatmap.html, which tolerates
# larger genome sets because values appear on hover and regions can be zoomed.
# Decoder also offers a tanglegram mode that compares the clustered metabolic
# output against a user-supplied Newick tree; it needs more than three genomes
# and matching names, and this pipeline does not use it.
#
# The Expander branch is optional and skips itself when the concatenated protein
# file, the HMM database or hmmsearch are missing.

log_info()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $*"; }
log_error() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2; }
log_warn()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $*"; }
log_debug() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] [DEBUG] $*"; }

trap 'log_error "Failed at line $LINENO with exit code $?"; exit 1' ERR

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${REPO_ROOT}/configs/server.env}"

if [[ ! -f "$CONFIG_FILE" ]]; then
    log_error "No configuration at $CONFIG_FILE"
    log_error "Copy configs/server.env.example to configs/server.env and edit it."
    exit 1
fi
# shellcheck source=/dev/null
source "$CONFIG_FILE"
log_info "Loaded configuration from $CONFIG_FILE"

OUTPUT_DIR="${PROJECT_DIR}/${METABOLIC_RECONSTRUCTION_SUBDIR}"
LOG_DIR="${PROJECT_DIR}/${LOGS_SUBDIR}"
TEMP_DIR="${PROJECT_DIR}/${TEMP_SUBDIR}"
TIMESTAMP=$(date +%Y%m%d)
MAIN_LOG="${LOG_DIR}/07_metabolic_reconstruction_${TIMESTAMP}.log"
SUMMARY_FILE="${OUTPUT_DIR}/07_metabolic_reconstruction_summ_${TIMESTAMP}_${SLURM_JOB_ID:-local}.tsv"

DECODER_OUTPUT_STATIC="${OUTPUT_DIR}/keggdecoder_static"
DECODER_OUTPUT_INTERACTIVE="${OUTPUT_DIR}/keggdecoder_interactive"
EXPANDER_OUTPUT="${OUTPUT_DIR}/keggexpander_output"
COMBINED_OUTPUT="${OUTPUT_DIR}/decoder_expander_combined"

CPU_CORES="${SLURM_CPUS_PER_TASK:-8}"

# The KofamKOALA submission date is part of the directory name, so find the
# table rather than hard-coding the date.
INPUT_FILE="${INPUT_FILE:-}"
if [[ -z "$INPUT_FILE" ]]; then
    INPUT_FILE=$(find "${PROJECT_DIR}/${KO_ASSIGNMENT_SUBDIR}" -maxdepth 2 \
                      -name "all_kofam_KOs.transformed.tsv" -type f | sort | tail -1)
fi

PROTEIN_FILE="${PROTEIN_FILE:-${PROJECT_DIR}/${GENE_PREDICTION_SUBDIR}/proteins/CMRA_all_proteins_final.faa}"

# KEGG-Expander ships inside the KEGGDecoder package rather than on PATH, so ask
# Python where the package lives instead of naming a conda prefix. Override
# KEGGDECODER_DIR to point somewhere else.
KEGGDECODER_DIR="${KEGGDECODER_DIR:-}"
if [[ -z "$KEGGDECODER_DIR" ]]; then
    KEGGDECODER_DIR=$(python -c \
        'import os, KEGGDecoder; print(os.path.dirname(KEGGDecoder.__file__))' 2>/dev/null || true)
fi

EXPANDER_HMM="${EXPANDER_HMM:-${KEGGDECODER_DIR}/HMM_Models/expander_dbv0.7.hmm}"
# The upstream repository calls this KEGG-expander.py, but the file that ships
# with KEGG-Decoder uses an underscore. Naming it with the hyphen finds nothing.
KEGG_EXPANDER_SCRIPT="${KEGG_EXPANDER_SCRIPT:-${KEGGDECODER_DIR}/KEGG_expander.py}"
KEGG_DECODE_AND_EXPANDER_SCRIPT="${KEGG_DECODE_AND_EXPANDER_SCRIPT:-${KEGGDECODER_DIR}/Decode_and_Expand.py}"

mkdir -p "$OUTPUT_DIR" "$DECODER_OUTPUT_STATIC" "$DECODER_OUTPUT_INTERACTIVE" \
         "$EXPANDER_OUTPUT" "$COMBINED_OUTPUT" "$LOG_DIR" "$TEMP_DIR"

log_info "KEGG-Decoder metabolic pathway reconstruction"
log_info ""
log_info "Configuration:"
log_info "  Input file:       $INPUT_FILE"
log_info "  Output directory: $OUTPUT_DIR"
log_info "  Log file:         $MAIN_LOG"
log_info "  CPUs:             $CPU_CORES"
log_info "  Host:             $(hostname)"
log_info "  SLURM job ID:     ${SLURM_JOB_ID:-NA}"
log_info ""

if ! command -v KEGG-decoder &> /dev/null; then
    log_error "KEGG-decoder not found in PATH"
    log_error "Install with: conda install -c bioconda kegg-decoder"
    exit 1
fi
log_info "Found KEGG-decoder"

if [[ -z "$INPUT_FILE" || ! -f "$INPUT_FILE" ]]; then
    log_error "No all_kofam_KOs.transformed.tsv found"
    log_error "Run scripts/pipeline/06_KO_assignment.sh to check the step 6 output."
    exit 1
fi

n_lines=$(wc -l < "$INPUT_FILE")
# Gene identifiers arrive as GENOME_XXX, so the genome name is the field before
# the first underscore.
n_genomes=$(cut -f1 "$INPUT_FILE" | cut -d'_' -f1 | sort -u | wc -l)
log_info "  Input lines:      $n_lines"
log_info "  Genomes detected: $n_genomes"
log_info ""

# ---- KEGG-Decoder, static ----

log_info "Step 1a: KEGG-Decoder, static"

decoder_static_list="${DECODER_OUTPUT_STATIC}/keggdecoder_output.list"
decoder_static_log="${DECODER_OUTPUT_STATIC}/keggdecoder.log"

if KEGG-decoder \
    --input "$INPUT_FILE" \
    --output "$decoder_static_list" \
    --vizoption "static" \
    > "$decoder_static_log" 2>&1; then

    log_info "  [OK] static mode finished"

    if [[ -f "$decoder_static_list" ]]; then
        n_pathways=$(head -1 "$decoder_static_list" | tr '\t' '\n' | wc -l)
        log_info "  Pathways analysed: $((n_pathways - 1))"
        log_info "  Matrix:            $decoder_static_list"
    fi

    svg_file="${decoder_static_list/.list/.svg}"
    if [[ -f "$svg_file" ]]; then
        log_info "  Heatmap (SVG):     $svg_file"
    else
        log_warn "  SVG heatmap not generated"
    fi
else
    exit_code=$?
    log_error "KEGG-Decoder static failed with exit code $exit_code"
    log_error "See $decoder_static_log"
    exit 1
fi

log_info ""

# ---- KEGG-Decoder, interactive ----

log_info "Step 1b: KEGG-Decoder, interactive"

decoder_interactive_list="${DECODER_OUTPUT_INTERACTIVE}/keggdecoder_output.list"
decoder_interactive_log="${DECODER_OUTPUT_INTERACTIVE}/keggdecoder.log"

if KEGG-decoder \
    --input "$INPUT_FILE" \
    --output "$decoder_interactive_list" \
    --vizoption "interactive" \
    > "$decoder_interactive_log" 2>&1; then

    log_info "  [OK] interactive mode finished"

    if [[ -f "$decoder_interactive_list" ]]; then
        n_pathways=$(head -1 "$decoder_interactive_list" | tr '\t' '\n' | wc -l)
        log_info "  Pathways analysed: $((n_pathways - 1))"
        log_info "  Matrix:            $decoder_interactive_list"
    fi

    html_file="${decoder_interactive_list/.list/.html}"
    if [[ -f "$html_file" ]]; then
        log_info "  Heatmap (HTML):    $html_file"
    else
        log_warn "  HTML heatmap not generated"
    fi
else
    exit_code=$?
    log_error "KEGG-Decoder interactive failed with exit code $exit_code"
    log_error "See $decoder_interactive_log"
    exit 1
fi

log_info ""

# ---- KEGG-Expander ----

run_expander=false

if [[ -f "$PROTEIN_FILE" ]] && [[ -f "$EXPANDER_HMM" ]]; then
    run_expander=true
    log_info "Step 2: KEGG-Expander"

    for tool in hmmsearch hmmpress; do
        if ! command -v "$tool" &> /dev/null; then
            log_warn "$tool not found, skipping KEGG-Expander"
            run_expander=false
        fi
    done

    if [[ ! -f "$KEGG_EXPANDER_SCRIPT" ]]; then
        log_warn "KEGG_expander.py not found at $KEGG_EXPANDER_SCRIPT"
        run_expander=false
    fi
else
    log_info "Skipping KEGG-Expander: protein file or HMM database missing"
fi

if [[ "$run_expander" == true ]]; then
    log_info ""
    log_info "Step 2a: preparing the HMM database"

    # hmmpress writes its index next to the .hmm, so the database is copied out
    # of the read-only package directory first.
    LOCAL_HMM_DIR="${PROJECT_DIR}/databases/hmm_models"
    mkdir -p "$LOCAL_HMM_DIR"
    LOCAL_HMM="${LOCAL_HMM_DIR}/expander_dbv0.7.hmm"

    if [[ ! -f "$LOCAL_HMM" ]] || [[ "$EXPANDER_HMM" -nt "$LOCAL_HMM" ]]; then
        log_info "  Copying HMM database into the project tree"
        if cp "$EXPANDER_HMM" "$LOCAL_HMM" 2>/dev/null; then
            log_info "  [OK] copied to $LOCAL_HMM"
        else
            log_error "Could not copy $EXPANDER_HMM to $LOCAL_HMM"
            run_expander=false
        fi
    else
        log_info "  Local HMM database is up to date"
    fi

    if [[ "$run_expander" == true ]]; then
        hmm_base="${LOCAL_HMM%.hmm}"
        needs_pressing=false

        if [[ ! -f "${hmm_base}.h3i" ]] || [[ ! -f "${hmm_base}.h3f" ]] || \
           [[ ! -f "${hmm_base}.h3m" ]] || [[ ! -f "${hmm_base}.h3p" ]]; then
            needs_pressing=true
            log_info "  HMM database not pressed, running hmmpress"
        elif [[ "$LOCAL_HMM" -nt "${hmm_base}.h3i" ]]; then
            needs_pressing=true
            log_info "  HMM database changed, re-running hmmpress"
        else
            log_info "  HMM database already pressed"
        fi

        if [[ "$needs_pressing" == true ]]; then
            hmmpress_log="${EXPANDER_OUTPUT}/hmmpress.log"

            if (cd "$LOCAL_HMM_DIR" && hmmpress -f "$LOCAL_HMM") > "$hmmpress_log" 2>&1; then
                log_info "  [OK] hmmpress finished"

                if [[ -f "${hmm_base}.h3i" ]] && [[ -f "${hmm_base}.h3f" ]] && \
                   [[ -f "${hmm_base}.h3m" ]] && [[ -f "${hmm_base}.h3p" ]]; then
                    log_debug "  Index files present in $LOCAL_HMM_DIR"
                else
                    log_warn "  hmmpress reported success but the index is incomplete"
                    for ext in h3i h3f h3m h3p; do
                        log_debug "    $ext: $([ -f "${hmm_base}.${ext}" ] && echo present || echo missing)"
                    done
                fi
            else
                exit_code=$?
                log_error "hmmpress failed with exit code $exit_code, see $hmmpress_log"
                run_expander=false
            fi
        fi

        EXPANDER_HMM="$LOCAL_HMM"
    fi
fi

if [[ "$run_expander" == true ]]; then
    log_info ""
    log_info "Step 2b: hmmsearch"

    expander_tbl="${EXPANDER_OUTPUT}/keggexpander.tbl"
    expander_log="${EXPANDER_OUTPUT}/hmmsearch.log"

    if hmmsearch \
        --tblout "$expander_tbl" \
        --cpu "$CPU_CORES" \
        -T 75 \
        "$EXPANDER_HMM" \
        "$PROTEIN_FILE" \
        > "$expander_log" 2>&1; then

        log_info "  [OK] hmmsearch finished"

        n_hits=$(grep -cv "^#" "$expander_tbl" 2>/dev/null || echo "0")
        log_info "  HMM hits: $n_hits"
        log_info "  Output:   $expander_tbl"

        if [[ $n_hits -eq 0 ]]; then
            log_warn "No HMM hits, skipping KEGG_expander.py"
            run_expander=false
        fi
    else
        exit_code=$?
        log_error "hmmsearch failed with exit code $exit_code, see $expander_log"
        tail -20 "$expander_log" 2>/dev/null | while read -r line; do
            log_error "  $line"
        done
        run_expander=false
    fi
fi

if [[ "$run_expander" == true ]]; then
    log_info ""
    log_info "Step 2c: KEGG_expander.py"

    expander_list="${EXPANDER_OUTPUT}/HMM_OUT.list"
    expander_py_log="${EXPANDER_OUTPUT}/keggexpander_py.log"

    # KEGG_expander.py writes its outputs into the current directory, so the
    # call runs from EXPANDER_OUTPUT in a subshell.
    if (cd "$EXPANDER_OUTPUT" && python "$KEGG_EXPANDER_SCRIPT" "$expander_tbl" "HMM_OUT.list") \
        > "$expander_py_log" 2>&1; then

        log_info "  [OK] KEGG_expander.py finished"

        if [[ -f "$expander_list" ]]; then
            n_funcs=$(tail -n +2 "$expander_list" 2>/dev/null | wc -l)
            log_info "  Expanded functions: $n_funcs"
            log_info "  Expander list:      $expander_list"
        else
            log_warn "  Expander list not generated at $expander_list"
            run_expander=false
        fi

        expander_heatmap="${EXPANDER_OUTPUT}/function_heatmap.svg"
        if [[ -f "$expander_heatmap" ]]; then
            log_info "  Heatmap (SVG):      $expander_heatmap"
        else
            log_warn "  Expander heatmap not generated"
        fi
    else
        exit_code=$?
        log_error "KEGG_expander.py failed with exit code $exit_code, see $expander_py_log"
        tail -20 "$expander_py_log" 2>/dev/null | while read -r line; do
            log_error "  $line"
        done
        run_expander=false
    fi
fi

log_info ""

# ---- combine both matrices ----

if [[ "$run_expander" == true ]] && [[ -f "${expander_list:-}" ]]; then
    log_info "Step 3: combining Decoder and Expander"

    if [[ -f "$KEGG_DECODE_AND_EXPANDER_SCRIPT" ]]; then
        combined_log="${COMBINED_OUTPUT}/combined.log"

        if (cd "$COMBINED_OUTPUT" && \
            python "$KEGG_DECODE_AND_EXPANDER_SCRIPT" "$decoder_static_list" "$expander_list") \
            > "$combined_log" 2>&1; then

            log_info "  [OK] combined heatmap generated"

            combined_heatmap="${COMBINED_OUTPUT}/decode-expand_heatmap.svg"
            if [[ -f "$combined_heatmap" ]]; then
                log_info "  Combined heatmap: $combined_heatmap"
            fi
        else
            exit_code=$?
            log_error "Decode_and_Expand.py failed with exit code $exit_code, see $combined_log"
            tail -20 "$combined_log" 2>/dev/null | while read -r line; do
                log_error "  $line"
            done
        fi
    else
        log_warn "Decode_and_Expand.py not found at $KEGG_DECODE_AND_EXPANDER_SCRIPT"
    fi
fi

log_info ""

# ---- summary ----

log_info "Summary"
log_info "  Input genomes: $n_genomes"
log_info ""
log_info "KEGG-Decoder:"
log_info "  static matrix:       $decoder_static_list"
log_info "  static heatmap:      ${decoder_static_list/.list/.svg}"
log_info "  interactive matrix:  $decoder_interactive_list"
log_info "  interactive heatmap: ${decoder_interactive_list/.list/.html}"
log_info ""

if [[ "$run_expander" == true ]]; then
    log_info "KEGG-Expander:"
    log_info "  HMM results:   $expander_tbl"
    log_info "  Expander list: $expander_list"
    log_info ""
fi

{
    echo -e "Analysis\tMode\tOutput_File\tFormat\tDescription"
    echo -e "KEGG-Decoder\tstatic\t${decoder_static_list}\tTSV\tPathway completeness matrix"
    echo -e "KEGG-Decoder\tstatic\t${decoder_static_list/.list/.svg}\tSVG\tStatic heatmap"
    echo -e "KEGG-Decoder\tinteractive\t${decoder_interactive_list}\tTSV\tPathway completeness matrix"
    echo -e "KEGG-Decoder\tinteractive\t${decoder_interactive_list/.list/.html}\tHTML\tInteractive heatmap"

    if [[ "$run_expander" == true ]] && [[ -f "${expander_list:-}" ]]; then
        echo -e "KEGG-Expander\t-\t${expander_list}\tTSV\tExpanded function matrix"
        echo -e "KEGG-Expander\t-\t${expander_tbl}\tTBL\thmmsearch output"
    fi

    if [[ -f "${COMBINED_OUTPUT}/decode-expand_heatmap.svg" ]]; then
        echo -e "Combined\t-\t${COMBINED_OUTPUT}/decode-expand_heatmap.svg\tSVG\tCombined Decoder and Expander heatmap"
    fi
} > "$SUMMARY_FILE"

log_info "Summary table:    $SUMMARY_FILE"
log_info "Output directory: $OUTPUT_DIR"
log_info "Log file:         $MAIN_LOG"

exit 0
