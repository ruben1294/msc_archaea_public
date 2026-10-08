#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 24. Runs anvi-compute-functional-enrichment-in-pan for every species,
# across each category (Depth, Station, ROMZ_flag) and each annotation source
# (KOfam, COG20_FUNCTION, COG20_CATEGORY, KEGG_Module).
#
# The test is Rao's score test over presence and absence. It stands in for the
# count-based methods like DESeq2 or CPM, which do not apply here because a
# pangenome records whether a gene cluster occurs in a genome, not how many
# reads mapped to it.
#
# A category with a single value inside a species is skipped, since the test
# needs at least two levels.

# --- Logging functions ---
log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $*"
}
log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $*" >&2
}
log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $*"
}

# --- Error handling ---
trap 'log_error "Script failed at line $LINENO with exit code $?"; exit 1' ERR

# --- Header ---
log_info "anvi'o9 functional enrichment pipeline"
log_info "Job ID: ${SLURM_JOB_ID:-N/A}"
log_info "Job name: ${SLURM_JOB_NAME:-interactive}"
log_info "Node: ${SLURMD_NODENAME:-$(hostname)}"
log_info "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
load_project_config

require_commands anvi-compute-functional-enrichment-in-pan

required_vars=(PROJECT_DIR RESULTS_DIR PANGENOME_SUBDIR)
for var in "${required_vars[@]}"; do
    if [[ -z "${!var:-}" ]]; then
        log_error "Required config variable '$var' is not set or is empty"
        exit 1
    fi
done

PANGENOME_DIR="$RESULTS_DIR/$PANGENOME_SUBDIR"
PANGENOME_WORK_DIR="$PANGENOME_DIR/pangenome"

if [[ ! -d "$PANGENOME_WORK_DIR" ]]; then
    log_error "Pangenome working directory does not exist: $PANGENOME_WORK_DIR"
    exit 1
fi

#### CONDA ACTIVATION (uncomment if your sbatch session does not pre-activate)
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate anvio-9

# Categories to test and annotation sources to use
CATEGORIES=(Depth Station ROMZ_flag)
ANNOTATION_SOURCES=(KOfam COG20_FUNCTION COG20_CATEGORY KEGG_Module)

# Helper: check that a per-species layers TSV has >= 2 unique values for a column
has_multiple_levels() {
    local layers_tsv="$1"
    local column="$2"
    [[ -f "$layers_tsv" ]] || return 1
    awk -F'\t' -v col="$column" '
        NR == 1 {
            for (i = 1; i <= NF; i++) if ($i == col) { idx = i; break }
            if (!idx) exit 1
            next
        }
        { vals[$idx] = 1 }
        END { exit (length(vals) >= 2 ? 0 : 1) }
    ' "$layers_tsv"
}

# ---- run enrichment for each species x category x annotation source ----
log_info "Running functional enrichment per species..."
for species_dir in "${PANGENOME_WORK_DIR}"/*/; do
    [[ -d "$species_dir" ]] || continue
    species_name=$(basename "$species_dir")
    pan_db="${species_dir}/pan_output/${species_name}-PAN.db"
    genomes_db="${species_dir}/${species_name}-GENOMES.db"
    species_layers_tsv="${species_dir}/misc_data_layers_${species_name}.tsv"
    enrichment_dir="${species_dir}/enrichment"

    if [[ ! -f "$pan_db" || ! -f "$genomes_db" ]]; then
        log_warn "  Inputs missing for ${species_name//_/ }; skipping"
        continue
    fi
    if [[ ! -f "$species_layers_tsv" ]]; then
        log_warn "  Layers TSV not found for ${species_name//_/ }; run 24 first; skipping"
        continue
    fi

    mkdir -p "$enrichment_dir"
    log_info "  Species: ${species_name//_/ }"

    for category in "${CATEGORIES[@]}"; do
        if ! has_multiple_levels "$species_layers_tsv" "$category"; then
            log_warn "    [$category] has < 2 levels for ${species_name//_/ }; skipping"
            continue
        fi
        for source in "${ANNOTATION_SOURCES[@]}"; do
            output="${enrichment_dir}/enrichment_${category}_${source}_${species_name}.tsv"
            if [[ -f "$output" ]]; then
                log_info "    [$category | $source] already exists; skipping"
                continue
            fi
            log_info "    [$category | $source] computing enrichment..."
            if anvi-compute-functional-enrichment-in-pan \
                -p "$pan_db" \
                -g "$genomes_db" \
                --category-variable "$category" \
                --annotation-source "$source" \
                -o "$output" 2>&1; then
                :
            else
                log_warn "      Enrichment failed for $category x $source on ${species_name//_/ }; continuing"
                rm -f "$output"
            fi
        done
    done
done
log_info "Enrichment runs complete"

# --- Done ---
log_info "Pipeline finished successfully"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
log_info "Outputs in: ${PANGENOME_WORK_DIR}/<species>/enrichment/"
log_info "Next, on the local machine: scripts/R/26 through 31"
