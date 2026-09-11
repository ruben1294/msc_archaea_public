#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 24. Builds one TSV per species carrying each MAG's Station, Depth,
# ROMZ_flag and Order from MAGs_metadata.tsv, and imports them into the
# corresponding PAN.db as layer misc-data.
#
# Step 25 cannot compute enrichment without these layers, since they are the
# categories it tests against. They also appear as coloured sidebars in
# anvi-display-pan.
#
# MAG names are sanitised here exactly as step 14 does it,
# ${mag_name//[^a-zA-Z0-9_]/_}. Any drift between the two leaves the layers
# keyed to names the pangenome does not recognise, and the import silently
# matches nothing.

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
log_info "anvi'o9 misc-data layers import pipeline"
log_info "Job ID: ${SLURM_JOB_ID:-N/A}"
log_info "Job name: ${SLURM_JOB_NAME:-interactive}"
log_info "Node: ${SLURMD_NODENAME:-$(hostname)}"
log_info "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
load_project_config

require_commands anvi-import-misc-data anvi-delete-misc-data

required_vars=(PROJECT_DIR PANGENOME_SUBDIR MAGS_METADATA)
for var in "${required_vars[@]}"; do
    if [[ -z "${!var:-}" ]]; then
        log_error "Required config variable '$var' is not set or is empty"
        exit 1
    fi
done

PANGENOME_DIR="$PROJECT_DIR/$PANGENOME_SUBDIR"
PANGENOME_WORK_DIR="$PANGENOME_DIR/pangenome"

if [[ ! -d "$PANGENOME_WORK_DIR" ]]; then
    log_error "Pangenome working directory does not exist: $PANGENOME_WORK_DIR"
    exit 1
fi
if [[ ! -f "$MAGS_METADATA" ]]; then
    log_error "Metadata file not found: $MAGS_METADATA"
    exit 1
fi

#### CONDA ACTIVATION (uncomment if your sbatch session does not pre-activate)
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate anvio-9

# Step 1: Generate the global misc-data layers TSV from MAGs_metadata.tsv.
# Columns of MAGs_metadata.tsv (1-indexed):
#   1=FullName, 2=NameWithoutExtension, 3=DisplayName, 4=Station, 5=Depth,
#   6=Assembler, 7=BinID, 8=Genus, 9=Species, 10=FullTaxonomy, 11=OriginalName
# We sanitize NameWithoutExtension to match the genome key used in
# external_genomes.tsv (script 14 step 9).
GLOBAL_LAYERS_TSV="${PANGENOME_DIR}/misc_data_layers_all_mags.tsv"
log_info "Step 1: Building global layers TSV at $GLOBAL_LAYERS_TSV"

awk -F'\t' '
BEGIN { OFS="\t"; print "item", "DisplayName", "Station", "Depth", "ROMZ_flag", "Order" }
NR > 1 {
    name = $2
    gsub(/[^a-zA-Z0-9_]/, "_", name)
    display = ($3 == "" ? "NA" : $3)
    station = ($4 == "" ? "NA" : $4)
    depth   = ($5 == "" ? "NA" : $5)
    romz    = (depth == "ROMZ" || depth == "OMZ") ? "TRUE" : "FALSE"
    order = "NA"
    n = split($10, parts, ";")
    for (i = 1; i <= n; i++) {
        if (parts[i] ~ /^o__/) {
            order = parts[i]
            sub(/^o__/, "", order)
            if (order == "") order = "NA"
            break
        }
    }
    print name, display, station, depth, romz, order
}' "$MAGS_METADATA" > "$GLOBAL_LAYERS_TSV"

log_info "Step 1 complete: $(($(wc -l < "$GLOBAL_LAYERS_TSV") - 1)) MAG entries written"

# Step 2: For each species, derive a per-species layers TSV restricted to
# the genomes present in its external_genomes.tsv, then import it into
# the corresponding PAN.db.
log_info "Step 2: Importing per-species layers into PAN.db..."
for species_dir in "${PANGENOME_WORK_DIR}"/*/; do
    [[ -d "$species_dir" ]] || continue
    species_name=$(basename "$species_dir")
    pan_db="${species_dir}/pan_output/${species_name}-PAN.db"
    external_genomes_file="${species_dir}/external_genomes.tsv"
    species_layers_tsv="${species_dir}/misc_data_layers_${species_name}.tsv"

    if [[ ! -f "$pan_db" ]]; then
        log_warn "  PAN.db not found for ${species_name//_/ }; skipping"
        continue
    fi
    if [[ ! -f "$external_genomes_file" ]]; then
        log_warn "  external_genomes.tsv not found for ${species_name//_/ }; skipping"
        continue
    fi

    # Subset global layers TSV to the genomes present in external_genomes.tsv
    awk -F'\t' -v global="$GLOBAL_LAYERS_TSV" '
        NR == FNR { if (FNR > 1) keep[$1] = 1; next }
        FNR == 1 { print; next }
        ($1 in keep) { print }
    ' "$external_genomes_file" "$GLOBAL_LAYERS_TSV" > "$species_layers_tsv"

    n_rows=$(($(wc -l < "$species_layers_tsv") - 1))
    if [[ "$n_rows" -le 0 ]]; then
        log_warn "  No layer rows produced for ${species_name//_/ }; skipping import"
        continue
    fi
    log_info "  Importing $n_rows layer rows into PAN.db for: ${species_name//_/ }"

    # Drop pre-existing layers misc-data tables to avoid stale conflicts
    for table in DisplayName Station Depth ROMZ_flag Order; do
        anvi-delete-misc-data \
            -p "$pan_db" \
            --target-data-table layers \
            --keys-to-remove "$table" 2>/dev/null || true
    done

    anvi-import-misc-data \
        "$species_layers_tsv" \
        -p "$pan_db" \
        --target-data-table layers
done
log_info "Step 2 complete"

# --- Done ---
log_info "Pipeline finished successfully"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
log_info "Next: run 25_pangenome_enrichment.sh to compute functional enrichment."
