#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 22. Refreshes the genomes storages so the KEGG annotations from step 21
# reach the pangenome, then runs anvi-summarize per species.
#
# The refresh is not optional: annotations added to a contigs database after the
# storage was built stay invisible to the pangenome until the storage is rebuilt.
#
# What this exports, gene_clusters_summary, is the matrix every downstream R
# analysis reads.

# --- Error handling ---
trap 'log_error "Script failed at line $LINENO with exit code $?"; exit 1' ERR

# --- Header ---
log_info "anvi'o9 pangenome summarize pipeline"
log_info "SLURM job information"
log_info "Job ID: ${SLURM_JOB_ID:-N/A}"
log_info "Job name: ${SLURM_JOB_NAME:-interactive}"
log_info "Node: ${SLURMD_NODENAME:-$(hostname)}"
log_info "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
load_project_config

require_commands anvi-gen-genomes-storage anvi-pan-genome anvi-script-add-default-collection anvi-summarize

# --- Validate required config variables ---
required_vars=(PROJECT_DIR RESULTS_DIR PANGENOME_SUBDIR)
for var in "${required_vars[@]}"; do
    if [[ -z "${!var:-}" ]]; then
        log_error "Required config variable '$var' is not set or is empty"
        exit 1
    fi
done

SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-52}"

# --- Define directory paths ---
PANGENOME_DIR="$RESULTS_DIR/$PANGENOME_SUBDIR"
ALL_MAGS_CONTIGS_DB="$PANGENOME_DIR/all_mags/contigs_db"
PANGENOME_WORK_DIR="$PANGENOME_DIR/pangenome"

# --- Validate inputs ---
if [[ ! -d "$PANGENOME_WORK_DIR" ]]; then
    log_error "Pangenome working directory does not exist: $PANGENOME_WORK_DIR"
    log_error "Run 20_pangenome_construction.sh first."
    exit 1
fi
if [[ ! -d "$ALL_MAGS_CONTIGS_DB" ]]; then
    log_error "All-MAGs contigs.db directory does not exist: $ALL_MAGS_CONTIGS_DB"
    exit 1
fi

#### CONDA ACTIVATION (uncomment if your sbatch session does not pre-activate)
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate anvio-9

# Step 1: Refresh species-level contigs.db copies from the KEGG-annotated
# originals in all_mags/contigs_db. anvi-pan-genome already linked them by
# path; what we need is for anvi-summarize to read KEGG functions from the
# refreshed contigs DBs through a regenerated genomes storage.
log_info "Step 1: Refreshing per-species contigs.db copies from KEGG-annotated originals..."
refresh_count=0
for species_dir in "${PANGENOME_WORK_DIR}"/*/; do
    [[ -d "$species_dir" ]] || continue
    species_name=$(basename "$species_dir")
    species_contigs_db_dir="${species_dir}/contigs_db"
    [[ -d "$species_contigs_db_dir" ]] || continue

    log_info "  Refreshing contigs.db for species: ${species_name//_/ }"
    for db in "${species_contigs_db_dir}"/*.db; do
        [[ -f "$db" ]] || continue
        mag_name=$(basename "$db")
        src="${ALL_MAGS_CONTIGS_DB}/${mag_name}"
        if [[ -f "$src" ]]; then
            cp -f "$src" "$db"
            ((++refresh_count))
        else
            log_warn "  Source contigs.db not found for $mag_name in $ALL_MAGS_CONTIGS_DB"
        fi
    done
done
log_info "Step 1 complete: $refresh_count contigs.db copies refreshed"

# Step 2: Regenerate genomes storage per species so KEGG annotations are
# accessible to anvi-summarize and to anvi-compute-functional-enrichment.
log_info "Step 2: Regenerating genomes storage per species..."
for species_dir in "${PANGENOME_WORK_DIR}"/*/; do
    [[ -d "$species_dir" ]] || continue
    species_name=$(basename "$species_dir")
    external_genomes_file="${species_dir}/external_genomes.tsv"
    genomes_db_out="${species_dir}/${species_name}-GENOMES.db"

    if [[ ! -f "$external_genomes_file" ]]; then
        log_warn "  external_genomes.tsv not found for ${species_name//_/ }; skipping"
        continue
    fi

    # Detect the actual gene caller used in the first contigs.db of this
    # species. anvi'o 9 default is pyrodigal_gv, NOT literal "prodigal".
    # Hardcoding the wrong name causes "None of your genomes seem to have
    # a gene call". Use anvi-export-gene-calls --list-gene-callers to read
    # the actual caller name from the DB.
    first_db=$(awk -F'\t' 'NR==2 {print $2; exit}' "$external_genomes_file")
    if [[ -z "$first_db" || ! -f "$first_db" ]]; then
        log_warn "  Could not locate a contigs.db to detect the gene caller; skipping ${species_name//_/ }"
        continue
    fi
    # anvi-export-gene-calls --list-gene-callers prints lines like
    #   * 'pyrodigal-gv' is a valid gene caller ...
    # We need the bare name (without surrounding quotes), e.g. pyrodigal-gv.
    raw=$(anvi-export-gene-calls -c "$first_db" --list-gene-callers 2>&1 \
        | grep -Ei "['\"](pyrodigal[-_]?gv|pyrodigal|prodigal|prokka|gtdb)['\"]" \
        | head -n 1)
    gene_caller=$(printf '%s\n' "$raw" \
        | sed -E "s/.*['\"]([^'\"]+)['\"].*/\1/")
    if [[ -z "$gene_caller" ]]; then
        log_warn "  Could not auto-detect gene caller for ${species_name//_/ }; relying on anvi'o default"
        anvi-gen-genomes-storage \
            -e "$external_genomes_file" \
            -o "$genomes_db_out"
    else
        log_info "  Regenerating genomes storage for: ${species_name//_/ } (gene-caller: $gene_caller)"
        anvi-gen-genomes-storage \
            -e "$external_genomes_file" \
            -o "$genomes_db_out" \
            --gene-caller "$gene_caller"
    fi
done
log_info "Step 2 complete"

# ---- Step 3: Add a default collection to each PAN.db (required by anvi-summarize) ----
log_info "Step 3: Ensuring default collection exists on each PAN.db..."
for species_dir in "${PANGENOME_WORK_DIR}"/*/; do
    [[ -d "$species_dir" ]] || continue
    species_name=$(basename "$species_dir")
    pan_db="${species_dir}/pan_output/${species_name}-PAN.db"
    genomes_db="${species_dir}/${species_name}-GENOMES.db"

    if [[ ! -f "$pan_db" ]]; then
        log_warn "  PAN.db not found for ${species_name//_/ }; skipping"
        continue
    fi
    if [[ ! -f "$genomes_db" ]]; then
        log_warn "  GENOMES.db not found for ${species_name//_/ }; skipping"
        continue
    fi

    sentinel="${species_dir}/.default_collection_done"
    if [[ -f "$sentinel" ]]; then
        log_info "  Skipping default collection for ${species_name//_/ } (already added)"
        continue
    fi
    log_info "  Adding default collection for: ${species_name//_/ }"
    anvi-script-add-default-collection \
        -p "$pan_db" \
        --collection-name default
    touch "$sentinel"
done
log_info "Step 3 complete"

# ---- Step 4: anvi-summarize per species ----
log_info "Step 4: Running anvi-summarize per species..."
for species_dir in "${PANGENOME_WORK_DIR}"/*/; do
    [[ -d "$species_dir" ]] || continue
    species_name=$(basename "$species_dir")
    pan_db="${species_dir}/pan_output/${species_name}-PAN.db"
    genomes_db="${species_dir}/${species_name}-GENOMES.db"
    summary_dir="${species_dir}/SUMMARY"

    if [[ ! -f "$pan_db" || ! -f "$genomes_db" ]]; then
        log_warn "  Inputs missing for ${species_name//_/ }; skipping"
        continue
    fi

    if [[ -d "$summary_dir" ]]; then
        log_info "  Removing previous SUMMARY for ${species_name//_/ }"
        rm -rf "$summary_dir"
    fi

    log_info "  Summarizing pangenome for: ${species_name//_/ }"
    anvi-summarize \
        -p "$pan_db" \
        -g "$genomes_db" \
        -C default \
        -o "$summary_dir"
done
log_info "Step 4 complete"

# --- Done ---
log_info "Pipeline finished successfully"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
log_info "Next: run 23_pangenome_misc_data_layers.sh to import metadata layers."
