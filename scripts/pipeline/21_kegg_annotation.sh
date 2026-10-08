#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 21. Runs anvi-run-kegg-kofams over every contigs database that step 20
# produced.
#
# The KEGG annotations are a prerequisite for the functional enrichment in step
# 25, and they are what lets the findings of this thesis be located inside the
# anvi'o pangenome: cbb3 as K00404 to K00407, HndC as K18006, and amoA as
# K10944.

# --- Error handling ---
trap 'log_error "Script failed at line $LINENO with exit code $?"; exit 1' ERR

# --- Header ---
log_info "anvi'o9 KEGG annotation pipeline"
log_info "SLURM job information"
log_info "Job ID: ${SLURM_JOB_ID:-N/A}"
log_info "Job name: ${SLURM_JOB_NAME:-interactive}"
log_info "Node: ${SLURMD_NODENAME:-$(hostname)}"
log_info "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
load_project_config

require_commands anvi-setup-kegg-data anvi-run-kegg-kofams anvi-setup-pfams anvi-run-pfams

# --- Validate required config variables ---
# KEGG_DATA_DIR and PFAM_DATA_DIR are where the anvi'o snapshots live, about
# 30 GB and 3 GB. They are downloaded on the first run and reused afterwards, so
# point them at a filesystem every compute node can read. On nodes without
# network access, download them beforehand with scripts/utils/setup_anvio9.sh.
required_vars=(PROJECT_DIR RESULTS_DIR PANGENOME_SUBDIR TEMP_SUBDIR KEGG_DATA_DIR PFAM_DATA_DIR)
for var in "${required_vars[@]}"; do
    if [[ -z "${!var:-}" ]]; then
        log_error "Required config variable '$var' is not set or is empty"
        exit 1
    fi
done

SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-52}"

# --- Define directory paths ---
PANGENOME_DIR="$RESULTS_DIR/$PANGENOME_SUBDIR"
TEMP_DIR="$PROJECT_DIR/$TEMP_SUBDIR"
ALL_MAGS="$PANGENOME_DIR/all_mags"
ALL_MAGS_CONTIGS_DB="$ALL_MAGS/contigs_db"

# --- Validate input directories ---
if [[ ! -d "$ALL_MAGS_CONTIGS_DB" ]]; then
    log_error "Required input directory does not exist: $ALL_MAGS_CONTIGS_DB"
    log_error "Run 20_pangenome_construction.sh first."
    exit 1
fi

mkdir -p "$TEMP_DIR" "$(dirname "$KEGG_DATA_DIR")" "$(dirname "$PFAM_DATA_DIR")"

#### CONDA ACTIVATION (uncomment if your sbatch session does not pre-activate)
# source "$(conda info --base)/etc/profile.d/conda.sh"
# conda activate anvio-9

# ---- Step 1: Setup KEGG data (one-time download) ----
log_info "Step 1: Setting up KEGG data snapshot..."
KEGG_SETUP_SENTINEL="${KEGG_DATA_DIR}/.kegg_setup_done"
if [[ -f "$KEGG_SETUP_SENTINEL" ]]; then
    log_info "  Skipping anvi-setup-kegg-data (already done at $KEGG_DATA_DIR)"
else
    log_info "  Downloading KEGG snapshot to: $KEGG_DATA_DIR"
    anvi-setup-kegg-data \
        --kegg-data-dir "$KEGG_DATA_DIR" \
        --just-do-it
    touch "$KEGG_SETUP_SENTINEL"
fi
log_info "Step 1 complete"

# ---- Step 2: Run anvi-run-kegg-kofams on every MAG contigs.db ----
log_info "Step 2: Running anvi-run-kegg-kofams on each contigs database..."
kegg_count=0
kegg_skip=0
for db in "${ALL_MAGS_CONTIGS_DB}"/*.db; do
    [[ -f "$db" ]] || continue
    sentinel="${db%.db}.kegg_done"
    if [[ -f "$sentinel" ]]; then
        log_info "  Skipping $(basename "$db") (KEGG already run)"
        ((++kegg_skip))
        continue
    fi
    ((++kegg_count))
    log_info "  Running KOfam on database $kegg_count: $(basename "$db")"
    anvi-run-kegg-kofams \
        -c "$db" \
        --kegg-data-dir "$KEGG_DATA_DIR" \
        --num-threads "$SLURM_CPUS_PER_TASK" \
        --just-do-it
    touch "$sentinel"
done
log_info "Step 2 complete: KOfam run on $kegg_count databases, $kegg_skip skipped"

# Step 3: Setup Pfam data (one-time download)
# Reinforces domain-level annotations alongside KEGG and COG. Pfam-A.hmm is
# the gold-standard protein domain database; useful for archaeal proteins
# that may lack KEGG orthologs but retain conserved Pfam domains.
log_info "Step 3: Setting up Pfam data snapshot..."
PFAM_SETUP_SENTINEL="${PFAM_DATA_DIR}/.pfam_setup_done"
if [[ -f "$PFAM_SETUP_SENTINEL" ]]; then
    log_info "  Skipping anvi-setup-pfams (already done at $PFAM_DATA_DIR)"
else
    # anvi'o refuses --reset on a non-default --pfam-data-dir (safety guard), so
    # a directory left by an interrupted download has to go first. That is left
    # to the user rather than done here, in case PFAM_DATA_DIR points at
    # something else.
    if [[ -e "$PFAM_DATA_DIR" ]]; then
        log_error "  $PFAM_DATA_DIR exists but holds no complete Pfam snapshot."
        log_error "  Remove it and run this step again, or run scripts/utils/setup_anvio9.sh."
        exit 1
    fi
    log_info "  Downloading Pfam snapshot to: $PFAM_DATA_DIR"
    anvi-setup-pfams \
        --pfam-data-dir "$PFAM_DATA_DIR"
    touch "$PFAM_SETUP_SENTINEL"
fi
log_info "Step 3 complete"

# ---- Step 4: Run anvi-run-pfams on every MAG contigs.db ----
log_info "Step 4: Running anvi-run-pfams on each contigs database..."
pfam_count=0
pfam_skip=0
for db in "${ALL_MAGS_CONTIGS_DB}"/*.db; do
    [[ -f "$db" ]] || continue
    sentinel="${db%.db}.pfam_done"
    if [[ -f "$sentinel" ]]; then
        log_info "  Skipping $(basename "$db") (Pfams already run)"
        ((++pfam_skip))
        continue
    fi
    ((++pfam_count))
    log_info "  Running Pfams on database $pfam_count: $(basename "$db")"
    anvi-run-pfams \
        -c "$db" \
        --pfam-data-dir "$PFAM_DATA_DIR" \
        --num-threads "$SLURM_CPUS_PER_TASK"
    touch "$sentinel"
done
log_info "Step 4 complete: Pfams run on $pfam_count databases, $pfam_skip skipped"

# --- Done ---
log_info "Pipeline finished successfully"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
log_info "Next: run 22_pangenome_summarize.sh to refresh genomes storage and summarize."
