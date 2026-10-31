#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# One-off setup of the anvi'o 9 databases that the pangenome steps read. Run it
# once, inside the anvio-9 environment of env/anvio9.yml, on a machine with
# network access:
#
#   conda activate anvio-9
#   bash scripts/utils/setup_anvio9.sh
#
# It installs the NCBI COGs that step 20 annotates with, in anvi'o's default
# location, and the KEGG and Pfam snapshots that step 21 uses, in KEGG_DATA_DIR
# and PFAM_DATA_DIR. Step 21 would download those two itself on its first run,
# but compute nodes often have no outbound network; the marker files written
# here tell it the snapshots are complete, so it goes straight to annotating.
#
# The pangenome steps ran on a second cluster because anvi'o 9 would not build on
# the first, where its datrie dependency failed to compile.

load_project_config
require_vars KEGG_DATA_DIR PFAM_DATA_DIR
require_commands anvi-setup-ncbi-cogs anvi-setup-kegg-data anvi-setup-pfams

THREADS="${THREADS:-4}"

log_info "NCBI COGs, in the anvi'o default location"
anvi-setup-ncbi-cogs --num-threads "$THREADS" --just-do-it

# Same marker files as step 21 checks.
KEGG_SETUP_SENTINEL="${KEGG_DATA_DIR}/.kegg_setup_done"
PFAM_SETUP_SENTINEL="${PFAM_DATA_DIR}/.pfam_setup_done"

if [[ -f "$KEGG_SETUP_SENTINEL" ]]; then
    log_info "KEGG snapshot already in $KEGG_DATA_DIR"
else
    log_info "KEGG snapshot, about 30 GB, into $KEGG_DATA_DIR"
    mkdir -p "$(dirname "$KEGG_DATA_DIR")"
    anvi-setup-kegg-data --kegg-data-dir "$KEGG_DATA_DIR" --just-do-it
    touch "$KEGG_SETUP_SENTINEL"
fi

if [[ -f "$PFAM_SETUP_SENTINEL" ]]; then
    log_info "Pfam snapshot already in $PFAM_DATA_DIR"
else
    # anvi'o refuses --reset outside its default Pfam location, so a directory
    # left by an interrupted download has to go before it can start again. That
    # is left to you rather than done here, in case PFAM_DATA_DIR holds
    # something else.
    if [[ -e "$PFAM_DATA_DIR" ]]; then
        log_error "$PFAM_DATA_DIR exists but holds no complete snapshot."
        log_error "Remove it and run this again."
        exit 1
    fi
    log_info "Pfam snapshot, about 3 GB, into $PFAM_DATA_DIR"
    mkdir -p "$(dirname "$PFAM_DATA_DIR")"
    anvi-setup-pfams --pfam-data-dir "$PFAM_DATA_DIR"
    touch "$PFAM_SETUP_SENTINEL"
fi

log_info "anvi'o databases ready for steps 20 and 21"
