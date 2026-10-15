#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 26. Opens each cleaned pangenome in the anvi'o interactive interface, one
# species at a time, so the figures can be laid out by hand.
#
# It comes last because it reads the *_modified.db copies that step 25 writes,
# not the originals.
#
# anvi-display-pan starts a server and blocks
# until you stop it, so each invocation below waits for you to close the browser
# before the next one starts. It runs on the local machine, not the cluster.
#
# There is no headless SVG export here. The layer ordering, colouring and gene
# cluster selection of the thesis figures were set by hand in the interface,
# and a scripted render would only produce the default view of the same data.

trap 'log_error "Failed at line $LINENO with exit code $?"; exit 1' ERR

# A laptop may have no configs/server.env, so the clone itself is the fallback
# root, the same one the R scripts use. MSC_ARCHAEA_DIR overrides it.
if [[ -f "${CONFIG_FILE:-${REPO_ROOT}/configs/server.env}" ]]; then
    load_project_config
fi
RESULTS_DIR="${RESULTS_DIR:-${MSC_ARCHAEA_DIR:-$REPO_ROOT}/results}"

# Where step 25 left the publication copies on this machine.
LOCAL_PANGENOME_DIR="${LOCAL_PANGENOME_DIR:-${RESULTS_DIR}/${PANGENOME_SUBDIR:-20_pangenome}/pangenome}"

# The five species with enough MAGs to build a pangenome from.
SPECIES=(
    Thalassarchaeum_sp002507125
    MGIIb_O1_sp002498525
    UBA8886_sp003193815
    UBA57_sp002494565
    UBA60_sp002503395
)

if ! command -v anvi-display-pan > /dev/null 2>&1; then
    log_error "anvi-display-pan not found. Activate the anvio-9 environment first:"
    log_error "  conda activate anvio-9"
    exit 1
fi

log_info "Pangenome directory: $LOCAL_PANGENOME_DIR"
log_info "Each pangenome opens in turn. Stop the server to move to the next one."

for species in "${SPECIES[@]}"; do
    pan_db="${LOCAL_PANGENOME_DIR}/${species}/pan_output/${species}-PAN_modified.db"
    genomes_db="${LOCAL_PANGENOME_DIR}/${species}/${species}-GENOMES_modified.db"

    if [[ ! -f "$pan_db" ]]; then
        log_warn "Skipping ${species}: no $pan_db"
        log_warn "  Run 25_pangenome_clean_for_publication.sh first."
        continue
    fi

    log_info "Opening $species"
    anvi-display-pan -p "$pan_db" -g "$genomes_db" -I localhost
done

log_info "All pangenomes visited."
