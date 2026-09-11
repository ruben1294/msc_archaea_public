#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 15. Opens each cleaned pangenome in the anvi'o interactive interface, one
# species at a time, so the figures can be laid out by hand.
#
# This is not a batch renderer. anvi-display-pan starts a server and blocks
# until you stop it, so each invocation below waits for you to close the browser
# before the next one starts. It runs on the local machine, not the cluster.
#
# There is no headless SVG export here, and that is deliberate. The layer
# ordering, colouring and gene cluster selection that the thesis figures use
# were chosen in the interface; a scripted render would not reproduce those
# choices, only a default view of the same data.
#
# It reads the *_modified.db copies that step 26 produces, not the originals.

trap 'log_error "Failed at line $LINENO with exit code $?"; exit 1' ERR

load_project_config
require_vars LOCAL_RESULTS_DIR

# Where step 26 left the publication copies on this machine.
LOCAL_PANGENOME_DIR="${LOCAL_PANGENOME_DIR:-${LOCAL_RESULTS_DIR}/14_pangenome/pangenome}"

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
        log_warn "  Run 26_pangenome_clean_for_publication.sh first."
        continue
    fi

    log_info "Opening $species"
    anvi-display-pan -p "$pan_db" -g "$genomes_db" -I localhost
done

log_info "All pangenomes visited."
