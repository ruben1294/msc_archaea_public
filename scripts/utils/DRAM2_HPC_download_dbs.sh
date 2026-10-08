#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Requests the DRAM2 databases that step 10 uses from the public DRAM2 Globus
# collection into DRAM2_DB_DIR on your own Globus endpoint. Needs the Globus CLI,
# logged in (globus login), and the UUID of the destination collection:
#
#   export GLOBUS_DEST_UUID=<your collection UUID>
#   bash scripts/utils/DRAM2_HPC_download_dbs.sh
#
# globus transfer only submits a task, which then runs on Globus's side, so a
# successful run here means every transfer was accepted, not that it finished.
# Follow them with `globus task list` and `globus task show <TASK_ID>`.
#
# KEGG is not among them: it is licensed, and the thesis used a licensed copy
# (see "KEGG licensing" in README.md). Without one, set use_kegg to false in the
# DRAM2 configuration.
#
# DRAM2_SOURCE_UUID overrides the source collection, in case the DRAM2 team
# moves it.

load_project_config
require_vars DRAM2_DB_DIR GLOBUS_DEST_UUID
require_commands globus

DRAM2_SOURCE_UUID="${DRAM2_SOURCE_UUID:-97ed64b9-dea0-4eb5-a7e0-0b50ac94e889}"
DATABASES=(db_descriptions camper dbcan fegenie kofam metals methyl pfam sulfur uniref)

LOG_DIR="${PROJECT_DIR}/${LOGS_SUBDIR:-logs}"
mkdir -p "$LOG_DIR"
LOG_FILE="${LOG_DIR}/DRAM2_download_$(date +%Y%m%d_%H%M%S).log"

log_info "DRAM2 databases: ${DRAM2_SOURCE_UUID} -> ${GLOBUS_DEST_UUID}:${DRAM2_DB_DIR}" | tee -a "$LOG_FILE"

SUBMITTED=()
FAILED=()

for db in "${DATABASES[@]}"; do
    log_info "Submitting transfer: $db" | tee -a "$LOG_FILE"
    # pipefail makes the if test globus's exit status, not tee's.
    if globus transfer \
        "${DRAM2_SOURCE_UUID}:/databases/${db}/" \
        "${GLOBUS_DEST_UUID}:${DRAM2_DB_DIR}/${db}/" \
        --recursive \
        --label "DRAM2_${db}_$(date +%Y%m%d)" 2>&1 | tee -a "$LOG_FILE"; then
        SUBMITTED+=("$db")
    else
        log_error "Transfer of $db was not accepted" 2>&1 | tee -a "$LOG_FILE"
        FAILED+=("$db")
    fi
done

log_info "Submitted ${#SUBMITTED[@]} of ${#DATABASES[@]}: ${SUBMITTED[*]:-none}" | tee -a "$LOG_FILE"
if [[ ${#FAILED[@]} -gt 0 ]]; then
    log_error "Not accepted: ${FAILED[*]}" 2>&1 | tee -a "$LOG_FILE"
fi
log_info "Follow the transfers with: globus task list" | tee -a "$LOG_FILE"
log_info "Log: $LOG_FILE"

[[ ${#FAILED[@]} -eq 0 ]]
