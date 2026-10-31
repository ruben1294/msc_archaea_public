# Helpers shared by the pipeline and utility scripts. Source it near the top of
# a script, before anything else runs:
#
#   source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"
#   load_project_config
#
# This file is meant to be sourced, so it has no shebang and sets no options of
# its own; the calling script owns set -euo pipefail.

# shellcheck shell=bash

_log_stamp() { date '+%Y-%m-%d %H:%M:%S'; }

log_info()  { echo "[$(_log_stamp)] [INFO] $*"; }
log_warn()  { echo "[$(_log_stamp)] [WARN] $*"; }
log_debug() { echo "[$(_log_stamp)] [DEBUG] $*"; }
log_error() { echo "[$(_log_stamp)] [ERROR] $*" >&2; }

# Repository root, derived from this file's own location rather than from the
# caller, so it holds no matter which directory the script runs from or how deep
# in scripts/ the caller is.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export REPO_ROOT

# Sources configs/server.env and exports what it defines, so that tools started
# from a step (Nextflow reading the DRAM2 config, the Python helpers) see the
# same paths. Override CONFIG_FILE to use another profile, for example one per
# machine when the work is split across two clusters.
load_project_config() {
    CONFIG_FILE="${CONFIG_FILE:-${REPO_ROOT}/configs/server.env}"

    if [[ ! -f "$CONFIG_FILE" ]]; then
        log_error "No configuration at $CONFIG_FILE"
        log_error "Copy configs/server.env.example to configs/server.env and edit it."
        return 1
    fi

    set -a
    # shellcheck source=/dev/null
    source "$CONFIG_FILE"
    set +a
    log_info "Loaded configuration from $CONFIG_FILE"
}

# Fails with a list of everything missing rather than one variable at a time, so
# every gap in a half-filled server.env shows up in a single run.
require_vars() {
    local missing=()
    local var
    for var in "$@"; do
        if [[ -z "${!var:-}" ]]; then
            missing+=("$var")
        fi
    done

    if [[ ${#missing[@]} -gt 0 ]]; then
        log_error "Unset or empty in ${CONFIG_FILE:-the configuration}: ${missing[*]}"
        return 1
    fi
}

# Checks that every external tool a script needs is reachable, and names all the
# missing ones at once. Without this a script gets partway through and stops on
# "command not found", after it has already written output.
require_commands() {
    local missing=()
    local cmd
    for cmd in "$@"; do
        command -v "$cmd" > /dev/null 2>&1 || missing+=("$cmd")
    done

    if [[ ${#missing[@]} -gt 0 ]]; then
        log_error "Not found on PATH: ${missing[*]}"
        log_error "Load the module or activate the conda environment that provides them."
        return 1
    fi
}

# Job context for the SLURM log. The defaults keep the script usable outside a
# scheduler, where none of these variables exist.
log_job_header() {
    log_info "$1"
    log_info "  Host:       ${SLURMD_NODENAME:-$(hostname)}"
    log_info "  Job ID:     ${SLURM_JOB_ID:-NA}"
    log_info "  Job name:   ${SLURM_JOB_NAME:-interactive}"
    log_info "  Start time: $(_log_stamp)"
}
