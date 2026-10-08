#!/bin/bash
#SBATCH --job-name=DRAM2_install
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.log

# One-off setup for step 10. Clones the DRAM2 fork at the commit that
# DRAM2_REVISION pins in configs/server.env, caches the Nextflow plugins the
# pipeline declares, and checks the container engine. It needs network access,
# so run it wherever there is some, from the repository root:
#
#   bash slurm/DRAM2_HPC_install.slurm
#
# or, on a cluster whose login nodes should not run it, as a job on a node that
# reaches the internet:
#
#   mkdir -p logs
#   sbatch --account=<account> --partition=<partition> \
#       --nodelist=<node with network access> slurm/DRAM2_HPC_install.slurm
#
# Compute nodes often have no outbound network, which is why everything that
# needs it happens here, once. Step 10 then runs with NXF_OFFLINE=true and finds
# the pipeline in DRAM2_REPO and the plugins in ~/.nextflow/plugins.
#
# The pipeline is cloned with git rather than fetched with nextflow pull because
# recent Nextflow releases keep pulled revisions under ~/.nextflow/assets/.repos,
# one directory per commit, and step 10 needs one fixed path to run from.

set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"

if [[ ! -f "$REPO_ROOT/scripts/lib/common.sh" ]]; then
    echo "Run this from the repository root, not from slurm/." >&2
    exit 1
fi
mkdir -p "$REPO_ROOT/logs"

# Cluster-specific module and conda setup, as in the other wrappers.
if [[ -f "$REPO_ROOT/slurm/cluster.env" ]]; then
    # shellcheck source=/dev/null
    source "$REPO_ROOT/slurm/cluster.env"
fi

if [[ -n "${CLUSTER_MODULE:-}" ]]; then
    module load "$CLUSTER_MODULE"
fi
if [[ -n "${CONDA_PROFILE:-}" && -f "$CONDA_PROFILE" ]]; then
    # shellcheck source=/dev/null
    source "$CONDA_PROFILE"
fi
conda activate "${DRAM2_CONDA_ENV:-DRAMV2}"

# shellcheck source=../scripts/lib/common.sh
source "$REPO_ROOT/scripts/lib/common.sh"
load_project_config
require_vars DRAM2_PIPELINE DRAM2_REVISION DRAM2_REPO
require_commands nextflow git

log_job_header "DRAM2 install"

url="https://github.com/${DRAM2_PIPELINE}.git"

if [[ -d "$DRAM2_REPO/.git" ]]; then
    log_info "Updating the existing checkout in $DRAM2_REPO"
    git -C "$DRAM2_REPO" fetch --quiet origin
elif [[ -e "$DRAM2_REPO" ]]; then
    log_error "$DRAM2_REPO exists but is not a git checkout. Move it aside first."
    exit 1
else
    log_info "Cloning $url into $DRAM2_REPO"
    mkdir -p "$(dirname "$DRAM2_REPO")"
    git clone --quiet "$url" "$DRAM2_REPO"
fi

if ! git -C "$DRAM2_REPO" cat-file -e "${DRAM2_REVISION}^{commit}" 2> /dev/null; then
    log_error "Commit $DRAM2_REVISION is not in $DRAM2_REPO."
    log_error "Its origin should be $url; check with: git -C $DRAM2_REPO remote -v"
    exit 1
fi

# A detached HEAD at the pinned commit, so that later pushes to the fork cannot
# change what step 10 runs.
git -C "$DRAM2_REPO" -c advice.detachedHead=false checkout --quiet "$DRAM2_REVISION"
log_info "DRAM2 checkout at $(git -C "$DRAM2_REPO" rev-parse HEAD)"

# Installing exactly the plugins the pipeline declares is what lets the offline
# runs of step 10 start without reaching the plugin registry.
mapfile -t plugins < <(grep -oE "id +['\"][^'\"]+@[^'\"]+['\"]" "$DRAM2_REPO/nextflow.config" \
    | grep -oE "[^'\"]+@[^'\"]+")

if [[ ${#plugins[@]} -gt 0 ]]; then
    plugin_list=$(IFS=,; echo "${plugins[*]}")
    log_info "Caching Nextflow plugins: $plugin_list"
    NXF_OFFLINE=false nextflow plugin install "$plugin_list"
else
    log_warn "No plugins declared in $DRAM2_REPO/nextflow.config"
fi

log_info "Nextflow version:"
nextflow -version

# The container images are pulled the first time each DRAM2 task runs, on
# whichever node runs it. If those nodes have no network either, set
# DRAM2_NODELIST in server.env to one that does, as the thesis runs did.
DRAM2_PROFILE="${DRAM2_PROFILE:-docker}"
case "$DRAM2_PROFILE" in
    docker)
        if command -v docker > /dev/null 2>&1 && docker info > /dev/null 2>&1; then
            log_info "docker responds on $(hostname)"
        else
            log_warn "docker is missing or not responding on $(hostname)"
        fi
        ;;
    singularity|apptainer|podman)
        if command -v "$DRAM2_PROFILE" > /dev/null 2>&1; then
            log_info "$DRAM2_PROFILE found on $(hostname)"
        else
            log_warn "$DRAM2_PROFILE not found on $(hostname)"
        fi
        ;;
    *)
        log_warn "Unrecognised DRAM2_PROFILE: $DRAM2_PROFILE"
        ;;
esac

log_info "Done. Step 10 can now run offline."
