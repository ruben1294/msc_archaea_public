#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 8. Runs the DRAM2 Nextflow workflow over the filtered MAGs and writes the
# distillate under 08_MAG_distill/: raw-annotations.tsv, traits.xlsx, product.tsv
# and the metabolism summaries. Almost every R script downstream reads one of
# those files, which makes this the pivot of the whole analysis.
#
# This supersedes the DRAM v1 pass in step 5. Both are kept because the two
# versions disagree on enough trait calls that the comparison is worth being
# able to reproduce.
#
# Nextflow itself needs almost nothing to run: it schedules the work and waits.
# The DRAM2 documentation suggests either launching it in the background with
# -bg, which survives an ssh logout and leaves its pid in .nextflow.pid, or
# putting it in a small job of one CPU and a gigabyte or two. The SLURM wrapper
# takes the second route, which is why it asks for so little while running for
# days.
#
# To stop a backgrounded run:
#   kill "$(cat .nextflow.pid)"
#
# DRAM2_REPO must point at a checkout of the DRAM2 pipeline; nextflow runs from
# inside it. NF_CONFIG_FILE selects the profile, and configs/ holds the two used
# here.

load_project_config
require_vars PROJECT_DIR FILTERED_SUBDIR QUALITY_SUBDIR TAXONOMY_SUBDIR
require_commands nextflow

log_job_header "DRAM2 MAG distillation"

OUT_DIR="${PROJECT_DIR}/${MAG_DISTILL_SUBDIR:-08_MAG_distill}"
RAW_DIR="${PROJECT_DIR}/${FILTERED_SUBDIR}"
QUALITY_DIR="${PROJECT_DIR}/${QUALITY_SUBDIR}"
TAXONOMY_DIR="${PROJECT_DIR}/${TAXONOMY_SUBDIR}"
export OUT_DIR RAW_DIR QUALITY_DIR TAXONOMY_DIR

NF_CONFIG_FILE="${NF_CONFIG_FILE:-${REPO_ROOT}/configs/20260113_CMRA_nf_DRAM2_HPC.config}"
DRAM2_REPO="${DRAM2_REPO:-${HOME}/.nextflow/assets/WrightonLabCSU/DRAM}"

if [[ ! -f "$NF_CONFIG_FILE" ]]; then
    log_error "Nextflow config not found: $NF_CONFIG_FILE"
    exit 1
fi

if [[ ! -d "$DRAM2_REPO" ]]; then
    log_error "DRAM2 checkout not found: $DRAM2_REPO"
    log_error "Set DRAM2_REPO, or fetch it with: nextflow pull WrightonLabCSU/DRAM"
    exit 1
fi

mkdir -p "$OUT_DIR"

log_info "Nextflow config: $NF_CONFIG_FILE"
log_info "DRAM2 checkout:  $DRAM2_REPO"
log_info "Output:          $OUT_DIR"
log_info "Databases:       ${DRAM2_DB_DIR:-set in the Nextflow config}"

log_info "Nextflow version:"
nextflow -version || log_warn "could not read the Nextflow version"

# The pipeline runs its tools in containers, so a broken container runtime shows
# up here rather than several hours in.
if command -v docker > /dev/null 2>&1; then
    docker info > /dev/null 2>&1 || log_warn "docker is on PATH but not responding"
    docker images | grep -q dram || log_warn "no DRAM docker image found locally"
else
    log_warn "docker not found; the config must select another container engine"
fi

# Offline, because the compute nodes have no outbound network and Nextflow would
# otherwise stall trying to check for updates.
export NXF_OFFLINE='true'

# -resume reuses the cached work directory, which matters because a full run
# takes days and a failure partway through is expensive to repeat.
exit_code=0
(cd "$DRAM2_REPO" && nextflow run . \
    --outdir "$OUT_DIR" \
    -c "$NF_CONFIG_FILE" \
    -profile docker \
    -resume) || exit_code=$?

if [[ $exit_code -eq 0 ]]; then
    log_info "Distillate written to $OUT_DIR"
else
    log_error "DRAM2 workflow failed with exit code $exit_code"
fi

log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
exit $exit_code
