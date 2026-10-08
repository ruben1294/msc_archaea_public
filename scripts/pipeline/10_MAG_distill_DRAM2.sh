#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 10. Runs the DRAM2 Nextflow workflow over the filtered MAGs and writes the
# distillate under results/10_MAG_distill/: raw-annotations.tsv, traits.xlsx,
# product.tsv and the metabolism summaries. Almost every R script downstream
# reads one of those files, which makes this the pivot of the whole analysis.
# It needs the GTDB-Tk summary of step 03 and the CheckM2 report of step 01,
# both of which DRAM2 carries into the distillate.
#
# This supersedes the DRAM v1 pass in step 07. Both are kept because the two
# versions disagree on enough trait calls that the comparison is worth being
# able to reproduce.
#
# The thesis ran a fork of DRAM2, ruben1294/DRAM, under Nextflow 25.10.2 with
# the docker profile. The fork carries fixes to the KEGG database formatting,
# the archaea-calibrated trait rules and the run configuration of this step;
# its README lists every change. All of that is DRAM2-derived and GPL-licensed,
# which is why it lives in the fork and not here. DRAM2_REVISION in server.env
# pins the commit; slurm/DRAM2_HPC_install.slurm fetches it into DRAM2_REPO,
# and this script refuses a checkout at any other commit unless
# DRAM2_ALLOW_OTHER_REVISION=true.
#
# Nextflow itself needs almost nothing to run: it schedules the work and waits.
# Either launch it in the background with -bg, which leaves its pid in
# .nextflow.pid, or put it in a small job of one CPU and a gigabyte or two. The
# SLURM wrapper takes the second route, which is why it asks for so little while
# running for days.
#
# To stop a backgrounded run:
#   kill "$(cat .nextflow.pid)"
#
# NF_CONFIG_FILE selects the configuration. The fork's conf/ holds the full
# run, gom_archaea_mags.config, and a shorter test run,
# gom_archaea_mags_test.config. Both read REPO_ROOT, RESULTS_DIR, DRAM2_DB_DIR
# and the DRAM2_* executor settings from the environment, which
# load_project_config exports.

load_project_config
require_vars PROJECT_DIR RESULTS_DIR FILTERED_SUBDIR QUALITY_SUBDIR TAXONOMY_SUBDIR \
    DRAM2_DB_DIR DRAM2_REPO
require_commands nextflow git

log_job_header "DRAM2 MAG distillation"

OUT_DIR="${RESULTS_DIR}/${MAG_DISTILL_SUBDIR:-10_MAG_distill}"
NF_CONFIG_FILE="${NF_CONFIG_FILE:-${DRAM2_REPO}/conf/gom_archaea_mags.config}"
DRAM2_PROFILE="${DRAM2_PROFILE:-docker}"

if [[ ! -f "$NF_CONFIG_FILE" ]]; then
    log_error "Nextflow config not found: $NF_CONFIG_FILE"
    exit 1
fi

if [[ ! -d "$DRAM2_REPO" ]]; then
    log_error "DRAM2 checkout not found: $DRAM2_REPO"
    log_error "Clone it with slurm/DRAM2_HPC_install.slurm, on a machine with network access."
    exit 1
fi

# A checkout at another commit produces a different distillate, so a mismatch
# stops the run unless it is asked for.
if [[ -n "${DRAM2_REVISION:-}" ]]; then
    actual_revision=$(git -C "$DRAM2_REPO" rev-parse HEAD 2>/dev/null || echo "unknown")
    if [[ "$actual_revision" != "$DRAM2_REVISION" ]]; then
        if [[ "${DRAM2_ALLOW_OTHER_REVISION:-false}" == true ]]; then
            log_warn "DRAM2 checkout is at $actual_revision, not the pinned $DRAM2_REVISION"
        else
            log_error "DRAM2 checkout is at $actual_revision, not the pinned $DRAM2_REVISION"
            log_error "Check it out with: git -C $DRAM2_REPO checkout $DRAM2_REVISION"
            log_error "or set DRAM2_ALLOW_OTHER_REVISION=true to run anyway."
            exit 1
        fi
    fi
fi

mkdir -p "$OUT_DIR"

log_info "Nextflow config: $NF_CONFIG_FILE"
log_info "DRAM2 checkout:  $DRAM2_REPO (${DRAM2_REVISION:-revision not pinned})"
log_info "Profile:         $DRAM2_PROFILE"
log_info "Output:          $OUT_DIR"
log_info "Databases:       $DRAM2_DB_DIR"

log_info "Nextflow version:"
nextflow -version || log_warn "could not read the Nextflow version"

# The pipeline runs its tools in containers. With the SLURM executor they start
# on the compute nodes rather than here, so a missing engine on this node is a
# warning, not an error.
case "$DRAM2_PROFILE" in
    docker)
        if command -v docker > /dev/null 2>&1; then
            docker info > /dev/null 2>&1 || log_warn "docker is on PATH but not responding"
        else
            log_warn "docker not found on this node; the nodes that run DRAM2 tasks need it"
        fi
        ;;
    singularity|apptainer|podman)
        command -v "$DRAM2_PROFILE" > /dev/null 2>&1 \
            || log_warn "$DRAM2_PROFILE not found on this node; the nodes that run DRAM2 tasks need it"
        ;;
esac

# Offline by default, because compute nodes often have no outbound network and
# Nextflow would otherwise stall checking for updates. Set NXF_OFFLINE=false
# where the network is available.
export NXF_OFFLINE="${NXF_OFFLINE:-true}"

# Every task leaves its intermediate files in the work directory, which grows
# large, so it sits under PROJECT_DIR, which can be scratch, rather than inside
# the pipeline checkout.
export NXF_WORK="${DRAM2_WORK_DIR:-${PROJECT_DIR}/work/DRAM2}"
log_info "Work directory:  $NXF_WORK"

# -resume reuses the cached work directory, which matters because a full run
# takes days and a failure partway through is expensive to repeat.
exit_code=0
(cd "$DRAM2_REPO" && nextflow run . \
    --outdir "$OUT_DIR" \
    -c "$NF_CONFIG_FILE" \
    -profile "$DRAM2_PROFILE" \
    -resume) || exit_code=$?

if [[ $exit_code -eq 0 ]]; then
    log_info "Distillate written to $OUT_DIR"
else
    log_error "DRAM2 workflow failed with exit code $exit_code"
fi

log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
exit $exit_code
