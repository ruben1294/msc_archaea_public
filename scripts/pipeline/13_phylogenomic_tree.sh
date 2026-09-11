#!/usr/bin/env bash
set -euo pipefail

# Step 13. Infers maximum likelihood trees with IQ-TREE 2 over each supermatrix
# from step 12. Three model strategies run against each of the three MAFFT
# branches, so nine trees in total, each under its own directory.
#
# The three strategies:
#
#   partitioned_merge  MFP+MERGE. ModelFinder picks a model per partition and
#                      merges partitions whose models are close enough, which
#                      cuts parameters without forcing one model on everything.
#   concatenated       MFP with no partition file. One model for the whole
#                      supermatrix, the baseline the other two are judged against.
#   partitioned_C60    LG+C60+F+G. A profile mixture model, which is the usual
#                      answer to compositional heterogeneity and the long-branch
#                      attraction it produces. Slowest of the three by far.
#
# On the partition flags, from the IQ-TREE documentation: -q makes all
# partitions share one set of branch lengths, -p gives each partition its own
# evolutionary rate over shared branch lengths, and -Q gives each partition its
# own branch lengths to absorb heterotachy. -p is the recommended default and
# the one used here; -q is unrealistic and -Q is very parameter rich. Running
# all three and comparing BIC is the way to settle it for a given dataset.
#
# Support values come from 1000 ultrafast bootstrap replicates (-B) and 1000
# SH-aLRT replicates (-alrt). Both need to be high for a branch to be trusted.

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

trap 'log_error "Failed at line $LINENO with exit code $?"; exit 1' ERR

log_job_header "IQ-TREE 2 phylogenomic inference"

load_project_config
require_vars PROJECT_DIR

CONCAT_DIR="${PROJECT_DIR}/${MARKER_CONCAT_SUBDIR:-12_marker_concatenation}"
PHYLO_DIR="${PROJECT_DIR}/${PHYLOGENY_SUBDIR:-13_phylogeny}"
CPUS="${SLURM_CPUS_PER_TASK:-24}"

mkdir -p "$PHYLO_DIR"

if ! command -v iqtree2 &> /dev/null; then
    log_error "iqtree2 not found. Load the module or activate its environment."
    exit 1
fi

log_info "IQ-TREE 2: $(iqtree2 --version | head -1)"
log_info "Output:    $PHYLO_DIR"
log_info "CPUs:      $CPUS"

METHODS=(genafpair localpair globalpair)
ANALYSES=(partitioned_merge concatenated partitioned_C60)

# Runs one strategy against one branch. The IQ-TREE arguments arrive as an
# array rather than a string, so paths with spaces survive and nothing goes
# through eval.
run_iqtree() {
    local method="$1" analysis_type="$2"
    shift 2
    local iqtree_opts=("$@")

    local supermatrix="${CONCAT_DIR}/${method}/output/FcC_supermatrix_${method}.fas"
    local output_dir="${PHYLO_DIR}/${method}/${analysis_type}"
    local prefix="${method}_${analysis_type}"

    log_info "$method / $analysis_type"

    if [[ ! -f "$supermatrix" ]]; then
        log_error "  Supermatrix not found: $supermatrix"
        return 1
    fi

    mkdir -p "$output_dir"
    log_info "  Supermatrix: $supermatrix"
    log_info "  Output:      $output_dir"

    if (cd "$output_dir" && iqtree2 "${iqtree_opts[@]}" --prefix "$prefix") \
        > "${output_dir}/${prefix}.iqtree2.log" 2>&1; then

        log_info "  [OK] finished at $(date '+%Y-%m-%d %H:%M:%S')"

        [[ -f "${output_dir}/${prefix}.treefile" ]] && log_info "    Tree:      ${prefix}.treefile"
        [[ -f "${output_dir}/${prefix}.contree" ]]  && log_info "    Consensus: ${prefix}.contree"

        if [[ -f "${output_dir}/${prefix}.iqtree" ]]; then
            local logL
            logL=$(grep "BEST SCORE FOUND" "${output_dir}/${prefix}.iqtree" | awk '{print $5}' || echo "NA")
            log_info "    Log-likelihood: ${logL:-NA}"
        fi
        return 0
    fi

    log_error "  [FAIL] see ${output_dir}/${prefix}.iqtree2.log"
    tail -20 "${output_dir}/${prefix}.iqtree2.log" | while read -r line; do
        log_error "    $line"
    done
    return 1
}

# Returns the IQ-TREE arguments for one strategy against one branch.
iqtree_opts_for() {
    local method="$1" analysis_type="$2"
    local supermatrix="${CONCAT_DIR}/${method}/output/FcC_supermatrix_${method}.fas"
    local partitions="${CONCAT_DIR}/${method}/output/FcC_supermatrix_${method}_partition.txt"

    case "$analysis_type" in
        partitioned_merge)
            printf '%s\n' -s "$supermatrix" -p "$partitions" -m MFP+MERGE ;;
        concatenated)
            printf '%s\n' -s "$supermatrix" -m MFP ;;
        partitioned_C60)
            printf '%s\n' -s "$supermatrix" -p "$partitions" -m LG+C60+F+G ;;
        *)
            log_error "Unknown analysis type: $analysis_type"
            return 1 ;;
    esac
    printf '%s\n' -B 1000 -alrt 1000 -T AUTO --threads-max "$CPUS"
}

declare -A analysis_status

for analysis_type in "${ANALYSES[@]}"; do
    log_info "Strategy: $analysis_type"

    for method in "${METHODS[@]}"; do
        mapfile -t opts < <(iqtree_opts_for "$method" "$analysis_type")

        if run_iqtree "$method" "$analysis_type" "${opts[@]}"; then
            analysis_status["${method}_${analysis_type}"]=0
        else
            analysis_status["${method}_${analysis_type}"]=1
        fi
    done
done

total_analyses=0
successful_analyses=0

log_info "Summary"
for method in "${METHODS[@]}"; do
    log_info "  $method"
    for analysis_type in "${ANALYSES[@]}"; do
        ((++total_analyses))
        if [[ ${analysis_status["${method}_${analysis_type}"]} -eq 0 ]]; then
            log_info "    [OK]   $analysis_type"
            ((++successful_analyses))
        else
            log_error "    [FAIL] $analysis_type"
        fi
    done
done

log_info "  $successful_analyses of $total_analyses analyses succeeded"
log_info "Output layout: ${PHYLO_DIR}/<branch>/<strategy>/"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"

if [[ $successful_analyses -ne $total_analyses ]]; then
    log_error "$((total_analyses - successful_analyses)) analyses failed."
    exit 1
fi
