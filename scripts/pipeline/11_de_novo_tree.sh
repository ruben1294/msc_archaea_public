#!/usr/bin/env bash
set -euo pipefail

# Step 11. Builds a de novo archaeal tree with gtdbtk de_novo_wf, placing the 54
# MAGs among the GTDB reference genomes. The result is the 7,022-leaf Newick
# that scripts/pipeline/curate_tree_and_itol.py later prunes for publication.
#
# The thesis ran this with GTDB-Tk 2.5.2 against GTDB R226, the same pairing as
# step 03.
#
# The outgroup is p__Euryarchaeota. Rooting on p__Thermoplasmatota, the
# alternative left in the code below, is not usable here because most of the
# ingroup sits inside that phylum.
#
# GTDB-Tk runs Prodigal itself and can die with "OSError: AF_UNIX path too long"
# when the temporary directory path is deep. Passing --genes with the proteins
# from step 05 avoids that, but GTDB-Tk then skips the ANI screen and
# classification steps, and this version rejected the flag
# ("unrecognized arguments"). The workaround in use is a short TMPDIR.

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

trap 'log_error "Failed at line $LINENO with exit code $?"; exit 1' ERR

log_job_header "GTDB-Tk de novo tree"

load_project_config
require_commands gtdbtk

# GTDB-Tk reads its reference data from GTDBTK_DATA_PATH.
if [[ -n "${GTDBTK_DB:-}" ]]; then
    export GTDBTK_DATA_PATH="$GTDBTK_DB"
fi
require_vars PROJECT_DIR RESULTS_DIR FILTERED_SUBDIR TAXONOMY_SUBDIR

RAW_DIR="${RESULTS_DIR}/${FILTERED_SUBDIR}"
OUT_DIR="${RESULTS_DIR}/${GTDBTK_DENOVO_TREE_SUBDIR:-11_de_novo_tree}"
TAX_FILE="${RESULTS_DIR}/${TAXONOMY_SUBDIR}/classify/tax_.ar53.summary.tsv"

PROT_MODEL="LG"                 # JTT, WAG or LG
PREFIX="54_archaea"
OUTGROUP_TAX="p__Euryarchaeota"
CPUS="${SLURM_CPUS_PER_TASK:-24}"

# Keep this path short. A deep TMPDIR is what triggers the AF_UNIX error above.
# mktemp rather than reusing an inherited TMPDIR: the cleanup below must only
# ever remove a directory this script created. SLURM wrappers export a TMPDIR
# inside the project tree.
TMPDIR="$(mktemp -d /tmp/gtdbtk_XXXXXXXX)"
export TMPDIR
trap 'rm -rf "${TMPDIR:?}"' EXIT

log_info "Input MAGs:  $RAW_DIR"
log_info "Output:      $OUT_DIR"
log_info "Taxonomy:    $TAX_FILE"
log_info "Outgroup:    $OUTGROUP_TAX"
log_info "Model:       $PROT_MODEL"
log_info "CPUs:        $CPUS"
log_info "TMPDIR:      $TMPDIR"

if [[ ! -d "$RAW_DIR" ]]; then
    log_error "Input directory not found: $RAW_DIR"
    exit 1
fi

if [[ ! -f "$TAX_FILE" ]]; then
    log_error "Classification file not found: $TAX_FILE"
    log_error "Run 03_taxonomic_assignment.sh first."
    exit 1
fi

mkdir -p "$OUT_DIR"

# Capture the exit code instead of letting set -e abort, so the closing message
# can tell success from failure.
exit_code=0
gtdbtk de_novo_wf \
    --genome_dir "$RAW_DIR" \
    --outgroup_taxon "$OUTGROUP_TAX" \
    --extension "fa" \
    --archaea \
    --out_dir "$OUT_DIR" \
    --prot_model "$PROT_MODEL" \
    --write_single_copy_genes \
    --prefix "$PREFIX" \
    --rnd_seed 42 \
    --gtdbtk_classification_file "$TAX_FILE" \
    --cpus "$CPUS" \
    --force \
    --tmpdir "$TMPDIR" \
    --keep_intermediates \
    --debug || exit_code=$?

if [[ $exit_code -eq 0 ]]; then
    log_info "Tree written under $OUT_DIR"
else
    log_error "gtdbtk de_novo_wf failed with exit code $exit_code"
fi

log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
exit $exit_code
