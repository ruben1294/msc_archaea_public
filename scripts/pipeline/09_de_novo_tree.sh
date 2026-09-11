#!/usr/bin/env bash
set -euo pipefail

# Step 9. Builds a de novo archaeal tree with gtdbtk de_novo_wf, placing the 54
# MAGs among the GTDB reference genomes. The result is the 7,022-leaf Newick
# that scripts/pipeline/curate_tree_and_itol.py later prunes for publication.
#
# The outgroup is p__Euryarchaeota. Rooting on p__Thermoplasmatota, the
# alternative left in the code below, is not usable here because most of the
# ingroup sits inside that phylum.
#
# GTDB-Tk runs Prodigal itself and can die with "OSError: AF_UNIX path too long"
# when the temporary directory path is deep. Passing --genes with the proteins
# from step 3 avoids that, but GTDB-Tk then skips the ANI screen and
# classification steps, and this version rejected the flag outright
# ("unrecognized arguments"). The workaround in use is a short TMPDIR.

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

trap 'log_error "Failed at line $LINENO with exit code $?"; exit 1' ERR

log_job_header "GTDB-Tk de novo tree"

load_project_config
require_vars PROJECT_DIR FILTERED_SUBDIR TAXONOMY_SUBDIR

RAW_DIR="${PROJECT_DIR}/${FILTERED_SUBDIR}"
OUT_DIR="${PROJECT_DIR}/${GTDBTK_DENOVO_TREE_SUBDIR:-09_de_novo_tree}"
TAX_FILE="${PROJECT_DIR}/${TAXONOMY_SUBDIR}/classify/tax_.ar53.summary.tsv"

PROT_MODEL="LG"                 # JTT, WAG or LG
PREFIX="54_archaea"
OUTGROUP_TAX="p__Euryarchaeota"
CPUS="${SLURM_CPUS_PER_TASK:-24}"

# Keep this path short. A deep TMPDIR is what triggers the AF_UNIX error above.
TMPDIR="${TMPDIR:-/tmp/gtdbtk_$$}"
mkdir -p "$TMPDIR"
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
    log_error "Run 02_taxonomic_assignment.sh first."
    exit 1
fi

mkdir -p "$OUT_DIR"

# set -e would abort before any $? check, so the exit code is captured here to
# let the closing message distinguish success from failure.
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
