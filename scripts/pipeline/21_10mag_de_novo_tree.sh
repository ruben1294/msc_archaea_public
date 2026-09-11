#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 21. A second, independent tree for the same 10 MAGs, this time through
# gtdbtk de_novo_wf over the ar53 markers rather than the hand-built supermatrix
# of steps 18 to 20. Agreement between the two is what makes the placement of
# the novel lineages credible.
#
# It symlinks the 10 MAGs into their own input directory so de_novo_wf sees only
# them, and runs once per outgroup listed below.

trap 'log_error "Script failed at line $LINENO with exit code $?"; exit 1' ERR

log_job_header "GTDB-Tk de novo tree, 10 MAGs without species assignment"
load_project_config

require_vars PROJECT_DIR FILTERED_SUBDIR TAXONOMY_SUBDIR

# GTDB-Tk writes some intermediates relative to the working directory.
cd "$PROJECT_DIR" || { log_error "Cannot enter $PROJECT_DIR"; exit 1; }

MAG_SOURCE="$PROJECT_DIR/$FILTERED_SUBDIR"
BASE_OUT_DIR="$PROJECT_DIR/21_10mag_de_novo_tree"
INPUT_DIR="${BASE_OUT_DIR}/00_input"
TAX_FILE="$PROJECT_DIR/$TAXONOMY_SUBDIR/classify/tax_.ar53.summary.tsv"
PROT_MODEL="LG" # Protein model for tree building. Options: JTT, WAG, and LG.
PREFIX="10mag_archaea"
SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-24}"

# Keep this short: GTDB-Tk runs Prodigal, which dies with "AF_UNIX path too
# long" when the temporary directory sits deep in a filesystem.
# mktemp rather than reusing an inherited TMPDIR: the cleanup below must only
# ever remove a directory this script created. SLURM wrappers export TMPDIR
# pointing into the project tree.
TMPDIR="$(mktemp -d /tmp/gtdbtk_XXXXXXXX)"
export TMPDIR
trap 'rm -rf "${TMPDIR:?}"' EXIT
# Methanobacteriota covers the 53 archaeal markers well and sits equidistant
# from both groups of unclassified MAGs, which is what makes it a fair root
# here. Both the class and the phylum run, since the rooting choice is exactly
# what a reviewer will question.
OUTGROUPS=(
    "c__Methanococci"
    "p__Methanobacteriota"
)

# The 10 target MAGs (no species-level assignment, MIMAG medium-quality or better)
TARGET_MAGS=(
    "A4.1000.megahit.bin.28"
    "A4.BTM.megahit.bin.15"
    "A4.BTM.megahit.bin.25"
    "A4.BTM.megahit.bin.33"
    "A4.BTM.megahit.bin.58"   # WARNING: 6.5% contamination
    "A6.BTM.megahit.bin.3"
    "A6.BTM.megahit.bin.9"
    "PO1.BTM.megahit.bin.42"
    "PO1.BTM.megahit.bin.7"
    "Y6.1000.megahit.bin.23"
)

log_debug "MAG source: $MAG_SOURCE"

# ---- Step 0, validate MAG inputs ----
log_info "Step 0: Validating target MAG files"
[[ -d "$MAG_SOURCE" ]] || { log_error "MAG source not found: $MAG_SOURCE"; exit 1; }

missing_mags=0
for mag in "${TARGET_MAGS[@]}"; do
    fa="${MAG_SOURCE}/${mag}.fa"
    if [[ ! -f "$fa" ]]; then
        log_error "  MAG not found: $fa"
        ((++missing_mags))
    else
        log_debug "  [OK] $mag"
    fi
done
[[ $missing_mags -eq 0 ]] || { log_error "$missing_mags MAG file(s) missing. Aborting."; exit 1; }
log_info "All 10 target MAGs found."

# ---- Step 1, build filtered input directory ----
log_info "Step 1: Building filtered input directory"
mkdir -p "${INPUT_DIR}"

log_info "Symlinking 10 target MAGs..."
for mag in "${TARGET_MAGS[@]}"; do
    src="${MAG_SOURCE}/${mag}.fa"
    dst="${INPUT_DIR}/${mag}.fa"
    if [[ ! -L "$dst" ]]; then
        ln -s "$src" "$dst"
        log_debug "  Linked: $mag.fa"
    else
        log_debug "  Already exists: $mag.fa"
    fi
done

total_genomes=$(find "${INPUT_DIR}" -maxdepth 1 -name "*.fa" | wc -l)
log_info "Total genomes in input dir: $total_genomes"

# ---- Step 2, GTDB-Tk de novo workflow (one run per outgroup) ----
log_info "Step 2: Running GTDB-Tk de novo workflow (${#OUTGROUPS[@]} outgroups)"
OVERALL_EXIT=0

for OUTGROUP_TAX in "${OUTGROUPS[@]}"; do
    # Derive a safe directory suffix from the taxon string (e.g. c__Methanococci to c__Methanococci)
    SUFFIX="${OUTGROUP_TAX//\//_}"   # replace any slashes, just in case
    OUT_DIR="${BASE_OUT_DIR}/${SUFFIX}"
    log_info "Outgroup: $OUTGROUP_TAX"
    log_info "Output:   $OUT_DIR"
    mkdir -p "$OUT_DIR" || { log_error "Failed to create $OUT_DIR"; exit 1; }

    # set -e would abort before $? could be read, so the status is captured
    # inline and the loop continues to the next outgroup.
    RUN_EXIT=0
    gtdbtk de_novo_wf \
        --genome_dir "$INPUT_DIR" \
        --outgroup_taxon "$OUTGROUP_TAX" \
        --extension "fa" \
        --archaea \
        --out_dir "$OUT_DIR" \
        --prot_model "$PROT_MODEL" \
        --write_single_copy_genes \
        --prefix "$PREFIX" \
        --rnd_seed 42 \
        --gtdbtk_classification_file "$TAX_FILE" \
        --cpus "${SLURM_CPUS_PER_TASK}" \
        --force \
        --tmpdir "$TMPDIR" \
        --keep_intermediates \
        --debug || RUN_EXIT=$?

    if [ "$RUN_EXIT" -eq 0 ]; then
        log_info "  [OK] Completed: $OUTGROUP_TAX"
    else
        log_error "  [FAIL] Failed (exit $RUN_EXIT): $OUTGROUP_TAX"
        OVERALL_EXIT=$RUN_EXIT
    fi
done

if [ $OVERALL_EXIT -eq 0 ]; then
    log_info "All de novo tree runs completed successfully"
    log_info "  Output structure:"
    for OUTGROUP_TAX in "${OUTGROUPS[@]}"; do
        log_info "    ${BASE_OUT_DIR}/${OUTGROUP_TAX}/"
    done
else
    log_error "One or more runs failed, check logs above"
fi
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
exit $OVERALL_EXIT
