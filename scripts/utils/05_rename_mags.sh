#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Replaces underscores with dots in MAG names, across the FASTA filenames and
# the two tables that reference them.
#
# The binner emitted names like A4_BTM_megahit_bin_25 while everything
# downstream expects A4.BTM.megahit.bin.25. Renaming the files alone would
# break the join, so quality_report.tsv and the GTDB-Tk summary are rewritten
# to match in the same pass.
#
# Only the first column of each table changes. In the taxonomy summary the
# header row is skipped, because substituting there would mangle the column
# names, which legitimately contain underscores.
#
# Both tables get a .backup alongside them before being touched.

load_project_config
require_vars PROJECT_DIR FILTERED_SUBDIR QUALITY_SUBDIR

log_info "Replacing underscores with dots in MAG names"

# ---- 1. the FASTA files ----

FILTERED_DIR="${PROJECT_DIR}/${FILTERED_SUBDIR}"

if [ ! -d "$FILTERED_DIR" ]; then
    log_error "Directory does not exist: $FILTERED_DIR"
    exit 1
fi

log_info "[1/3] Renaming .fa files in ${FILTERED_SUBDIR}/"

count=0
shopt -s nullglob
# The list is taken before any renaming, so the glob does not pick up names
# this loop has already changed.
files=("$FILTERED_DIR"/*.fa)
shopt -u nullglob

for file in "${files[@]}"; do
    [ -f "$file" ] || continue
    dir=$(dirname "$file")
    base=$(basename "$file")
    newbase="${base//_/.}"

    if [ "$base" != "$newbase" ]; then
        mv "$file" "${dir}/${newbase}"
        log_info "  $base to $newbase"
        ((++count))
    fi
done

log_info "  Renamed: $count"

# ---- 2. quality_report.tsv ----

QUALITY_FILE="${PROJECT_DIR}/${QUALITY_SUBDIR}/quality_report.tsv"

if [ ! -f "$QUALITY_FILE" ]; then
    log_error "File does not exist: $QUALITY_FILE"
    exit 1
fi

log_info "[2/3] Rewriting quality_report.tsv"
cp "$QUALITY_FILE" "${QUALITY_FILE}.backup"

# CheckM2 writes no header row here, so every line is a record.
awk 'BEGIN {FS=OFS="\t"} {gsub(/_/, ".", $1); print}' "$QUALITY_FILE" > "${QUALITY_FILE}.tmp"
mv "${QUALITY_FILE}.tmp" "$QUALITY_FILE"
log_info "  Backup: ${QUALITY_FILE}.backup"

# ---- 3. the GTDB-Tk summary ----

TAXONOMY_FILE="${PROJECT_DIR}/${TAXONOMY_SUBDIR:-03_taxonomy_assigned}/classify/tax_.ar53.summary.tsv"

if [ ! -f "$TAXONOMY_FILE" ]; then
    log_error "File does not exist: $TAXONOMY_FILE"
    exit 1
fi

log_info "[3/3] Rewriting tax_.ar53.summary.tsv"
cp "$TAXONOMY_FILE" "${TAXONOMY_FILE}.backup"

awk 'BEGIN {FS=OFS="\t"} NR==1 {print; next} {gsub(/_/, ".", $1); print}' \
    "$TAXONOMY_FILE" > "${TAXONOMY_FILE}.tmp"
mv "${TAXONOMY_FILE}.tmp" "$TAXONOMY_FILE"
log_info "  Backup: ${TAXONOMY_FILE}.backup"

log_info "Done."
log_info "  $count .fa files renamed under ${FILTERED_SUBDIR}/"
log_info "  ${QUALITY_SUBDIR}/quality_report.tsv rewritten"
log_info "  classify/tax_.ar53.summary.tsv rewritten"
log_info "To undo the table edits:"
log_info "  mv ${QUALITY_FILE}.backup ${QUALITY_FILE}"
log_info "  mv ${TAXONOMY_FILE}.backup ${TAXONOMY_FILE}"
