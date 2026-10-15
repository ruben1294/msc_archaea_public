#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Replaces underscores with dots in MAG names: A4_BTM_megahit_bin.25 becomes
# A4.BTM.megahit.bin.25, the form data/MAGs_metadata.tsv calls
# NameWithoutExtension and every step from 02 onward uses.
#
# Run it once, right after utils/01_filter_mags.py and before step 02, so that
# GUNC, GTDB-Tk and everything after them see the final names. It renames the
# filtered FASTA files and rewrites the Name column of the CheckM2 report to
# match, keeping a .backup of the report. Running it again changes nothing.
#
# If a GTDB-Tk summary already exists, because classification ran before the
# rename, its first column is rewritten the same way. The header row is left
# alone there, since the column names contain underscores of their own.
#
# The R scripts join on NameWithoutExtension after turning underscores into
# dots themselves, so they read either form of these tables.

load_project_config
require_vars RESULTS_DIR FILTERED_SUBDIR QUALITY_SUBDIR

log_info "Replacing underscores with dots in MAG names"

# ---- 1. the FASTA files ----

FILTERED_DIR="${RESULTS_DIR}/${FILTERED_SUBDIR}"

if [ ! -d "$FILTERED_DIR" ]; then
    log_error "Directory does not exist: $FILTERED_DIR"
    log_error "Run utils/01_filter_mags.py first."
    exit 1
fi

log_info "[1/3] Renaming .fa files in ${FILTERED_DIR}/"

count=0
shopt -s nullglob
# The list is taken before any renaming, so the glob does not pick up names
# this loop has already changed.
files=("$FILTERED_DIR"/*.fa)
shopt -u nullglob

for file in "${files[@]}"; do
    [ -f "$file" ] || continue
    dir=$(dirname "$file")
    base=$(basename "$file" .fa)
    newbase="${base//_/.}.fa"

    if [ "${base}.fa" != "$newbase" ]; then
        mv "$file" "${dir}/${newbase}"
        log_info "  ${base}.fa to $newbase"
        ((++count))
    fi
done

log_info "  Renamed: $count"

# ---- 2. quality_report.tsv ----

QUALITY_FILE="${RESULTS_DIR}/${QUALITY_SUBDIR}/quality_report.tsv"

if [ ! -f "$QUALITY_FILE" ]; then
    log_error "File does not exist: $QUALITY_FILE"
    exit 1
fi

log_info "[2/3] Rewriting quality_report.tsv"
[ -f "${QUALITY_FILE}.backup" ] || cp "$QUALITY_FILE" "${QUALITY_FILE}.backup"

# The header row is "Name ...", which has no underscore, so every row can go
# through the same substitution.
awk 'BEGIN {FS=OFS="\t"} {gsub(/_/, ".", $1); print}' "$QUALITY_FILE" > "${QUALITY_FILE}.tmp"
mv "${QUALITY_FILE}.tmp" "$QUALITY_FILE"
log_info "  Backup of the original: ${QUALITY_FILE}.backup"

# ---- 3. the GTDB-Tk summary, only if classification already ran ----

TAXONOMY_FILE="${RESULTS_DIR}/${TAXONOMY_SUBDIR:-03_taxonomy_assigned}/classify/tax_.ar53.summary.tsv"

if [ -f "$TAXONOMY_FILE" ]; then
    log_info "[3/3] Rewriting tax_.ar53.summary.tsv"
    [ -f "${TAXONOMY_FILE}.backup" ] || cp "$TAXONOMY_FILE" "${TAXONOMY_FILE}.backup"

    awk 'BEGIN {FS=OFS="\t"} NR==1 {print; next} {gsub(/_/, ".", $1); print}' \
        "$TAXONOMY_FILE" > "${TAXONOMY_FILE}.tmp"
    mv "${TAXONOMY_FILE}.tmp" "$TAXONOMY_FILE"
    log_info "  Backup of the original: ${TAXONOMY_FILE}.backup"
else
    log_info "[3/3] No GTDB-Tk summary yet, nothing to rewrite"
fi

log_info "Done. Step 02 can run on ${FILTERED_DIR}/."
