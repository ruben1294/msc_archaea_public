#!/usr/bin/env bash
set -euo pipefail

# Step 2. Screens the 54 MAGs for chimerism with GUNC 1.1 against two reference
# databases, ProGenomes3 and GTDB r214. It closes the quality block that step 01
# opens. CheckM2 estimates completeness and contamination from the genes a MAG
# carries; GUNC tests whether those genes come from more than one lineage. It
# runs on the MIMAG-filtered set that utils/01_filter_mags.py leaves in
# results/01_mag_quality/filtered/, after utils/05_rename_mags.sh has given the
# files their final names, and needs nothing from any later step.
#
# GUNC assigns single-copy genes to clades and derives a clade separation score
# (CSS) from how far the assignments spread. Genes spread across several
# distinct clades point to a chimeric assembly. GUNC calls a genome chimeric at
# CSS above 0.45, with pass.GUNC False.
#
# ProGenomes3 gives broad prokaryotic coverage; GTDB r214 carries the deeper
# archaeal taxonomy and is the one to trust for this dataset.
#
# GUNC flags genomes and this step removes none: all 54 go on to step 03. In
# this dataset 22 of the 54 fail against GTDB r214 and 24 against ProGenomes3.
# Many belong to archaeal lineages the references represent poorly, where a high
# CSS comes with a low reference representation score and reflects gaps in the
# database more than contamination. MAG_quality_assessment.R reports both flags
# per MAG. The contig-level output (--contig_taxonomy_output) is kept so the
# contigs behind a flag can be checked one by one against the step 10
# annotations.
#
# Databases come from utils/download_gunc_1.1_db.sh; set GUNC_DB in server.env
# to where that put them. Pass --clean to discard previous output and start
# over; without it, completed passes are detected and skipped.

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

trap 'log_error "Failed at line $LINENO with exit code $?"; exit 1' ERR

log_job_header "GUNC 1.1 chimerism screening"

load_project_config
require_vars PROJECT_DIR RESULTS_DIR FILTERED_SUBDIR GUNC_DB TEMP_SUBDIR
require_commands gunc

CLEAN_RUN=false
for arg in "$@"; do
    [[ "$arg" == "--clean" ]] && CLEAN_RUN=true
done

THREADS="${SLURM_CPUS_PER_TASK:-24}"
MAG_EXT="fa"

MAG_DIR="$RESULTS_DIR/$FILTERED_SUBDIR"
OUT_DIR="$RESULTS_DIR/${CHIMERISM_SUBDIR:-02_chimerism}"
GUNC_OUT="$OUT_DIR/gunc_output"
PLOT_DIR="$OUT_DIR/plots"
TEMP_DIR="$PROJECT_DIR/$TEMP_SUBDIR"

# Label, GUNC database file, output subdirectory name.
DB_LABELS=("ProGenomes3" "GTDB r214")
DB_FILES=("${GUNC_DB}/gunc_db_progenomes3.dmnd" "${GUNC_DB}/gunc_db_gtdb214.dmnd")
DB_SLUGS=("progenomes_3" "gtdb_214")
# The maxCSS filename GUNC writes is keyed to the database name it was given.
DB_OUTNAMES=("GUNC.progenomes_3.maxCSS_level.tsv" "GUNC.gtdb_214.maxCSS_level.tsv")
# The summaries MAG_quality_assessment.R reads, one per database.
SUMMARY_NAMES=("gunc_summary_progenomes3.tsv" "gunc_summary_gtdb214.tsv")

log_info "MAG input directory : $MAG_DIR"
log_info "Output directory    : $OUT_DIR"
log_info "GUNC database dir   : $GUNC_DB"
log_info "Threads             : $THREADS"

if [[ ! -d "$MAG_DIR" ]]; then
    log_error "MAG directory not found: $MAG_DIR"
    exit 1
fi

N_MAGS=$(find "$MAG_DIR" -maxdepth 1 -name "*.${MAG_EXT}" | wc -l)
log_info "MAGs found: $N_MAGS"

if [[ "$N_MAGS" -eq 0 ]]; then
    log_error "No .${MAG_EXT} files in $MAG_DIR"
    exit 1
fi

# Both databases are checked before any work starts, so a missing one stops the
# step here instead of leaving half the screen undone behind a clean exit.
missing_dbs=()
for db_file in "${DB_FILES[@]}"; do
    [[ -f "$db_file" ]] || missing_dbs+=("$db_file")
done
if [[ ${#missing_dbs[@]} -gt 0 ]]; then
    log_error "GUNC database not found: ${missing_dbs[*]}"
    log_error "Run utils/download_gunc_1.1_db.sh, then set GUNC_DB in $CONFIG_FILE"
    exit 1
fi

if [[ "$CLEAN_RUN" == true ]]; then
    log_info "--clean given, removing $OUT_DIR"
    rm -rf "${OUT_DIR:?}"
fi

mkdir -p "$TEMP_DIR" "$OUT_DIR"
for slug in "${DB_SLUGS[@]}"; do
    mkdir -p "${GUNC_OUT}/${slug}" "${PLOT_DIR}/${slug}"
done

# ---- step 1: run GUNC against each database ----
#
# gunc run calls Prodigal internally, maps the genes to the reference with
# DIAMOND, and reports per genome:
#
#   CSS                    clade separation score, chimeric above 0.45
#   contamination_portion  fraction of genes outside the dominant clade
#   pass.GUNC              False marks a genome as chimeric
#
# --sensitive uses DIAMOND's sensitive mode, which recovers more hits for
# archaeal MAGs that map sparsely against the references. --detailed_output
# and --contig_taxonomy_output produce the per-gene and per-contig files that
# the plots and the contig-level cross-check both need.

for i in "${!DB_LABELS[@]}"; do
    label="${DB_LABELS[$i]}"
    db_file="${DB_FILES[$i]}"
    slug="${DB_SLUGS[$i]}"
    out_dir="${GUNC_OUT}/${slug}"
    main_output="${out_dir}/${DB_OUTNAMES[$i]}"

    log_info "Step 1: GUNC against $label"

    if [[ -f "$main_output" ]]; then
        log_info "  Output exists, skipping: $main_output"
        continue
    fi

    gunc run \
        --input_dir "$MAG_DIR" \
        --file_suffix ".${MAG_EXT}" \
        --out_dir "$out_dir" \
        --db_file "$db_file" \
        --threads "$THREADS" \
        --temp_dir "$TEMP_DIR" \
        --sensitive \
        --detailed_output \
        --contig_taxonomy_output

    log_info "  Finished: $main_output"
done

# ---- step 2: summarise ----
#
# Columns of the maxCSS output worth reading:
#
#   genome                 file basename without extension
#   n_genes_called         genes Prodigal predicted
#   n_genes_mapped         genes with a DIAMOND hit
#   n_contigs              contig count; poorer MAGs have more
#   taxonomic_level        rank at which CSS peaks
#   CSS                    clade separation score, 0 to 1
#   contamination_portion  genes outside the dominant clade
#   pass.GUNC              False marks the genome chimeric

summarise_gunc_output() {
    local db_label="$1" main_output="$2" summary_out="$3"

    if [[ ! -f "$main_output" ]]; then
        log_warn "  [$db_label] no output to summarise"
        return 1
    fi

    local n_total n_pass n_fail
    n_total=$(( $(wc -l < "$main_output") - 1 ))
    n_pass=$(awk -F'\t' 'NR>1 && $NF=="True"' "$main_output" | wc -l)
    n_fail=$(( n_total - n_pass ))

    log_info "  [$db_label] evaluated: $n_total, pass: $n_pass, chimeric: $n_fail"

    if [[ "$n_fail" -gt 0 ]]; then
        log_warn "  [$db_label] chimeric genomes:"
        awk -F'\t' 'NR>1 && $NF=="False" {printf "    %-40s CSS=%-6s contamination=%-6s level=%s\n", $1, $(NF-5), $(NF-4), $5}' \
            "$main_output" | while read -r line; do
            log_warn "$line"
        done
    fi

    # Chimeric genomes first, then descending CSS, so the worst cases are at the
    # top of the file. The sort keys are located by column name because GUNC has
    # moved them between versions.
    local pass_col css_col
    pass_col=$(head -1 "$main_output" | tr '\t' '\n' | grep -n "^pass.GUNC$" | cut -d: -f1)
    css_col=$(head -1 "$main_output" | tr '\t' '\n' | grep -n "^clade_separation_score$" | cut -d: -f1)

    {
        head -1 "$main_output"
        tail -n +2 "$main_output" | sort -t$'\t' -k"$pass_col" -k"$css_col" -gr
    } > "$summary_out"

    log_info "  [$db_label] summary: $summary_out"
}

log_info "Step 2: summarising"

declare -A summary_files
n_summarised=0
for i in "${!DB_LABELS[@]}"; do
    slug="${DB_SLUGS[$i]}"
    main_output="${GUNC_OUT}/${slug}/${DB_OUTNAMES[$i]}"
    summary_out="${OUT_DIR}/${SUMMARY_NAMES[$i]}"
    summary_files[$slug]="$summary_out"
    if summarise_gunc_output "${DB_LABELS[$i]}" "$main_output" "$summary_out"; then
        n_summarised=$((n_summarised + 1))
    fi
done

if [[ "$n_summarised" -lt "${#DB_LABELS[@]}" ]]; then
    log_error "GUNC left no output for $(( ${#DB_LABELS[@]} - n_summarised )) of ${#DB_LABELS[@]} databases"
    exit 1
fi

# ---- step 3: per-genome diagnostic plots ----
#
# gunc plot draws one stacked bar per genome showing how gene assignments
# distribute across taxonomic ranks, which identifies the contigs
# carrying the contaminating signal. It takes a single .diamond.*.out file from
# the diamond_output/ directory that --detailed_output produced.

log_info "Step 3: diagnostic plots"

for i in "${!DB_LABELS[@]}"; do
    label="${DB_LABELS[$i]}"
    slug="${DB_SLUGS[$i]}"
    main_output="${GUNC_OUT}/${slug}/${DB_OUTNAMES[$i]}"
    plot_dir="${PLOT_DIR}/${slug}"

    [[ -f "$main_output" ]] || continue

    n_genomes=$(( $(wc -l < "$main_output") - 1 ))
    n_plots=$(find "$plot_dir" -maxdepth 1 -name "*.png" | wc -l)

    if [[ "$n_plots" -ge "$n_genomes" ]]; then
        log_info "  [$label] $n_plots plots already present, skipping"
        continue
    fi

    log_info "  [$label] plotting $n_genomes genomes"
    shopt -s nullglob
    for diamond_file in "${GUNC_OUT}/${slug}/diamond_output/"*.diamond.*.out; do
        gunc plot -d "$diamond_file" -o "$plot_dir"
    done
    shopt -u nullglob
    log_info "  [$label] plots: $plot_dir"
done

log_info "Output files:"
for i in "${!DB_LABELS[@]}"; do
    slug="${DB_SLUGS[$i]}"
    log_info "  ${DB_LABELS[$i]}"
    log_info "    raw      ${GUNC_OUT}/${slug}/"
    log_info "    maxCSS   ${GUNC_OUT}/${slug}/${DB_OUTNAMES[$i]}"
    log_info "    summary  ${summary_files[$slug]}"
    log_info "    plots    ${PLOT_DIR}/${slug}/"
done

log_info "Chimeric when CSS > 0.45 and pass.GUNC is False"
log_info "Genomes evaluated: $N_MAGS"
log_info "These are flags, not a filter: every MAG continues to step 03."
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
