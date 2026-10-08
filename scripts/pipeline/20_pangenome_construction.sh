#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 20. Builds the anvi'o 9 pangenomes. It reformats every MAG, turns each
# into a contigs database, annotates it, and groups the MAGs by species into one
# pangenome per species.
#
# Only species with at least four MAGs get a pangenome; below that the core and
# accessory split is not meaningful. The gene clusters this produces are what
# steps 21 to 25 annotate, summarise and clean, and what the R scripts under
# scripts/R/ read for the pangenome figures.

# --- Error handling ---
trap 'log_error "Script failed at line $LINENO with exit code $?"; exit 1' ERR

# --- Header ---
log_info "EASY anvi'o9 pangenome construction pipeline"
log_info "SLURM job information"
log_info "Job ID: ${SLURM_JOB_ID:-N/A}"
log_info "Job name: ${SLURM_JOB_NAME:-interactive}"
log_info "Node: ${SLURMD_NODENAME:-$(hostname)}"
log_info "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
load_project_config

require_commands anvi-script-reformat-fasta anvi-gen-contigs-database anvi-run-hmms anvi-run-ncbi-cogs anvi-gen-genomes-storage anvi-pan-genome

# --- Validate required config variables ---
required_vars=(PROJECT_DIR RESULTS_DIR FILTERED_SUBDIR PANGENOME_SUBDIR TEMP_SUBDIR MAGS_METADATA)
for var in "${required_vars[@]}"; do
    if [[ -z "${!var:-}" ]]; then
        log_error "Required config variable '$var' is not set or is empty"
        exit 1
    fi
done

SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-52}"

# --- Define directory paths ---
ORIGINAL_INPUT="$RESULTS_DIR/$FILTERED_SUBDIR"
PANGENOME_DIR="$RESULTS_DIR/$PANGENOME_SUBDIR"
TEMP_DIR="$PROJECT_DIR/$TEMP_SUBDIR"
ALL_MAGS="$PANGENOME_DIR/all_mags"
ALL_MAGS_REFORMATTED="$ALL_MAGS/reformatted"
ALL_MAGS_CONTIGS_DB="$ALL_MAGS/contigs_db"
ALL_MAGS_GENE_REPORTS="$ALL_MAGS/gene_calling_reports"
# Step 5 variables (commented out while step is disabled)
#PRODIGAL_DIR="$RESULTS_DIR/$MAG_DISTILL_SUBDIR/ANNOTATE/PRODIGAL"
#NCBI_COGS_OUTPUT_DIR="$ALL_MAGS/COGS_from_genes"
#CLEAN_PRODIGAL_DIR="${PANGENOME_DIR}/genes_clean"

# --- Validate input directories ---
if [[ ! -d "$ORIGINAL_INPUT" ]]; then
    log_error "Input directory does not exist: $ORIGINAL_INPUT"
    exit 1
fi

if [[ ! -f "$MAGS_METADATA" ]]; then
    log_error "Metadata file not found: $MAGS_METADATA"
    exit 1
fi

# --- Create output directories ---
mkdir -p "$PANGENOME_DIR" \
    "$ALL_MAGS_REFORMATTED" "$ALL_MAGS_CONTIGS_DB" "$ALL_MAGS_GENE_REPORTS" "$TEMP_DIR"
log_info "Output directory: $PANGENOME_DIR"

# ---- part 1: process every MAG ----

# --- Step 1: Reformat MAG FASTA files for anvi'o9 ---
log_info "Step 1: Reformatting MAG FASTA files..."
mag_count=0
mag_skip=0
for mag_file in "${ORIGINAL_INPUT}"/*.fa; do
    [[ -f "$mag_file" ]] || continue
    mag=$(basename "$mag_file" .fa)
    output="${ALL_MAGS_REFORMATTED}/${mag}_fixed.fa"
    if [[ -f "$output" ]]; then
        log_info "  Skipping MAG $mag (output already exists)"
        ((++mag_skip))
        continue
    fi
    ((++mag_count))
    log_info "  Reformatting MAG $mag_count: $mag"
    anvi-script-reformat-fasta "$mag_file" \
        -o "$output" \
        -l 1000 \
        --simplify-names \
        --report-file "${ALL_MAGS_REFORMATTED}/${mag}_deflines.txt"
done
log_info "Step 1 complete: $mag_count MAGs reformatted, $mag_skip skipped"

# --- Step 2: Generate contigs databases ---
log_info "Step 2: Generating contigs databases..."
db_count=0
db_skip=0
for mag in "${ALL_MAGS_REFORMATTED}"/*_fixed.fa; do
    [[ -f "$mag" ]] || continue
    sample_name=$(basename "$mag" _fixed.fa)
    output="${ALL_MAGS_CONTIGS_DB}/${sample_name}.db"
    if [[ -f "$output" ]]; then
        log_info "  Skipping $sample_name (contigs database already exists)"
        ((++db_skip))
        continue
    fi
    ((++db_count))
    log_info "  Generating contigs database $db_count: $sample_name"
    anvi-gen-contigs-database \
        -f "$mag" \
        -o "$output" \
        -n "$sample_name MAG" \
        --num-threads "$SLURM_CPUS_PER_TASK" \
        --full-gene-calling-report "${ALL_MAGS_GENE_REPORTS}/${sample_name}_genes.tsv"
done
log_info "Step 2 complete: $db_count contigs databases generated, $db_skip skipped"

# --- Step 3: Run HMM searches ---
log_info "Step 3: Running HMM searches..."
hmm_count=0
hmm_skip=0
for db in "${ALL_MAGS_CONTIGS_DB}"/*.db; do
    [[ -f "$db" ]] || continue
    sentinel="${db%.db}.hmms_done"
    if [[ -f "$sentinel" ]]; then
        log_info "  Skipping $(basename "$db") (HMMs already run)"
        ((++hmm_skip))
        continue
    fi
    ((++hmm_count))
    log_info "  Running HMMs on database $hmm_count: $(basename "$db")"
    anvi-run-hmms \
        -c "$db" \
        --also-scan-trnas \
        --num-threads "$SLURM_CPUS_PER_TASK" \
        -I "Archaea_76,Ribosomal_RNA_16S" \
        --just-do-it
    touch "$sentinel"
done
log_info "Step 3 complete: HMMs run on $hmm_count databases, $hmm_skip skipped"

# --- Step 4: Run NCBI COG annotation on contigs databases ---
# 2026/02/17: This function has two input types: a contig database or a FASTA file.
# For now, I don't know exactly the gene calling output from anvi-gen-contigs-database,
# so I will run it on the contigs databases.
log_info "Step 4: Running NCBI COG annotation on contigs databases..."
cog_count=0
cog_skip=0
for db in "${ALL_MAGS_CONTIGS_DB}"/*.db; do
    [[ -f "$db" ]] || continue
    sentinel="${db%.db}.cogs_done"
    if [[ -f "$sentinel" ]]; then
        log_info "  Skipping $(basename "$db") (NCBI COGs already run)"
        ((++cog_skip))
        continue
    fi
    ((++cog_count))
    log_info "  Running NCBI COGs on database $cog_count: $(basename "$db")"
    anvi-run-ncbi-cogs \
        -c "$db" \
        --num-threads "$SLURM_CPUS_PER_TASK" \
        --temporary-dir-path "$TEMP_DIR"
    touch "$sentinel"
done
log_info "Step 4 complete: NCBI COGs run on $cog_count contigs databases, $cog_skip skipped"

# --- Step 5: Run NCBI COG annotation on Prodigal gene calls from DRAM2 ---
log_info "Step 5: Running NCBI COG annotation on Prodigal gene calls..."
log_info "Skipping step 5 for now..."

#for faa in "${PRODIGAL_DIR}"/*.faa; do
#    [[ -f "$faa" ]] || continue
#    base=$(basename "$faa")
#    clean_output="${CLEAN_PRODIGAL_DIR}/${base}"
#    if [[ ! -f "$clean_output" ]]; then
        # Replace stop codons with X so anvi'o accepts the file
#        sed 's/\*/X/g' "$faa" > "$clean_output"
#    fi
#done

#faa_count=0
#faa_skip=0
#for faa in "${CLEAN_PRODIGAL_DIR}"/*.faa; do
#    [[ -f "$faa" ]] || continue
#    output="${NCBI_COGS_OUTPUT_DIR}/$(basename "$faa" .faa).tsv"
#    if [[ -f "$output" ]]; then
#        log_info "  Skipping $(basename "$faa") (NCBI COGs output already exists)"
#        ((++faa_skip))
#        continue
#    fi
#    ((++faa_count))
#    log_info "  Running NCBI COGs on gene file $faa_count: $(basename "$faa")"
#    anvi-run-ncbi-cogs \
#        -f "$faa" \
#        -o "$output" \
#        --num-threads "$SLURM_CPUS_PER_TASK" \
#        --temporary-dir-path "$TEMP_DIR"
#done
#log_info "Step 5 complete: NCBI COGs run on $faa_count gene files, $faa_skip skipped"

# ---- part 2: pangenome construction ----
log_info "Starting pangenomics workflow..."

# --- Cleanup: remove output files from previous runs ---
SPECIES_MAGS_FILE="${PANGENOME_DIR}/MAGs_pangenome.tsv"
CONTIGS_BINS_FILE="${PANGENOME_DIR}/contigs_to_bins.tsv"
for f in "$SPECIES_MAGS_FILE" "$CONTIGS_BINS_FILE"; do
    if [[ -f "$f" ]]; then
        log_info "  Removing previous run file: $(basename "$f")"
        rm "$f"
    fi
done

# --- Step 6: Filter species with >= 4 MAGs from metadata ---
log_info "Step 6: Counting MAGs per species from MAGs_metadata.tsv..."

awk -F'\t' '
NR > 1 {
    if ($9 == "" || $9 == "N/A" || $9 == "NA") next
    name = $1; sub(/\.fa$/, "", name)
    ++count[$9]
    mag_names[$9] = (length(mag_names[$9]) > 0 ? mag_names[$9] SUBSEP name : name)
}
END {
    for (sp in count) {
        if (count[sp] >= 4) {
            n = split(mag_names[sp], arr, SUBSEP)
            for (i = 1; i <= n; i++) {
                print arr[i] "\t" sp
            }
        }
    }
}' "$MAGS_METADATA" > "$SPECIES_MAGS_FILE"

species_mag_count=$(wc -l < "$SPECIES_MAGS_FILE")
log_info "Step 6 complete: $species_mag_count MAGs belonging to species with >= 4 representatives"

if [[ "$species_mag_count" -eq 0 ]]; then
    log_warn "No species found with >= 4 MAGs. Pangenomics workflow has nothing to process."
    log_info "Pipeline finished (no pangenome to build)"
    exit 0
fi

# --- Step 7: Extract contig headers from filtered MAG bins ---
log_info "Step 7: Extracting contig headers from filtered MAG bins..."
: > "$CONTIGS_BINS_FILE"

while IFS=$'\t' read -r mag species; do
    bin_file="${ALL_MAGS_REFORMATTED}/${mag}_fixed.fa"
    if [[ -f "$bin_file" ]]; then
        { grep '^>' "$bin_file" || true; } | while IFS= read -r header; do
            printf '%s\t%s\n' "$header" "${mag}_fixed.fa"
        done >> "$CONTIGS_BINS_FILE"
    else
        log_warn "No bin file found for MAG: $mag in $ALL_MAGS_REFORMATTED"
    fi
done < "$SPECIES_MAGS_FILE"

log_info "Step 7 complete: $(wc -l < "$CONTIGS_BINS_FILE") contig-to-bin entries written"

# --- Step 8: Organize filtered MAGs by species ---
log_info "Step 8: Organizing filtered MAGs by species..."
PANGENOME_WORK_DIR="$PANGENOME_DIR/pangenome"
mkdir -p "$PANGENOME_WORK_DIR"

while IFS=$'\t' read -r mag species; do
    species_dir="${PANGENOME_WORK_DIR}/${species//[^a-zA-Z0-9_]/_}"
    mkdir -p "${species_dir}/reformatted" "${species_dir}/contigs_db" "${species_dir}/gene_calling_reports"

    # Reformatted FASTA and deflines
    if [[ -f "${ALL_MAGS_REFORMATTED}/${mag}_fixed.fa" ]]; then
        cp "${ALL_MAGS_REFORMATTED}/${mag}_fixed.fa" "${species_dir}/reformatted/"
    else
        log_warn "Reformatted FASTA not found for MAG: $mag"
    fi
    if [[ -f "${ALL_MAGS_REFORMATTED}/${mag}_deflines.txt" ]]; then
        cp "${ALL_MAGS_REFORMATTED}/${mag}_deflines.txt" "${species_dir}/reformatted/"
    fi
    # Contigs database
    if [[ -f "${ALL_MAGS_CONTIGS_DB}/${mag}.db" ]]; then
        cp "${ALL_MAGS_CONTIGS_DB}/${mag}.db" "${species_dir}/contigs_db/"
    else
        log_warn "Contigs database not found for MAG: $mag"
    fi
    # Gene calling reports
    if [[ -f "${ALL_MAGS_GENE_REPORTS}/${mag}_genes.tsv" ]]; then
        cp "${ALL_MAGS_GENE_REPORTS}/${mag}_genes.tsv" "${species_dir}/gene_calling_reports/"
    fi
done < "$SPECIES_MAGS_FILE"

log_info "Step 8 complete: MAGs organized into $(find "$PANGENOME_WORK_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l) species folders"

# --- Step 9: Build external_genomes.tsv and run anvi-gen-genomes-storage per species ---
log_info "Step 9: Building genomes storage per species..."
for species_dir in "${PANGENOME_WORK_DIR}"/*/; do
    [[ -d "$species_dir" ]] || continue
    species_name=$(basename "$species_dir")

    external_genomes_file="${species_dir}/external_genomes.tsv"
    printf 'name\tcontigs_db_path\n' > "$external_genomes_file"
    for db in "${species_dir}/contigs_db"/*.db; do
        [[ -f "$db" ]] || continue
        mag_name=$(basename "$db" .db)
        sanitized_name="${mag_name//[^a-zA-Z0-9_]/_}"
        printf '%s\t%s\n' "$sanitized_name" "$db"
    done >> "$external_genomes_file"

    genomes_db_out="${species_dir}/${species_name}-GENOMES.db"
    if [[ -f "$genomes_db_out" ]]; then
        log_info "  Skipping anvi-gen-genomes-storage for ${species_name//_/ } (output already exists)"
        continue
    fi
    log_info "  Running anvi-gen-genomes-storage for species: ${species_name//_/ }"
    anvi-gen-genomes-storage \
        -e "$external_genomes_file" \
        -o "$genomes_db_out"
done
log_info "Step 9 complete"

# --- Step 10: Run anvi-pan-genome per species ---
log_info "Step 10: Running pangenome analysis per species..."
for species_dir in "${PANGENOME_WORK_DIR}"/*/; do
    [[ -d "$species_dir" ]] || continue
    species_name=$(basename "$species_dir")
    genomes_db="${species_dir}/${species_name}-GENOMES.db"

    if [[ ! -f "$genomes_db" ]]; then
        log_warn "Genomes storage not found for ${species_name//_/ }, skipping pangenome analysis"
        continue
    fi

    if [[ -d "${species_dir}/pan_output" ]]; then
        log_info "  Skipping anvi-pan-genome for ${species_name//_/ } (output directory already exists)"
        continue
    fi
    log_info "  Running anvi-pan-genome for species: ${species_name//_/ }"
    anvi-pan-genome \
        --genomes-storage "$genomes_db" \
        --project-name "$species_name" \
        --output-dir "${species_dir}/pan_output" \
        --num-threads "$SLURM_CPUS_PER_TASK"
done
log_info "Step 10 complete"

# --- Done ---
log_info "Pipeline finished successfully"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
