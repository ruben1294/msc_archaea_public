#!/usr/bin/env bash

# Script: 03_post_KO_assign_concat_KOs.sh
# Description: Concatenate all KO assignment files, filter best hits for duplicates
# Date: 2025/12/30
# Last modified: 2026/01/02

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

load_project_config
require_vars RESULTS_DIR KO_ASSIGNMENT_SUBDIR

# ---- Configuration ----

# The KofamKOALA submission date is part of the directory name, so the most
# recent one is taken unless a path is given on the command line.
INPUT_DIR="${1:-}"
if [[ -z "$INPUT_DIR" ]]; then
    INPUT_DIR=$(find "${RESULTS_DIR}/${KO_ASSIGNMENT_SUBDIR}" -maxdepth 1 -type d -name "*kofam_web" | sort | tail -1)
fi
if [[ -z "$INPUT_DIR" || ! -d "$INPUT_DIR" ]]; then
    log_error "No KofamKOALA output directory found. Pass one as the first argument."
    exit 1
fi

# Output filenames
RAW_FILE="${INPUT_DIR}/all_kofam_KOs.raw.tsv"
TRANSFORMED_FILE="${INPUT_DIR}/all_kofam_KOs.transformed.tsv"
PARSED_FILE="${INPUT_DIR}/all_kofam_KOs.parsed.tsv"
BY_MAG_FILE="${INPUT_DIR}/all_kofam_KOs.byMAG.tsv"
MAG_STATS_FILE="${INPUT_DIR}/all_kofam_KOs.MAG_stats.tsv"

# Log file
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${INPUT_DIR}/03_post_KO_assign_${TIMESTAMP}.log"

# ---- Functions ----

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg" | tee -a "$LOG_FILE" >&2
}

error() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*"
    echo "$msg" | tee -a "$LOG_FILE" >&2
    exit 1
}

# ---- Validation ----

log "=========================================="
log "KO Assignment Post-Processing Pipeline"
log "=========================================="
log "Script: 03_post_KO_assign_concat_KOs.sh"
log "Working directory: $(pwd)"
log "Input directory: $INPUT_DIR"
log "Log file: $LOG_FILE"
log ""

[[ -d "$INPUT_DIR" ]] || error "Input directory does not exist: $INPUT_DIR"

# Check for result_ko_part*.txt files (primary source)
shopt -s nullglob
input_ko_files=("$INPUT_DIR"/result_ko_part*.txt)
[[ ${#input_ko_files[@]} -gt 0 ]] || error "No result_ko_part*.txt files found"

# Check for result_part*.txt files (for scoring info)
input_part_files=("$INPUT_DIR"/result_part*.txt)
[[ ${#input_part_files[@]} -gt 0 ]] || error "No result_part*.txt files found"

log "Found ${#input_ko_files[@]} result_ko_part*.txt files (primary source)"
log "Found ${#input_part_files[@]} result_part*.txt files (for scoring)"
log ""

# ---- Cleanup ----

log "Cleaning up previous output files..."
rm -f "$INPUT_DIR"/all_kofam_*.tsv || true
rm -f "$INPUT_DIR"/all_kofam_*.tmp || true
log "Cleanup completed"
log ""

# ---- Step 1: Extract score information from result_part*.txt files ----

log "=========================================="
log "Step 1: Extracting score information"
log "=========================================="

SCORES_FILE="${INPUT_DIR}/all_scores.tmp"

# Extract: gene_name, KO, score, e-value from result_part*.txt
for result_file in "${input_part_files[@]}"; do
    log "  Processing $(basename "$result_file") for scores..."
    grep "^\*" "$result_file" | awk '{print $2 "\t" $3 "\t" $5 "\t" $6}'
done > "$SCORES_FILE"

n_scores=$(wc -l < "$SCORES_FILE")
log "[OK] Extracted $n_scores score entries"
log ""

# ---- Step 2: Build best-hit lookup table (for proteins with multiple KOs) ----

log "=========================================="
log "Step 2: Building best-hit lookup table"
log "=========================================="

BEST_HITS_FILE="${INPUT_DIR}/best_hits.tmp"
DISCARDED_FILE="${INPUT_DIR}/all_kofam_KOs.discarded_hits.tsv"

# Find best hit for each protein (highest score, then lowest e-value)
awk -F'\t' -v discarded="$DISCARDED_FILE" '
{
    prot = $1
    ko = $2
    score = $3 + 0  # Convert to number
    evalue = $4 + 0
    
    # Check if this is a better hit
    if (!(prot in best_score) ||
        score > best_score[prot] ||
        (score == best_score[prot] && evalue < best_evalue[prot])) {
        
        # Save previous best as discarded (if exists)
        if (prot in best_ko) {
            print prot "\t" best_ko[prot] "\t" best_score[prot] "\t" best_evalue[prot] >> discarded
        }
        
        # Update best hit
        best_ko[prot] = ko
        best_score[prot] = score
        best_evalue[prot] = evalue
    }
    else {
        # This is not the best hit
        print $0 >> discarded
    }
}
END {
    # Output best hits: protein\tKO
    for (p in best_ko) {
        print p "\t" best_ko[p]
    }
}' "$SCORES_FILE" > "$BEST_HITS_FILE"

n_best=$(wc -l < "$BEST_HITS_FILE")
n_discarded=$(wc -l < "$DISCARDED_FILE" 2>/dev/null || echo 0)
log "[OK] Best hits selected: $n_best"
log "[OK] Discarded duplicate hits: $n_discarded"
log ""

# ---- Step 3: Concatenate result_ko_part*.txt files (ALL proteins) ----

log "=========================================="
log "Step 3: Concatenating all proteins"
log "=========================================="

TEMP_RAW="${INPUT_DIR}/all_kofam_KOs.temp.tsv"

# Concatenate all result_ko_part*.txt files
cat "${input_ko_files[@]}" > "$TEMP_RAW"

n_total=$(wc -l < "$TEMP_RAW")
log "[OK] Total proteins from result_ko_part*.txt: $n_total"
log ""

# ---- Step 4: Replace multiple KO assignments with best hit ----

log "=========================================="
log "Step 4: Filtering by best hits"
log "=========================================="

awk -F'\t' -v OFS='\t' -v best_file="$BEST_HITS_FILE" '
BEGIN {
    # Load best hits into associative array
    while ((getline < best_file) > 0) {
        best_ko[$1] = $2
    }
    close(best_file)
}
{
    prot_id = $1
    ko = $2
    
    # If protein has a best hit in our lookup table, use it
    if (prot_id in best_ko) {
        ko = best_ko[prot_id]
    }
    
    # Output: protein_id\tKO (or protein_id alone if no KO)
    if (ko != "") {
        print prot_id, ko
    } else {
        print prot_id
    }
}' "$TEMP_RAW" > "$RAW_FILE"

rm "$TEMP_RAW" "$SCORES_FILE" "$BEST_HITS_FILE"

n_raw=$(wc -l < "$RAW_FILE")
n_with_ko=$(awk 'NF==2' "$RAW_FILE" | wc -l)
n_without_ko=$((n_raw - n_with_ko))

log "[OK] RAW file created: $RAW_FILE"
log "  Total proteins (deduplicated): $n_raw"
log "  With KO: $n_with_ko"
log "  Without KO: $n_without_ko"
log ""

# ---- Step 5: Transform protein IDs (replace _ with . and | with _) ----

log "=========================================="
log "Step 5: Transforming protein IDs"
log "=========================================="

awk -F'\t' -v OFS='\t' '
{
    prot_id = $1
    ko = $2
    
    # Replace underscores with dots
    gsub(/_/, ".", prot_id)
    
    # Replace pipes with underscores
    gsub(/\|/, "_", prot_id)
    
    # Print transformed ID and KO
    if (ko != "") {
        print prot_id, ko
    } else {
        print prot_id
    }
}' "$RAW_FILE" > "$TRANSFORMED_FILE"

n_transformed=$(wc -l < "$TRANSFORMED_FILE")
log "[OK] Transformed file created: $TRANSFORMED_FILE ($n_transformed lines)"
log ""

# ---- Step 6: Parse MAG names and genes ----

log "=========================================="
log "Step 6: Parsing MAG names and genes"
log "=========================================="

awk -F'\t' -v OFS='\t' '
{
    full_id = $1
    ko = $2
    
    # Split by first underscore to get MAG and gene
    idx = index(full_id, "_")
    if (idx > 0) {
        mag_name = substr(full_id, 1, idx-1)
        gene_name = substr(full_id, idx+1)
        
        # Print MAG, gene, KO
        if (ko != "") {
            print mag_name, gene_name, ko
        } else {
            print mag_name, gene_name, ""
        }
    }
}' "$TRANSFORMED_FILE" > "$PARSED_FILE"

n_parsed=$(wc -l < "$PARSED_FILE")
log "[OK] Parsed file created: $PARSED_FILE ($n_parsed lines)"
log ""

# ---- Step 7: Generate unique MAG-KO assignments ----

log "=========================================="
log "Step 7: Generating MAG-KO assignments"
log "=========================================="

awk -F'\t' 'NF >= 3 && $3 != "" {print $1 "\t" $3}' "$PARSED_FILE" | \
    sort -u > "$BY_MAG_FILE"

n_mag_ko=$(wc -l < "$BY_MAG_FILE")
log "[OK] MAG-KO assignments: $BY_MAG_FILE ($n_mag_ko unique pairs)"
log ""

# ---- Step 8: Generate per-MAG statistics (with a multi-KO column) ----

log "=========================================="
log "Step 8: Calculating per-MAG statistics"
log "=========================================="

# Create a mapping file: original_MAG -> transformed_MAG
MAPPING_FILE="${INPUT_DIR}/mag_name_mapping.tmp"

awk -F'\t' '{
    # Get original protein ID from RAW_FILE
    orig_id = $1
    
    # Store original MAG name (before pipe)
    split(orig_id, parts, "|")
    orig_mag = parts[1]
    
    # Mark that we have seen this MAG
    seen[orig_mag] = 1
}
END {
    for (mag in seen) {
        # Transform MAG name
        transformed = mag
        gsub(/_/, ".", transformed)
        gsub(/\|/, "_", transformed)
        print mag "\t" transformed
    }
}' "$RAW_FILE" > "$MAPPING_FILE"

# Get proteins with multiple KOs from discarded file
MULTI_KO_FILE="${INPUT_DIR}/multi_ko_proteins.tmp"

awk -F'\t' '{
    # Extract protein ID from discarded hits
    prot_id = $1
    
    # Count unique proteins (those that had multiple KOs)
    seen[prot_id] = 1
}
END {
    for (prot in seen) {
        print prot
    }
}' "$DISCARDED_FILE" > "$MULTI_KO_FILE"

# Calculate statistics using PARSED_FILE (which has deduplicated data)
awk -F'\t' -v OFS='\t' -v mapping_file="$MAPPING_FILE" -v multi_ko_file="$MULTI_KO_FILE" '
BEGIN {
    # Load mapping: transformed_mag -> original_mag
    while ((getline < mapping_file) > 0) {
        trans_to_orig[$2] = $1
    }
    close(mapping_file)
    
    # Load proteins with multiple KOs
    while ((getline < multi_ko_file) > 0) {
        # Transform protein ID to match PARSED_FILE format
        prot_id = $1
        gsub(/_/, ".", prot_id)
        gsub(/\|/, "_", prot_id)
        multi_ko_prots[prot_id] = 1
    }
    close(multi_ko_file)
}
{
    mag_transformed = $1
    gene = $2
    ko = $3
    
    # Reconstruct full protein ID
    full_id = mag_transformed "_" gene
    
    # Get original MAG name
    mag_original = trans_to_orig[mag_transformed]
    if (mag_original == "") {
        mag_original = mag_transformed  # Fallback
    }
    
    # Count total proteins per MAG (using transformed name as key)
    total[mag_transformed]++
    orig_name[mag_transformed] = mag_original
    
    # Count annotated proteins (with KO)
    if (ko != "") {
        annotated[mag_transformed]++
    }
    
    # Count proteins with multiple KOs
    if (full_id in multi_ko_prots) {
        multi_ko_count[mag_transformed]++
    }
}
END {
    # Store data for sorting
    n = 0
    for (mag in total) {
        data[n++] = mag "\t" orig_name[mag] "\t" total[mag] "\t" ((mag in annotated) ? annotated[mag] : 0) "\t" ((mag in multi_ko_count) ? multi_ko_count[mag] : 0)
    }
    
    # Sort by original MAG name (second column)
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            split(data[i], a, "\t")
            split(data[j], b, "\t")
            if (a[2] > b[2]) {  # Compare original names
                tmp = data[i]
                data[i] = data[j]
                data[j] = tmp
            }
        }
    }
    
    # Print header first
    print "MAG_Name_Original", "MAG_Name_Transformed", "Total_Proteins", "Annotated_Proteins", "Unannotated_Proteins", "Proteins_With_Multiple_KOs", "Percent_Annotated"
    
    # Print sorted data
    for (i = 0; i < n; i++) {
        split(data[i], fields, "\t")
        mag_trans = fields[1]
        mag_orig = fields[2]
        total_prots = fields[3]
        annot_prots = fields[4]
        multi_ko = fields[5]
        unannot_prots = total_prots - annot_prots
        pct_annot = (total_prots > 0) ? (annot_prots / total_prots) * 100 : 0
        
        printf "%s\t%s\t%d\t%d\t%d\t%d\t%.2f\n", mag_orig, mag_trans, total_prots, annot_prots, unannot_prots, multi_ko, pct_annot
    }
}' "$PARSED_FILE" > "$MAG_STATS_FILE"

rm "$MAPPING_FILE" "$MULTI_KO_FILE"

n_mags=$(tail -n +2 "$MAG_STATS_FILE" | wc -l)
log "[OK] MAG statistics file created: $MAG_STATS_FILE"
log "  Number of MAGs: $n_mags"

# Show top 5 MAGs by annotation percentage
log ""
log "Top 5 MAGs by annotation percentage:"
tail -n +2 "$MAG_STATS_FILE" | sort -t$'\t' -k7,7nr | head -5 | \
    awk -F'\t' '{printf "  %s: %.2f%% (%d/%d, %d multi-KO)\n", $1, $7, $4, $3, $6}'
log ""

# ---- Summary ----

log "=========================================="
log "Pipeline Summary"
log "=========================================="
log "Input files:"
log "  result_ko_part*.txt: ${#input_ko_files[@]} files"
log "  result_part*.txt: ${#input_part_files[@]} files"
log ""
log "Processing results:"
log "  Total proteins (deduplicated): $n_raw"
log "  Proteins with KO: $n_with_ko ($(awk "BEGIN {printf \"%.1f\", ($n_with_ko/$n_raw)*100}")%)"
log "  Proteins without KO: $n_without_ko ($(awk "BEGIN {printf \"%.1f\", ($n_without_ko/$n_raw)*100}")%)"
log "  Discarded duplicate hits: $n_discarded"
log "  Unique MAG-KO pairs: $n_mag_ko"
log "  Total MAGs: $n_mags"
log ""
log "Output files:"
log "  1. Raw (deduplicated): $RAW_FILE"
log "  2. Transformed IDs: $TRANSFORMED_FILE"
log "  3. Parsed (MAG/gene/KO): $PARSED_FILE"
log "  4. MAG-KO pairs: $BY_MAG_FILE"
log "  5. MAG statistics: $MAG_STATS_FILE"
log "  6. Discarded hits: $DISCARDED_FILE"
log "  7. Log file: $LOG_FILE"
log ""
log "=========================================="
log "Processing completed successfully!"
log "=========================================="
