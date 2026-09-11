#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 18, first of four that build a focused tree for the 10 MAGs GTDB-Tk left
# without a species. The 54-MAG tree from steps 09 to 13 places them in a
# 7,022-leaf context that is too coarse to resolve where exactly they sit, so
# this branch rebuilds a small tree around them with hand-picked references.
#
# It downloads archaeal reference genomes and bacterial outgroups from NCBI into
# data/reference_genomes/, combines them with the 10 MAGs from
# 01_raw_files/filtered/, and runs GTDB-Tk identify and align over the ar53
# marker set.
#
# The per-marker HMM-aligned .faa files it leaves in
# 02_align/intermediate_results/markers/ are what step 19 picks up.

trap 'log_error "Script failed at line $LINENO with exit code $?"; exit 1' ERR

log_info "10-MAG focused phylogenomic tree, GTDB-Tk identify + align"
log_info "SLURM Job Information"
log_info "  Job ID:    ${SLURM_JOB_ID:-N/A}"
log_info "  Job Name:  ${SLURM_JOB_NAME:-interactive}"
log_info "  Node:      ${SLURMD_NODENAME:-$(hostname)}"
log_info "  Start:     $(date '+%Y-%m-%d %H:%M:%S')"
# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
load_project_config

# Paths
WORK_DIR="${PROJECT_DIR}/18_10mag_phylogeny"
MAG_SOURCE="${PROJECT_DIR}/01_raw_files/filtered"   # 10 MAGs .fa
REF_SOURCE="${PROJECT_DIR}/data/reference_genomes"  # archaeal refs + outgroups (.fna)
COMBINED_DIR="${WORK_DIR}/00_input"                          # combined input for GTDB-Tk
IDENTIFY_DIR="${WORK_DIR}/01_identify"
ALIGN_DIR="${WORK_DIR}/02_align"

PREFIX="10mag_refs"
SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-16}"

# Prodigal, which GTDB-Tk calls internally, fails with "AF_UNIX path too long"
# when the temporary directory sits deep in a filesystem. mktemp keeps it short
# and keeps it ours, so an inherited TMPDIR is never touched.
TMPDIR="$(mktemp -d /tmp/gtdbtk_XXXXXXXX)"
export TMPDIR
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
log_info "GTDB-Tk temp dir: ${TMPDIR}"

# The 10 target MAGs (no-species assignment, MIMAG medium-quality or better)
# Note: A4.BTM.megahit.bin.58 has 6.5% contamination, included but flagged.
TARGET_MAGS=(
    "A6.BTM.megahit.bin.3"
    "A4.BTM.megahit.bin.15"
    "PO1.BTM.megahit.bin.42"
    "Y6.1000.megahit.bin.23"
    "A4.1000.megahit.bin.28"
    "A4.BTM.megahit.bin.25"
    "A4.BTM.megahit.bin.33"
    "PO1.BTM.megahit.bin.7"
    "A4.BTM.megahit.bin.58"   # WARNING: 6.5% contamination
    "A6.BTM.megahit.bin.9"
)

# ---------------------------------------------------------------------------
# Reference genome list
# ---------------------------------------------------------------------------
# Accessions marked as closest hits come straight from the GTDB-Tk classify_wf
# output. The rest are established archaeal references chosen for context.
# Comment out a line to drop that genome from the tree.

declare -A REF_LABELS    # accession to short label for log messages
declare -A REF_SPECIES   # accession to official display name used in final tree

# --- Thalassarchaeaceae (Poseidoniales / MGIIb-Q1), close GTDB-Tk hits ---
REFS_THALASSO=(
    "GCA_012964335.1"   # Closest ref: A6.3, A4.15, PO1.42  (ANI 93%, sublineage A)
    "GCA_002505495.1"   # Closest ref: Y6.23, A4.28          (ANI 94-95%, sublineage B)
    "GCA_963580235.1"   # Closest placement ref: A4.28       (via closest_placement_reference)
)
REF_LABELS["GCA_012964335.1"]="Thalasso_SublineajeA_ref"
REF_LABELS["GCA_002505495.1"]="Thalasso_SublineajeB_ref1"
REF_LABELS["GCA_963580235.1"]="Thalasso_SublineajeB_ref2"
REF_SPECIES["GCA_012964335.1"]="Thalassarchaeaceae_sp._GCA_012964335.1"
REF_SPECIES["GCA_002505495.1"]="Thalassarchaeaceae_sp._GCA_002505495.1"
REF_SPECIES["GCA_963580235.1"]="Thalassarchaeaceae_sp._GCA_963580235.1"

# --- Thalassarchaeaceae (UBA36), close GTDB-Tk hits ---
REFS_UBA36=(
    "GCA_036153975.1"   # Closest ref: A4.58, A6.9  (ANI 93%, UBA36)
    "GCA_036336595.1"   # UBA36 additional representative
)
REF_LABELS["GCA_036153975.1"]="UBA36_ref1"
REF_LABELS["GCA_036336595.1"]="UBA36_ref2"
REF_SPECIES["GCA_036153975.1"]="UBA36_sp._GCA_036153975.1"
REF_SPECIES["GCA_036336595.1"]="UBA36_sp._GCA_036336595.1"

# --- Nitrosopumilaceae (Nitrososphaerales / Nitrosopelagicus), close GTDB-Tk hits ---
# U25 (GCA_000402075.1) was suppressed by NCBI, so CN25 (GCF_000812185.2)
# stands in for it: a closed 1.23 Mb genome on the single circular contig
# CP007026, and a formally described type strain.
# SCM1 (GCF_000018465.1) joins them as the AOA type strain from the same
# family, which sharpens resolution inside this clade.
REFS_NITROSOPELAGICUS=(
    "GCA_013390375.1"   # Candidatus Nitrosopelagicus sp. (GTDB ref for A4.25, A4.33, PO1.7)
    "GCF_000812185.2"   # Candidatus Nitrosopelagicus brevis CN25, closed genome (REPLACES U25)
    "GCF_000018465.1"   # Nitrosopumilus maritimus SCM1, type strain, closed genome (Nitrosopumilaceae)
)
REF_LABELS["GCA_013390375.1"]="Nitrosopelagicus_ref1"
REF_LABELS["GCF_000812185.2"]="Nitrosopelagicus_brevis_CN25"
REF_LABELS["GCF_000018465.1"]="Nitrosopumilus_maritimus_SCM1"
REF_SPECIES["GCA_013390375.1"]="Ca._Nitrosopelagicus_sp._GCA_013390375.1"
REF_SPECIES["GCF_000812185.2"]="Ca._Nitrosopelagicus_brevis_CN25"
REF_SPECIES["GCF_000018465.1"]="Nitrosopumilus_maritimus_SCM1"

# --- Broader archaeal context (domain-level context, GTDB-R226) ---
# Haloarcula marismortui and Methanothermobacter thermautotrophicus were taken
# out: they weighed the tree down without resolving the clades of interest.
# Candidatus Lokiarchaeum sp. GC14_75 (GCA_900232055.1) came out for a harder
# reason. Its branch length reached 1.386 in the first tree and the long-branch
# attraction that followed distorted the topology throughout.
REFS_BROAD_ARCHAEA=(
    # Thermoproteota, Sulfolobales (hyperthermophilic Crenarchaeota)
    "GCA_000012285.1"   # Sulfolobus acidocaldarius DSM 639

    # Thermoproteota, Nitrososphaerales (broader genus Nitrososphaera)
    "GCA_000259595.1"   # Candidatus Nitrososphaera viennensis EN76

    # Methanobacteriota_B, Methanococcales
    "GCA_000091665.1"   # Methanocaldococcus jannaschii DSM 2661
)
REF_LABELS["GCA_000012285.1"]="Sulfolobus_acidocaldarius"
REF_LABELS["GCA_000259595.1"]="Nitrososphaera_viennensis"
REF_LABELS["GCA_000091665.1"]="Methanocaldococcus_jannaschii"
REF_SPECIES["GCA_000012285.1"]="Sulfolobus_acidocaldarius_DSM_639"
REF_SPECIES["GCA_000259595.1"]="Ca._Nitrososphaera_viennensis_EN76"
REF_SPECIES["GCA_000091665.1"]="Methanocaldococcus_jannaschii_DSM_2661"

# --- Bacterial outgroups (NOT DPANN; placed as external root) ---
# Using GCF_ (RefSeq) accessions, closed reference genomes, type strains
REFS_OUTGROUPS=(
    "GCF_000008545.1"   # Thermotoga maritima MSB8        (Thermotogota)
    "GCF_000008625.1"   # Aquifex aeolicus VF5            (Aquificota)
)
REF_LABELS["GCF_000008545.1"]="Thermotoga_maritima_outgroup"
REF_LABELS["GCF_000008625.1"]="Aquifex_aeolicus_outgroup"
REF_SPECIES["GCF_000008545.1"]="Thermotoga_maritima_MSB8"
REF_SPECIES["GCF_000008625.1"]="Aquifex_aeolicus_VF5"

# Combined list (all groups)
ALL_REFS=(
    "${REFS_THALASSO[@]}"
    "${REFS_UBA36[@]}"
    "${REFS_NITROSOPELAGICUS[@]}"
    "${REFS_BROAD_ARCHAEA[@]}"
    "${REFS_OUTGROUPS[@]}"
)

# ---------------------------------------------------------------------------
# Download function
# ---------------------------------------------------------------------------
# Uses ncbi-datasets-cli if available; falls back to NCBI REST API via wget.
download_genome() {
    local accession="$1"
    local outdir="$2"
    local outfile="${outdir}/${accession}.fna"

    if [[ -f "$outfile" ]]; then
        log_info "  [OK] Already downloaded: $accession"
        return 0
    fi

    local label="${REF_LABELS[$accession]:-$accession}"
    local tmpdir="${outdir}/tmp_${accession}"
    mkdir -p "$tmpdir"

    log_info "  Downloading: $accession  ($label)"

    if command -v datasets &>/dev/null; then
        # ---- ncbi-datasets-cli path ----
        datasets download genome accession "$accession" \
            --include genome \
            --no-progressbar \
            --filename "${tmpdir}/download.zip" 2>/dev/null || {
            log_warn "    datasets CLI failed for $accession; trying wget fallback..."
            rm -rf "$tmpdir"
            mkdir -p "$tmpdir"
            _wget_download "$accession" "$tmpdir" "$outfile" || return 1
            return 0
        }
    else
        # ---- wget path ----
        _wget_download "$accession" "$tmpdir" "$outfile" || return 1
        return 0
    fi

    # Extract .fna from the zip (datasets CLI zip layout)
    unzip -q "${tmpdir}/download.zip" -d "${tmpdir}/extracted" 2>/dev/null || {
        log_error "    Failed to unzip download for $accession"
        rm -rf "$tmpdir"
        return 1
    }

    local fna_file
    fna_file=$(find "${tmpdir}/extracted" -name "*.fna" | head -1)

    # Some NCBI zips contain .fna.gz instead of plain .fna
    if [[ -z "$fna_file" ]]; then
        local fna_gz
        fna_gz=$(find "${tmpdir}/extracted" -name "*.fna.gz" | head -1)
        if [[ -n "$fna_gz" ]]; then
            gunzip -c "$fna_gz" > "$outfile"
            log_info "    [OK] Saved (from .fna.gz): $(basename "$outfile")  ($(du -sh "$outfile" | cut -f1))"
            rm -rf "$tmpdir"
            return 0
        fi
    fi

    if [[ -n "$fna_file" && -f "$fna_file" ]]; then
        mv "$fna_file" "$outfile"
        log_info "    [OK] Saved: $(basename "$outfile")  ($(du -sh "$outfile" | cut -f1))"
    else
        log_error "    No .fna found in zip for $accession"
        rm -rf "$tmpdir"
        return 1
    fi

    rm -rf "$tmpdir"
}

_wget_download() {
    local accession="$1"
    local tmpdir="$2"
    local outfile="$3"
    local api_url="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${accession}/download?include_annotation_type=GENOME_FASTA"

    wget -q --tries=3 --timeout=120 \
        "$api_url" \
        -O "${tmpdir}/download.zip" || {
        log_error "    wget failed for $accession"
        rm -rf "$tmpdir"
        return 1
    }

    unzip -q "${tmpdir}/download.zip" -d "${tmpdir}/extracted" 2>/dev/null || {
        log_error "    Failed to unzip wget download for $accession"
        rm -rf "$tmpdir"
        return 1
    }

    local fna_file
    fna_file=$(find "${tmpdir}/extracted" -name "*.fna" | head -1)

    # Some NCBI zips contain .fna.gz instead of plain .fna
    if [[ -z "$fna_file" ]]; then
        local fna_gz
        fna_gz=$(find "${tmpdir}/extracted" -name "*.fna.gz" | head -1)
        if [[ -n "$fna_gz" ]]; then
            gunzip -c "$fna_gz" > "$outfile"
            log_info "    [OK] Saved via wget (from .fna.gz): $(basename "$outfile")  ($(du -sh "$outfile" | cut -f1))"
            rm -rf "$tmpdir"
            return 0
        fi
    fi

    if [[ -n "$fna_file" && -f "$fna_file" ]]; then
        mv "$fna_file" "$outfile"
        log_info "    [OK] Saved via wget: $(basename "$outfile")  ($(du -sh "$outfile" | cut -f1))"
        rm -rf "$tmpdir"
        return 0
    else
        log_error "    No .fna in wget download for $accession"
        rm -rf "$tmpdir"
        return 1
    fi
}

# ---------------------------------------------------------------------------
# Step 0, validate MAG inputs
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Step 1, download reference genomes
# ---------------------------------------------------------------------------
log_info "Step 1: Downloading reference genomes (${#ALL_REFS[@]} accessions)"
mkdir -p "${REF_SOURCE}"

if command -v datasets &>/dev/null; then
    log_info "Download tool: ncbi-datasets-cli ($(datasets version 2>/dev/null | head -1))"
else
    log_info "Download tool: wget + NCBI REST API  (datasets CLI not found in PATH)"
fi

n_ok=0
n_fail=0

log_info "--- Thalassarchaeaceae (Poseidoniales / MGIIb-Q1 close refs) ---"
for acc in "${REFS_THALASSO[@]}"; do
    if download_genome "$acc" "$REF_SOURCE"; then ((++n_ok)); else ((++n_fail)); fi
done

log_info "--- Thalassarchaeaceae (UBA36 close refs) ---"
for acc in "${REFS_UBA36[@]}"; do
    if download_genome "$acc" "$REF_SOURCE"; then ((++n_ok)); else ((++n_fail)); fi
done

log_info "--- Nitrosopumilaceae (Nitrosopelagicus close refs) ---"
for acc in "${REFS_NITROSOPELAGICUS[@]}"; do
    if download_genome "$acc" "$REF_SOURCE"; then ((++n_ok)); else ((++n_fail)); fi
done

log_info "--- Broader archaeal context (1-2 per phylum) ---"
for acc in "${REFS_BROAD_ARCHAEA[@]}"; do
    if download_genome "$acc" "$REF_SOURCE"; then ((++n_ok)); else ((++n_fail)); fi
done

log_info "--- Bacterial outgroups ---"
for acc in "${REFS_OUTGROUPS[@]}"; do
    if download_genome "$acc" "$REF_SOURCE"; then ((++n_ok)); else ((++n_fail)); fi
done

log_info "Download summary:"
log_info "  Downloaded / already present: $n_ok"
log_info "  Failed:                       $n_fail"
[[ $n_fail -eq 0 ]] || log_warn "Some genomes failed to download, check internet access and accessions."

ref_count=$(find "$REF_SOURCE" -maxdepth 1 -name "*.fna" | wc -l)
log_info "Reference genomes in $REF_SOURCE: $ref_count"
[[ $ref_count -gt 0 ]] || { log_error "No reference genomes available. Aborting."; exit 1; }

# Write reference taxonomy mapping for use by script 20 (tree label renaming)
# Format: tree_name (REF_<accession>) TAB official_species_name
REF_TAXONOMY_FILE="${WORK_DIR}/reference_taxonomy.tsv"
{
    echo -e "tree_name\tspecies"
    for acc in "${ALL_REFS[@]}"; do
        species="${REF_SPECIES[$acc]:-REF_${acc}}"
        echo -e "REF_${acc}\t${species}"
    done
} > "$REF_TAXONOMY_FILE"
log_info "Reference taxonomy file written: $REF_TAXONOMY_FILE"

# ---------------------------------------------------------------------------
# Step 2, build combined input directory
# ---------------------------------------------------------------------------
log_info "Step 2: Building combined input directory"
mkdir -p "${COMBINED_DIR}"

# Symlink 10 target MAGs (.fa)
log_info "Symlinking 10 target MAGs..."
for mag in "${TARGET_MAGS[@]}"; do
    src="${MAG_SOURCE}/${mag}.fa"
    dst="${COMBINED_DIR}/${mag}.fa"
    if [[ ! -L "$dst" ]]; then
        ln -s "$src" "$dst"
        log_debug "  Linked: $mag.fa"
    else
        log_debug "  Already exists: $mag.fa"
    fi
done

# Symlink reference genomes, ONLY accessions in ALL_REFS, prefixed with REF_ to
# avoid GTDB-Tk ID collision (genomes are keyed by filename stem; GCA/GCF accessions
# match internal DB entries and cause "same id" errors without the prefix).
#
# IMPORTANT: iterate over ALL_REFS explicitly instead of globbing REF_SOURCE/*.fna.
# Globbing would include stale .fna files from previous runs (removed accessions,
# NCBI-suppressed genomes, duplicate GCA/GCF versions) and add them to the tree.
# Clearing old symlinks first ensures COMBINED_DIR contains exactly the current list.
log_info "Symlinking reference genomes (prefix REF_, normalise to .fa)..."
find "${COMBINED_DIR}" -maxdepth 1 -name "REF_*.fa" -type l -delete 2>/dev/null || true
log_info "  Stale REF_ symlinks cleared"

linked_refs=0
for acc in "${ALL_REFS[@]}"; do
    ref_file="${REF_SOURCE}/${acc}.fna"
    if [[ ! -f "$ref_file" ]]; then
        log_error "  Reference genome not found after download step: $ref_file"
        exit 1
    fi
    dst="${COMBINED_DIR}/REF_${acc}.fa"
    ln -sf "$ref_file" "$dst"
    ((++linked_refs))
    log_debug "  Linked: ${acc}.fna to REF_${acc}.fa"
done
log_info "  Linked $linked_refs reference genomes"

total_genomes=$(find "${COMBINED_DIR}" -maxdepth 1 -name "*.fa" | wc -l)
log_info "Total genomes in combined dir: $total_genomes (10 MAGs + ${linked_refs} references)"

# ---------------------------------------------------------------------------
# Step 3, GTDB-Tk identify
# ---------------------------------------------------------------------------
log_info "Step 3: GTDB-Tk identify (ar53 markers)"
mkdir -p "${IDENTIFY_DIR}"

log_info "Running gtdbtk identify..."
log_info "  Input:  ${COMBINED_DIR}"
log_info "  Output: ${IDENTIFY_DIR}"
log_info "  CPUs:   ${SLURM_CPUS_PER_TASK}"

gtdbtk identify \
    --genome_dir "${COMBINED_DIR}" \
    --out_dir "${IDENTIFY_DIR}" \
    --extension "fa" \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --prefix "${PREFIX}" \
    --tmpdir "${TMPDIR}"

log_info "GTDB-Tk identify completed."

# ---------------------------------------------------------------------------
# Step 4, GTDB-Tk align
# ---------------------------------------------------------------------------
log_info "Step 4: GTDB-Tk align (ar53 markers)"
mkdir -p "${ALIGN_DIR}"

log_info "Running gtdbtk align..."
log_info "  Identify dir: ${IDENTIFY_DIR}"
log_info "  Output:       ${ALIGN_DIR}"

gtdbtk align \
    --identify_dir "${IDENTIFY_DIR}" \
    --out_dir "${ALIGN_DIR}" \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --prefix "${PREFIX}" \
    --tmpdir "${TMPDIR}"

log_info "GTDB-Tk align completed."

# Report output files
MSA_CONCAT="${ALIGN_DIR}/${PREFIX}.ar53.msa.fasta.gz"
MARKERS_DIR="${ALIGN_DIR}/intermediate_results/markers"

if [[ -f "$MSA_CONCAT" ]]; then
    ntaxa=$(zcat "$MSA_CONCAT" | grep -c "^>" || true)
    log_info "  Concatenated MSA: $MSA_CONCAT ($ntaxa taxa)"
fi

if [[ -d "$MARKERS_DIR" ]]; then
    nmarkers=$(find "$MARKERS_DIR" -name "*.faa" | wc -l)
    log_info "  Per-marker files: $MARKERS_DIR/ ($nmarkers .faa files)"
    log_info "  These will be used by script 19 for trimAl + FASconCAT-G"
fi

# ---------------------------------------------------------------------------
# Cleanup temp
# ---------------------------------------------------------------------------
rm -rf "${TMPDIR}"
log_info "Temp dir removed: ${TMPDIR}"
log_info "Pipeline completed successfully"
log_info "  Output structure:"
log_info "    ${WORK_DIR}/"
log_info "    ├── 00_input/    (symlinked genomes)"
log_info "    ├── 01_identify/ (GTDB-Tk identify output)"
log_info "    └── 02_align/    (GTDB-Tk align output)"
log_info "  Per-marker HMM-aligned files:"
log_info "    ${MARKERS_DIR}/*.faa"
log_info "  Next step: sbatch 19_10mag_phylo_trimal_fasconcat.slurm"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
