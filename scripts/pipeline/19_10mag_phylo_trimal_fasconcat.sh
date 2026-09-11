#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 19. Rebuilds per-marker alignments from the GTDB-Tk identify output of
# step 18, then trims and concatenates them into the supermatrix and partition
# file that step 20 hands to IQ-TREE.
#
# Four passes:
#   1. Pull the raw protein sequence for each ar53 marker out of the identify
#      intermediates: _pfam_tophit.tsv, _tigrfam_tophit.tsv and _protein.faa.
#   2. Align each marker on its own with MAFFT --localpair --maxiterate 1000
#      (L-INS-i), the setting used in the published archaeal phylogenomics this
#      follows: Dombrowski et al. 2020 Nat Commun 10.1038/s41467-020-16835-7,
#      Seitz et al. 2019 Nat Commun 10.1038/s41467-019-09364-x, Spang et al.
#      2018 Nat Microbiol 10.1038/s41564-018-0176-9, Farag et al. 2023 ISME J
#      10.1038/s41396-023-01501-1.
#   3. Trim each alignment with trimAl -automated1.
#   4. Concatenate with FASconCAT-G (-s -l), which also writes the partitions.
#
# Why rebuild rather than reuse what GTDB-Tk align already produced: that step
# emits one concatenated MSA with column masking applied, 8,062 columns against
# 13,540 summed raw HMM lengths (Parks et al. 2022 Nat Biotechnol). The masking
# is not recorded per marker, so the concatenated MSA cannot be split back into
# per-marker alignments. Starting from the identify intermediates, which are
# pre-alignment and pre-masking, avoids the problem entirely.
#
# Input:  18_10mag_phylogeny/01_identify/identify/intermediate_results/
#           marker_genes/{genome_id}/
# Output: 03_extracted_markers/  raw unaligned per-marker .faa
#         03_mafft_aligned/      MAFFT-aligned per-marker .faa
#         03_trimal/             trimmed per-marker .fa and stats
#         04_fasconcat/          supermatrix.fas and partition.txt

trap 'log_error "Script failed at line $LINENO with exit code $?"; exit 1' ERR

log_info "10-MAG phylogenomics, per-marker MAFFT L-INS-i + trimAl + FASconCAT-G"
log_info "SLURM Job Information"
log_info "  Job ID:    ${SLURM_JOB_ID:-N/A}"
log_info "  Job Name:  ${SLURM_JOB_NAME:-interactive}"
log_info "  Node:      ${SLURMD_NODENAME:-$(hostname)}"
log_info "  Start:     $(date '+%Y-%m-%d %H:%M:%S')"
# ---- configuration ----
load_project_config

WORK_DIR="${PROJECT_DIR}/18_10mag_phylogeny"
PREFIX="10mag_refs"

# Identify intermediates (input, from script 18)
MARKERS_DIR="${WORK_DIR}/01_identify/identify/intermediate_results/marker_genes"

# Output directories
EXTRACT_DIR="${WORK_DIR}/03_extracted_markers"  # raw (unaligned) per-marker .faa
MAFFT_DIR="${WORK_DIR}/03_mafft_aligned"        # MAFFT L-INS-i aligned per-marker .faa
TRIMAL_DIR="${WORK_DIR}/03_trimal"              # trimAl-trimmed per-marker .fa
FASCONCAT_DIR="${WORK_DIR}/04_fasconcat"        # supermatrix + partition file

# Minimum number of genomes that must have a marker to include it in the tree
# (with 22 taxa total, 5 is a conservative lower bound)
MIN_SEQUENCES=5

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# ---- tool resolution ----

# MAFFT
# MAFFT_PATH in server.env names a specific build; otherwise take what is on PATH.
MAFFT_LINK="${MAFFT_PATH:-}"
if [[ -x "${MAFFT_LINK}" ]]; then
    MAFFT_EXEC="$MAFFT_LINK"
elif command -v mafft >/dev/null 2>&1; then
    MAFFT_EXEC="$(command -v mafft)"
else
    log_error "MAFFT not found. Set MAFFT_PATH in $CONFIG_FILE or put mafft on PATH."
    exit 1
fi
log_info "MAFFT:       $MAFFT_EXEC  ($("$MAFFT_EXEC" --version 2>&1 | head -1 || true))"

# trimAl
TRIMAL_LINK="${TRIMAL_PATH:-}"
if [[ -x "${TRIMAL_LINK}" ]]; then
    TRIMAL_EXEC="$TRIMAL_LINK"
elif command -v trimal >/dev/null 2>&1; then
    TRIMAL_EXEC="$(command -v trimal)"
else
    log_error "trimAl not found. Set TRIMAL_PATH in $CONFIG_FILE or put trimal on PATH."
    exit 1
fi
log_info "trimAl:      $TRIMAL_EXEC  ($(${TRIMAL_EXEC} --version 2>&1 | head -1 || true))"

# FASconCAT-G
# utils/install_FASconCAT-G.sh writes this file next to the installed script.
FASCONCAT_CONFIG="${PROJECT_DIR}/${MARKER_CONCAT_SUBDIR:-12_marker_concatenation}/fasconcat_config.env"
if [[ -z "${FASCONCAT_EXEC:-}" && -f "$FASCONCAT_CONFIG" ]]; then
    # shellcheck source=/dev/null
    source "$FASCONCAT_CONFIG"
fi
[[ -f "${FASCONCAT_EXEC:-}" ]] || {
    log_error "FASconCAT-G not found. Run utils/install_FASconCAT-G.sh first."
    exit 1
}
log_info "FASconCAT-G: $FASCONCAT_EXEC"

log_info "Threads:     $THREADS"

# ---- Step 0, validate inputs ----
log_info "Step 0: Validating inputs"
[[ -d "$MARKERS_DIR" ]] || {
    log_error "GTDB-Tk identify intermediates not found: $MARKERS_DIR"
    log_error "Run script 18 (gtdbtk identify) first."
    exit 1
}

n_genomes=$(find "$MARKERS_DIR" -mindepth 1 -maxdepth 1 -type d | wc -l)
log_info "  [OK] Identify intermediates: $MARKERS_DIR"
log_info "  [OK] Genome directories found: $n_genomes"
[[ "$n_genomes" -gt 0 ]] || { log_error "No genome directories in $MARKERS_DIR"; exit 1; }

mkdir -p "$EXTRACT_DIR" "$MAFFT_DIR"

# ---- Step 1, extract raw per-marker protein sequences from identify intermediates ----
log_info "Step 1: Extracting per-marker sequences from GTDB-Tk identify intermediates"
log_info "  Input (markers dir): $MARKERS_DIR"
log_info "  Output (raw .faa):   $EXTRACT_DIR"
log_info "  Min genomes/marker:  $MIN_SEQUENCES"

# Python prints informational messages to stderr (visible in logs),
# and only the integer count to stdout (captured by bash).
n_markers_extracted=$(python3 - "$MARKERS_DIR" "$EXTRACT_DIR" "$MIN_SEQUENCES" <<'PYEOF'
import sys, os
from collections import defaultdict

markers_dir = sys.argv[1]
output_dir  = sys.argv[2]
min_seqs    = int(sys.argv[3])

os.makedirs(output_dir, exist_ok=True)

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, flush=True, **kwargs)

# marker_id to {genome_id: (sequence, bitscore)}
# If multiple genes in the same genome match the same marker, keep the best score.
marker_seqs = defaultdict(dict)

genome_dirs = sorted(
    entry.path for entry in os.scandir(markers_dir) if entry.is_dir()
)
eprint(f"  Genomes found: {len(genome_dirs)}")

for genome_dir in genome_dirs:
    genome_id = os.path.basename(genome_dir)

    # ---- Read protein sequences indexed by gene_id ----
    prot_file = os.path.join(genome_dir, f"{genome_id}_protein.faa")
    if not os.path.exists(prot_file):
        eprint(f"  WARN: protein.faa missing for {genome_id}, skipping")
        continue

    proteins = {}
    cur_id, cur_seq = None, []
    with open(prot_file) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if cur_id:
                    proteins[cur_id] = ''.join(cur_seq)
                cur_id = line[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id:
            proteins[cur_id] = ''.join(cur_seq)

    def parse_tophit(filepath, strip_pfam_version=False):
        """
        Parse a GTDB-Tk _pfam_tophit.tsv or _tigrfam_tophit.tsv file.

        File format (tab-separated):
            Gene Id    Top hits (Family id,e-value,bitscore)
            k141_123   PF00466.21,7.5e-12,41.9

        Returns list of (gene_id, marker_id, bitscore).
        """
        hits = []
        if not os.path.exists(filepath):
            return hits
        with open(filepath) as fh:
            next(fh)  # skip header line
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                gene_id   = parts[0]
                hit_field = parts[1]  # e.g. "PF00466.21,7.5e-12,41.9"
                hit_parts = hit_field.split(',')
                if len(hit_parts) < 3:
                    continue
                marker_raw = hit_parts[0]
                try:
                    bitscore = float(hit_parts[2])
                except ValueError:
                    bitscore = 0.0
                # Strip Pfam version number: PF00466.21 to PF00466
                if strip_pfam_version and '.' in marker_raw:
                    marker_id = marker_raw.split('.')[0]
                else:
                    marker_id = marker_raw
                hits.append((gene_id, marker_id, bitscore))
        return hits

    pfam_hits = parse_tophit(
        os.path.join(genome_dir, f"{genome_id}_pfam_tophit.tsv"),
        strip_pfam_version=True
    )
    tigr_hits = parse_tophit(
        os.path.join(genome_dir, f"{genome_id}_tigrfam_tophit.tsv"),
        strip_pfam_version=False
    )

    for gene_id, marker_id, bitscore in pfam_hits + tigr_hits:
        if gene_id not in proteins:
            continue
        # Keep best-scoring gene when multiple genes hit the same marker in one genome
        if genome_id not in marker_seqs[marker_id] or \
                bitscore > marker_seqs[marker_id][genome_id][1]:
            marker_seqs[marker_id][genome_id] = (proteins[gene_id], bitscore)

eprint(f"  Unique markers detected across all genomes: {len(marker_seqs)}")

# ---- Write per-marker .faa files ----
n_written = n_skipped = 0
for marker_id in sorted(marker_seqs):
    genome_hits = marker_seqs[marker_id]
    if len(genome_hits) < min_seqs:
        eprint(f"  SKIP {marker_id}: {len(genome_hits)} genomes (< {min_seqs})")
        n_skipped += 1
        continue
    out_file = os.path.join(output_dir, f"{marker_id}.faa")
    with open(out_file, 'w') as out:
        for genome_id in sorted(genome_hits):
            seq = genome_hits[genome_id][0]
            out.write(f">{genome_id}\n{seq}\n")
    n_written += 1

eprint(f"  Written: {n_written} marker files  |  Skipped (< {min_seqs} genomes): {n_skipped}")
print(n_written)   # ← only this integer goes to stdout to captured by bash $()
PYEOF
)

# Guard against unexpected Python output
if [[ ! "$n_markers_extracted" =~ ^[0-9]+$ ]]; then
    log_error "Unexpected output from extraction step: '$n_markers_extracted'"
    exit 1
fi

log_info "  Per-marker files extracted: $n_markers_extracted"
[[ "$n_markers_extracted" -gt 0 ]] || { log_error "No marker files produced. Aborting."; exit 1; }

# Step 2, MAFFT L-INS-i per-marker alignment
#
# Method: --localpair --maxiterate 1000 (L-INS-i)
# Justification: standard approach in published de novo archaeal phylogenomics
#   Dombrowski et al. 2020 Nat Commun doi:10.1038/s41467-020-16835-7
#   Seitz et al. 2019 Nat Commun doi:10.1038/s41467-019-09364-x
#   Spang et al. 2018 Nat Microbiol doi:10.1038/s41564-018-0176-9
#   Farag et al. 2023 ISME J      doi:10.1038/s41396-023-01501-1
log_info "Step 2: MAFFT L-INS-i alignment (${n_markers_extracted} markers)"
mafft_ok=0
mafft_fail=0

for marker_file in "${EXTRACT_DIR}"/*.faa; do
    [[ -f "$marker_file" ]] || continue
    marker_name=$(basename "$marker_file" .faa)
    out_file="${MAFFT_DIR}/${marker_name}_aligned.faa"

    seq_count=$(grep -c "^>" "$marker_file" || echo "0")

    if [[ "$seq_count" -lt 2 ]]; then
        # Only one sequence: copy as-is (nothing to align)
        cp "$marker_file" "$out_file"
        log_debug "  COPY  $marker_name (${seq_count} seq, no alignment needed)"
        ((++mafft_ok))
        continue
    fi

    start_t=$(date +%s)
    if "$MAFFT_EXEC" \
            --localpair \
            --maxiterate 1000 \
            --thread "$THREADS" \
            "$marker_file" > "$out_file" 2>/dev/null; then
        elapsed=$(( $(date +%s) - start_t ))
        log_debug "  [OK] $marker_name  (${seq_count} seqs, ${elapsed}s)"
        ((++mafft_ok))
    else
        log_error "  [FAIL] MAFFT failed for: $marker_name"
        ((++mafft_fail))
    fi
done

log_info "  MAFFT summary: $mafft_ok OK, $mafft_fail failed"
[[ $mafft_fail -eq 0 ]] || { log_error "Some MAFFT alignments failed."; exit 1; }
[[ $mafft_ok   -gt 0 ]] || { log_error "No MAFFT alignments produced.";  exit 1; }

# ---- Step 3, trimAl -automated1 on each MAFFT-aligned marker ----
log_info "Step 3: trimAl -automated1 (${mafft_ok} alignments)"
mkdir -p "${TRIMAL_DIR}/trimmed"
mkdir -p "${TRIMAL_DIR}/reports"

DATE_TAG=$(date '+%Y%m%d')
STATS_FILE="${TRIMAL_DIR}/${DATE_TAG}_trim_statistics.tsv"
echo -e "marker\tseqs_before\tseqs_after\tlen_before\tlen_after\tpct_retained" > "$STATS_FILE"

total_ok=0
total_fail=0

for marker_file in "${MAFFT_DIR}"/*_aligned.faa; do
    [[ -f "$marker_file" ]] || continue
    marker_name=$(basename "$marker_file" _aligned.faa)
    out_file="${TRIMAL_DIR}/trimmed/${marker_name}_trimmed.fa"
    html_file="${TRIMAL_DIR}/reports/${marker_name}_trimal.html"

    seqs_before=$(grep -c "^>" "$marker_file" || echo "0")
    len_before=$(awk '/^>/{next}{print length($0); exit}' "$marker_file")

    if "$TRIMAL_EXEC" \
            -in  "$marker_file" \
            -out "$out_file" \
            -fasta \
            -automated1 \
            -htmlout "$html_file" 2>/dev/null; then

        seqs_after=$(grep -c "^>" "$out_file" || echo "0")
        len_after=$(awk '/^>/{next}{print length($0); exit}' "$out_file")
        pct=$(awk "BEGIN{printf \"%.1f\", ($len_before > 0 ? ($len_after / $len_before) * 100 : 0)}")
        echo -e "${marker_name}\t${seqs_before}\t${seqs_after}\t${len_before}\t${len_after}\t${pct}" >> "$STATS_FILE"
        log_debug "  [OK] $marker_name: ${len_before} to ${len_after} cols (${pct}% retained)"
        ((++total_ok))
    else
        log_error "  [FAIL] trimAl failed for: $marker_name"
        ((++total_fail))
    fi
done

log_info "trimAl summary:"
log_info "  OK:     $total_ok markers"
log_info "  Failed: $total_fail markers"
log_info "  Stats:  $STATS_FILE"
[[ $total_fail -eq 0 ]] || { log_error "Some markers failed trimming."; exit 1; }
[[ $total_ok   -gt 0 ]] || { log_error "No markers trimmed. Aborting.";  exit 1; }

# ---- Step 4, FASconCAT-G: concatenate trimmed markers + generate partition file ----
log_info "Step 4: FASconCAT-G concatenation (-s -l)"
mkdir -p "${FASCONCAT_DIR}/input"
mkdir -p "${FASCONCAT_DIR}/output"
mkdir -p "${FASCONCAT_DIR}/logs"

FASCONCAT_LOG="${FASCONCAT_DIR}/logs/fasconcat.log"

# FASconCAT-G requires all input .fas files in its working directory.
# Copy trimmed markers there.
log_info "Staging trimmed markers in FASconCAT-G output dir..."
n_copied=0
for trimmed_file in "${TRIMAL_DIR}/trimmed"/*_trimmed.fa; do
    [[ -f "$trimmed_file" ]] || continue
    base=$(basename "$trimmed_file" _trimmed.fa)
    cp "$trimmed_file" "${FASCONCAT_DIR}/output/${base}.fas"
    ((++n_copied))
done
log_info "  $n_copied marker files staged"
[[ $n_copied -gt 0 ]] || { log_error "No trimmed markers to concatenate. Aborting."; exit 1; }

log_info "Running FASconCAT-G..."
cd "${FASCONCAT_DIR}/output"

if LC_ALL=C perl "$FASCONCAT_EXEC" -s -l > "$FASCONCAT_LOG" 2>&1; then
    log_info "  [OK] FASconCAT-G completed"

    # Globs rather than ls: the partition file also starts with
    # FcC_supermatrix, so the .fas suffix is what separates the two.
    shopt -s nullglob
    supermatrix_candidates=(FcC_supermatrix*.fas)
    partition_candidates=(FcC_supermatrix*_partition.txt)
    shopt -u nullglob
    supermatrix_raw="${supermatrix_candidates[0]:-}"
    partition_raw="${partition_candidates[0]:-}"

    if [[ -n "$supermatrix_raw" && -f "$supermatrix_raw" ]]; then
        ntaxa=$(grep -c "^>" "$supermatrix_raw")
        ncols=$(awk '/^>/{next}{print length($0); exit}' "$supermatrix_raw")
        mv "$supermatrix_raw" "${PREFIX}_supermatrix.fas"
        log_info "  Supermatrix: ${PREFIX}_supermatrix.fas  ($ntaxa taxa × $ncols positions)"
    else
        log_error "  Supermatrix file not found. Check: $FASCONCAT_LOG"
        exit 1
    fi

    if [[ -n "$partition_raw" && -f "$partition_raw" ]]; then
        n_parts=$(grep -c "." "$partition_raw" || echo "0")
        mv "$partition_raw" "${PREFIX}_partition.txt"
        log_info "  Partition file: ${PREFIX}_partition.txt  ($n_parts partitions)"
    else
        log_warn "  Partition file not found (FASconCAT-G may not have generated it)"
    fi

    # Remove staged input .fas copies (supermatrix and partition already renamed)
    find "${FASCONCAT_DIR}/output" -maxdepth 1 -name "*.fas" \
        ! -name "${PREFIX}_supermatrix.fas" -delete 2>/dev/null || true

else
    log_error "  [FAIL] FASconCAT-G failed. Check: $FASCONCAT_LOG"
    tail -30 "$FASCONCAT_LOG" | while IFS= read -r line; do log_error "    $line"; done
    exit 1
fi

# ---- summary ----
log_info "Pipeline completed successfully"
log_info "  Output structure:"
log_info "    ${WORK_DIR}/"
log_info "    ├── 03_extracted_markers/  (raw unaligned per-marker .faa)"
log_info "    ├── 03_mafft_aligned/      (MAFFT L-INS-i aligned .faa)"
log_info "    ├── 03_trimal/"
log_info "    │   ├── trimmed/           (*_trimmed.fa)"
log_info "    │   ├── reports/           (*_trimal.html)"
log_info "    │   └── ${DATE_TAG}_trim_statistics.tsv"
log_info "    └── 04_fasconcat/"
log_info "        ├── output/            (${PREFIX}_supermatrix.fas + ${PREFIX}_partition.txt)"
log_info "        └── logs/"
log_info ""
log_info "  Supermatrix (input for IQ-TREE2, script 20):"
log_info "    ${FASCONCAT_DIR}/output/${PREFIX}_supermatrix.fas"
log_info "  Partition file (for partitioned tree, script 20 Step C):"
log_info "    ${FASCONCAT_DIR}/output/${PREFIX}_partition.txt"
log_info ""
log_info "  Next step: sbatch 20_10mag_phylo_iqtree.slurm"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
