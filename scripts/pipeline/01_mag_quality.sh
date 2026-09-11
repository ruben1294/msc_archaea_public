#!/usr/bin/env bash
set -euo pipefail

# Step 1. Runs CheckM2 over the raw MAGs and writes quality_report.tsv under
# 02_mags_quality/. Before that it strips duplicate contig IDs, which CheckM2
# rejects outright and which show up when MAGs binned from separate assemblies
# of the same station carry contigs with the same name.
#
# utils/01_filter_mags.py reads the report afterwards and applies the MIMAG
# thresholds set in configs/server.env.

VERSION="1.0.0"

print_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Assess the quality of metagenome-assembled genomes with CheckM2.

OPTIONS:
    -c, --cpu-cores NUM          CPU cores (default: \$SLURM_CPUS_PER_TASK, or 4)
    -p, --project-dir PATH       Project root (default: PROJECT_DIR from server.env)
    -w, --workdir PATH           Working directory (default: \$PROJECT_DIR)
    -r, --raw-dir PATH           Input MAGs (default: \$PROJECT_DIR/01_raw_files)
    -t, --temp-dir PATH          Scratch space (default: \$WORKDIR/temp)
    -l, --logdir PATH            Log directory (default: \$WORKDIR/logs)
    -h, --help                   Show this message and exit
    -v, --version                Show the version and exit

ENVIRONMENT:
    CONFIG_FILE                  Alternative path to server.env
    SLURM_CPUS_PER_TASK          Sets CPU_CORES when --cpu-cores is absent

Any of the path options can also arrive as an environment variable of the same
name in upper case, which is how the SLURM wrappers in slurm/ pass them.

EXAMPLES:
    bash $0
    bash $0 --cpu-cores 8 --workdir /scratch/archaea
    RAW_DIR=/data/mags bash $0
EOF
}

print_version() {
    echo "01_mag_quality.sh ${VERSION}"
}

parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)        print_help; exit 0 ;;
            -v|--version)     print_version; exit 0 ;;
            -c|--cpu-cores)   CPU_CORES="$2"; shift 2 ;;
            -p|--project-dir) PROJECT_DIR="$2"; shift 2 ;;
            -w|--workdir)     WORKDIR="$2"; shift 2 ;;
            -r|--raw-dir)     RAW_DIR="$2"; shift 2 ;;
            -t|--temp-dir)    TEMP_DIR="$2"; shift 2 ;;
            -l|--logdir)      LOGDIR="$2"; shift 2 ;;
            *)
                echo "Unknown option '$1'. Run with --help for usage." >&2
                exit 1
                ;;
        esac
    done
}

# --help and --version answer before the configuration check, so someone can
# read the usage on a machine that has no server.env yet.
for arg in "$@"; do
    case $arg in
        -h|--help)    print_help; exit 0 ;;
        -v|--version) print_version; exit 0 ;;
    esac
done

# The config sits two levels up from scripts/pipeline/, so the script finds it
# whatever directory it runs from.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CONFIG_FILE="${CONFIG_FILE:-${REPO_ROOT}/configs/server.env}"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "No configuration at $CONFIG_FILE" >&2
    echo "Copy configs/server.env.example to configs/server.env and edit it." >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$CONFIG_FILE"
echo "Loaded configuration from $CONFIG_FILE"

CPU_CORES="${SLURM_CPUS_PER_TASK:-4}"

parse_arguments "$@"

# Resolved after parsing so that --project-dir propagates into the paths derived
# from it.
PROJECT_DIR="${PROJECT_DIR:?PROJECT_DIR is unset; check $CONFIG_FILE}"
WORKDIR="${WORKDIR:-$PROJECT_DIR}"
RAW_DIR="${RAW_DIR:-$PROJECT_DIR/${RAWFILES_SUBDIR:-01_raw_files}}"
TEMP_DIR="${TEMP_DIR:-$WORKDIR/${TEMP_SUBDIR:-temp}}"
TS="${TS:-$(date +%Y%m%d_%H%M%S)}"
LOGDIR="${LOGDIR:-$WORKDIR/${LOGS_SUBDIR:-logs}}"
mkdir -p "$LOGDIR"
LOGFILE="${LOGFILE:-$LOGDIR/01_mag_quality_${TS}.log}"

exec > >(tee -a "$LOGFILE") 2>&1
echo "== CheckM2 quality assessment =="
echo "Log file:     $LOGFILE"
echo "Host:         $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID:-NA}"
echo "CPUs:         $CPU_CORES"
echo "Working dir:  $WORKDIR"
echo "Input MAGs:   $RAW_DIR"

D02_QUALITY="$WORKDIR/${QUALITY_SUBDIR:-02_mags_quality}"
D03_SUMMARY="$WORKDIR/03_summary"

for d in "$D02_QUALITY" "$D03_SUMMARY" "$TEMP_DIR"; do
    mkdir -p "$d"
done

if ! compgen -G "$RAW_DIR/*.fa" > /dev/null; then
    echo "No .fa files in $RAW_DIR" >&2
    exit 1
fi

# ---- duplicate contig IDs ----

DUPLICATES_FOUND=0
for mag in "$RAW_DIR"/*.fa; do
    dups=$(grep "^>" "$mag" | sort | uniq -d)
    if [ -n "$dups" ]; then
        echo "Duplicate contig IDs in $(basename "$mag"):"
        echo "$dups"
        DUPLICATES_FOUND=1
    fi
done

if [ "$DUPLICATES_FOUND" -eq 1 ]; then
    echo ""
    echo "Removing duplicates, keeping the first occurrence of each contig."

    for mag in "$RAW_DIR"/*.fa; do
        mag_basename=$(basename "$mag")
        mag_temp="${TEMP_DIR}/${mag_basename}.tmp"

        # The in_duplicate flag carries across the sequence lines that follow a
        # repeated header, so the whole record is dropped rather than its header
        # alone.
        awk '
          /^>/ {
            contig_id=$1
            if (contig_id in seen) {
              in_duplicate=1
              next
            }
            seen[contig_id]=1
            in_duplicate=0
            print
            next
          }
          !in_duplicate { print }
        ' "$mag" > "$mag_temp"

        original_count=$(grep -c "^>" "$mag")
        deduplicated_count=$(grep -c "^>" "$mag_temp")
        duplicates_removed=$((original_count - deduplicated_count))

        if [ "$duplicates_removed" -gt 0 ]; then
            mv "$mag_temp" "$mag"
            echo "  [OK] $mag_basename: dropped $duplicates_removed contig(s), ${original_count} to ${deduplicated_count}"
        else
            rm -f "$mag_temp"
            echo "  [OK] $mag_basename: nothing to drop"
        fi
    done

    echo "Deduplication finished."
else
    echo "No duplicate contigs. Skipping deduplication."
fi
echo ""

# ---- CheckM2 ----

checkm2 predict \
    --threads "$CPU_CORES" \
    --input "$RAW_DIR" \
    --output-directory "$D02_QUALITY" \
    --extension fa \
    --force

echo "Quality report: $D02_QUALITY/quality_report.tsv"
