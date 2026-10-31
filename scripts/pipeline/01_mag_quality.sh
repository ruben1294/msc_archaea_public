#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 1. Runs CheckM2 over the input MAGs in data/mags/ and writes
# quality_report.tsv under results/01_mag_quality/. Before that it strips
# duplicate contig IDs, which CheckM2 rejects and which show up when
# MAGs binned from separate assemblies of the same station contain contigs with
# the same name.
#
# Two helpers finish the quality block before step 02:
#   utils/01_filter_mags.py  keeps the MAGs of MIMAG medium quality or better,
#                            copying them into results/01_mag_quality/filtered/
#   utils/05_rename_mags.sh  gives those copies their final dotted names, which
#                            every later step and data/MAGs_metadata.tsv use

VERSION="1.0.0"

print_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Assess the quality of metagenome-assembled genomes with CheckM2.

OPTIONS:
    -c, --cpu-cores NUM          CPU cores (default: \$SLURM_CPUS_PER_TASK, or 4)
    -p, --project-dir PATH       Project root (default: PROJECT_DIR from server.env)
    -w, --workdir PATH           Results root (default: \$PROJECT_DIR/results)
    -r, --raw-dir PATH           Input MAGs (default: \$PROJECT_DIR/data/mags)
    -t, --temp-dir PATH          Scratch space (default: \$PROJECT_DIR/temp)
    -l, --logdir PATH            Log directory (default: \$PROJECT_DIR/logs)
    -h, --help                   Show this message and exit
    -v, --version                Show the version and exit

ENVIRONMENT:
    CONFIG_FILE                  Alternative path to server.env
    SLURM_CPUS_PER_TASK          Sets CPU_CORES when --cpu-cores is absent

Any of the path options can also arrive as an environment variable of the same
name in upper case, which is how the SLURM wrappers in slurm/ pass them.

EXAMPLES:
    bash $0
    bash $0 --cpu-cores 8 --project-dir /scratch/archaea
    RAW_DIR=/elsewhere/mags bash $0
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
            -p|--project-dir) OPT_PROJECT_DIR="$2"; shift 2 ;;
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

load_project_config
require_commands checkm2

CPU_CORES="${SLURM_CPUS_PER_TASK:-4}"

parse_arguments "$@"

# Resolved after parsing so that --project-dir moves every default derived from
# it; --workdir and --raw-dir then override one each.
if [[ -n "${OPT_PROJECT_DIR:-}" ]]; then
    PROJECT_DIR="$OPT_PROJECT_DIR"
    RESULTS_DIR="$PROJECT_DIR/results"
    MAGS_DIR="$PROJECT_DIR/data/mags"
fi
PROJECT_DIR="${PROJECT_DIR:?PROJECT_DIR is unset; check $CONFIG_FILE}"
WORKDIR="${WORKDIR:-${RESULTS_DIR:-$PROJECT_DIR/results}}"
RAW_DIR="${RAW_DIR:-${MAGS_DIR:-$PROJECT_DIR/data/mags}}"
TEMP_DIR="${TEMP_DIR:-$PROJECT_DIR/${TEMP_SUBDIR:-temp}}"
TS="${TS:-$(date +%Y%m%d_%H%M%S)}"
LOGDIR="${LOGDIR:-$PROJECT_DIR/${LOGS_SUBDIR:-logs}}"
mkdir -p "$LOGDIR"
LOGFILE="${LOGFILE:-$LOGDIR/01_mag_quality_${TS}.log}"

exec > >(tee -a "$LOGFILE") 2>&1
echo "== CheckM2 quality assessment =="
echo "Log file:     $LOGFILE"
echo "Host:         $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID:-NA}"
echo "CPUs:         $CPU_CORES"
echo "Results root: $WORKDIR"
echo "Input MAGs:   $RAW_DIR"

QUALITY_DIR="$WORKDIR/${QUALITY_SUBDIR:-01_mag_quality}"

for d in "$QUALITY_DIR" "$TEMP_DIR"; do
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

        # The in_duplicate flag persists across the sequence lines that follow a
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
    --output-directory "$QUALITY_DIR" \
    --extension fa \
    --force

echo "Quality report: $QUALITY_DIR/quality_report.tsv"
echo "Next: python3 scripts/utils/01_filter_mags.py, then bash scripts/utils/05_rename_mags.sh"
