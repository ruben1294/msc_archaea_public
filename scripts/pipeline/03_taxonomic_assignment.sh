#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 3. Classifies the filtered MAGs against GTDB with gtdbtk classify_wf and
# writes the summary tables under results/03_taxonomy_assigned/. The DRAM2
# config, the tree scripts and the R scripts all read
# classify/tax_.ar53.summary.tsv from there, so the prefix defaults to "tax_",
# the one in that file name; change it and they stop finding the table.
#
# The thesis ran this with GTDB-Tk 2.5.2 against GTDB R226.
#
# min_af 0.5 is the GTDB-Tk default. min_perc_aa stays at 10 because several
# DPANN and MGIII bins are fragmented enough to fall under the stricter value
# and would otherwise come back unclassified.

VERSION="1.0.0"

print_help() {
    cat <<EOF
Usage: bash $0 [options]

Assign GTDB taxonomy to MAGs with GTDB-Tk.

Options:
  -h, --help                Show this help and exit
  -v, --version             Show the version and exit
  -p, --project-dir DIR     Project root (default: PROJECT_DIR from server.env)
  -r, --raw-dir DIR         Input MAGs (default: \$PROJECT_DIR/results/01_mag_quality/filtered)
  -w, --workdir DIR         Results root (default: \$PROJECT_DIR/results)
  -x, --extension EXT       MAG file extension (default: fa)
      --min-perc-aa N       Minimum percent of amino acids (default: 10)
      --prefix PREFIX       Output prefix (default: tax_, what every consumer reads)
      --min-af FLOAT        Minimum alignment fraction for pplacer (default: 0.5)
  -c, --cpu-cores N         CPU cores (default: \$SLURM_CPUS_PER_TASK, or 1)
      --logfile FILE        Log file (default: computed from workdir and timestamp)
      --ts TIMESTAMP        Timestamp used in output names

ENVIRONMENT:
  CONFIG_FILE               Alternative path to server.env
  GTDBTK_DATA_PATH          GTDB reference data; set from GTDBTK_DB in server.env

Example:
  bash $0 --raw-dir ./mags --extension fna -c 8
EOF
}

for arg in "$@"; do
    case $arg in
        -h|--help)    print_help; exit 0 ;;
        -v|--version) echo "03_taxonomic_assignment.sh ${VERSION}"; exit 0 ;;
    esac
done

load_project_config
require_commands gtdbtk

EXTENSION="${EXTENSION:-fa}"
MIN_PERC_AA="${MIN_PERC_AA:-10}"
PREFIX="${PREFIX:-tax_}"
MIN_AF="${MIN_AF:-0.5}"
CPU_CORES="${SLURM_CPUS_PER_TASK:-1}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -p|--project-dir) OPT_PROJECT_DIR="$2"; shift 2 ;;
        -r|--raw-dir)     RAW_DIR="$2"; shift 2 ;;
        -w|--workdir)     WORKDIR="$2"; shift 2 ;;
        -x|--extension)   EXTENSION="$2"; shift 2 ;;
        --min-perc-aa)    MIN_PERC_AA="$2"; shift 2 ;;
        --prefix)         PREFIX="$2"; shift 2 ;;
        --min-af)         MIN_AF="$2"; shift 2 ;;
        -c|--cpu-cores)   CPU_CORES="$2"; shift 2 ;;
        --logfile)        LOGFILE="$2"; shift 2 ;;
        --ts)             TS="$2"; shift 2 ;;
        --)               shift; break ;;
        *) echo "Unknown option: $1" >&2; print_help; exit 1 ;;
    esac
done

# --project-dir moves every default derived from it; --workdir and --raw-dir
# then override one each.
if [[ -n "${OPT_PROJECT_DIR:-}" ]]; then
    PROJECT_DIR="$OPT_PROJECT_DIR"
    RESULTS_DIR="$PROJECT_DIR/results"
fi
PROJECT_DIR="${PROJECT_DIR:?PROJECT_DIR is unset; check $CONFIG_FILE}"
WORKDIR="${WORKDIR:-${RESULTS_DIR:-$PROJECT_DIR/results}}"
# Classification runs on the MIMAG-filtered, renamed set, not on data/mags.
RAW_DIR="${RAW_DIR:-$WORKDIR/${FILTERED_SUBDIR:-01_mag_quality/filtered}}"
TEMP_DIR="${TEMP_DIR:-$PROJECT_DIR/${TEMP_SUBDIR:-temp}}"
TS="${TS:-$(date +%Y%m%d_%H%M%S)}"
LOGDIR="${LOGDIR:-$PROJECT_DIR/${LOGS_SUBDIR:-logs}}"

# GTDB-Tk reads its reference data from GTDBTK_DATA_PATH.
if [[ -n "${GTDBTK_DB:-}" ]]; then
    export GTDBTK_DATA_PATH="$GTDBTK_DB"
fi
mkdir -p "$LOGDIR"
LOGFILE="${LOGFILE:-$LOGDIR/03_taxonomic_assignment_${TS}.log}"

exec > >(tee -a "$LOGFILE") 2>&1
echo "== GTDB-Tk taxonomic assignment =="
echo "Log file:     $LOGFILE"
echo "Host:         $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID:-NA}"
echo "CPUs:         $CPU_CORES"
echo "Results root: $WORKDIR"
echo "Input MAGs:   $RAW_DIR"

D03_TAXONOMY="$WORKDIR/${TAXONOMY_SUBDIR:-03_taxonomy_assigned}"
mkdir -p "$D03_TAXONOMY" "$TEMP_DIR"

if [ ! -d "$RAW_DIR" ]; then
    echo "Input directory does not exist: $RAW_DIR" >&2
    exit 1
fi

gtdbtk classify_wf \
    --genome_dir "$RAW_DIR" \
    --out_dir "$D03_TAXONOMY" \
    --extension "$EXTENSION" \
    --min_perc_aa "$MIN_PERC_AA" \
    --prefix "$PREFIX" \
    --cpus "$CPU_CORES" \
    --pplacer_cpus "$CPU_CORES" \
    --force \
    --min_af "$MIN_AF" \
    --tmpdir "$TEMP_DIR" \
    --debug

echo "Summary tables: $D03_TAXONOMY/classify/"
