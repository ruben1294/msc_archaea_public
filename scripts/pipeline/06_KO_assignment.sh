#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 6 has no automation. KO assignment ran through the KofamKOALA web server
# rather than a local kofam_scan install, so this file documents the manual step
# and checks that its output landed where step 7 expects it.
#
# The procedure:
#
#   1. utils/02_concat_protein_files.sh merges the per-MAG .faa files from step 3
#      into one FASTA, prefixing every header with its MAG name so the origin
#      survives the round trip.
#   2. Upload that FASTA to https://www.genome.jp/tools/kofamkoala/ and request
#      the detail-tsv output. The server splits large submissions, so this comes
#      back as several files.
#   3. Drop the results into $PROJECT_DIR/06_KO_assignment/<date>_kofam_web/.
#   4. utils/03_post_KO_assign_concat_KOs.sh concatenates them, keeps the best
#      hit per gene and writes all_kofam_KOs.transformed.tsv, which is the file
#      07_metabolic_reconstruction.sh reads.
#
# Running this script verifies step 4 produced that file.

load_project_config

KO_DIR="${PROJECT_DIR}/${KO_ASSIGNMENT_SUBDIR}"

if [[ ! -d "$KO_DIR" ]]; then
    echo "No KO assignment directory at $KO_DIR" >&2
    echo "Work through the manual steps in the header of this script first." >&2
    exit 1
fi

# The submission date is part of the directory name, so glob rather than assume.
mapfile -t transformed < <(find "$KO_DIR" -maxdepth 2 -name "all_kofam_KOs.transformed.tsv" -type f | sort)

if [[ ${#transformed[@]} -eq 0 ]]; then
    echo "No all_kofam_KOs.transformed.tsv under $KO_DIR" >&2
    echo "Run utils/03_post_KO_assign_concat_KOs.sh over the KofamKOALA output." >&2
    exit 1
fi

for f in "${transformed[@]}"; do
    n_lines=$(wc -l < "$f")
    echo "[OK] $f ($n_lines lines)"
done

echo "Step 7 can now run."
