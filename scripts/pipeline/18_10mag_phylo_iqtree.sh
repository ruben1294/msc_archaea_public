#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 18. Infers three trees over the supermatrix from step 17, covering the 10
# target MAGs, the archaeal references and the bacterial outgroups.
#
#   A. Guide tree, LG+F+G, with 1000 ultrafast bootstrap and 1000 SH-aLRT
#      replicates. Fast, and its only job is to seed B.
#   B. Final PMSF tree, LG+C60+F+G, using the guide tree through -ft. The
#      posterior mean site frequency approximation makes the C60 profile
#      mixture tractable at this size while keeping its resistance to
#      long-branch attraction, which is the artefact that matters when placing
#      divergent archaea against bacterial outgroups.
#   C. Partitioned tree, ModelFinder per partition (-p -m MFP), as an
#      independent check on B rather than a replacement for it.
#   D. Relabels the tips, from NameWithoutExtension to DisplayName for the MAGs
#      and to species names for the references, so the figure is readable.
#
# Input:  16_10mag_phylogeny/04_fasconcat/output/
# Output: 16_10mag_phylogeny/05_iqtree/

trap 'log_error "Script failed at line $LINENO with exit code $?"; exit 1' ERR

log_info "10-MAG focused phylogenomic tree, IQ-TREE2 PMSF (LG+C60+F+G)"
log_info "SLURM Job Information"
log_info "  Job ID:    ${SLURM_JOB_ID:-N/A}"
log_info "  Job Name:  ${SLURM_JOB_NAME:-interactive}"
log_info "  Node:      ${SLURMD_NODENAME:-$(hostname)}"
log_info "  Start:     $(date '+%Y-%m-%d %H:%M:%S')"
# ---- configuration ----
load_project_config

WORK_DIR="${RESULTS_DIR}/${FOCUSED_PHYLOGENY_SUBDIR:-16_10mag_phylogeny}"
FASCONCAT_DIR="${WORK_DIR}/04_fasconcat"
IQTREE_DIR="${WORK_DIR}/05_iqtree"
PREFIX="10mag_refs"

SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-32}"

# Verify IQ-TREE2 is available
if ! command -v iqtree2 &>/dev/null; then
    log_error "iqtree2 not found in PATH. Load the iqtree conda env."
    exit 1
fi
IQTREE_VERSION=$(iqtree2 --version 2>&1 | head -1)
log_info "IQ-TREE2: $IQTREE_VERSION"

# Input: supermatrix + partition file from FASconCAT-G (script 19)
SUPERMATRIX="${FASCONCAT_DIR}/output/${PREFIX}_supermatrix.fas"
PARTITION_FILE="${FASCONCAT_DIR}/output/${PREFIX}_partition.txt"

# ---- Step 0, validate inputs ----
log_info "Step 0: Validating inputs"
[[ -f "$SUPERMATRIX" ]] || {
    log_error "Supermatrix not found: $SUPERMATRIX"
    log_error "Run script 19 (per-marker trimAl + FASconCAT-G) first."
    exit 1
}
[[ -f "$PARTITION_FILE" ]] || {
    log_error "Partition file not found: $PARTITION_FILE"
    log_error "Run script 19 (per-marker trimAl + FASconCAT-G) first."
    exit 1
}

ntaxa=$(grep -c "^>" "$SUPERMATRIX")
ncols=$(awk '/^>/{next}{print length($0);exit}' "$SUPERMATRIX")
n_parts=$(grep -c "." "$PARTITION_FILE" || echo "0")
log_info "Supermatrix:     $ntaxa taxa × $ncols positions"
log_info "Partition file:  $n_parts partitions"

mkdir -p "${IQTREE_DIR}/guide_tree"
mkdir -p "${IQTREE_DIR}/pmsf_final"
mkdir -p "${IQTREE_DIR}/partition_tree"

# ---- Step A, Guide tree: LG+F+G ----
log_info "Step A: Guide tree (LG+F+G) with UFBoot2 + SH-aLRT"
GUIDE_DIR="${IQTREE_DIR}/guide_tree"

log_info "  Output dir:  $GUIDE_DIR"
log_info "  Model:       LG+F+G"
log_info "  Supports:    UFBoot2 -B 1000 | SH-aLRT --alrt 1000"
log_info "  Threads:     $SLURM_CPUS_PER_TASK"
log_info "  Start:       $(date '+%Y-%m-%d %H:%M:%S')"

cd "$GUIDE_DIR"

iqtree2 \
    -s "$SUPERMATRIX" \
    -m LG+F+G \
    -B 1000 \
    --alrt 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    --seed 42 \
    --prefix 10mags_guide \
    --redo \
    > "${GUIDE_DIR}/10mags_guide.iqtree2.log" 2>&1

GUIDE_TREEFILE="${GUIDE_DIR}/10mags_guide.treefile"

if [[ -f "$GUIDE_TREEFILE" ]]; then
    log_info "  [OK] Guide tree: $GUIDE_TREEFILE"
    logL_guide=$(grep "BEST SCORE FOUND" "${GUIDE_DIR}/10mags_guide.iqtree" | awk '{print $5}' || echo "N/A")
    log_info "  Log-likelihood: $logL_guide"
else
    log_error "  [FAIL] Guide tree not produced. Check: ${GUIDE_DIR}/guide.iqtree2.log"
    exit 1
fi

log_info "  Guide tree done: $(date '+%Y-%m-%d %H:%M:%S')"

# ---- Step B, PMSF final tree: LG+C60+F+G ----
log_info "Step B: PMSF final tree (LG+C60+F+G) with UFBoot2 + SH-aLRT"
log_info "  Using guide tree: $GUIDE_TREEFILE"
log_info "  The -ft flag computes posterior mean site frequencies from the guide"
log_info "  tree under the C60 mixture model (PMSF approach)."
log_info "  Model:    LG+C60+F+G"
log_info "  Supports: UFBoot2 -B 1000 | SH-aLRT --alrt 1000"
log_info "  Threads:  $SLURM_CPUS_PER_TASK"
log_info "  Start:    $(date '+%Y-%m-%d %H:%M:%S')"

PMSF_DIR="${IQTREE_DIR}/pmsf_final"
cd "$PMSF_DIR"

iqtree2 \
    -s "$SUPERMATRIX" \
    -m LG+C60+F+G \
    -ft "$GUIDE_TREEFILE" \
    -B 1000 \
    --alrt 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    --seed 42 \
    --prefix 10mags_pmsf \
    --redo \
    > "${PMSF_DIR}/10mags_pmsf.iqtree2.log" 2>&1

PMSF_TREEFILE="${PMSF_DIR}/10mags_pmsf.treefile"
PMSF_CONTREE="${PMSF_DIR}/10mags_pmsf.contree"

if [[ -f "$PMSF_TREEFILE" ]]; then
    log_info "  [OK] PMSF tree: $PMSF_TREEFILE"
    if [[ -f "$PMSF_CONTREE" ]]; then
        log_info "  [OK] Consensus tree: $PMSF_CONTREE"
    fi
    logL_pmsf=$(grep "BEST SCORE FOUND" "${PMSF_DIR}/10mags_pmsf.iqtree" | awk '{print $5}' || echo "N/A")
    log_info "  Log-likelihood: $logL_pmsf"
else
    log_error "  [FAIL] PMSF tree not produced. Check: ${PMSF_DIR}/pmsf.iqtree2.log"
    exit 1
fi

log_info "  PMSF tree done: $(date '+%Y-%m-%d %H:%M:%S')"

# ---- Step C, Partition tree: ModelFinder per partition + UFBoot2 + SH-aLRT ----
log_info "Step C: Partition tree (MFP per partition) with UFBoot2 + SH-aLRT"
log_info "  Partition file: $PARTITION_FILE"
log_info "  Model:          -p (per-partition ModelFinder)"
log_info "  Supports:       UFBoot2 -B 1000 | SH-aLRT --alrt 1000"
log_info "  Threads:        $SLURM_CPUS_PER_TASK"
log_info "  Start:          $(date '+%Y-%m-%d %H:%M:%S')"

PART_DIR="${IQTREE_DIR}/partition_tree"
cd "$PART_DIR"

iqtree2 \
    -s "$SUPERMATRIX" \
    -p "$PARTITION_FILE" \
    -m MFP \
    -B 1000 \
    --alrt 1000 \
    -T "$SLURM_CPUS_PER_TASK" \
    --seed 42 \
    --prefix 10mags_partition \
    --redo \
    > "${PART_DIR}/10mags_partition.iqtree2.log" 2>&1

PART_TREEFILE="${PART_DIR}/10mags_partition.treefile"
PART_CONTREE="${PART_DIR}/10mags_partition.contree"

if [[ -f "$PART_TREEFILE" ]]; then
    log_info "  [OK] Partition tree: $PART_TREEFILE"
    if [[ -f "$PART_CONTREE" ]]; then
        log_info "  [OK] Consensus tree: $PART_CONTREE"
    fi
    logL_part=$(grep "BEST SCORE FOUND" "${PART_DIR}/10mags_partition.iqtree" | awk '{print $5}' || echo "N/A")
    log_info "  Log-likelihood: $logL_part"
else
    log_error "  [FAIL] Partition tree not produced. Check: ${PART_DIR}/10mags_partition.iqtree2.log"
    exit 1
fi

log_info "  Partition tree done: $(date '+%Y-%m-%d %H:%M:%S')"

# ---- Step D, rename taxa: NameWithoutExtension to DisplayName ----
log_info "Step D: Renaming taxa (NameWithoutExtension to DisplayName)"
METADATA="${MAGS_METADATA}"
REF_TAXONOMY="${WORK_DIR}/reference_taxonomy.tsv"
PMSF_RENAMED="${PMSF_DIR}/10mags_pmsf.renamed.treefile"
PMSF_CONTREE_RENAMED="${PMSF_DIR}/10mags_pmsf.renamed.contree"
PART_RENAMED="${PART_DIR}/10mags_partition.renamed.treefile"
PART_CONTREE_RENAMED="${PART_DIR}/10mags_partition.renamed.contree"

[[ -f "$METADATA" ]]     || log_warn "MAG metadata not found: $METADATA"
[[ -f "$REF_TAXONOMY" ]] || log_warn "Reference taxonomy not found: $REF_TAXONOMY, run script 18 first"

if [[ -f "$METADATA" || -f "$REF_TAXONOMY" ]]; then
    python3 - "$METADATA" "$REF_TAXONOMY" "$PMSF_TREEFILE" "$PMSF_RENAMED" <<'PYEOF'
import re, sys

metadata_file, ref_taxonomy_file, tree_in, tree_out = sys.argv[1:5]

def build_mapping(metadata_file, ref_taxonomy_file):
    mapping = {}

    # MAG metadata: NameWithoutExtension (col 2) to DisplayName (col 3)
    try:
        with open(metadata_file) as f:
            next(f)  # skip header
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 3 and parts[1] and parts[2]:
                    mapping[parts[1]] = parts[2].replace(' ', '_')
    except FileNotFoundError:
        pass

    # Reference taxonomy: tree_name (col 1) to species (col 2)
    try:
        with open(ref_taxonomy_file) as f:
            next(f)  # skip header
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 2 and parts[0] and parts[1]:
                    mapping[parts[0]] = parts[1]
    except FileNotFoundError:
        pass

    return mapping

mapping = build_mapping(metadata_file, ref_taxonomy_file)

with open(tree_in) as f:
    tree = f.read()

n_renamed = 0
# Sort longest-first to avoid partial matches (e.g. bin.1 before bin.15)
for old in sorted(mapping, key=len, reverse=True):
    new = mapping[old]
    # Newick: taxon names follow '(' or ',' and precede ':', ')', ',', or ';'
    replaced, count = re.subn(
        r'(?<=[,(])' + re.escape(old) + r'(?=[:,);])',
        new, tree
    )
    if count:
        tree = replaced
        n_renamed += count

with open(tree_out, 'w') as f:
    f.write(tree)

print(f"Renamed {n_renamed} taxon occurrence(s)")
PYEOF

    log_info "  [OK] Renamed PMSF tree: $PMSF_RENAMED"

    # Shared renaming helper as a bash function calling python3
    rename_tree() {
        local tree_in="$1" tree_out="$2"
        python3 - "$METADATA" "$REF_TAXONOMY" "$tree_in" "$tree_out" <<'PYEOF'
import re, sys

metadata_file, ref_taxonomy_file, tree_in, tree_out = sys.argv[1:5]

mapping = {}
try:
    with open(metadata_file) as f:
        next(f)
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3 and parts[1] and parts[2]:
                mapping[parts[1]] = parts[2].replace(' ', '_')
except FileNotFoundError:
    pass

try:
    with open(ref_taxonomy_file) as f:
        next(f)
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[0] and parts[1]:
                mapping[parts[0]] = parts[1]
except FileNotFoundError:
    pass

with open(tree_in) as f:
    tree = f.read()

for old in sorted(mapping, key=len, reverse=True):
    new = mapping[old]
    tree = re.sub(r'(?<=[,(])' + re.escape(old) + r'(?=[:,);])', new, tree)

with open(tree_out, 'w') as f:
    f.write(tree)
PYEOF
    }

    [[ -f "$PMSF_CONTREE"  ]] && rename_tree "$PMSF_CONTREE"  "$PMSF_CONTREE_RENAMED"  && log_info "  [OK] Renamed PMSF contree: $PMSF_CONTREE_RENAMED"
    [[ -f "$PART_TREEFILE" ]] && rename_tree "$PART_TREEFILE" "$PART_RENAMED"           && log_info "  [OK] Renamed partition tree: $PART_RENAMED"
    [[ -f "$PART_CONTREE"  ]] && rename_tree "$PART_CONTREE"  "$PART_CONTREE_RENAMED"   && log_info "  [OK] Renamed partition contree: $PART_CONTREE_RENAMED"
fi

# ---- summary ----
log_info "Pipeline completed successfully"
log_info "  Output structure:"
log_info "    ${IQTREE_DIR}/"
log_info "    ├── guide_tree/"
log_info "    │   └── guide.treefile          (guide tree, LG+F+G)"
log_info "    ├── pmsf_final/"
log_info "    │   ├── pmsf.treefile           (raw PMSF tree)"
log_info "    │   ├── pmsf.renamed.treefile   (PMSF tree, display labels) ← MAIN"
log_info "    │   └── pmsf.renamed.contree    (PMSF consensus, display labels)"
log_info "    └── partition_tree/"
log_info "        ├── partition.treefile      (raw partition tree)"
log_info "        ├── partition.renamed.treefile  (partition tree, display labels)"
log_info "        └── partition.renamed.contree   (partition consensus, display labels)"
log_info ""
log_info "  MAIN TREE (for figures):"
log_info "    $PMSF_RENAMED"
log_info "  ROBUSTNESS CHECK tree:"
log_info "    $PART_RENAMED"
log_info ""
log_info "  Branch support format: SH-aLRT/UFBoot2 (e.g., 95/98)"
log_info "    Cutoffs: SH-aLRT ≥ 80 AND UFBoot2 ≥ 95"
log_info ""
log_info "  Post-tree checks (CLAUDE.md):"
log_info "    1. Outgroup (Thermotoga/Aquifex) forms external clade to Archaea"
log_info "    2. Thermoplasmatota and Thermoproteota form two separate clades"
log_info "    3. MGIIb-Q1 and UBA36 are sister clades (same family)"
log_info "    4. Sublineages A and B of MGIIb-Q1 resolve as distinct clades"
log_info "    5. No branch length > 1.0 (LBA check)"
log_info ""
log_info "  Next step: visualise tree (FigTree / R ggtree) and root on"
log_info "             bacterial outgroup (Thermotoga / Aquifex branch)"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
