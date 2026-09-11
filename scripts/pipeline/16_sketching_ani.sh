#!/usr/bin/env bash
set -euo pipefail

# Step 16. Computes ANI over the 54 MAGs with skani, in four passes:
#
#   1. All-vs-all pairwise ANI across the 54 MAGs, for phylogenomic divergence.
#   2. The 10 species-unassigned MAGs against all 54, for intra-study position.
#   3. All 54 against GTDB representatives, to confirm the GTDB-Tk assignments.
#   4. The 10 unassigned MAGs against GTDB, for novel-species detection.
#
# The 10 MAGs with no GTDB-Tk species (Species = N/A in MAGs_metadata.tsv) fall
# in MGIIb-Q1 (5), Nitrosopelagicus (3) and UBA36 (2). The question for each is
# whether ANI stays below 95% against every GTDB reference in its genus, how far
# it sits from its described relatives, and where it lands among the other 53.
#
# Steps 3 and 4 need a pre-built skani GTDB database; pass --skip-gtdb to run
# only the intra-study passes. To build one:
#   wget http://faust.compbio.cs.cmu.edu/skani-files/skani_gtdb_r226-v0.3.tar.gz
#   tar -zxvf skani_gtdb_r226-v0.3.tar.gz
# then point GTDB_SKETCH_DIR at the extracted path.

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

trap 'log_error "Failed at line $LINENO with exit code $?"; exit 1' ERR

log_job_header "skani ANI over 54 archaeal MAGs"

# This step ran on the second cluster, so its profile is server_ixachi.env.
# Override CONFIG_FILE to use the other one.
CONFIG_FILE="${CONFIG_FILE:-${REPO_ROOT}/configs/server_ixachi.env}"
load_project_config

SKIP_GTDB=false
for arg in "$@"; do
    [[ "$arg" == "--skip-gtdb" ]] && SKIP_GTDB=true
done

require_vars PROJECT_DIR FILTERED_SUBDIR MAGS_METADATA
if [[ "$SKIP_GTDB" == false ]]; then
    require_vars GTDB_SKETCH_DIR
fi

THREADS="${SLURM_CPUS_PER_TASK:-24}"

MAG_DIR="$PROJECT_DIR/$FILTERED_SUBDIR"
OUT_DIR="$PROJECT_DIR/16_sketching_ani"
SKETCH_DIR="$OUT_DIR/sketches"                    # unified skani DB of the 54 MAGs
ALL_MAG_DIR="$OUT_DIR/all_mags"                   # all-vs-all triangle results
NOVEL_MAG_DIR="$OUT_DIR/novel_mags"               # results for the 10 unassigned MAGs
GTDB_DIR="$OUT_DIR/gtdb_search"                   # GTDB search results
GTDB_SKETCH_LIST="$OUT_DIR/gtdb_sketch_list.txt"  # list of individual GTDB .sketch files
MAG_EXT="fa"

# ANI parameters, applied to the search and dist passes. triangle accepts only
# some of them.
#
# --min-af 0.1, under the skani default of 0.15:
#   Archaeal MAGs from metagenomes are fragmented, and a lower aligned-fraction
#   floor keeps partially overlapping or distant genomes from being dropped
#   without a trace. It matters most for the novel MAGs, whose closest GTDB
#   reference may share little synteny.
#
# -s 75, under the skani default of 80:
#   Novel archaeal genera can sit at 78 to 82% ANI from their closest GTDB
#   representative. Dropping the screen by five points, as skani recommends,
#   keeps those borderline hits alive long enough for the accurate estimator
#   to run on them.
#
# --slow (-c 30):
#   Best accuracy on fragmented genomes with N50 in the low kb. It trades speed
#   for dependable aligned-fraction estimates, which is what the species and
#   genus boundary calls rest on.
MIN_AF="0.1"
ANI_SCREEN="75"

# The 10 MAGs whose Species column reads N/A in MAGs_metadata.tsv.
NOVEL_MAGS=(
    "A4.1000.megahit.bin.28"   # MGIIb-Q1
    "A4.BTM.megahit.bin.15"    # MGIIb-Q1
    "A4.BTM.megahit.bin.25"    # Nitrosopelagicus
    "A4.BTM.megahit.bin.33"    # Nitrosopelagicus
    "A4.BTM.megahit.bin.58"    # UBA36
    "A6.BTM.megahit.bin.3"     # MGIIb-Q1
    "A6.BTM.megahit.bin.9"     # UBA36
    "PO1.BTM.megahit.bin.42"   # MGIIb-Q1
    "PO1.BTM.megahit.bin.7"    # Nitrosopelagicus
    "Y6.1000.megahit.bin.23"   # MGIIb-Q1
)

log_info "MAG input directory : $MAG_DIR"
log_info "Output directory    : $OUT_DIR"
log_info "GTDB sketch dir     : ${GTDB_SKETCH_DIR:-not set}"
log_info "Threads             : $THREADS"
log_info "Novel MAGs          : ${#NOVEL_MAGS[@]}"

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

# Every novel MAG is checked up front, so a typo in the list above surfaces
# before any of the long passes start.
for mag in "${NOVEL_MAGS[@]}"; do
    if [[ ! -f "${MAG_DIR}/${mag}.${MAG_EXT}" ]]; then
        log_error "Novel MAG file not found: ${MAG_DIR}/${mag}.${MAG_EXT}"
        exit 1
    fi
done
log_info "All ${#NOVEL_MAGS[@]} novel MAG files present"

if [[ "$SKIP_GTDB" == false && ! -d "$GTDB_SKETCH_DIR" ]]; then
    log_error "GTDB sketch directory not found: $GTDB_SKETCH_DIR"
    log_error "Build it as described at the top of this script, or pass --skip-gtdb."
    exit 1
fi

mkdir -p "$ALL_MAG_DIR" "$NOVEL_MAG_DIR" "$GTDB_DIR"
# skani sketch insists on creating this directory itself.
rm -rf "${SKETCH_DIR:?}"

# ---- step 1: sketch the 54 MAGs into one database ----
#
# A compact indexed database in the skani 0.3 unified format, one directory
# rather than individual .sketch files. Steps 3 and 4 query it, which keeps
# their memory use flat. Step 2 goes straight from the FASTA files instead,
# since triangle sketches as it goes.

log_info "Step 1: sketching 54 MAGs into a unified database"

skani sketch \
    "$MAG_DIR"/*.${MAG_EXT} \
    -o "$SKETCH_DIR" \
    -t "$THREADS"
log_info "  Database: $SKETCH_DIR"

# ---- step 2: all-vs-all ANI triangle ----
#
# 54 * 53 / 2 = 1431 unique pairwise values. What they give:
#   (a) an ANI distance matrix for clustering, complementary to the
#       marker-gene phylogenomics of steps 10 to 13;
#   (b) same-species pairs inside the study, at ANI >= 95% with AF >= 50%,
#       which is how several Thalassarchaeum sp002507125 MAGs from different
#       stations turn out to be the same species;
#   (c) between-genus divergence, expected around 75 to 90% ANI.
#
# triangle rather than dist because dist computes all n^2 comparisons, which is
# redundant for a symmetric matrix. triangle does the unique pairs and, with
# --full-matrix, writes the square phylip form that R reads through as.dist.
#
# --min-af is not available on triangle, so -s 75 does the work of keeping
# distant inter-genus pairs in the output. Pairs under the screen are simply
# absent from the sparse output and count as 0 downstream.

log_info "Step 2: all-vs-all ANI triangle"

FULL_MATRIX="$ALL_MAG_DIR/ani_triangle_full.tsv"
SPARSE_MATRIX="$ALL_MAG_DIR/ani_triangle_sparse.tsv"

if [[ -f "$FULL_MATRIX" ]]; then
    log_info "  Full matrix exists, skipping: $FULL_MATRIX"
else
    skani triangle \
        "$MAG_DIR"/*.${MAG_EXT} \
        --full-matrix \
        --slow \
        -s "$ANI_SCREEN" \
        -t "$THREADS" \
        -o "$FULL_MATRIX"
    log_info "  Full matrix: $FULL_MATRIX"
fi

if [[ -f "$SPARSE_MATRIX" ]]; then
    log_info "  Sparse matrix exists, skipping: $SPARSE_MATRIX"
else
    skani triangle \
        "$MAG_DIR"/*.${MAG_EXT} \
        --sparse \
        --slow \
        -s "$ANI_SCREEN" \
        -t "$THREADS" \
        -o "$SPARSE_MATRIX"
    log_info "  Sparse matrix: $SPARSE_MATRIX"
fi

# ---- step 3: the 10 novel MAGs against all 54 ----
#
# Queries the unassigned MAGs against the unified database to find the closest
# MAG inside the study for each one, and to settle whether novel MAGs of the
# same genus (the five MGIIb-Q1, say) sit above 95% ANI from each other, which
# would make them one novel species spread across stations rather than several
# independent lineages.

log_info "Step 3: 10 novel MAGs against all 54"

NOVEL_VS_ALL="$NOVEL_MAG_DIR/novel10_vs_all54.tsv"

novel_fasta_paths=()
for mag in "${NOVEL_MAGS[@]}"; do
    novel_fasta_paths+=("${MAG_DIR}/${mag}.${MAG_EXT}")
done

if [[ -f "$NOVEL_VS_ALL" ]]; then
    log_info "  Output exists, skipping: $NOVEL_VS_ALL"
else
    skani search \
        "${novel_fasta_paths[@]}" \
        -d "$SKETCH_DIR" \
        --min-af "$MIN_AF" \
        -s "$ANI_SCREEN" \
        -t "$THREADS" \
        -o "$NOVEL_VS_ALL"

    log_info "  Saved: $NOVEL_VS_ALL"
    log_info "  Hits: $(( $(wc -l < "$NOVEL_VS_ALL") - 1 ))"
fi

# skani ships the GTDB database in either of two shapes: consolidated, a single
# directory holding markers.bin, or individual, the GTDB-Tk style hierarchy of
# .sketch files. The flag that selects the reference differs between them.
GTDB_REF_ARG=()
if [[ "$SKIP_GTDB" == false ]]; then
    if [[ -f "$GTDB_SKETCH_DIR/markers.bin" ]]; then
        log_info "GTDB database: consolidated format"
        GTDB_REF_ARG=(-d "$GTDB_SKETCH_DIR")
    else
        find "$GTDB_SKETCH_DIR" -name "*.sketch" | sort > "$GTDB_SKETCH_LIST"
        N_GTDB_SKETCHES=$(wc -l < "$GTDB_SKETCH_LIST")

        if [[ "$N_GTDB_SKETCHES" -gt 0 ]]; then
            log_info "GTDB database: individual format, $N_GTDB_SKETCHES sketches"
            GTDB_REF_ARG=(--rl "$GTDB_SKETCH_LIST")
        else
            log_warn "Neither markers.bin nor .sketch files under $GTDB_SKETCH_DIR"
            log_warn "Build the database as described at the top of this script."
            log_warn "Skipping steps 4 and 5."
            SKIP_GTDB=true
        fi
    fi
fi

if [[ "$SKIP_GTDB" == true ]]; then
    log_warn "Steps 4 and 5 skipped."
else
    # ---- step 4: all 54 MAGs against GTDB ----
    #
    # Confirms the GTDB-Tk assignments: an assigned MAG should come back at
    # ANI >= 95% with AF >= 50% against its designated reference. For the 10
    # unassigned ones it gives the closest reference and the distance to it,
    # which cross-references against the classify_wf columns
    # closest_placement_ani, closest_placement_af and red_value.
    #
    # skani search columns: Ref_file, Query_file, ANI, Align_fraction_ref,
    # Align_fraction_query, Ref_name, Query_name.

    log_info "Step 4: all 54 MAGs against GTDB"

    ALL_VS_GTDB="$GTDB_DIR/all54_vs_gtdb.tsv"

    if [[ -f "$ALL_VS_GTDB" ]]; then
        log_info "  Output exists, skipping: $ALL_VS_GTDB"
    else
        skani search \
            "$MAG_DIR"/*.${MAG_EXT} \
            "${GTDB_REF_ARG[@]}" \
            --min-af "$MIN_AF" \
            -s "$ANI_SCREEN" \
            -t "$THREADS" \
            -o "$ALL_VS_GTDB"

        log_info "  Saved: $ALL_VS_GTDB"
        log_info "  Hits: $(( $(wc -l < "$ALL_VS_GTDB") - 1 ))"
    fi

    # ---- step 5: the 10 novel MAGs against GTDB ----
    #
    # The same search restricted to the MAGs classify_wf left at s__, written
    # separately so the downstream analysis does not have to filter a large
    # file. How to read it:
    #
    #   Novel species  ANI below 95% against every GTDB reference in the same
    #                  genus, corroborated by a GTDB-Tk red_value above the
    #                  genus threshold, around 0.45 for Archaea.
    #   Novel genus    no hit at ANI >= 80% to any genome in the same family.
    #                  Unlikely here, since GTDB-Tk did assign genera, but the
    #                  -s 75 screen is low enough to catch it.
    #   Divergence     best-hit ANI and aligned fraction per genus, which is
    #                  what quantifies how far MGIIb-Q1, Nitrosopelagicus and
    #                  UBA36 sit from their nearest described species.
    #
    # skani is tuned for ANI at or above 82%. Hits near or below that need
    # cross-validation against the GTDB-Tk placement metrics.

    log_info "Step 5: 10 novel MAGs against GTDB"

    NOVEL_VS_GTDB="$GTDB_DIR/novel10_vs_gtdb.tsv"

    if [[ -f "$NOVEL_VS_GTDB" ]]; then
        log_info "  Output exists, skipping: $NOVEL_VS_GTDB"
    else
        skani search \
            "${novel_fasta_paths[@]}" \
            "${GTDB_REF_ARG[@]}" \
            --min-af "$MIN_AF" \
            -s "$ANI_SCREEN" \
            -t "$THREADS" \
            -o "$NOVEL_VS_GTDB"

        log_info "  Saved: $NOVEL_VS_GTDB"
        log_info "  Hits: $(( $(wc -l < "$NOVEL_VS_GTDB") - 1 ))"
    fi
fi

log_info "Output files:"
log_info "  All 54, triangle full    $FULL_MATRIX"
log_info "  All 54, triangle sparse  $SPARSE_MATRIX"
log_info "  Novel 10 vs 54 MAGs      $NOVEL_VS_ALL"
if [[ "$SKIP_GTDB" == false ]]; then
    log_info "  All 54 vs GTDB           $ALL_VS_GTDB"
    log_info "  Novel 10 vs GTDB         $NOVEL_VS_GTDB"
else
    log_warn "  GTDB passes skipped"
fi
log_info "Same species when ANI >= 95% and aligned fraction >= 50%"
log_info "Novel species when no GTDB hit reaches ANI 95%"
log_info "Aligned fraction floor: ${MIN_AF} (search passes only)"
log_info "ANI screen: ${ANI_SCREEN}%"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
