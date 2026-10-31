#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Step 25. Writes a publication-ready copy of each pangenome,
# <species>-PAN_modified.db, alongside a copy of <species>-GENOMES.db. The
# originals are never touched.
#
# What the copy changes:
#   - drops the selected items and layers misc-data keys that were working
#     state rather than results;
#   - translates the parameter names anvi'o displays into Spanish, since the
#     thesis figures are in Spanish;
#   - leaves genome_name keys alone, because step 23 already imported
#     DisplayName as a layer and anvi-display-pan labels genomes from it.
#
# Pfam, COG24_PATHWAY and COG24_FUNCTION are deleted from
# GENOMES_modified.db.gene_function_calls, where they are stored;
# they are not in PAN.db. EVERYTHING stays, because it is synthetic, the union
# of the remaining functional sources, and removing it would mean removing all
# of them.
#
# This one runs locally: it is sqlite work and small
# TSV edits, nothing that needs a scheduler.

trap 'log_error "Script failed at line $LINENO with exit code $?"; exit 1' ERR

log_info "anvi'o 9 pangenome cleanup for publication"
log_info "Host: $(hostname)"
log_info "Start time: $(date '+%Y-%m-%d %H:%M:%S')"

require_commands anvi-export-misc-data anvi-import-misc-data anvi-delete-misc-data sqlite3

# ---- local paths ----
# A laptop may have no configs/server.env, so the clone itself is the fallback
# root, the same one the R scripts use. MSC_ARCHAEA_DIR overrides it.
if [[ -f "${CONFIG_FILE:-${REPO_ROOT}/configs/server.env}" ]]; then
    load_project_config
fi
RESULTS_DIR="${RESULTS_DIR:-${MSC_ARCHAEA_DIR:-$REPO_ROOT}/results}"
PANGENOME_DIR="${RESULTS_DIR}/${PANGENOME_SUBDIR:-20_pangenome}"
PANGENOME_WORK_DIR="${PANGENOME_DIR}/pangenome"
MAGS_METADATA="${MAGS_METADATA:-${REPO_ROOT}/data/MAGs_metadata.tsv}"

if [[ ! -d "$PANGENOME_WORK_DIR" ]]; then
    log_error "Pangenome working directory does not exist: $PANGENOME_WORK_DIR"
    exit 1
fi
if [[ ! -f "$MAGS_METADATA" ]]; then
    log_error "Metadata file not found: $MAGS_METADATA"
    exit 1
fi

# Keys to delete from items / layers misc-data tables
ITEMS_KEYS_TO_REMOVE=(
    "EVERYTHING"
    "max_num_paralogs"
    "combined_homogeneity_index"
    "Pfam"
    "KOfam"
    "KEGG_Module"
    "KEGG_Class"
    "COG24_PATHWAY"
    "COG24_FUNCTION"
)
LAYERS_KEYS_TO_REMOVE=(
    "num_genes_per_kb"
    "gc_content"
)

# English -> Spanish renames (real anvi'o v9 pangenome keys)
# Item keys, capitalised for display
declare -A RENAME_ITEMS=(
    [num_genomes_gene_cluster_has_hits]="Prevalencia_genómica"
    [num_genes_in_gene_cluster]="Genes_en_cluster"
    [geometric_homogeneity_index]="Homogeneidad_geometrica"
    [functional_homogeneity_index]="Homogeneidad_funcional"
    [SCG]="Clusteres_SCG"
    [COG24_CATEGORY]="Categoria_COG24"
    [AAI_min]="AAI_minimo"
    [AAI_max]="AAI_maximo"
    [AAI_avg]="AAI_promedio"
)
# Layers
declare -A RENAME_LAYERS=(
    [num_gene_clusters]="Clusteres_genicos_totales"
    [singleton_gene_clusters]="Clusteres_singleton"
    [percent_redundancy]="Redundancia"
    [percent_completion]="Completitud"
    [total_length]="Longitud_del_MAG"
    [num_genes]="Numero_genes"
    [avg_gene_length]="Longitud_promedio_gen"
)

# Function annotation sources to delete from GENOMES_modified.db.gene_function_calls
GENE_FUNCTION_SOURCES_TO_REMOVE=(
    "Pfam"
    "COG24_PATHWAY"
    "COG24_FUNCTION"
)

# ---- per-species cleanup loop ----
log_info "Starting per-species cleanup..."
for species_dir in "${PANGENOME_WORK_DIR}"/*/; do
    [[ -d "$species_dir" ]] || continue
    species_name=$(basename "$species_dir")
    pan_db="${species_dir}/pan_output/${species_name}-PAN.db"
    pan_mod="${species_dir}/pan_output/${species_name}-PAN_modified.db"
    gen_db="${species_dir}/${species_name}-GENOMES.db"
    gen_mod="${species_dir}/${species_name}-GENOMES_modified.db"
    external_genomes_file="${species_dir}/external_genomes.tsv"

    if [[ ! -f "$pan_db" ]]; then
        log_warn "  PAN.db not found for ${species_name//_/ }; skipping"
        continue
    fi
    if [[ ! -f "$gen_db" ]]; then
        log_warn "  GENOMES.db not found for ${species_name//_/ }; skipping"
        continue
    fi

    log_info "  Species: ${species_name//_/ }"

# ---- step 0: working copies, with integrity checks ----
    log_info "    Step 0: copying PAN.db -> PAN_modified.db; GENOMES.db -> GENOMES_modified.db"
    cp -f "$pan_db" "$pan_mod"
    cp -f "$gen_db" "$gen_mod"
    for db in "$pan_mod" "$gen_mod"; do
        integrity=$(sqlite3 "$db" "PRAGMA integrity_check;")
        if [[ "$integrity" != "ok" ]]; then
            log_error "    Integrity check failed for $db: $integrity"
            continue 2
        fi
    done
    log_info "    Originals intact: $(basename "$pan_db"), $(basename "$gen_db")"
    log_info "    Working copies:   $(basename "$pan_mod"), $(basename "$gen_mod")"

    # Snapshot of the keys before anything is removed
    log_info "    Initial item_additional_data keys:"
    sqlite3 "$pan_mod" "SELECT DISTINCT data_key FROM item_additional_data;" \
        | sed 's/^/      /'
    log_info "    Initial layer_additional_data keys:"
    sqlite3 "$pan_mod" "SELECT DISTINCT data_key FROM layer_additional_data;" \
        | sed 's/^/      /'

# ---- task 1: drop items keys ----
    log_info "    Task 1a: deleting items keys from PAN_modified.db.item_additional_data"
    for key in "${ITEMS_KEYS_TO_REMOVE[@]}"; do
        if anvi-delete-misc-data \
            -p "$pan_mod" \
            --target-data-table items \
            --keys-to-remove "$key" >/dev/null 2>&1; then
            log_info "      removed items key: $key"
        else
            log_info "      items key not present (skipped): $key"
        fi
    done

    log_info "    Task 1b: deleting function annotation sources from GENOMES_modified.db.gene_function_calls"
    for src in "${GENE_FUNCTION_SOURCES_TO_REMOVE[@]}"; do
        deleted=$(sqlite3 "$gen_mod" "DELETE FROM gene_function_calls WHERE source = '$src'; SELECT changes();")
        log_info "      removed source '$src': $deleted rows"
    done
    # Sync GENOMES.self.gene_function_sources: drop the removed sources from the comma-separated list
    python3 - <<PYEOF
import sqlite3
remove = set("""${GENE_FUNCTION_SOURCES_TO_REMOVE[*]}""".split())
con = sqlite3.connect("${gen_mod}")
cur = con.cursor()
cur.execute("SELECT value FROM self WHERE key='gene_function_sources'")
row = cur.fetchone()
if row and row[0]:
    sources = [s for s in row[0].split(",") if s and s not in remove]
    new_val = ",".join(sources)
    cur.execute("UPDATE self SET value=? WHERE key='gene_function_sources'", (new_val,))
    con.commit()
    print(f"      GEN.self.gene_function_sources updated -> {new_val}")
con.close()
PYEOF
    sqlite3 "$gen_mod" "VACUUM;" >/dev/null 2>&1 || true

# ---- task 2: drop layers keys ----
    log_info "    Task 2: deleting layers keys"
    for key in "${LAYERS_KEYS_TO_REMOVE[@]}"; do
        if anvi-delete-misc-data \
            -p "$pan_mod" \
            --target-data-table layers \
            --keys-to-remove "$key" >/dev/null 2>&1; then
            log_info "      removed layers key: $key"
        else
            log_info "      layers key not present (skipped): $key"
        fi
    done

    # Task 3: genome_name stays as it is.
    # Step 23 already imported DisplayName as a layer, and anvi-display-pan
    # labels genomes from it. The internal genome_name stays as the sanitised
    # NameWithoutExtension so it keeps matching GENOMES.db and every internal
    # table in PAN.db.
    log_info "    Task 3: skipped (genome_name kept; DisplayName layer already provides relabeling)"

# ---- task 4: translate the parameter names, items and layers ----
    log_info "    Task 4: translating parameter names to Spanish"
    items_tsv_en="${species_dir}/items_data.tsv"
    items_tsv_es="${species_dir}/items_data_esp.tsv"
    layers_tsv_en="${species_dir}/layers_data.tsv"
    layers_tsv_es="${species_dir}/layers_data_esp.tsv"

    anvi-export-misc-data \
        -p "$pan_mod" \
        --target-data-table items \
        -o "$items_tsv_en"

    anvi-export-misc-data \
        -p "$pan_mod" \
        --target-data-table layers \
        -o "$layers_tsv_en"

    # The output TSV keeps the ID column and the renamed columns, nothing else.
    # Keeping the untouched columns would collide on re-import with the
    # original keys, which are still in the database at this point.
    python3 - <<PYEOF
import pandas as pd

rename_items = {
    "num_genomes_gene_cluster_has_hits": "Prevalencia_genómica",
    "num_genes_in_gene_cluster": "Genes_en_cluster",
    "geometric_homogeneity_index": "Homogeneidad_geometrica",
    "functional_homogeneity_index": "Homogeneidad_funcional",
    "SCG": "Clusteres_SCG",
    "COG24_CATEGORY": "Categoria_COG24",
    "AAI_min": "AAI_minimo",
    "AAI_max": "AAI_maximo",
    "AAI_avg": "AAI_promedio",
}
rename_layers = {
    "num_gene_clusters": "Clusteres_genicos_totales",
    "singleton_gene_clusters": "Clusteres_singleton",
    "percent_redundancy": "Redundancia",
    "percent_completion": "Completitud",
    "total_length": "Longitud_del_MAG",
    "num_genes": "Numero_genes",
    "avg_gene_length": "Longitud_promedio_gen",
}

def rewrite(path_in, path_out, rename_map):
    df = pd.read_csv(path_in, sep="\t")
    id_col = df.columns[0]
    present = {k: v for k, v in rename_map.items() if k in df.columns}
    keep_cols = [id_col] + list(present.keys())
    out = df[keep_cols].rename(columns=present)
    out.to_csv(path_out, sep="\t", index=False)
    print(f"  {path_in.split('/')[-1]}: kept {id_col} + renamed {len(present)} cols -> {list(present.values())}")

rewrite("${items_tsv_en}", "${items_tsv_es}", rename_items)
rewrite("${layers_tsv_en}", "${layers_tsv_es}", rename_layers)
PYEOF

    # Import the Spanish versions
    anvi-import-misc-data \
        "$items_tsv_es" \
        -p "$pan_mod" \
        --target-data-table items

    anvi-import-misc-data \
        "$layers_tsv_es" \
        -p "$pan_mod" \
        --target-data-table layers

    # Drop the original English keys, skipping any that were never there.
    for key in "${!RENAME_ITEMS[@]}"; do
        if anvi-delete-misc-data \
            -p "$pan_mod" \
            --target-data-table items \
            --keys-to-remove "$key" >/dev/null 2>&1; then
            log_info "      removed legacy items key: $key"
        fi
    done
    for key in "${!RENAME_LAYERS[@]}"; do
        if anvi-delete-misc-data \
            -p "$pan_mod" \
            --target-data-table layers \
            --keys-to-remove "$key" >/dev/null 2>&1; then
            log_info "      removed legacy layers key: $key"
        fi
    done

    # Task 5: import only Station and Depth, translated, as layers.
    # DisplayName, Order and ROMZ_flag go: the short genome_name replaces
    # DisplayName in task 6 below.
    # Depth translations:
    #   OMZ / ROMZ -> ZMOR
    #   1000       -> "1000 m"
    #   BTM        -> Fondo
    log_info "    Task 5: building & importing layers (Estacion, Profundidad ES)"
    species_layers_tsv="${species_dir}/misc_data_layers_${species_name}.tsv"

    if [[ ! -f "$external_genomes_file" ]]; then
        log_warn "      external_genomes.tsv not found; skipping Tarea 5"
    else
        awk -F'\t' '
            NR == FNR {
                if (FNR > 1) keep[$1] = 1
                next
            }
            FNR == 1 {
                print "item\tEstacion\tProfundidad"
                next
            }
            {
                name = $2
                gsub(/[^a-zA-Z0-9_]/, "_", name)
                if (!(name in keep)) next
                station = ($4 == "" ? "NA" : $4)
                depth_orig = ($5 == "" ? "NA" : $5)
                if (depth_orig == "OMZ" || depth_orig == "ROMZ")      depth = "ZMOR"
                else if (depth_orig == "1000")                         depth = "1000 m"
                else if (depth_orig == "BTM")                          depth = "Fondo"
                else                                                    depth = depth_orig
                print name "\t" station "\t" depth
            }
        ' "$external_genomes_file" "$MAGS_METADATA" > "$species_layers_tsv"

        n_rows=$(($(wc -l < "$species_layers_tsv") - 1))
        if [[ "$n_rows" -le 0 ]]; then
            log_warn "      No metadata rows produced; skipping import"
        else
            log_info "      Importing $n_rows metadata rows"
            # Removing the keys first makes a re-run idempotent, and also clears
            # the DisplayName, Order and ROMZ_flag layers earlier versions imported.
            for key in Station Depth Estacion Profundidad DisplayName Order ROMZ_flag; do
                anvi-delete-misc-data \
                    -p "$pan_mod" \
                    --target-data-table layers \
                    --keys-to-remove "$key" >/dev/null 2>&1 || true
            done
            anvi-import-misc-data \
                "$species_layers_tsv" \
                -p "$pan_mod" \
                --target-data-table layers
        fi
    fi

    # Task 6: rename genome_name to PanShortName (Station_Depth_BinID),
    # cascading through PAN_modified.db and GENOMES_modified.db.
    # Tables affected:
    #   PAN: gene_clusters.genome_name, gene_cluster_presence_absence.layer,
    #        gene_cluster_frequencies.layer, layer_additional_data.item_name,
    #        layer_orders.data (newick string).
    #   GEN: genome_info.genome_name, gene_function_calls.genome_name,
    #        gene_info.genome_name, self.value (where key='genome_names').
    # genomes_storage_hash stays identical across both databases;
    # anvi'o only checks that the two agree; it does not recompute them.
    log_info "    Task 6: renaming genome_name -> PanShortName (cascade across PAN+GENOMES)"
    python3 - <<PYEOF
import sqlite3, csv, re, sys

def sanitize(s):
    return re.sub(r"[^a-zA-Z0-9_]", "_", s or "")

# The genomes this pangenome covers
keep = set()
with open("${external_genomes_file}") as f:
    next(f)
    for line in f:
        keep.add(line.split("\t")[0])

# Old to new name mapping, read from MAGs_metadata.tsv
mapping = {}
with open("${MAGS_METADATA}") as f:
    rdr = csv.DictReader(f, delimiter="\t")
    if "PanShortName" not in rdr.fieldnames:
        print("ERROR: MAGs_metadata.tsv no tiene la columna PanShortName")
        sys.exit(1)
    for row in rdr:
        old = sanitize(row["NameWithoutExtension"])
        if old not in keep:
            continue
        new = sanitize(row["PanShortName"])
        if old == new:
            continue
        mapping[old] = new

if not mapping:
    print("      no genomes to rename for this species")
    sys.exit(0)

print(f"      Renaming {len(mapping)} genomes:")
for o, n in sorted(mapping.items()):
    print(f"        {o}  ->  {n}")

# Sanity check: the target names have to be unique
if len(set(mapping.values())) != len(mapping):
    print("ERROR: PanShortName values colisionan dentro de esta especie")
    sys.exit(1)

# --- Update PAN_modified.db ---
con = sqlite3.connect("${pan_mod}")
cur = con.cursor()
for old, new in mapping.items():
    cur.execute("UPDATE gene_clusters SET genome_name = ? WHERE genome_name = ?", (new, old))
    cur.execute("UPDATE layer_additional_data SET item_name = ? WHERE item_name = ?", (new, old))
    cur.execute("UPDATE gene_cluster_presence_absence SET layer = ? WHERE layer = ?", (new, old))
    cur.execute("UPDATE gene_cluster_frequencies SET layer = ? WHERE layer = ?", (new, old))

# layer_orders is (data_key, data_type, data_value); data_value holds
# the newick, and any other string that contains genome_names
cur.execute("SELECT data_key, data_value FROM layer_orders")
rows = cur.fetchall()
for data_key, data_value in rows:
    if not data_value:
        continue
    new_value = data_value
    for old, new in mapping.items():
        new_value = re.sub(rf"\b{re.escape(old)}\b", new, new_value)
    if new_value != data_value:
        cur.execute("UPDATE layer_orders SET data_value = ? WHERE data_key = ?", (new_value, data_key))

# PAN.self: update the values holding genome names (external_genome_names,
# internal_genome_names, default_item_order, available_item_orders, items_ordered)
cur.execute("SELECT key, value FROM self")
for k, v in cur.fetchall():
    if v and isinstance(v, str):
        new_v = v
        for old, new in mapping.items():
            new_v = re.sub(rf"\b{re.escape(old)}\b", new, new_v)
        if new_v != v:
            cur.execute("UPDATE self SET value = ? WHERE key = ?", (new_v, k))
            print(f"      updated PAN.self[{k}]")

con.commit()
con.close()

# --- Update GENOMES_modified.db ---
con = sqlite3.connect("${gen_mod}")
cur = con.cursor()
for old, new in mapping.items():
    cur.execute("UPDATE genome_info SET genome_name = ? WHERE genome_name = ?", (new, old))
    cur.execute("UPDATE gene_function_calls SET genome_name = ? WHERE genome_name = ?", (new, old))
    cur.execute("UPDATE gene_info SET genome_name = ? WHERE genome_name = ?", (new, old))

# self table: a possible 'genome_names' entry (comma-separated list)
cur.execute("SELECT key, value FROM self")
for k, v in cur.fetchall():
    if v and isinstance(v, str):
        new_v = v
        for old, new in mapping.items():
            new_v = re.sub(rf"\b{re.escape(old)}\b", new, new_v)
        if new_v != v:
            cur.execute("UPDATE self SET value = ? WHERE key = ?", (new_v, k))
            print(f"      updated GEN.self[{k}]")

con.commit()
con.close()
print("      Cascade rename complete.")
PYEOF

# ---- final snapshot of the keys ----
    log_info "    Final item_additional_data keys:"
    sqlite3 "$pan_mod" "SELECT DISTINCT data_key FROM item_additional_data;" \
        | sed 's/^/      /'
    log_info "    Final layer_additional_data keys:"
    sqlite3 "$pan_mod" "SELECT DISTINCT data_key FROM layer_additional_data;" \
        | sed 's/^/      /'

    log_info "  Done: ${species_name//_/ }"
done
log_info "Pipeline finished successfully"
log_info "End time: $(date '+%Y-%m-%d %H:%M:%S')"
log_info "Per-species outputs:"
log_info "  <species>/pan_output/<species>-PAN_modified.db   (publication-ready)"
log_info "  <species>/<species>-GENOMES_modified.db          (without Pfam/COG24_PATHWAY/COG24_FUNCTION)"
log_info "  <species>/items_data.tsv                (intermediate, English)"
log_info "  <species>/items_data_esp.tsv            (intermediate, Spanish)"
log_info "  <species>/layers_data.tsv               (intermediate, English)"
log_info "  <species>/layers_data_esp.tsv           (intermediate, Spanish)"
log_info "  <species>/misc_data_layers_<species>.tsv (Tarea 5: Estacion + Profundidad ES)"
log_info "Inspect: anvi-display-pan -p <species>-PAN_modified.db -g <species>-GENOMES_modified.db"
log_info "Note: layer colors will fall back to greys (anvi'o uses hardcoded names);"
log_info "      use 'Save State' in the GUI after customizing colors to persist them."
