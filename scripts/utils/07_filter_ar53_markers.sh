#!/bin/bash

# Drops short protein sequences from the ar53 marker files, then drops the
# MAGs left carrying fewer than half the markers.

set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

load_project_config

AR53_DIR="${PROJECT_DIR}/${GTDBTK_DENOVO_TREE_SUBDIR:-09_de_novo_tree}/identify/intermediate_results/single_copy_fasta/ar53"

OUTPUT_DIR="${PROJECT_DIR}/${MARKER_ALIGNMENT_SUBDIR:-10_marker_alignment}/ar53_filtered"

# Scratch space for the intermediate files
TEMP_DIR="${OUTPUT_DIR}/.temp_processing"

# Configuración
TOTAL_MARKERS=53
MIN_MARKER_PERCENT=50  # Porcentaje mínimo de marcadores requerido
MIN_MARKERS=$(awk -v total=$TOTAL_MARKERS -v percent=$MIN_MARKER_PERCENT 'BEGIN {print int((total * percent / 100) + 0.5)}')

mkdir -p "$OUTPUT_DIR"
mkdir -p "$TEMP_DIR"

# Timestamped log file
LOG_FILE="${OUTPUT_DIR}/filtrado_$(date +%Y%m%d_%H%M%S).log"
MAG_REPORT="${OUTPUT_DIR}/mag_completeness_$(date +%Y%m%d_%H%M%S).tsv"

log_msg() {
    echo "$1" | tee -a "$LOG_FILE"
}

# Iniciar log
log_msg "========================================"
log_msg "Filtrado de marcadores AR53"
log_msg "Fecha: $(date '+%Y-%m-%d %H:%M:%S')"
log_msg "========================================"
log_msg "Directorio de entrada: $AR53_DIR"
log_msg "Directorio de salida: $OUTPUT_DIR"
log_msg "Directorio temporal: $TEMP_DIR"
log_msg "Archivo de log: $LOG_FILE"
log_msg ""
log_msg "Configuración:"
log_msg "  Total de marcadores AR53: $TOTAL_MARKERS"
log_msg "  Mínimo porcentaje requerido: ${MIN_MARKER_PERCENT}%"
log_msg "  Mínimo marcadores requeridos: $MIN_MARKERS"
log_msg ""

# ============================================================
# PASO 1: Filtrar secuencias cortas
# ============================================================

log_msg "========================================"
log_msg "PASO 1: Filtrado por longitud de secuencia"
log_msg "========================================"

processed=0
total_seqs_kept_length=0
total_seqs_removed_length=0
total_files=$(ls "$AR53_DIR"/*.fa 2>/dev/null | wc -l)

log_msg "Total de archivos a procesar: $total_files"
log_msg ""

# Procesar cada archivo .fa
for fasta_file in "$AR53_DIR"/*.fa; do
    [ -f "$fasta_file" ] || continue
    
    basename_file=$(basename "$fasta_file")
    temp_output_file="$TEMP_DIR/$basename_file"
    
    log_msg "Procesando: $basename_file"
    
    # Calcular longitud promedio usando awk
    avg_length=$(awk '
    BEGIN { total_length=0; seq_count=0; current_seq="" }
    /^>/ {
        if (current_seq != "") {
            total_length += length(current_seq)
            ++seq_count
        }
        current_seq = ""
        next
    }
    { current_seq = current_seq $0 }
    END {
        if (current_seq != "") {
            total_length += length(current_seq)
            ++seq_count
        }
        if (seq_count > 0) {
            print int(total_length / seq_count)
        } else {
            print 0
        }
    }
    ' "$fasta_file")
    
    # A sequence shorter than 30% of the marker mean carries too little signal
    # to align usefully, so that is the cut.
    min_length=$(awk -v avg="$avg_length" 'BEGIN { print int(avg * 0.30) }')
    
    log_msg "  Longitud promedio: $avg_length AA"
    log_msg "  Umbral mínimo (30%): $min_length AA"
    
    temp_stats="${TEMP_DIR}/.stats_${basename_file}.tmp"
    
    # Filtrar secuencias y escribir archivo temporal
    awk -v min_len="$min_length" -v stats_file="$temp_stats" '
    BEGIN { current_header=""; current_seq=""; kept=0; removed=0 }
    /^>/ {
        # Procesar secuencia anterior
        if (current_header != "") {
            seq_length = length(current_seq)
            if (seq_length > min_len) {
                print current_header
                print current_seq
                ++kept
            } else {
                ++removed
            }
        }
        # Guardar nuevo header
        current_header = $0
        current_seq = ""
        next
    }
    { current_seq = current_seq $0 }
    END {
        # Procesar última secuencia
        if (current_header != "") {
            seq_length = length(current_seq)
            if (seq_length > min_len) {
                print current_header
                print current_seq
                ++kept
            } else {
                ++removed
            }
        }
        # Escribir estadísticas al archivo temporal
        print kept "," removed > stats_file
    }
    ' "$fasta_file" > "$temp_output_file"
    
    # awk ran in a subshell, so its counters come back through the temp file.
    if [ -f "$temp_stats" ]; then
        kept=$(cut -d',' -f1 "$temp_stats")
        removed=$(cut -d',' -f2 "$temp_stats")
        rm -f "$temp_stats"
    else
        kept=0
        removed=0
    fi
    
    log_msg "  Secuencias conservadas: $kept"
    log_msg "  Secuencias eliminadas: $removed"
    
    # Actualizar contadores totales
    total_seqs_kept_length=$((total_seqs_kept_length + kept))
    total_seqs_removed_length=$((total_seqs_removed_length + removed))
    processed=$((processed + 1))
    
    log_msg ""
done

log_msg "Paso 1 completado:"
log_msg "  Archivos procesados: $processed"
log_msg "  Secuencias conservadas: $total_seqs_kept_length"
log_msg "  Secuencias eliminadas por longitud: $total_seqs_removed_length"
log_msg ""

# ============================================================
# ---- step 2: marker completeness per MAG ----
# ============================================================

log_msg "========================================"
log_msg "PASO 2: Análisis de completitud de marcadores"
log_msg "========================================"

ALL_MAGS_FILE="${TEMP_DIR}/all_mags.txt"
: > "$ALL_MAGS_FILE"

for temp_file in "$TEMP_DIR"/*.fa; do
    [ -f "$temp_file" ] || continue
    grep "^>" "$temp_file" | sed 's/^>//' >> "$ALL_MAGS_FILE"
done

sort "$ALL_MAGS_FILE" | uniq > "${TEMP_DIR}/unique_mags.txt"
total_mags=$(wc -l < "${TEMP_DIR}/unique_mags.txt")

log_msg "MAGs found: $total_mags"
log_msg ""
log_msg "Computing marker completeness per MAG..."

echo -e "MAG\tMarkers_Present\tMarkers_Total\tCompleteness_%\tStatus" > "$MAG_REPORT"

MAGS_TO_KEEP="${TEMP_DIR}/mags_to_keep.txt"
MAGS_TO_REMOVE="${TEMP_DIR}/mags_to_remove.txt"
: > "$MAGS_TO_KEEP"
: > "$MAGS_TO_REMOVE"

mags_kept=0
mags_removed=0

while IFS= read -r mag_name; do
    # Contar en cuántos archivos aparece este MAG
    marker_count=0
    
    for temp_file in "$TEMP_DIR"/*.fa; do
        [ -f "$temp_file" ] || continue
        if grep -q "^>${mag_name}$" "$temp_file"; then
            ((++marker_count))
        fi
    done
    
    # Calcular porcentaje de completitud
    completeness=$(awk -v present=$marker_count -v total=$TOTAL_MARKERS 'BEGIN {printf "%.2f", (present/total)*100}')
    
    # Decidir si mantener o eliminar
    if [ "$marker_count" -ge "$MIN_MARKERS" ]; then
        echo "$mag_name" >> "$MAGS_TO_KEEP"
        echo -e "${mag_name}\t${marker_count}\t${TOTAL_MARKERS}\t${completeness}\tKEPT" >> "$MAG_REPORT"
        ((++mags_kept))
    else
        echo "$mag_name" >> "$MAGS_TO_REMOVE"
        echo -e "${mag_name}\t${marker_count}\t${TOTAL_MARKERS}\t${completeness}\tREMOVED" >> "$MAG_REPORT"
        ((++mags_removed))
    fi
    
done < "${TEMP_DIR}/unique_mags.txt"

log_msg "Análisis completado:"
log_msg "  MAGs con ≥${MIN_MARKER_PERCENT}% marcadores: $mags_kept (conservados)"
log_msg "  MAGs con <${MIN_MARKER_PERCENT}% marcadores: $mags_removed (eliminados)"
log_msg "  Reporte de completitud: $MAG_REPORT"
log_msg ""

# ============================================================
# PASO 3: Filtrar archivos removiendo MAGs incompletos
# ============================================================

log_msg "========================================"
log_msg "PASO 3: Generando archivos finales"
log_msg "========================================"

total_seqs_kept_final=0
total_seqs_removed_mag=0
files_with_seqs=0
empty_files=0

for temp_file in "$TEMP_DIR"/*.fa; do
    [ -f "$temp_file" ] || continue
    
    basename_file=$(basename "$temp_file")
    final_output="${OUTPUT_DIR}/${basename_file}"
    
    # Keep only the MAGs that cleared the completeness threshold.
    awk '
    BEGIN { 
        current_header=""; 
        current_seq=""; 
        keep_current=0;
        kept_seqs=0;
        removed_seqs=0;
    }
    /^>/ {
        # Procesar secuencia anterior
        if (current_header != "" && keep_current == 1) {
            print current_header
            print current_seq
            ++kept_seqs
        } else if (current_header != "") {
            ++removed_seqs
        }
        
        mag_name = substr($0, 2)
        keep_current = 0
        
        # Buscar en lista de MAGs a mantener
        while ((getline line < mags_keep_file) > 0) {
            if (line == mag_name) {
                keep_current = 1
                break
            }
        }
        close(mags_keep_file)
        
        current_header = $0
        current_seq = ""
        next
    }
    { current_seq = current_seq $0 }
    END {
        # Procesar última secuencia
        if (current_header != "" && keep_current == 1) {
            print current_header
            print current_seq
            ++kept_seqs
        } else if (current_header != "") {
            ++removed_seqs
        }
        print kept_seqs "," removed_seqs > stats_file
    }
    ' mags_keep_file="$MAGS_TO_KEEP" stats_file="${TEMP_DIR}/.final_stats_${basename_file}.tmp" \
      "$temp_file" > "$final_output"
    
    # Leer estadísticas
    stats_file="${TEMP_DIR}/.final_stats_${basename_file}.tmp"
    if [ -f "$stats_file" ]; then
        kept_final=$(cut -d',' -f1 "$stats_file")
        removed_final=$(cut -d',' -f2 "$stats_file")
        rm -f "$stats_file"
    else
        kept_final=0
        removed_final=0
    fi
    
    # A marker that lost every sequence would break the alignment step.
    if [ "$kept_final" -eq 0 ]; then
        rm -f "$final_output"
        ((++empty_files))
        log_msg "  ${basename_file}: sin secuencias después de filtrado (archivo no creado)"
    else
        ((++files_with_seqs))
        log_msg "  ${basename_file}: $kept_final secuencias conservadas, $removed_final removidas por MAG incompleto"
    fi
    
    total_seqs_kept_final=$((total_seqs_kept_final + kept_final))
    total_seqs_removed_mag=$((total_seqs_removed_mag + removed_final))
done

log_msg ""
log_msg "Paso 3 completado:"
log_msg "  Archivos con datos: $files_with_seqs"
log_msg "  Archivos vacíos (no creados): $empty_files"
log_msg "  Secuencias conservadas: $total_seqs_kept_final"
log_msg "  Secuencias removidas por MAG incompleto: $total_seqs_removed_mag"
log_msg ""

# ============================================================
# LIMPIEZA Y RESUMEN FINAL
# ============================================================

log_msg "========================================"
log_msg "Limpiando archivos temporales..."
log_msg "========================================"

rm -rf "$TEMP_DIR"
log_msg "Directorio temporal eliminado"
log_msg ""

# Calcular estadísticas finales
total_seqs_original=$((total_seqs_kept_length + total_seqs_removed_length))
total_seqs_removed=$((total_seqs_removed_length + total_seqs_removed_mag))

if [ $total_seqs_original -gt 0 ]; then
    final_retention=$(awk -v kept=$total_seqs_kept_final -v total=$total_seqs_original 'BEGIN {printf "%.2f%%", (kept/total)*100}')
else
    final_retention="N/A"
fi

log_msg "========================================"
log_msg "RESUMEN FINAL"
log_msg "========================================"
log_msg "Secuencias originales: $total_seqs_original"
log_msg ""
log_msg "Filtrado por longitud:"
log_msg "  Eliminadas: $total_seqs_removed_length"
log_msg ""
log_msg "Filtrado por completitud de MAGs:"
log_msg "  MAGs analizados: $total_mags"
log_msg "  MAGs conservados (≥${MIN_MARKER_PERCENT}%): $mags_kept"
log_msg "  MAGs eliminados (<${MIN_MARKER_PERCENT}%): $mags_removed"
log_msg "  Secuencias eliminadas: $total_seqs_removed_mag"
log_msg ""
log_msg "Resultado final:"
log_msg "  Secuencias conservadas: $total_seqs_kept_final"
log_msg "  Secuencias eliminadas (total): $total_seqs_removed"
log_msg "  Tasa de retención: $final_retention"
log_msg ""
log_msg "Archivos marcadores finales: $files_with_seqs / $total_files"
log_msg ""
log_msg "Resultados en: $OUTPUT_DIR"
log_msg "Reporte de MAGs: $MAG_REPORT"
log_msg "Log guardado en: $LOG_FILE"
log_msg "========================================"
