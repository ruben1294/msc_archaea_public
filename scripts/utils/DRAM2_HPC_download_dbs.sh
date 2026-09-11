#!/bin/bash
set -e
set -u

export DRAM2_UUID="97ed64b9-dea0-4eb5-a7e0-0b50ac94e889"
export CICESE_HPC_UUID="262feade-f044-11f0-9f57-0affd850fd3b"
# Override DRAM2_DB_DIR to put the databases somewhere else. They are large.
export DRAM2_dbs="${DRAM2_DB_DIR:-${PROJECT_DIR}/databases/DRAM2}"

# Timestamped log file
LOG_FILE="DRAM2_download_$(date +%Y%m%d_%H%M%S).log"
SUMMARY_FILE="DRAM2_download_summary_$(date +%Y%m%d_%H%M%S).txt"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Iniciar log
log "=== Iniciando descarga de bases de datos DRAM2 ==="
log "DRAM2_UUID: $DRAM2_UUID"
log "CICESE_HPC_UUID: $CICESE_HPC_UUID"
log "Directorio destino: $DRAM2_dbs"
log ""

# Per-database outcome, reported at the end
declare -a DATABASES=("db_descriptions" "camper" "dbcan" "fegenie" "kofam" "metals" "methyl" "pfam" "sulfur" "uniref")
declare -a SUCCESS=()
declare -a FAILED=()

transfer_db() {
    local db_name=$1
    log "----------------------------------------"
    log "Iniciando transferencia: $db_name"
    
    if globus transfer \
        "$DRAM2_UUID:/databases/$db_name/" \
        "$CICESE_HPC_UUID:$DRAM2_dbs/$db_name/" \
        --recursive \
        --label "DRAM2_${db_name}_$(date +%Y%m%d)" 2>&1 | tee -a "$LOG_FILE"; then
        
        log "[OK] Transferencia enviada exitosamente: $db_name"
        SUCCESS+=("$db_name")
    else
        log "[FAIL] ERROR: Falló la transferencia de $db_name"
        FAILED+=("$db_name")
    fi
    log ""
}

# Transferir cada base de datos
for db in "${DATABASES[@]}"; do
    transfer_db "$db"
done

# Generar resumen
log "=== RESUMEN DE TRANSFERENCIAS ==="
log "Total de bases de datos: ${#DATABASES[@]}"
log "Transferencias exitosas: ${#SUCCESS[@]}"
log "Transferencias fallidas: ${#FAILED[@]}"
log ""

# Detalles de éxitos
if [ ${#SUCCESS[@]} -gt 0 ]; then
    log "Bases de datos transferidas exitosamente:"
    for db in "${SUCCESS[@]}"; do
        log "  [OK] $db"
    done
    log ""
fi

# Detalles de fallos
if [ ${#FAILED[@]} -gt 0 ]; then
    log "Bases de datos que fallaron:"
    for db in "${FAILED[@]}"; do
        log "  [FAIL] $db"
    done
    log ""
    log "ADVERTENCIA: Algunas transferencias fallaron. Revisa el log para más detalles."
else
    log "Todas las transferencias fueron enviadas exitosamente."
fi

log ""
log "Para monitorear el progreso de las transferencias, ejecuta:"
log "  globus task list"
log ""
log "Para ver detalles de una transferencia específica:"
log "  globus task show <TASK_ID>"
log ""
log "Archivo de log guardado en: $LOG_FILE"

# Human-readable summary alongside the log
cat > "$SUMMARY_FILE" << EOF
=================================================
RESUMEN DE DESCARGA DE BASES DE DATOS DRAM2
=================================================
Fecha: $(date '+%Y-%m-%d %H:%M:%S')
Directorio destino: $DRAM2_dbs

ESTADÍSTICAS:
-------------
Total de bases de datos: ${#DATABASES[@]}
Transferencias exitosas: ${#SUCCESS[@]}
Transferencias fallidas: ${#FAILED[@]}

BASES DE DATOS EXITOSAS:
-------------------------
EOF

if [ ${#SUCCESS[@]} -gt 0 ]; then
    for db in "${SUCCESS[@]}"; do
        echo "[OK] $db" >> "$SUMMARY_FILE"
    done
else
    echo "(ninguna)" >> "$SUMMARY_FILE"
fi

cat >> "$SUMMARY_FILE" << EOF

BASES DE DATOS FALLIDAS:
------------------------
EOF

if [ ${#FAILED[@]} -gt 0 ]; then
    for db in "${FAILED[@]}"; do
        echo "[FAIL] $db" >> "$SUMMARY_FILE"
    done
else
    echo "(ninguna)" >> "$SUMMARY_FILE"
fi

cat >> "$SUMMARY_FILE" << EOF

COMANDOS ÚTILES:
----------------
Monitorear todas las transferencias:
  globus task list

Ver detalles de una transferencia:
  globus task show <TASK_ID>

Cancelar una transferencia:
  globus task cancel <TASK_ID>

=================================================
EOF

log "Resumen guardado en: $SUMMARY_FILE"

# Mostrar resumen en consola
echo ""
cat "$SUMMARY_FILE"

# Exit code basado en éxito/fallo
if [ ${#FAILED[@]} -gt 0 ]; then
    exit 1
else
    exit 0
fi
