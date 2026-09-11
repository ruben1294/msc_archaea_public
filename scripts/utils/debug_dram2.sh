#!/usr/bin/env bash
# The variables inside the ssh strings below expand on this side on purpose: the
# remote shell has no idea which container or job is under inspection.
# shellcheck disable=SC2029

# Collects everything useful about a stalled DRAM2 Nextflow job on SLURM: the
# job accounting, which node it landed on, what is running there, the container
# it started, the Nextflow work directory and the command that was actually
# issued.
#
# Written while chasing an mmseqs2 step that hung without failing. It reports
# and never changes anything, so it is safe to run against a live job.
#
# Usage: bash debug_dram2.sh <parent_job_id> [stuck_subjob_id]

set -uo pipefail

PARENT_JOB="${1:-}"
STUCK_JOB="${2:-$PARENT_JOB}"

if [[ -z "$PARENT_JOB" ]]; then
    echo "Usage: bash $0 <parent_job_id> [stuck_subjob_id]" >&2
    echo "Find the parent id with: squeue -u \$USER" >&2
    exit 1
fi

echo "DRAM2 and Nextflow diagnostics"
echo "Parent job: $PARENT_JOB   Subjob under inspection: $STUCK_JOB"
echo ""

# 1. SLURM job accounting
echo "=== 1. SLURM job accounting ==="
echo ""
sacct -j "$PARENT_JOB","$STUCK_JOB" --format=JobID,JobName,State,Elapsed,AllocCPUS,MaxRSS,MaxVMSize,NodeList
echo ""

echo "=== All subjobs of $PARENT_JOB ==="
sacct -j "$PARENT_JOB" --format=JobID,JobName,State,Elapsed,AllocCPUS,MaxRSS,NodeList
echo ""

echo "=== Detail for $STUCK_JOB ==="
scontrol show job "$STUCK_JOB"
echo ""

# 2. execution node
echo "=== 2. execution node ==="
NODO=$(squeue -j "$STUCK_JOB" -h -o %N)
if [ -z "$NODO" ]; then
    echo "Could not read the node from squeue; $STUCK_JOB may have finished."
    NODO=$(sacct -j "$STUCK_JOB" -n -o NodeList | head -1)
    echo "Intentando obtener nodo de sacct: $NODO"
fi
echo "Nodo: $NODO"
echo ""

if [ -z "$NODO" ]; then
    echo "ERROR CRÍTICO: No se pudo identificar el nodo. Abortando."
    exit 1
fi

# 3. processes on that node
echo "=== 3. processes on that node ==="
echo "--- Procesos relacionados con mmseqs, docker, nextflow ---"
ssh "$NODO" "ps aux | grep -E '(mmseqs|docker|nextflow|slurm)' | grep -v grep" || echo "No matching processes"
echo ""

echo "--- Árbol de procesos del usuario ---"
ssh "$NODO" "pstree -p $USER | head -50" || echo "pstree unavailable"
echo ""

echo "--- Top 20 procesos por CPU ---"
ssh "$NODO" "ps aux --sort=-%cpu | head -20"
echo ""

# 4. container detection
echo "=== 4. container detection ==="

echo "--- All running Docker containers ---"
ssh "$NODO" "docker ps" 2>&1
echo ""

echo "--- Listando containeres Docker (formato compacto) ---"
CONTAINERS=$(ssh "$NODO" "docker ps --format 'ID={{.ID}} Image={{.Image}} Names={{.Names}}'")
echo "$CONTAINERS"
echo ""

# The container name is not stable, so try several ways to find it.
echo "--- Intentando detectar container de DRAM2/MMSeqs2 ---"

# By image name
CONTAINER_ID=$(ssh "$NODO" "docker ps --filter ancestor=microbiomedata/nmdc_taxa_classify --format '{{.ID}}'" 2>/dev/null | head -1)

# Fallback: match on a pattern in the image name
if [ -z "$CONTAINER_ID" ]; then
    CONTAINER_ID=$(ssh "$NODO" "docker ps --format '{{.ID}} {{.Image}}' | grep -i 'taxa\|dram\|mmseqs' | head -1 | awk '{print \$1}'")
fi

# Last resort: the user's most recent container
if [ -z "$CONTAINER_ID" ]; then
    CONTAINER_ID=$(ssh "$NODO" "docker ps --format '{{.ID}}'" | head -1)
fi

if [ -z "$CONTAINER_ID" ]; then
    echo "No running Docker container found."
    echo "Checking for Singularity instead..."
    ssh "$NODO" "ps aux | grep singularity | grep -v grep"
    echo ""
    USING_DOCKER=false
else
    echo "Container detectado: $CONTAINER_ID"
    USING_DOCKER=true
    
    echo ""
    echo "--- Container detail ---"
    ssh "$NODO" "docker inspect $CONTAINER_ID --format 'Image: {{.Config.Image}}
Name: {{.Name}}
Status: {{.State.Status}}
CPUs: {{.HostConfig.NanoCpus}}
Memory: {{.HostConfig.Memory}}'"
    echo ""
fi

# 5. environment and the mmseqs2 command
if [ "$USING_DOCKER" = true ]; then
    echo "=== 5. VARIABLES DE ENTORNO EN CONTENEDOR ==="
    ssh "$NODO" "docker exec $CONTAINER_ID env | grep -E '(THREAD|CPU|OMP|MMSEQS|NUMEXPR|SLURM)' | sort"
    echo ""
    
    echo "=== COMANDO MMSEQS2 EJECUTÁNDOSE ==="
    ssh "$NODO" "docker exec $CONTAINER_ID ps aux"
    echo ""
    
    echo "=== Container logs, last 100 lines ==="
    ssh "$NODO" "docker logs $CONTAINER_ID --tail 100"
    echo ""
else
    echo "=== 5. Skipping container inspection, no Docker found ==="
    echo ""
fi

# 6. Nextflow work directory
echo "=== 6. Nextflow work directory ==="

NF_WORK="${NF_WORK:-${HOME}/.nextflow/assets/WrightonLabCSU/DRAM/work}"

echo "Looking for directorios de trabajo recientes en: $NF_WORK"

# The five most recent work directories
echo "--- Últimos 5 directorios de trabajo modificados ---"
find "$NF_WORK" -type f -name ".command.sh" -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -5 | while read -r timestamp path; do
    dir=$(dirname "$path")
    date=$(date -d "@${timestamp%.*}" '+%Y-%m-%d %H:%M:%S')
    echo "$date - $dir"
done
echo ""

# Newest one
LAST_WORK=$(find "$NF_WORK" -type f -name ".command.sh" -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -1 | cut -d' ' -f2- | xargs dirname)

if [ -z "$LAST_WORK" ]; then
    echo "No Nextflow work directory found"
    echo ""
else
    echo "Most recent work directory: $LAST_WORK"
    echo ""
    
    echo "--- Contenido del directorio ---"
    ls -lh "$LAST_WORK"/
    echo ""
    
    echo "=== COMANDO EJECUTADO (.command.sh) ==="
    cat "$LAST_WORK"/.command.sh
    echo ""
    
    echo "=== stdout (.command.log), last 100 lines ==="
    if [ -f "$LAST_WORK/.command.log" ]; then
        tail -100 "$LAST_WORK"/.command.log
    else
        echo ".command.log is missing or empty"
    fi
    echo ""
    
    echo "=== STDERR (.command.err) ==="
    if [ -f "$LAST_WORK/.command.err" ]; then
        cat "$LAST_WORK"/.command.err
    else
        echo ".command.err is missing or empty"
    fi
    echo ""
    
    echo "=== Run script (.command.run), first 50 lines ==="
    if [ -f "$LAST_WORK/.command.run" ]; then
        head -50 "$LAST_WORK"/.command.run
    else
        echo "Archivo .command.run no existe"
    fi
    echo ""
fi

# 7. active Nextflow configuration
echo "=== 7. active Nextflow configuration ==="

echo "--- Configuration de proceso (threads, CPUs) ---"
grep -n -E "cpus|threads|memory|THREAD" *.config configs/*.config 2>/dev/null | grep -v "^Binary"
echo ""

echo "--- Configuration de Docker ---"
grep -n -E "docker\.|runOptions" *.config configs/*.config 2>/dev/null | grep -v "^Binary"
echo ""

# 8. the docker run command Nextflow issued
echo "=== 8. the docker run command Nextflow issued ==="
ssh "$NODO" "ps aux | grep 'docker run' | grep -v grep"
echo ""

# 9. resource usage
echo "=== 9. resource usage ==="
sstat -j "$STUCK_JOB" --format=JobID,MaxRSS,AveCPU,AveVMSize,NTasks 2>/dev/null || echo "$STUCK_JOB is not running, or sstat is unavailable"
echo ""

echo "================================================================================"
echo "                          FIN DEL DEBUGGING"
echo "================================================================================"
