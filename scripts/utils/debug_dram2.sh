#!/usr/bin/env bash
# The variables inside the ssh strings below expand locally: the remote shell
# does not know which container or job is under inspection.
# shellcheck disable=SC2029

# Collects what is needed to diagnose a stalled DRAM2 Nextflow job on SLURM:
# the job accounting, the node it landed on, what is running there, the
# container it started, the Nextflow work directory and the command Nextflow
# issued.
#
# Written while chasing an mmseqs2 step that hung without failing. It reports
# and never changes anything, so it is safe to run against a live job. It needs
# ssh access to the compute nodes, which not every cluster grants.
#
# Usage: bash scripts/utils/debug_dram2.sh <parent_job_id> [stuck_subjob_id]

set -uo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

PARENT_JOB="${1:-}"
STUCK_JOB="${2:-$PARENT_JOB}"

if [[ -z "$PARENT_JOB" ]]; then
    echo "Usage: bash $0 <parent_job_id> [stuck_subjob_id]" >&2
    echo "Find the parent id with: squeue -u \$USER" >&2
    exit 1
fi

# The configuration only supplies the work directory and the DRAM2 checkout, so
# the script still runs without it.
load_project_config > /dev/null 2>&1 || true

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
NODE=$(squeue -j "$STUCK_JOB" -h -o %N)
if [[ -z "$NODE" ]]; then
    echo "Could not read the node from squeue; $STUCK_JOB may have finished."
    NODE=$(sacct -j "$STUCK_JOB" -n -o NodeList | head -1 | xargs)
    echo "Node according to sacct: $NODE"
fi
echo "Node: $NODE"
echo ""

if [[ -z "$NODE" ]]; then
    echo "ERROR: could not identify the node. Stopping." >&2
    exit 1
fi

# 3. processes on that node
echo "=== 3. processes on that node ==="
echo "--- mmseqs, docker, nextflow and slurm processes ---"
ssh "$NODE" "ps aux | grep -E '(mmseqs|docker|nextflow|slurm)' | grep -v grep" || echo "No matching processes"
echo ""

echo "--- Process tree for $USER ---"
ssh "$NODE" "pstree -p $USER | head -50" || echo "pstree unavailable"
echo ""

echo "--- Top 20 processes by CPU ---"
ssh "$NODE" "ps aux --sort=-%cpu | head -20"
echo ""

# 4. container detection
echo "=== 4. container detection ==="

echo "--- All running Docker containers ---"
ssh "$NODE" "docker ps" 2>&1
echo ""

echo "--- Running containers, compact ---"
CONTAINERS=$(ssh "$NODE" "docker ps --format 'ID={{.ID}} Image={{.Image}} Names={{.Names}}'")
echo "$CONTAINERS"
echo ""

# The container name is not stable, so try several ways to find it.
echo "--- Looking for the DRAM2 or mmseqs2 container ---"

# By image name
CONTAINER_ID=$(ssh "$NODE" "docker ps --filter ancestor=microbiomedata/nmdc_taxa_classify --format '{{.ID}}'" 2>/dev/null | head -1)

# Fallback: match on a pattern in the image name
if [[ -z "$CONTAINER_ID" ]]; then
    CONTAINER_ID=$(ssh "$NODE" "docker ps --format '{{.ID}} {{.Image}}' | grep -i 'taxa\|dram\|mmseqs' | head -1 | awk '{print \$1}'")
fi

# Last resort: the most recent container
if [[ -z "$CONTAINER_ID" ]]; then
    CONTAINER_ID=$(ssh "$NODE" "docker ps --format '{{.ID}}'" | head -1)
fi

if [[ -z "$CONTAINER_ID" ]]; then
    echo "No running Docker container found."
    echo "Checking for Singularity instead..."
    ssh "$NODE" "ps aux | grep singularity | grep -v grep"
    echo ""
    USING_DOCKER=false
else
    echo "Container found: $CONTAINER_ID"
    USING_DOCKER=true

    echo ""
    echo "--- Container detail ---"
    ssh "$NODE" "docker inspect $CONTAINER_ID --format 'Image: {{.Config.Image}}
Name: {{.Name}}
Status: {{.State.Status}}
CPUs: {{.HostConfig.NanoCpus}}
Memory: {{.HostConfig.Memory}}'"
    echo ""
fi

# 5. environment and the mmseqs2 command
if [[ "$USING_DOCKER" == true ]]; then
    echo "=== 5. environment inside the container ==="
    ssh "$NODE" "docker exec $CONTAINER_ID env | grep -E '(THREAD|CPU|OMP|MMSEQS|NUMEXPR|SLURM)' | sort"
    echo ""

    echo "=== Processes inside the container ==="
    ssh "$NODE" "docker exec $CONTAINER_ID ps aux"
    echo ""

    echo "=== Container logs, last 100 lines ==="
    ssh "$NODE" "docker logs $CONTAINER_ID --tail 100"
    echo ""
else
    echo "=== 5. Skipping container inspection, no Docker found ==="
    echo ""
fi

# 6. Nextflow work directory, wherever step 10 put it
echo "=== 6. Nextflow work directory ==="

NF_WORK="${NF_WORK:-${DRAM2_WORK_DIR:-${PROJECT_DIR:-$REPO_ROOT}/work/DRAM2}}"

echo "Looking for recent task directories in: $NF_WORK"

echo "--- The five most recently modified task directories ---"
find "$NF_WORK" -type f -name ".command.sh" -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -5 | while read -r timestamp path; do
    dir=$(dirname "$path")
    date=$(date -d "@${timestamp%.*}" '+%Y-%m-%d %H:%M:%S')
    echo "$date - $dir"
done
echo ""

# Newest one
LAST_WORK=$(find "$NF_WORK" -type f -name ".command.sh" -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -1 | cut -d' ' -f2- | xargs -r dirname)

if [[ -z "$LAST_WORK" ]]; then
    echo "No Nextflow work directory found"
    echo ""
else
    echo "Most recent task directory: $LAST_WORK"
    echo ""

    echo "--- Contents ---"
    ls -lh "$LAST_WORK"/
    echo ""

    echo "=== Command issued (.command.sh) ==="
    cat "$LAST_WORK"/.command.sh
    echo ""

    echo "=== stdout (.command.log), last 100 lines ==="
    if [[ -f "$LAST_WORK/.command.log" ]]; then
        tail -100 "$LAST_WORK"/.command.log
    else
        echo ".command.log is missing or empty"
    fi
    echo ""

    echo "=== stderr (.command.err) ==="
    if [[ -f "$LAST_WORK/.command.err" ]]; then
        cat "$LAST_WORK"/.command.err
    else
        echo ".command.err is missing or empty"
    fi
    echo ""

    echo "=== Run script (.command.run), first 50 lines ==="
    if [[ -f "$LAST_WORK/.command.run" ]]; then
        head -50 "$LAST_WORK"/.command.run
    else
        echo ".command.run does not exist"
    fi
    echo ""
fi

# 7. Nextflow configuration in use
echo "=== 7. Nextflow configuration in use ==="

config_files=()
if [[ -n "${DRAM2_REPO:-}" && -d "$DRAM2_REPO" ]]; then
    config_files+=("$DRAM2_REPO/nextflow.config" "$DRAM2_REPO"/conf/gom_archaea_*.config)
fi

echo "--- Process settings (threads, CPUs, memory) ---"
grep -n -E "cpus|threads|memory|THREAD" "${config_files[@]}" 2>/dev/null
echo ""

echo "--- Docker settings ---"
grep -n -E "docker\.|runOptions" "${config_files[@]}" 2>/dev/null
echo ""

# 8. the docker run command Nextflow issued
echo "=== 8. the docker run command Nextflow issued ==="
ssh "$NODE" "ps aux | grep 'docker run' | grep -v grep"
echo ""

# 9. resource usage
echo "=== 9. resource usage ==="
sstat -j "$STUCK_JOB" --format=JobID,MaxRSS,AveCPU,AveVMSize,NTasks 2>/dev/null || echo "$STUCK_JOB is not running, or sstat is unavailable"
echo ""

echo "Diagnostics finished."
