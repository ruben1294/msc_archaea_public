#!/bin/bash
#SBATCH --job-name=DRAM2_install
#SBATCH --account=metagenomica
#SBATCH --partition=cicese
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nodelist=nodo5
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.log

# One-off setup for step 8: fetches the DRAM2 Nextflow pipeline and checks that
# Nextflow and Docker both work. Submit from the repository root:
#
#   mkdir -p logs && sbatch slurm/DRAM2_HPC_install.slurm
#
# The --nodelist line is not incidental. On this cluster only nodo5 reaches the
# internet; the compute nodes do not. So the pull, and the plugin download it
# triggers, have to happen here, once. Every later DRAM2 run sets NXF_OFFLINE
# and works from what this leaves in ~/.nextflow.
#
# The plugins do not need installing separately. nextflow pull fetches them, and
# declaring them in the Nextflow config is what makes the offline runs find them.
#
# Account, partition and nodelist are specific to the CICESE cluster.

set -euo pipefail

REPO_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "$REPO_ROOT/logs"

echo "DRAM2 install and Docker check"
echo "  Job ID: ${SLURM_JOB_ID:-NA}"
echo "  Node:   $(hostname)"
echo "  Start:  $(date)"

if [[ -f "$REPO_ROOT/slurm/cluster.env" ]]; then
    # shellcheck source=/dev/null
    source "$REPO_ROOT/slurm/cluster.env"
fi

if [[ -n "${CLUSTER_MODULE:-}" ]]; then
    module load "$CLUSTER_MODULE"
fi
if [[ -n "${CONDA_PROFILE:-}" && -f "$CONDA_PROFILE" ]]; then
    # shellcheck source=/dev/null
    source "$CONDA_PROFILE"
fi
conda activate "${DRAM2_CONDA_ENV:-DRAMV2}"

if ! command -v nextflow > /dev/null 2>&1; then
    echo "[FAIL] nextflow not found" >&2
    exit 1
fi
echo "[OK] nextflow"
nextflow -version

if ! command -v docker > /dev/null 2>&1; then
    echo "[FAIL] docker not found on $(hostname)" >&2
    exit 1
fi
echo "[OK] docker"
docker --version

echo "Pulling the DRAM2 pipeline..."
nextflow pull WrightonLabCSU/DRAM -r dev
nextflow info WrightonLabCSU/DRAM

echo "Docker check:"
docker info
docker images
docker pull hello-world
docker run hello-world

echo "Nextflow plugins now available offline:"
ls ~/.nextflow/plugins/

echo "Done. $(date)"
