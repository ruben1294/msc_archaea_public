#!/bin/bash
#SBATCH --job-name=23_pangenome_misc_data_layers
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Metadata layers. Submit from the repository root, naming your account and
# partition:
#
#   mkdir -p logs
#   sbatch --account=<account> --partition=<partition> slurm/23_pangenome_misc_data_layers.slurm
#
# Exporting SBATCH_ACCOUNT and SBATCH_PARTITION once does the same for every
# wrapper. The module and conda setup come from slurm/cluster.env, and
# pipeline/23_pangenome_misc_data_layers.sh runs on any machine with its tools
# on PATH. The resource requests are the ones the thesis runs used; adjust them
# to your data.

set -euo pipefail

# sbatch runs a copy of this file, so the script cannot locate itself.
# REPO_ROOT comes from SLURM_SUBMIT_DIR, the directory sbatch was called from,
# which has to be the repository root.
REPO_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "$REPO_ROOT/logs"

SCRIPT_PATH="$REPO_ROOT/scripts/pipeline/23_pangenome_misc_data_layers.sh"

if [[ ! -f "$SCRIPT_PATH" ]]; then
    echo "Script not found: $SCRIPT_PATH" >&2
    echo "Submit from the repository root, not from slurm/." >&2
    exit 1
fi

# Cluster-specific module and conda setup. Copy slurm/cluster.env.example to
# slurm/cluster.env and edit it; without it, conda has to be on PATH already.
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
conda activate anvio-9

echo "Running $SCRIPT_PATH"

# Capture the exit code instead of letting set -e abort, so the closing
# message can report it.
exit_code=0
bash "$SCRIPT_PATH" || exit_code=$?

if [[ $exit_code -eq 0 ]]; then
    echo "[OK] Metadata layers finished"
else
    echo "[FAIL] Metadata layers exited $exit_code" >&2
fi

exit $exit_code
