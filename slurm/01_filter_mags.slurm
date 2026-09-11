#!/bin/bash
#SBATCH --job-name=01_filter_mags
#SBATCH --account=metagenomica
#SBATCH --partition=d30
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# MIMAG filtering. Submit from the repository root:
#
#   mkdir -p logs && sbatch slurm/01_filter_mags.slurm
#
# The account, partition and module lines above and below are specific to the
# CICESE cluster. Edit them for yours; utils/01_filter_mags.py itself needs
# none of them and runs anywhere the tools are on PATH.
#
# Resource requests are the ones these runs actually used, so they are a
# starting point rather than a recommendation.

set -euo pipefail

# sbatch copies this file before running it, so its original path is not
# recoverable at run time. SLURM_SUBMIT_DIR is wherever sbatch was invoked,
# which is why submitting from the repository root matters.
REPO_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "$REPO_ROOT/logs"

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
conda activate checkm2

SCRIPT_PATH="$REPO_ROOT/scripts/utils/01_filter_mags.py"

if [[ ! -f "$SCRIPT_PATH" ]]; then
    echo "Script not found: $SCRIPT_PATH" >&2
    echo "Submit from the repository root, not from slurm/." >&2
    exit 1
fi

echo "Running $SCRIPT_PATH"

# set -e would abort before an exit code could be read, so it is captured here
# to let the closing message say what happened.
exit_code=0
python3 "$SCRIPT_PATH" || exit_code=$?

if [[ $exit_code -eq 0 ]]; then
    echo "[OK] MIMAG filtering finished"
else
    echo "[FAIL] MIMAG filtering exited $exit_code" >&2
fi

exit $exit_code
