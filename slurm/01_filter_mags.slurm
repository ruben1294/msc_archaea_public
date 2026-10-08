#!/bin/bash
#SBATCH --job-name=01_filter_mags
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# MIMAG filtering, between step 01 and utils/05_rename_mags.sh. Submit from the
# repository root, naming your account and partition:
#
#   mkdir -p logs
#   sbatch --account=<account> --partition=<partition> slurm/01_filter_mags.slurm
#
# Exporting SBATCH_ACCOUNT and SBATCH_PARTITION once does the same for every
# wrapper. Nothing here is tied to one cluster: the module and conda setup come
# from slurm/cluster.env, and utils/01_filter_mags.py itself runs anywhere its
# tools are on PATH. Resource requests are the ones these runs actually used, a
# starting point rather than a recommendation.

set -euo pipefail

# sbatch copies this file before running it, so its original path is not
# recoverable at run time. SLURM_SUBMIT_DIR is wherever sbatch was invoked,
# which is why submitting from the repository root matters.
REPO_ROOT="${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "$REPO_ROOT/logs"

SCRIPT_PATH="$REPO_ROOT/scripts/utils/01_filter_mags.py"

if [[ ! -f "$SCRIPT_PATH" ]]; then
    echo "Script not found: $SCRIPT_PATH" >&2
    echo "Submit from the repository root, not from slurm/." >&2
    exit 1
fi

# The filter takes its paths and MIMAG thresholds from configs/server.env.
# shellcheck source=/dev/null
source "$REPO_ROOT/scripts/lib/common.sh"
load_project_config

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
