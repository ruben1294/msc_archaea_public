#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Downloads the two GUNC 1.1 reference databases that step 02 runs against,
# ProGenomes3 and GTDB r214, into GUNC_DB. Run it inside the GUNC_v1.1
# environment of env/gunc.yml:
#
#   conda activate GUNC_v1.1
#   bash scripts/utils/download_gunc_1.1_db.sh
#
# GUNC's own downloader cannot resume, and the files are large enough that a
# dropped connection is likely, so this swaps in one that continues a partial
# file with an HTTP Range request. GUNC then checks each file against the md5
# sum the server publishes and deletes it on a mismatch, so a resumed download
# that went wrong is caught rather than used.
#
# A database already in place is skipped, which makes the script safe to rerun
# after an interruption.

load_project_config
require_vars GUNC_DB
require_commands python3

if ! python3 -c "import gunc.gunc_database" 2> /dev/null; then
    log_error "GUNC is not importable from python3. Activate the GUNC_v1.1 environment first."
    exit 1
fi

mkdir -p "$GUNC_DB"

download_gunc_db() {
    local db_name=$1
    local dmnd=$2

    if [[ -s "${GUNC_DB}/${dmnd}" ]]; then
        log_info "${dmnd} already in ${GUNC_DB}, skipping"
        return 0
    fi

    log_info "Downloading GUNC database ${db_name} into ${GUNC_DB}"

    python3 - "$db_name" "$GUNC_DB" <<'EOF'
import os
import sys
import time

import requests
import gunc.gunc_database as gd

MAX_RETRIES = 20
WAIT_SECONDS = 300  # between attempts, so the server does not block the address


def resumable_download(url, dest):
    for attempt in range(1, MAX_RETRIES + 1):
        start = os.path.getsize(dest) if os.path.exists(dest) else 0
        headers = {"Range": f"bytes={start}-"} if start else {}
        print(f"  Attempt {attempt}/{MAX_RETRIES}: {url}", flush=True)
        if start:
            print(f"  Resuming from byte {start:,}", flush=True)
        try:
            with requests.get(url, headers=headers, stream=True, timeout=120) as r:
                if start and r.status_code == 416:
                    # Nothing past what is on disk; the md5 check decides
                    # whether that is the whole file.
                    return
                r.raise_for_status()
                # A 200 to a ranged request means the server ignored the range
                # and is sending the whole file, so start over, never append.
                mode = "ab" if start and r.status_code == 206 else "wb"
                with open(dest, mode) as f:
                    for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                        f.write(chunk)
            return
        except requests.exceptions.RequestException as e:
            print(f"  Failed ({type(e).__name__}): {e}", flush=True)
            if attempt < MAX_RETRIES:
                print(f"  Waiting {WAIT_SECONDS}s before the next attempt", flush=True)
                time.sleep(WAIT_SECONDS)
    sys.exit(f"All {MAX_RETRIES} attempts failed for {url}")


gd.download_file = resumable_download
gd.get_db(sys.argv[2], sys.argv[1])
EOF

    log_info "Done: ${GUNC_DB}/${dmnd}"
}

log_info "GUNC databases into ${GUNC_DB}"

download_gunc_db "progenomes_3" "gunc_db_progenomes3.dmnd"
download_gunc_db "gtdb_214" "gunc_db_gtdb214.dmnd"

ls -lh "${GUNC_DB}"/*.dmnd
