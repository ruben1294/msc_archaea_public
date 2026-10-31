#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/../lib/common.sh"

# Downloads the pre-built skani sketches of GTDB R226 that step 04 searches the
# MAGs against, and unpacks them under GTDB_SKETCH_DIR. Needs python3 with
# requests, which most conda environments include:
#
#   bash scripts/utils/download_skani_gtdb_r226.sh
#
# The archive is large and the server slow, so a dropped connection resumes
# with an HTTP Range request instead of starting over. GTDB_SKETCH_DIR in
# configs/server.env is the directory the sketches end up in; the archive is
# kept next to it, in its parent directory, until the extraction succeeds.

load_project_config
require_vars GTDB_SKETCH_DIR
require_commands python3 tar

URL="https://faust.compbio.cs.cmu.edu/skani-files/skani_gtdb_r226-v0.3.tar.gz"
PARENT_DIR="$(dirname "$GTDB_SKETCH_DIR")"
TAR_FILE="${PARENT_DIR}/skani_gtdb_r226-v0.3.tar.gz"

if [[ -d "$GTDB_SKETCH_DIR" && -n "$(ls -A "$GTDB_SKETCH_DIR" 2> /dev/null)" ]]; then
    log_info "Sketches already in $GTDB_SKETCH_DIR, nothing to do"
    exit 0
fi

mkdir -p "$PARENT_DIR"
log_info "Downloading $URL"
log_info "  to $TAR_FILE"

python3 - "$URL" "$TAR_FILE" <<'EOF'
import os
import sys
import time

import requests

url, dest = sys.argv[1], sys.argv[2]
MAX_RETRIES = 20
WAIT_SECONDS = 300

for attempt in range(1, MAX_RETRIES + 1):
    start = os.path.getsize(dest) if os.path.exists(dest) else 0
    headers = {"Range": f"bytes={start}-"} if start else {}
    print(f"  Attempt {attempt}/{MAX_RETRIES}", flush=True)
    if start:
        print(f"  Resuming from byte {start:,}", flush=True)
    try:
        with requests.get(url, headers=headers, stream=True, timeout=120) as r:
            if start and r.status_code == 416:
                print("  Already complete.", flush=True)
                sys.exit(0)
            r.raise_for_status()
            # A 200 to a ranged request means the range was ignored: restart.
            mode = "ab" if start and r.status_code == 206 else "wb"
            with open(dest, mode) as f:
                for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                    f.write(chunk)
        print("  Download complete.", flush=True)
        sys.exit(0)
    except requests.exceptions.RequestException as e:
        print(f"  Failed ({type(e).__name__}): {e}", flush=True)
        if attempt < MAX_RETRIES:
            print(f"  Waiting {WAIT_SECONDS}s before the next attempt", flush=True)
            time.sleep(WAIT_SECONDS)

sys.exit(f"All {MAX_RETRIES} attempts failed for {url}")
EOF

log_info "Downloaded $(du -sh "$TAR_FILE" | cut -f1); testing the archive"
# A truncated or corrupted download fails here, before anything is extracted.
tar -tzf "$TAR_FILE" > /dev/null

log_info "Extracting"
tar -xzf "$TAR_FILE" -C "$PARENT_DIR"

if [[ ! -d "$GTDB_SKETCH_DIR" ]]; then
    log_error "The archive did not unpack into $GTDB_SKETCH_DIR. It holds:"
    tar -tzf "$TAR_FILE" | head -5 >&2
    log_error "Set GTDB_SKETCH_DIR in configs/server.env to that directory."
    exit 1
fi

rm -f "$TAR_FILE"
log_info "Sketches ready in $GTDB_SKETCH_DIR"
