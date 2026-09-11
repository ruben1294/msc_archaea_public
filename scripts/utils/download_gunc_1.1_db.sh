#!/bin/bash
# Downloads GUNC v1.1 databases with resume support.
# Monkey-patches GUNC's built-in downloader (no resume) with a Python/requests
# implementation that uses HTTP Range headers to continue partial downloads.

# Activate the environment that provides gunc before running this.
conda activate GUNC_v1.1

# Set GUNC_DB in configs/server.env to this same path afterwards.
DB_DIR="${GUNC_DB:-${PROJECT_DIR}/databases/GUNC_v1.1.0}"
mkdir -p "$DB_DIR"

download_gunc_db() {
    local DB_NAME=$1
    echo "$(date '+%T') : Downloading GUNC DB: ${DB_NAME}"

    python - "$DB_NAME" "$DB_DIR" <<'EOF'
import sys
import os
import gunc.gunc_database as gd

db_name = sys.argv[1]
db_path = sys.argv[2]

def resumable_wget(url, dest):
    import requests
    import urllib3
    import time
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

    # Try HTTPS first, then HTTP fallback
    urls_to_try = [url, url.replace("https://", "http://", 1)]

    start_byte = os.path.getsize(dest) if os.path.exists(dest) else 0
    headers = {"Range": f"bytes={start_byte}-"} if start_byte > 0 else {}

    max_retries = 20
    wait_seconds = 300  # 5 min between retries (avoid server IP block)

    for attempt in range(1, max_retries + 1):
        for try_url in urls_to_try:
            print(f"  Attempt {attempt}/{max_retries}  URL: {try_url}")
            if start_byte > 0:
                print(f"  Resuming from byte {start_byte:,}")
            try:
                with requests.get(try_url, headers=headers, stream=True,
                                  verify=False, timeout=120) as r:
                    r.raise_for_status()
                    mode = "ab" if start_byte > 0 else "wb"
                    with open(dest, mode) as f:
                        for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                            f.write(chunk)
                print(f"  Download complete.")
                return  # success
            except Exception as e:
                print(f"  Failed ({type(e).__name__}): {e}")
                # Update start_byte in case partial progress was made
                start_byte = os.path.getsize(dest) if os.path.exists(dest) else 0
                headers = {"Range": f"bytes={start_byte}-"} if start_byte > 0 else {}

        if attempt < max_retries:
            print(f"  Waiting {wait_seconds}s before next attempt...")
            time.sleep(wait_seconds)

    raise RuntimeError(f"All {max_retries} download attempts failed for {url}")

gd.download_file = resumable_wget
gd.get_db(db_path, db_name)
EOF

    local EXIT=$?
    if [[ $EXIT -ne 0 ]]; then
        echo "$(date '+%T') : ERROR downloading ${DB_NAME} (exit ${EXIT})"
        return 1
    fi
    echo "$(date '+%T') : Done: ${DB_NAME}"
}

echo "$(date '+%T') : START $(date '+%Y-%m-%d')"

download_gunc_db "progenomes_3"
download_gunc_db "gtdb_214"

echo ""
echo "$(date '+%T') : Final files:"
ls -lh "${DB_DIR}"/*.dmnd 2>/dev/null || echo "  No .dmnd files found."
