#!/bin/bash
# Downloads the pre-built skani GTDB-R226 sketch database with resume support.
# Output goes under $GTDB_SKETCH_DIR, or the project databases directory.

URL="http://faust.compbio.cs.cmu.edu/skani-files/skani_gtdb_r226-v0.3.tar.gz"
DB_DIR="${GTDB_SKETCH_DIR:-${PROJECT_DIR}/databases/GTDB/r226}"
TAR_FILE="${DB_DIR}/skani_gtdb_r226-v0.3.tar.gz"

mkdir -p "$DB_DIR"

echo "$(date '+%T') : START $(date '+%Y-%m-%d')"
echo "$(date '+%T') : Destination: $DB_DIR"

# ── Download with resume support ─────────────────────────────────────────────
python3 - "$URL" "$TAR_FILE" <<'EOF'
import sys, os, time, requests, urllib3

url, dest = sys.argv[1], sys.argv[2]
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

max_retries = 20
wait_seconds = 300  # 5 min between retries

for attempt in range(1, max_retries + 1):
    start_byte = os.path.getsize(dest) if os.path.exists(dest) else 0
    headers = {"Range": f"bytes={start_byte}-"} if start_byte > 0 else {}
    print(f"  Attempt {attempt}/{max_retries}", flush=True)
    if start_byte > 0:
        print(f"  Resuming from byte {start_byte:,}", flush=True)
    try:
        with requests.get(url, headers=headers, stream=True,
                          verify=False, timeout=120) as r:
            r.raise_for_status()
            mode = "ab" if start_byte > 0 else "wb"
            with open(dest, mode) as f:
                for chunk in r.iter_content(chunk_size=8 * 1024 * 1024):
                    f.write(chunk)
        print("  Download complete.", flush=True)
        sys.exit(0)
    except Exception as e:
        print(f"  Failed ({type(e).__name__}): {e}", flush=True)
        if attempt < max_retries:
            print(f"  Waiting {wait_seconds}s before next attempt...", flush=True)
            time.sleep(wait_seconds)

print(f"ERROR: all {max_retries} attempts failed.", flush=True)
sys.exit(1)
EOF

if [[ $? -ne 0 ]]; then
    echo "$(date '+%T') : ERROR: download failed."
    exit 1
fi

echo "$(date '+%T') : Download complete: $(du -sh "$TAR_FILE" | cut -f1)"

# ── Extract ───────────────────────────────────────────────────────────────────
echo "$(date '+%T') : Extracting..."
tar -zxvf "$TAR_FILE" -C "$DB_DIR"

echo "$(date '+%T') : Extracted contents:"
ls -lh "$DB_DIR"/skani_gtdb_r226-v0.3/ 2>/dev/null || ls -lh "$DB_DIR"/

echo ""
echo "$(date '+%T') : Done."
echo "$(date '+%T') : Update GTDB_SKETCH_DIR in configs/server.env to:"
echo "  ${DB_DIR}/skani_gtdb_r226-v0.3"
