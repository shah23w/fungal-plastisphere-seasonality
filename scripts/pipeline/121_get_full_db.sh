#!/bin/bash
# Download the uncapped Kraken2 PlusPF database, build 2026-06-26. Classification
# throughout uses the pre-built 16 GB capped index of the same build date. The
# uncapped index is needed only to supply the genome set for the fungal-only
# index built in 161_fungal_db_fetch.sh.
#
# Download only. Verification and extraction run as a separate batch job.
set -u
ARIA=$HOME/miniconda3/envs/assembly_env/bin/aria2c
D=$HOME/fpar/db
U=https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20260626.tar.gz
MU=https://genome-idx.s3.amazonaws.com/kraken/pluspf_20260626/pluspf.md5
F=k2_pluspf_20260626.tar.gz

[ -x "$ARIA" ] || { echo "FATAL: no aria2c at $ARIA"; exit 1; }
mkdir -p "$D"
echo "started: $(date)"
df -h "$D" | tail -1

curl -sL "$MU" -o "$D/pluspf.md5" && { echo "--- official md5 ---"; cat "$D/pluspf.md5"; }

for try in 1 2 3; do
  echo "--- download attempt $try/3  $(date +%H:%M)"
  "$ARIA" -x 8 -s 8 -k 1M -c --file-allocation=none \
          --max-tries=5 --retry-wait=10 --timeout=60 \
          --console-log-level=warn --summary-interval=0 \
          -d "$D" -o "$F" "$U" && break
  echo "    attempt $try failed"
done

SZ=$(stat -c%s "$D/$F" 2>/dev/null || echo 0)
echo "on disk: $SZ bytes ($(echo "scale=1; $SZ/1000000000" | bc 2>/dev/null || echo ?) GB)"
EXP=84800000000
if [ "$SZ" -gt 80000000000 ]; then echo "size plausible, verify and extract next"
else echo "SIZE TOO SMALL (expected ~84.8 GB) -- rerun this script, aria2 -c resumes"; fi
echo "finished: $(date)"
