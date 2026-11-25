#!/bin/bash
# Per-run worker for the archived corpus. Fetches one run, checks gzip integrity
# and pairing, trims with fastp, classifies with Kraken2 against the capped
# PlusPF index at confidence 0.1 with three minimum hit groups, writes the
# per-run summary and deletes the reads.
#
# Launched as a job array, where the array index selects one line of the manifest.
# Each run writes its own result file rather than appending to a shared table,
# so parallel elements cannot race. Re-running skips any run that already has a
# completion marker.
set -eu
IDX=${LSB_JOBINDEX:?FATAL: not running as part of a job array}
M=$HOME/fpar/data/corpus625_manifest.tsv
EM=$HOME/fpar/data/ena_manifest.tsv
BIN=$HOME/miniconda3/envs/assembly_env/bin
DB=$HOME/fpar/db/k2pf16
O=$HOME/fpar/corpus
R=$HOME/fpar/results/corpus625
T=${LSB_DJOB_NUMPROC:-8}
mkdir -p "$O" "$R"

LINE=$(sed -n "$((IDX+1))p" "$M")
[ -n "$LINE" ] || { echo "FATAL: no manifest line at index $IDX"; exit 1; }
RUN=$(echo "$LINE" | cut -f1)
STUDY=$(echo "$LINE" | cut -f2)
ROLE=$(echo "$LINE" | cut -f3)

DONE="$R/${RUN}.done"
if [ -s "$DONE" ]; then
  echo "SKIP $RUN -- already done"
  exit 0
fi
rm -f "$R/${RUN}.FAILED"

echo "host: $(hostname)  run: $RUN  study: $STUDY  role: $ROLE  started: $(date)"

MLINE=$(awk -F'\t' -v r="$RUN" '$1==r{print; exit}' "$EM")
if [ -z "$MLINE" ]; then
  echo "$RUN not found in ena_manifest.tsv" > "$R/${RUN}.FAILED"
  exit 1
fi
URL1=$(echo "$MLINE" | cut -f3)
URL2=$(echo "$MLINE" | cut -f4)

D="$O/$RUN"
mkdir -p "$D"
T0=$(date +%s)

dl_one() {
  local url=$1 out=$2
  for try in 1 2 3; do
    "$BIN/aria2c" -x 4 -s 4 -c --file-allocation=none \
        --max-tries=5 --retry-wait=10 --timeout=60 \
        --console-log-level=warn --summary-interval=0 \
        -d "$(dirname "$out")" -o "$(basename "$out")" "$url" && return 0
    echo "  download attempt $try failed for $url"
  done
  return 1
}
F1="$D/${RUN}_1.fastq.gz"; F2="$D/${RUN}_2.fastq.gz"
if ! { [ -s "$F1" ] || dl_one "$URL1" "$F1"; }; then
  echo "download failed: R1" > "$R/${RUN}.FAILED"; rm -rf "$D"; exit 1
fi
if ! { [ -s "$F2" ] || dl_one "$URL2" "$F2"; }; then
  echo "download failed: R2" > "$R/${RUN}.FAILED"; rm -rf "$D"; exit 1
fi
T1=$(date +%s); DL_SEC=$((T1-T0))

if ! gzip -t "$F1" 2>/dev/null || ! gzip -t "$F2" 2>/dev/null; then
  echo "gzip integrity check failed" > "$R/${RUN}.FAILED"; rm -rf "$D"; exit 1
fi
N1=$(( $(zcat "$F1" | wc -l) / 4 ))
N2=$(( $(zcat "$F2" | wc -l) / 4 ))
if [ "$N1" != "$N2" ]; then
  echo "pairing mismatch N1=$N1 N2=$N2" > "$R/${RUN}.FAILED"; rm -rf "$D"; exit 1
fi

T0=$(date +%s)
"$BIN/fastp" -i "$F1" -I "$F2" \
      -o "$D/c1.fq.gz" -O "$D/c2.fq.gz" \
      -w $((T>16?16:T)) -q 20 -l 50 --detect_adapter_for_pe \
      -j "$D/fastp.json" -h "$D/fastp.html" 2> "$D/fastp.log" \
  || { echo "fastp failed" > "$R/${RUN}.FAILED"; rm -rf "$D"; exit 1; }
T1=$(date +%s); FASTP_SEC=$((T1-T0))

"$BIN/kraken2" --help 2>&1 | grep -q -- "--minimum-hit-groups" \
  || { echo "kraken2 lacks --minimum-hit-groups"; exit 1; }

T0=$(date +%s)
"$BIN/kraken2" --db "$DB" --threads "$T" --paired --gzip-compressed \
        --confidence 0.1 --minimum-hit-groups 3 \
        --report-minimizer-data \
        --report "$R/${RUN}.kreport" --output /dev/null \
        "$D/c1.fq.gz" "$D/c2.fq.gz" 2> "$R/${RUN}.kraken2.log" \
  || { echo "kraken2 failed" > "$R/${RUN}.FAILED"; rm -rf "$D"; exit 1; }
T1=$(date +%s); KRAKEN_SEC=$((T1-T0))

declare -A G=( [Aspergillus]=5052 [Candida]=5475 [Cryptococcus]=5206 \
               [Fusarium]=5506 [Trichosporon]=5552 [Exophiala]=5583 \
               [Scedosporium]=41687 [Lomentospora]=1549750 )
K="$R/${RUN}.kreport"
get() { awk -v t="$1" '$7==t{print $2; exit}' "$K"; }
CLS=$(get 1); CLS=${CLS:-0}
UNC=$(get 0); UNC=${UNC:-0}
FUN=$(get 4751); FUN=${FUN:-0}
MZC=$(awk '$7==4751{print $4; exit}' "$K"); MZC=${MZC:-0}
row="$RUN\t$STUDY\t$ROLE\t$((CLS+UNC))\t$CLS\t$FUN"
for g in Aspergillus Candida Cryptococcus Fusarium Trichosporon Exophiala Scedosporium Lomentospora; do
  v=$(get "${G[$g]}"); row="$row\t${v:-0}"
done
printf "$row\t$MZC\n" > "$R/${RUN}.summary.tsv"

rm -rf "$D"

TOTAL_SEC=$(( DL_SEC + FASTP_SEC + KRAKEN_SEC ))
printf "run=%s dl=%ss fastp=%ss kraken2=%ss total=%ss finished=%s\n" \
  "$RUN" "$DL_SEC" "$FASTP_SEC" "$KRAKEN_SEC" "$TOTAL_SEC" "$(date)" > "$DONE"
echo "DONE $RUN  total=${TOTAL_SEC}s"
