#!/bin/bash
# Classify the spiked libraries with the same settings used for the corpus, so
# recovery is measured through the pipeline that produced every other result.
set -u
K2=$HOME/miniconda3/envs/assembly_env/bin/kraken2
DB=$HOME/fpar/db/k2pf16
W=$HOME/fpar/work_spike/libs
OUT=$HOME/fpar/results/spikein_kreports
mkdir -p "$OUT"
i=${LSB_JOBINDEX:-1}
f=$(ls -1 "$W"/*_L*_1.fq.gz | sed -n "${i}p")
[ -n "$f" ] || { echo "no file for index $i"; exit 1; }
tag=$(basename "$f"); tag=${tag%_1.fq.gz}
r2="${f%_1.fq.gz}_2.fq.gz"
# validity guard, not an existence guard: a kreport without a root line is junk
if [ -s "$OUT/$tag.kreport" ] && grep -qP "\tR\t" "$OUT/$tag.kreport" 2>/dev/null; then
  echo "$tag already classified and valid, skipping"; exit 0
fi
echo "[$(date +%H:%M)] $tag start"
"$K2" --db "$DB" --threads 4 --paired --gzip-compressed \
  --confidence 0.1 --minimum-hit-groups 3 --report-minimizer-data \
  --report "$OUT/$tag.kreport" --output /dev/null "$f" "$r2" 2> "$OUT/$tag.log"
rc=$?
echo "[$(date +%H:%M)] $tag kraken2 exit $rc"
grep -qP "\tR\t" "$OUT/$tag.kreport" && echo "$tag kreport OK" || echo "$tag kreport INVALID"
tail -3 "$OUT/$tag.log"
