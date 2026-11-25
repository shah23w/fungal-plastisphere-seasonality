#!/bin/bash
# Classify the nine fungus-free mock libraries under five configurations, which
# separate the contribution of reference scope from that of the confidence
# threshold. Every fungal read in these libraries is a false positive.
#
#   cell 1  fungal-only index, confidence 0.0, min hit groups 2
#   cell 2  fungal-only index, confidence 0.1, min hit groups 2
#   cell 3  capped PlusPF,     confidence 0.0, min hit groups 2
#   cell 4  capped PlusPF,     confidence 0.1, min hit groups 2
#   cell 5  capped PlusPF,     confidence 0.1, min hit groups 3
#
# Cell 5 is the setting used for every other analysis in the paper.
set -u
K2=$HOME/miniconda3/envs/assembly_env/bin/kraken2
DBF=$HOME/fpar/db/k2fungi
DBC=$HOME/fpar/db/k2pf16
M=$HOME/fpar/mock
OUT=$HOME/fpar/results/mock2x2
mkdir -p "$OUT"
RUNS=(SRR11075447 ERR10009752 SRR21461810 DRR510841 DRR188001 SRR14759986 SRR7762337 ERR3727302 ERR11260221)
i=${LSB_JOBINDEX:-1}
li=$(( (i-1) % 9 )); ci=$(( (i-1) / 9 + 1 ))
run=${RUNS[$li]}
case $ci in
  1) db=$DBF; conf=0.0; mhg=2 ;;
  2) db=$DBF; conf=0.1; mhg=2 ;;
  3) db=$DBC; conf=0.0; mhg=2 ;;
  4) db=$DBC; conf=0.1; mhg=2 ;;
  5) db=$DBC; conf=0.1; mhg=3 ;;
  *) echo "bad cell $ci"; exit 1 ;;
esac
tag="${run}.cell${ci}"
r1="$M/$run/c_1.fq.gz"; r2="$M/$run/c_2.fq.gz"
[ -s "$r1" ] && [ -s "$r2" ] || { echo "MISSING reads for $run"; exit 1; }
if [ -s "$OUT/$tag.kreport" ] && grep -qP "\tR\t" "$OUT/$tag.kreport" 2>/dev/null; then
  echo "$tag already valid, skipping"; exit 0
fi
echo "[$(date +%H:%M)] $tag  db=$(basename $db) conf=$conf mhg=$mhg"
"$K2" --db "$db" --threads 4 --paired --gzip-compressed \
  --confidence "$conf" --minimum-hit-groups "$mhg" --report-minimizer-data \
  --report "$OUT/$tag.kreport" --output /dev/null "$r1" "$r2" 2> "$OUT/$tag.log"
echo "[$(date +%H:%M)] $tag exit $?"
grep -qP "\tR\t" "$OUT/$tag.kreport" && echo "$tag kreport OK" || echo "$tag kreport INVALID"
