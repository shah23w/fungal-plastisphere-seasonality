#!/bin/bash
# Build the fungus-free mock community used to measure the false-positive floor.
# Total depth is fixed at 20,000,000 read pairs for every mixture and fungal
# reads replace background rather than adding to it, so the zero control sits at
# exactly matched depth. The same background reads are used in all mixtures, so
# they differ only in the spike. seqkit head is deterministic and preserves
# pairing.
set -eu
BIN=$HOME/miniconda3/envs/assembly_env/bin
B=$HOME/fpar/reads/mock
F=$HOME/fpar/reads/recall
O=$HOME/fpar/reads/mix
G=$HOME/fpar/results/123_ground_truth.tsv
mkdir -p "$O"
echo "host: $(hostname)  started: $(date)"
command -v "$BIN/seqkit" >/dev/null || { echo "FATAL: no seqkit at $BIN"; exit 1; }

TOTAL=20000000
BACT="SRR11075447 ERR10009752 SRR21461810 DRR510841"
EUK="DRR188001 SRR14759986 SRR7762337 ERR3727302"
HUM="ERR11260221"
FUN="DRR289668 DRR163506 ERR12798583 DRR093944 DRR1004064 DRR829880 DRR584958"

# take n pairs from run $1 in dir $2, append to mixture prefix $3
take() {
  local run=$1 dir=$2 pfx=$3 n=$4
  [ "$n" -le 0 ] && return 0
  for r in 1 2; do
    "$BIN/seqkit" head -n "$n" "$dir/${run}_${r}.fastq.gz" 2>/dev/null | gzip -c >> "${pfx}_${r}.fastq.gz"
  done
}

printf "mixture\tsource\trun\tpairs_requested\n" > "$G"

for spec in "mixA 200000" "mixB 20000" "mixC 2000" "mixD 200" "mixE 20" "mixF 0"; do
  M=${spec%% *}; NF=${spec##* }
  P="$O/$M"
  rm -f "${P}_1.fastq.gz" "${P}_2.fastq.gz"
  NB=$(( TOTAL - NF ))
  nb_each=$(( NB * 85 / 100 / 4 ))
  ne_each=$(( NB * 12 / 100 / 4 ))
  nh=$(( NB - 4*nb_each - 4*ne_each ))
  echo "--- $M : fungal=$NF background=$NB (bact 4x$nb_each, euk 4x$ne_each, human $nh)  $(date +%H:%M)"
  for r in $BACT; do take "$r" "$B" "$P" "$nb_each"; printf "%s\tbacteria\t%s\t%s\n" "$M" "$r" "$nb_each" >> "$G"; done
  for r in $EUK;  do take "$r" "$B" "$P" "$ne_each"; printf "%s\teukaryote\t%s\t%s\n" "$M" "$r" "$ne_each" >> "$G"; done
  take "$HUM" "$B" "$P" "$nh"; printf "%s\thuman\t%s\t%s\n" "$M" "$HUM" "$nh" >> "$G"
  if [ "$NF" -gt 0 ]; then
    nf_each=$(( NF / 7 )); rem=$(( NF - 7*nf_each ))
    i=0
    for r in $FUN; do
      n=$nf_each; i=$((i+1)); [ "$i" -le "$rem" ] && n=$((n+1))
      take "$r" "$F" "$P" "$n"; printf "%s\tfungus\t%s\t%s\n" "$M" "$r" "$n" >> "$G"
    done
  fi
done

echo
echo "=== verification ==="
printf "%-7s %12s %12s %10s %12s %s\n" "mix" "pairs_R1" "pairs_R2" "paired?" "fungal_true" "fungal_fraction"
VOK=1
for M in mixA mixB mixC mixD mixE mixF; do
  P="$O/$M"
  n1=$(( $("$BIN/seqkit" stats -T "${P}_1.fastq.gz" | awk 'NR==2{print $4}') ))
  n2=$(( $("$BIN/seqkit" stats -T "${P}_2.fastq.gz" | awk 'NR==2{print $4}') ))
  # ground truth: count reads whose header accession is one of the 7 fungal runs
  ft=$(zcat "${P}_1.fastq.gz" | awk 'NR%4==1' | cut -c2- | cut -d. -f1 \
        | grep -cE '^(DRR289668|DRR163506|ERR12798583|DRR093944|DRR1004064|DRR829880|DRR584958)$' || true)
  pr="OK"; [ "$n1" != "$n2" ] && { pr="BROKEN"; VOK=0; }
  [ "$n1" != "$TOTAL" ] && VOK=0
  printf "%-7s %12d %12d %10s %12d %15.3e\n" "$M" "$n1" "$n2" "$pr" "$ft" "$(echo "scale=12; $ft/$n1" | bc)"
  printf "%s\tOBSERVED_TOTAL\t-\t%s\n%s\tOBSERVED_FUNGAL\t-\t%s\n" "$M" "$n1" "$M" "$ft" >> "$G"
done
echo
du -sh "$O"
echo "finished: $(date)"
[ "$VOK" == "1" ] && echo "MIXTURES BUILT AND VERIFIED" || echo "*** VERIFICATION FAILED -- do not use these mixtures ***"
