#!/bin/bash
# Build the spiked libraries for the recovery experiment. Reads from three fungal
# donor genomes are simulated with wgsim at seven known read-pair counts and
# added to two deep, fungus-poor backgrounds. All three donors are added at each
# level because they are separable by taxon in the Kraken2 report, so one
# classification measures recovery for all three. Seeds are a deterministic
# function of background, level and donor.
set -u
ENV=$HOME/miniconda3/envs/assembly_env/bin
SEQKIT=$ENV/seqkit
WGSIM=$ENV/wgsim
LIB=$HOME/fpar/db/k2fungi/library/fungi
CORPUS=$HOME/fpar/corpus
W=$HOME/fpar/work_spike
mkdir -p "$W/donors" "$W/libs" "$W/sim"

NSUB=10000000
BGS="SRR30635593 SRR35440449"
LEVELS="2 10 50 250 1250 6250 31250"
DNAMES="Pichia Yarrowia Fusarium"
Pichia_FA=GCF_003054445.1_ASM305444v1_genomic.fna.gz
Yarrowia_FA=GCF_001761485.1_ASM176148v1_genomic.fna.gz
Fusarium_FA=GCF_000240135.3_ASM24013v3_genomic.fna.gz

for n in $DNAMES; do
  eval "src=\$${n}_FA"
  if [ ! -s "$W/donors/$n.fa" ]; then
    gunzip -c "$LIB/$src" > "$W/donors/$n.fa" || exit 1
  fi
  echo "donor $n  $(grep -c '^>' "$W/donors/$n.fa") sequences  $(du -h "$W/donors/$n.fa" | cut -f1)"
done

first=$(echo $BGS | awk '{print $1}')
RL=$(zcat "$CORPUS/$first/${first}_1.fastq.gz" | head -2 | tail -1 | tr -d '\n' | wc -c)
echo "background read length detected: $RL bp"
[ "$RL" -ge 50 ] && [ "$RL" -le 301 ] || { echo "implausible read length; stop"; exit 1; }

for bg in $BGS; do
  for r in 1 2; do
    out="$W/libs/${bg}_sub_${r}.fq.gz"
    if ! gzip -t "$out" 2>/dev/null; then
      rm -f "$out"
      "$SEQKIT" head -n "$NSUB" "$CORPUS/$bg/${bg}_${r}.fastq.gz" -o "$out" 2>/dev/null
      gzip -t "$out" 2>/dev/null || { echo "FATAL: subsample $out failed integrity"; exit 1; }
    fi
  done
  echo "background $bg subsampled: $(du -sh "$W/libs/${bg}_sub_1.fq.gz" | cut -f1)"
  c1=$(( $(zcat "$W/libs/${bg}_sub_1.fq.gz" | wc -l) / 4 ))
  c2=$(( $(zcat "$W/libs/${bg}_sub_2.fq.gz" | wc -l) / 4 ))
  echo "   subsample pairs: R1 $c1  R2 $c2  (want $NSUB)"
  [ "$c1" -eq "$c2" ] && [ "$c1" -eq "$NSUB" ] || { echo "FATAL: subsample count wrong"; exit 1; }
done

TRUTH=$HOME/fpar/results/166_spikein_truth.tsv
echo -e "library\tbackground\tlevel_pairs_per_donor\tdonor\tinjected_pairs\tbackground_pairs\texpected_rpm_per_donor\texpected_rpm_total" > "$TRUTH"

seed=1000
for bg in $BGS; do
  for L in $LEVELS; do
    tag="${bg}_L${L}"
    o1="$W/libs/${tag}_1.fq.gz"; o2="$W/libs/${tag}_2.fq.gz"
    if gzip -t "$o1" 2>/dev/null && gzip -t "$o2" 2>/dev/null; then echo "$tag exists, skip"; else
      cp "$W/libs/${bg}_sub_1.fq.gz" "$o1.tmp.gz"
      cp "$W/libs/${bg}_sub_2.fq.gz" "$o2.tmp.gz"
      for n in $DNAMES; do
        bi=0; for q in $BGS; do [ "$q" = "$bg" ] && break; bi=$((bi+1)); done
        di=0; for q in $DNAMES; do [ "$q" = "$n" ] && break; di=$((di+1)); done
        seed=$((7000 + bi*1000 + L + di))
        "$WGSIM" -N "$L" -1 "$RL" -2 "$RL" -e 0.005 -r 0 -R 0 -X 0 -S "$seed" \
                 -d 350 -s 50 "$W/donors/$n.fa" \
                 "$W/sim/${tag}_${n}_1.fq" "$W/sim/${tag}_${n}_2.fq" > /dev/null 2>&1
        got=$(( $(wc -l < "$W/sim/${tag}_${n}_1.fq") / 4 ))
        [ "$got" -eq "$L" ] || echo "   WARN $tag $n: wgsim made $got pairs, wanted $L"
        gzip -c "$W/sim/${tag}_${n}_1.fq" >> "$o1.tmp.gz"
        gzip -c "$W/sim/${tag}_${n}_2.fq" >> "$o2.tmp.gz"
        rm -f "$W/sim/${tag}_${n}_1.fq" "$W/sim/${tag}_${n}_2.fq"
      done
      mv "$o1.tmp.gz" "$o1"; mv "$o2.tmp.gz" "$o2"
      echo "built $tag  $(du -h "$o1" | cut -f1)"
    fi
    tot=$((NSUB + 3*L))
    for n in $DNAMES; do
      rpm=$(awk -v l="$L" -v t="$tot" 'BEGIN{printf "%.4f", l/(t/1000000)}')
      rpmt=$(awk -v l="$L" -v t="$tot" 'BEGIN{printf "%.4f", 3*l/(t/1000000)}')
      echo -e "${tag}\t${bg}\t${L}\t${n}\t${L}\t${NSUB}\t${rpm}\t${rpmt}" >> "$TRUTH"
    done
  done
done

echo
echo "=== verification ==="
n=$(ls -1 "$W/libs"/*_L*_1.fq.gz 2>/dev/null | wc -l)
echo "spiked libraries built: $n  (expect 14)"
bad=0
for f in "$W/libs"/*_L*_[12].fq.gz; do gzip -t "$f" 2>/dev/null || { echo "CORRUPT $f"; bad=$((bad+1)); }; done
echo "gzip integrity failures: $bad"
echo "pair-count check on one library:"
t=$(ls -1 "$W/libs"/*_L31250_1.fq.gz | head -1)
t2=${t%_1.fq.gz}_2.fq.gz
a=$(( $(zcat "$t" | wc -l) / 4 )); b=$(( $(zcat "$t2" | wc -l) / 4 ))
echo "   $(basename "$t"): R1 $a pairs, R2 $b pairs, expect $((NSUB + 3*31250))"
[ "$a" -eq "$b" ] && echo "   R1/R2 balanced: OK" || echo "   R1/R2 MISMATCH"
echo "truth table rows: $(( $(wc -l < "$TRUTH") - 1 ))  (expect 42)"
echo "total size: $(du -sh "$W/libs" | cut -f1)"
[ "$n" -eq 14 ] && [ "$bad" -eq 0 ] && [ "$a" -eq "$b" ] && echo "BUILD OK" || echo "BUILD INCOMPLETE"
