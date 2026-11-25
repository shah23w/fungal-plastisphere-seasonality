#!/bin/bash
# Build the fungal-only reference index from the fungal genome assemblies listed
# in the PlusPF library report of the same build date, reusing that build's
# taxonomy and seqid2taxid. The fungal-only index therefore differs from PlusPF
# only in which non-fungal libraries are present, which isolates reference scope
# from every other classifier setting.
set -u
shopt -s nullglob
ARIA=$HOME/miniconda3/envs/assembly_env/bin/aria2c
SRC=$HOME/fpar/db/k2pf_full
DB=$HOME/fpar/db/k2fungi
mkdir -p "$DB/taxonomy" "$DB/library/fungi"

awk -F'\t' '$1=="fungi" && $3!=""{print $3}' "$SRC/library_report.tsv" \
  | sort -u > "$DB/fungal_urls.txt"
want=$(wc -l < "$DB/fungal_urls.txt")
echo "distinct fungal genome files to fetch: $want"

cp -n "$SRC/names.dmp" "$SRC/nodes.dmp" "$DB/taxonomy/" 2>/dev/null
cd "$DB/library/fungi" || exit 1

if [ -x "$ARIA" ]; then
  "$ARIA" -i "$DB/fungal_urls.txt" -j 3 -x 2 -s 2 --continue=true \
          --retry-wait=5 --max-tries=8 --auto-file-renaming=false \
          --console-log-level=warn --summary-interval=0
  echo "aria2 exit $?"
else
  wget -c -nv -i "$DB/fungal_urls.txt"
  echo "wget exit $?"
fi

got=( *.fna.gz )
n=${#got[@]}
echo "downloaded $n / $want files   total $(du -sh . | cut -f1)"
bad=0
for f in "${got[@]}"; do gzip -t "$f" 2>/dev/null || { echo "CORRUPT: $f"; bad=$((bad+1)); }; done
echo "gzip integrity failures: $bad"
if [ "$n" -eq "$want" ] && [ "$bad" -eq 0 ]; then echo "FETCH OK"; else echo "FETCH INCOMPLETE -- rerun, --continue resumes"; fi
