#!/bin/bash

# ==============================================================================
# FBX Analysis 6 (LOCAL) — TE annotation of fbxl1 intron 4
# for Rubén reply memo §10 item 3.
#
# Question: what TE families sit in the ~5 kb B73-lineage insertion in
# intron 4 of Zm00001eb375600 (fbxl1)? Intron-4 span is chr9:17925728-
# 17933692 (7,965 bp) — derived from the T001 transcript exon boundaries
# in Zea_mays.gtf.
#
# Uses the MaizeGDB EDTA-based B73 v5 TE annotation. Downloads once into
# data/external/ (gitignored) if not present, then extracts the intron-4
# overlap into data/FBX_intron4_TEs.gff3 (tracked).
# ==============================================================================

set -euo pipefail

TE_URL="https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz"
TE_FILE="data/external/Zm-B73-REFERENCE-NAM-5.0.TE.gff3.gz"
OUT_GFF="data/FBX_intron4_TEs.gff3"

# fbxl1 intron 4 coords (T001 transcript, chr9, + strand)
INTRON_CHR="chr9"
INTRON_START=17925728
INTRON_END=17933692

mkdir -p data/external

if [ ! -f "$TE_FILE" ]; then
  echo "Downloading B73 v5 TE annotation (~39 MB)..."
  curl -s -L -o "$TE_FILE" "$TE_URL"
fi

echo "Intersecting $INTRON_CHR:$INTRON_START-$INTRON_END against TE annotation..."
zcat "$TE_FILE" \
  | awk -F'\t' -v chr="$INTRON_CHR" -v s=$INTRON_START -v e=$INTRON_END \
        '$1==chr && $4 <= e && $5 >= s' \
  > "$OUT_GFF"

n=$(wc -l < "$OUT_GFF")
echo "Wrote $OUT_GFF ($n records)"
echo ""
echo "=== Feature type summary ==="
awk -F'\t' '{print $3}' "$OUT_GFF" | sort | uniq -c
echo ""
echo "=== TE hits sorted by size (bp) ==="
awk -F'\t' '{
  len = $5 - $4 + 1
  # pull Name= and Classification= from the attributes
  name = ""; class = ""
  n_attr = split($9, aa, ";")
  for (i=1; i<=n_attr; i++) {
    if (aa[i] ~ /^Name=/)            { name  = substr(aa[i], 6) }
    if (aa[i] ~ /^Classification=/)  { class = substr(aa[i], 16) }
  }
  printf "%7d  %-30s  %-25s  %s\n", len, $3, class, name
}' "$OUT_GFF" | sort -rn
