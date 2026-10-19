#!/usr/bin/env bash
#BSUB -J demux[1-4]
#BSUB -n 2
#BSUB -R "span[hosts=1] rusage[mem=4GB]"
#BSUB -W 12:00
#BSUB -o logs/demux.%J.%I.out
#BSUB -e logs/demux.%J.%I.err
# Split each BRB-seq pool into per-sample cDNA (R2) FASTQs by exact R1 barcode match.
# Needs barcodes.tsv: pool<TAB>barcode<TAB>sample_id   (pool = pool_1 ... pool_4)
# Output: $OUT/demux/<sample_id>.fastq.gz  -> samplesheet rows with fq2 = NA (single-end)
set -euo pipefail; source config.sh
IN=/rsstu/users/r/rrellan/CERCA-Cold/PTxB73xBRBseq/trimmed
BARCODES=/rsstu/users/r/rrellan/CERCA-Cold/PTxB73xBRBseq/data/barcodes.tsv
BC_LEN=14                                   # EDIT to your kit's barcode length
pool=pool_$LSB_JOBINDEX
o=$OUT/demux/$pool; mkdir -p "$o"

paste <(zcat "$IN/${pool}_R1.fastq.gz" | paste - - - -) \
      <(zcat "$IN/${pool}_R2.fastq.gz" | paste - - - -) |
awk -F'\t' -v L="$BC_LEN" -v P="$pool" -v O="$o" '
  NR==FNR { if ($1==P) bc[$2]=$3; next }
  { b=substr($2,1,L); if (b in bc) { n[b]++; print $5"\n"$6"\n"$7"\n"$8 > (O"/"bc[b]".fastq") } else miss++ }
  END { for (b in bc) printf "%s\t%s\t%d\n", bc[b], b, n[b] > (O"/counts.tsv")
        printf "unmatched\t-\t%d\n", miss >> (O"/counts.tsv") }' \
  "$BARCODES" -

for f in "$o"/*.fastq; do gzip "$f"; done
