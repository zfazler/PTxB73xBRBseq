#!/bin/bash

# ==============================================================================
# FBX Analysis 3 (HPC) — per-base coverage across fbxl1 terminal exon
# for the per-base coverage-profile test in Rubén's reply memo §5.
#
# Runs `samtools depth -a` on every sample BAM across a 660 bp window that
# spans upstream of the CDS window through past the poly(A) site:
#
#   REGION            chr9:17933600-17934260   (660 bp on + strand)
#   CDS window        chr9:17933693-17933902   (210 bp)
#   3'UTR window      chr9:17933903-17934180   (278 bp)
#   B73 poly(A) site  chr9:17934180
#
# samtools depth called ONCE with all BAMs -> single tab-sep matrix
# (chr, pos, depth_sample1, depth_sample2, ...). No per-sample intermediates.
#
# Writes directly INTO the cloned repo (hannah/BZeaBRBseq/data/) so that
# `git add data/FBX_depth_matrix.tsv && git commit && git push` from HPC
# ships the output to local via git (no WinSCP).
#
# Output:
#   $repoDir/data/FBX_depth_matrix.tsv    (660 rows) x (2 + N samples) columns
# ==============================================================================

set -euo pipefail

baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
repoDir="$baseDir/BZeaBRBseq"
alignDir="$baseDir/alignments"
outDir="$repoDir/data"
mkdir -p "$outDir"

REGION="chr9:17933600-17934260"

# gather all sample BAMs
mapfile -t BAMS < <(ls "$alignDir"/*_Aligned.sortedByCoord.out.bam | sort)
n=${#BAMS[@]}
if [ "$n" -eq 0 ]; then
  echo "ERROR: no BAMs found in $alignDir" >&2
  exit 1
fi
echo "Depth over $REGION for $n BAMs..."

# bump open-file limit; 384 BAMs + stderr/stdout is fine under 4096
ulimit -n 4096 || true

# ---- index any BAMs missing .bai (samtools depth -r needs an index) ----
missing_idx=()
for b in "${BAMS[@]}"; do
  if [ ! -f "${b}.bai" ] && [ ! -f "${b%.bam}.bai" ]; then
    missing_idx+=("$b")
  fi
done
if [ "${#missing_idx[@]}" -gt 0 ]; then
  echo "Indexing ${#missing_idx[@]} BAMs (missing .bai)..."
  # xargs -P for parallelism; samtools index is fast on sorted BAMs
  printf '%s\n' "${missing_idx[@]}" \
    | xargs -n 1 -P 4 -I{} samtools index -@ 1 "{}"
  echo "Indexing done."
else
  echo "All BAMs already indexed."
fi

out="$outDir/FBX_depth_matrix.tsv"
{
  # header row: chr, pos, then one column per sample
  printf "chr\tpos"
  for b in "${BAMS[@]}"; do
    n=$(basename "$b" _Aligned.sortedByCoord.out.bam)
    printf "\t%s" "$n"
  done
  printf "\n"
  samtools depth -a -r "$REGION" "${BAMS[@]}"
} > "$out"

rows=$(wc -l < "$out")
cols=$(head -n 1 "$out" | awk -F'\t' '{print NF}')
echo "Done. Wrote $out ($rows rows x $cols cols)"
