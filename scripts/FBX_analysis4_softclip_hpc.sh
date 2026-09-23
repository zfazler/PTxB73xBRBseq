#!/bin/bash

# ==============================================================================
# FBX Analysis 4 (HPC) — soft-clip and start-position audit in the CDS window
# for Rubén's reply memo §6.2:
#
# > Displaced UTR reads cluster against the 3' edge of the window and carry
# > soft clips; genuine CDS reads distribute evenly across it.
#
# For every primary mapped read overlapping the CDS window
# (chr9:17933693-17933902), extract:
#   - leftmost aligned 5' position
#   - leading soft-clip length (from CIGAR ^NS)
#   - trailing soft-clip length (from CIGAR NS$)
#   - read length
#
# Writes ONE TSV directly into hannah/BZeaBRBseq/data/ so a subsequent
# `git add data/FBX_softclip_reads.tsv && git push` from HPC ships the output
# to local via git (no WinSCP).
#
# Output:
#   $repoDir/data/FBX_softclip_reads.tsv
#     sample_id  pos  softL  softR  read_len
# ==============================================================================

set -euo pipefail

baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
repoDir="$baseDir/BZeaBRBseq"
alignDir="$baseDir/alignments"
outFile="$repoDir/data/FBX_softclip_reads.tsv"

REGION="chr9:17933693-17933902"

mapfile -t BAMS < <(ls "$alignDir"/*_Aligned.sortedByCoord.out.bam | sort)
n=${#BAMS[@]}
if [ "$n" -eq 0 ]; then
  echo "ERROR: no BAMs in $alignDir" >&2
  exit 1
fi
echo "Processing $n BAMs over CDS window $REGION..."

# header
printf "sample_id\tpos\tsoftL\tsoftR\tread_len\n" > "$outFile"

for b in "${BAMS[@]}"; do
  sid=$(basename "$b" _Aligned.sortedByCoord.out.bam)
  # -F 260 = drop unmapped (0x4) + secondary (0x100)
  samtools view -F 260 "$b" "$REGION" | awk -v sid="$sid" '
    {
      pos      = $4
      cigar    = $6
      read_len = length($10)
      softL = 0; softR = 0
      # leading soft clip (^[0-9]+S)
      if (match(cigar, /^[0-9]+S/)) {
        softL = substr(cigar, 1, RLENGTH - 1) + 0
      }
      # trailing soft clip (NS at end)
      n = length(cigar)
      if (substr(cigar, n, 1) == "S") {
        i = n - 1
        while (i > 0 && substr(cigar, i, 1) ~ /[0-9]/) i--
        softR = substr(cigar, i + 1, n - i - 1) + 0
      }
      print sid "\t" pos "\t" softL "\t" softR "\t" read_len
    }' >> "$outFile"
done

rows=$(wc -l < "$outFile")
echo "Done. Wrote $outFile ($rows rows including header)"
