#!/bin/bash

# ==============================================================================
# 06b -- UMI deduplication on aligned BAMs
#
# Alithea's own Mercurius workflow (STARsolo) uses UMI-collapsed counts as
# the default. Our pipeline uses BRBseqTools + STAR + featureCounts, which
# preserves the 14 nt UMI in the read name but never uses it. This stage
# runs `umi_tools dedup` on each BAM so we can produce a parallel UMI-
# collapsed count matrix (see 06c_featureCounts_UMI.R).
#
# BRBseqTools writes UMIs at the end of the read name after the separator
# used below (default ':', overridable via UMI_SEPARATOR env var). The
# script prints a sample of read names first so you can eyeball that the
# separator is right for your data.
#
# Inputs:  hannah/alignments/*_Aligned.sortedByCoord.out.bam  (+ .bai)
# Outputs: hannah/alignments_dedup/*_dedup.bam                (+ .bai)
#          plus per-sample dedup logs & stats
#
# Cost: ~30-90 min for 384 BAMs at -P 8 parallelism.
# ==============================================================================

set -euo pipefail

baseDir="/rsstu/users/r/rrellan/sara/RNA_Sequencing_raw/BZea_CLY23D1/NVS205B_RellanAlvarez/hannah"
alignDir="$baseDir/alignments"
dedupDir="$baseDir/alignments_dedup"

UMI_SEPARATOR="${UMI_SEPARATOR:-:}"   # BRBseqTools default
PARALLEL="${PARALLEL:-8}"

module load conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /usr/local/usrapps/maize/hdpil/hdpil

# --- verify tools ---------------------------------------------------------
if ! command -v umi_tools >/dev/null; then
  echo "ERROR: umi_tools not on PATH." >&2
  echo "Install into the env with:" >&2
  echo "  conda install -n hdpil -c bioconda -c conda-forge umi_tools" >&2
  exit 1
fi
command -v samtools >/dev/null || { echo "ERROR: samtools missing"; exit 1; }

mkdir -p "$dedupDir"

# --- confirm UMI location in read names -----------------------------------
sample_bam=$(ls "$alignDir"/*_Aligned.sortedByCoord.out.bam | head -1)
echo "=== Sample read names from $(basename "$sample_bam") ==="
samtools view "$sample_bam" | head -3 | awk '{print $1}'
echo "UMI separator we'll use with umi_tools: '$UMI_SEPARATOR'"
echo "  (override with UMI_SEPARATOR=... if the character between the read"
echo "   ID and the UMI above is different)"
echo ""

# --- collect BAMs to process ---------------------------------------------
mapfile -t BAMS < <(ls "$alignDir"/*_Aligned.sortedByCoord.out.bam | sort)
n=${#BAMS[@]}
echo "$n BAMs to dedup, up to $PARALLEL in parallel"

# --- ensure .bai exists for every BAM (umi_tools dedup needs it) ---------
missing_idx=()
for b in "${BAMS[@]}"; do
  [ ! -f "${b}.bai" ] && [ ! -f "${b%.bam}.bai" ] && missing_idx+=("$b")
done
if [ "${#missing_idx[@]}" -gt 0 ]; then
  echo "Indexing ${#missing_idx[@]} BAMs..."
  printf '%s\n' "${missing_idx[@]}" \
    | xargs -n 1 -P "$PARALLEL" -I{} samtools index -@ 1 "{}"
fi

# --- one dedup task per BAM, parallel via xargs --------------------------
dedup_one() {
  local bam="$1"
  local sid; sid=$(basename "$bam" _Aligned.sortedByCoord.out.bam)
  local out="$dedupDir/${sid}_dedup.bam"
  local log="$dedupDir/${sid}_dedup.log"

  [ -f "$out" ] && [ -f "${out}.bai" ] && { echo "skip $sid (already done)"; return; }

  umi_tools dedup \
    --stdin="$bam" \
    --stdout="$out" \
    --umi-separator="$UMI_SEPARATOR" \
    --log="$log" \
    --output-stats="$dedupDir/${sid}_dedup_stats" \
    >/dev/null 2>&1

  samtools index "$out"
  echo "done $sid"
}
export -f dedup_one
export dedupDir UMI_SEPARATOR

printf '%s\n' "${BAMS[@]}" \
  | xargs -n 1 -P "$PARALLEL" -I{} bash -c 'dedup_one "$@"' _ {}

n_out=$(ls "$dedupDir"/*_dedup.bam 2>/dev/null | wc -l)
echo ""
echo "Dedup complete: $n_out / $n dedup BAMs at $dedupDir"
echo "Next: Rscript ../scripts/06c_featureCounts_UMI.R"
